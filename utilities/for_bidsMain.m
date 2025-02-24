function for_bidsMain(dataDir, bidsDir)
% FOR_BIDSMAIN This function converts cannon data to BIDS format
%
%   Input
%       dataDir: Directory with raw data
%       bidsDir: Directory for BIDS data
%
%   Output
%       None


% Add data directory and jsonlab library
addpath(genpath(dataDir));
addpath(genpath('JSONLAB'));

% Extract folders and data filnenames
whichFolders = dir(fullfile(dataDir));
allFilenames = whichFolders(~ismember({whichFolders.name}, {'.','..', '.DS_Store'}));

% Cycle over filenames
for j = 1:length(allFilenames)

    % Go into subject folder and extract data files
    whichFiles = dir(fullfile(dataDir, allFilenames(j).name));
    subjFilenames = whichFiles(~ismember({whichFiles.name},{'.','..', '.DS_Store'}));

    % Cycle over subject filenames
    for i = 1:length(subjFilenames)

        % Select current subject
        currentFilename = subjFilenames(i).name;

        % Load data
        blockData = load(currentFilename);
        blockData = blockData.taskData;

        % Recode 360 outcome to 0
        blockData.outcome(blockData.outcome == 360) = 0;

        % Add subject number
        blockData.ID = repmat(j, length(blockData.actRew), 1);

        % Combine data sets of the subject
        if i == 1
            subData = blockData;
        else
            subData = catStruct(subData, blockData);
        end
    end

    % Run BIDS conversion for current subject
    % ---------------------------------------

    % Data directory with variable subject string
    if j<10
        sub_string = 'sub_0';
    else
        sub_string = 'sub_';
    end

    % Full data directory
    data_dir = fullfile(bidsDir, [sub_string num2str(j)], 'behav');

    % Check if we need a new directory
    if exist(data_dir, 'dir') == false
        mkdir(data_dir)
    end

    % Extract variables of interest
    ID = subData.ID;
    block = subData.block;
    x_t = subData.outcome;
    b_t = subData.pred;
    delta_t = subData.predErr;
    a_t = subData.UP;
    e_t = subData.estErr;
    mu_t = subData.distMean;
    c_t = subData.cp;
    tac = subData.TAC;
    r_t = subData.hit;
    kappa_t = subData.concentration;
    v_t = subData.catchTrial;
    RT = subData.RT;
    initRT = subData.initiationRTs;

    % Compute new-block index
    new_block = [true; diff(subData.block) ~= 0];

    % Combine variables of interest in events file
    events_t = table(ID, block, new_block, x_t, b_t, delta_t, a_t, e_t, mu_t, c_t, tac, r_t, kappa_t, v_t, RT, initRT);

    % Events.csv and .tsv filenames
    events_csv = fullfile(data_dir, [[sub_string, num2str(j)] '_task-cannon_behav.csv']);
    events_tsv = fullfile(data_dir, [[sub_string, num2str(j)] '_task-cannon_behav.tsv']);

    % Check if files exists
    if exist(events_tsv, 'file')

        % If so, load load file for comparison
        existingBidsFile = readtable(events_tsv, 'FileType', 'text', 'Delimiter', '\t');
        createNewFile = false;

        % Check if current and new file are equal
        isEqual = safeSave(existingBidsFile, events_t);

    else
        createNewFile = true;
    end


    if createNewFile

        % Create new file
        writetable(events_t, events_csv,'Delimiter', '\t'); % as .csv
        copyfile(events_csv, events_tsv); % as .tsv
        delete(events_csv); % delete .csv

        % Inform user that new file is being created
        fprintf('%s%i: BIDS file did not exist. Creating it now.\n', sub_string, j);

    elseif isEqual == false

        % Create file with different name to avoid accidental overwriting
        events_csv = fullfile(data_dir, [[sub_string, num2str(j)] '_task-cannon_behav_new.csv']);
        events_tsv = fullfile(data_dir, [[sub_string, num2str(j)] '_task-cannon_behav_new.tsv']);
        writetable(events_t, events_csv,'Delimiter', '\t'); % as .csv
        copyfile(events_csv, events_tsv); % as .tsv
        delete(events_csv); % delete .csv

        % Inform user that new file with new name is being created
        fprintf('%s%i: BIDS file exists but incompatible with new one! Saving with suffix "_new".\n', sub_string, j)

    else

        % Inform user that nothing was updated
        fprintf('%s%i: BIDS file exists and is consistent. No need to update.\n', sub_string, j);

    end

    % Add metadata
    % ------------

    % Events description
    metadata_beh = for_bidsEventsDescr();

    % Save .json metadata file
    savejson('', metadata_beh, fullfile(data_dir, [[sub_string, num2str(j)] '_task-cannon_behav.json']));

end

% Add supplementary information function later
% for_bidsSupplnf(bidsDir)

    function isEqual = safeSave(existingTable, newTable)
        %SAFESAVE This function compares two tables to ensure that the
        % new table is saved safely
        %
        %   The function compares an existing and new structure and checks
        %   if they are equal.
        %
        %   Input
        %       existingTable: Current structure
        %       newTable: New structure
        %
        %   Output
        %       isEqual: Binary variable indicating equality


        % Define a tolerance
        tolerance = 1e-6;

        % First, check if tables have the same size and variable names
        if isequal(size(existingTable), size(newTable)) && isequal(existingTable.Properties.VariableNames, newTable.Properties.VariableNames)

            % Assume tables contain only numeric data; handle non-numeric columns separately if needed
            numericColumns = varfun(@isnumeric, existingTable, 'OutputFormat', 'uniform');
            isEqual = true;

            % Loop through each numeric column to compare within tolerance
            for c = find(numericColumns)

                % Get the data from the corresponding columns
                col1 = existingTable{:, c};
                col2 = newTable{:, c};

                % Compare the columns element-wise, accounting for NaNs and tolerance
                isEqualColumn = isequaln(col1, col2) || all(abs(col1 - col2) < tolerance | (isnan(col1) & isnan(col2)));

                if ~isEqualColumn
                    isEqual = false;
                    break;
                end
            end
        else
            isEqual = false;
        end
    end
end


