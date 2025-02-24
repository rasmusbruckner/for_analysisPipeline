function allSubBehavData = for_preprocessing(bidsDir)
%FOR_PREPROCESSING This function preprocesses FOR data based on BIDS data
%
%   Input
%       bidsDir: BIDS directory path
%
%   Output
%       allSubBehavData: Struct with all subject data


% Extract folders and data filnenames
whichFolders = dir(fullfile(bidsDir));
allFilenames = whichFolders(contains({whichFolders.name}, 'sub_'));

% Cycle ove filenames
for j = 1:length(allFilenames)

    % Construct the full path to the current folder
    currentPath = fullfile(allFilenames(j).folder, allFilenames(j).name, 'behav');

    % Get all files in the folder and find the first .tsv file
    fileNames = {dir(currentPath).name};
    tsvFile = fullfile(currentPath, fileNames{contains(fileNames, '.tsv')});

    % Read the .tsv file
    behavData = readtable(tsvFile, 'FileType', 'text', 'Delimiter', '\t');
    behavData = table2struct(behavData, 'ToScalar', true);

    % Recode variables to radians
    behavData.x_t = deg2rad(behavData.x_t);
    behavData.b_t = deg2rad(behavData.b_t);
    behavData.mu_t = deg2rad(behavData.mu_t);
    behavData.delta_t = deg2rad(behavData.delta_t);
    behavData.e_t = deg2rad(behavData.e_t);
    behavData.a_t = deg2rad(behavData.a_t);

    % Identify new blocks
    newB = find(behavData.new_block(1:end));

    % Compute circular PE and UP in radians to check task-based value
    PE = circ_dist(behavData.x_t, behavData.b_t);
    UP = [circ_dist(behavData.b_t(2:end), behavData.b_t(1:end-1)); nan];
    UP(newB(2:end)-1) = nan; % ensure that update is nan on last trial

    % Check prediction error
    if any(PE - behavData.delta_t > 1.e-5)
        error('Issue with prediction error');
    end

    % Realign update transformed to radians in order to match index of prediction error
    a_t = nan(length(behavData.a_t), 1);
    a_t(1:end-1) = behavData.a_t(2:end);
    a_t(newB(2:end)-1) = nan; % ensure that update is nan on last trial
    behavData.a_t = a_t;

    % Check update
    if any(UP - a_t > 1.e-5)
        error('Issue with update');
    end

    % Compute estimation error manually for checking
    estErr = abs(circ_dist(behavData.mu_t, behavData.b_t));

    % Check estimation error
    if any(abs(estErr) - abs(behavData.e_t) > 1.e-5)
        error('Issue with estimation error');
    end

    % Translate van Mises concentration into Gaussian standard deviation
    vmVar = 1./behavData.kappa_t;
    gaussStd = vmVar.^.5;

    % Add standard deviation to behavioral structure for modeling
    behavData.sigma = gaussStd;

    % Add sigma, hit, and catch-trial-dummy variables for regression
    behavData.sigma_dummy = nan(length(behavData.kappa_t),1);
    behavData.sigma_dummy(behavData.kappa_t == 8) = 1;
    behavData.sigma_dummy(behavData.kappa_t == 16) = 0;
    behavData.hit_dummy = nan(length(behavData.kappa_t),1);
    behavData.hit_dummy(behavData.r_t == 1) = 1;
    behavData.hit_dummy(behavData.r_t == 0) = 0;
    behavData.visible_dummy = nan(length(behavData.v_t), 1);
    behavData.visible_dummy(behavData.v_t == 1) = 1;
    behavData.visible_dummy(behavData.v_t == 0) = 0;

    % Combine data sets of all subjects
    if j == 1
        allSubBehavData = behavData;
    else
        allSubBehavData = catStruct(allSubBehavData, behavData);
    end
end

% Get the parent directory
[parentDir, ~, ~] = fileparts(bidsDir);

% Save newly created struct
save([parentDir filesep 'allSubBehavData.mat'], '-struct', 'allSubBehavData');

end