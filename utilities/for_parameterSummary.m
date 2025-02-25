function for_parameterSummary(parameters, behavLabels, gridSize)
% FOR_PARAMETERSUMMARY This function creates a simple plot showing
% parameter values
%
%   Input
%       parameters: All parameters
%       behavLabels: Labels for plot
%       gridSize: Plot grid
%
%   Output
%       None

% Translate to array
parameters = table2array(parameters);

% Variability of coefficients for jittering
varNormCoeffs = parameters(:,1:end)./repmat((std(parameters(:,1:end))), size(parameters(:,1:end), 1), 1);

% Number of subjects
nSubj = size(varNormCoeffs, 1);

% Create figure
figure();

% Get smart jitter values
xJit = smartJitter(varNormCoeffs,.01, 1);

% Marker size
markerSize = 2;

% Cycle over coefficients
for i = 1:size(varNormCoeffs, 2)

    % Create subplot
    subplot(gridSize(1), gridSize(2), i);
    hold on

    % Barplot
    bar(i, [mean(parameters(:,i))], 'FaceColor', 'b', 'facealpha', 0.5)

    % Jitter
    plot(ones(nSubj, 1).*i+xJit(:,i), parameters(:,i), 'o',...
        'markerSize', markerSize, 'markerFaceColor','b',...
        'markerEdgeColor', 'k', 'lineWidth', 1);

    % Add y-label
    ylabel(behavLabels{i})

end
end

