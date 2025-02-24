function for_recoverySummary(trueParams, estParams, behavLabels, gridSize)
% FOR_RECOVERYSUMMARY This function creates a simple plot showing
% parameter values
%
%   Input
%       trueParams: True parameter values
%       estParams: Estimated parameter values
%       behavLabels: Labels for plot
%       gridSize: Plot grid
%
%   Output
%       None

% Create figure
figure();

% Cycle over parameters
for i = 1:size(estParams.parameters, 2)

    % Create subplot
    subplot(gridSize(1), gridSize(2), i);
    hold on

    % Extract parameter values
    trueParamValue = trueParams.(behavLabels{i});
    estParamValue = estParams.parameters.(behavLabels{i});

    % Plot parameters
    plot(trueParamValue, estParamValue, 'o')
    r = corr(trueParamValue, estParamValue, 'type','Spearman');
    title([behavLabels{i} ': r=' num2str(round(r, 2))])

    % Add axis labels
    xlabel('True parameter')
    ylabel('Estimated parameter')

end
end