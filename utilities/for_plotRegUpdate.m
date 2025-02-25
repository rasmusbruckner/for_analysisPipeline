function for_plotRegUpdate(allSubBehavData,samples, ID)
%FOR_PLOTREGUPDATE This function compares updates of participants and
% predictions of the regression model
%
%   Input
%       allSubBehavData: All participant data
%       samples: Samples from the regression model
%       ID: If given, only specified subject will be plotted
%
%   Output
%       None

% Check if ID is provided
if ~exist('ID', 'var') || isempty(ID)
    realUp = allSubBehavData.a_t;
    predUp = samples.a_t;
else
    realUp = allSubBehavData.a_t(allSubBehavData.ID == ID);
    predUp = samples.a_t(allSubBehavData.ID == ID);
end

% Create figure
figure()
hold on

% Plot histograms
histogram(realUp, "facecolor", "k", "facealpha", 0.8, 'BinWidth', 0.25, 'Normalization', 'pdf')
histogram(predUp, "facecolor", "b", "facealpha", 0.5, 'BinWidth', 0.25, 'Normalization', 'pdf')

% Estimate and plot distributions
[f, xi] = ksdensity(realUp); % Kernel Density Estimation
plot(xi, f, 'r', 'LineWidth', 2, "color", "k"); % Overlay KDE curve
[f, xi] = ksdensity(predUp); % Kernel Density Estimation
plot(xi, f, 'r', 'LineWidth', 2, "color", "b"); % Overlay KDE curve

% Add legend
legend(["Participants", "Model predictions"])

end