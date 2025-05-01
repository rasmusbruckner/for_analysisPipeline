% Parameter recover regression model
%
% 1. Simulate data for recovery
% 2. Estimate regression model
% 3. Plot correlations

% Number of random starting points for regression estimation
n_sp = 15;

% Number of simulations for recovery
n_subj = 100;

% Identify parent directory of this config script
parentDirectory = fileparts(mfilename('fullpath'));
cd(parentDirectory)
addpath(genpath(parentDirectory));

% -----------------------------
% 1. Simulate data for recovery
% -----------------------------

% Initialize regression variables
reg_vars = ForRegVars();
reg_vars.n_subj = n_subj;
reg_vars.n_sp = n_sp;
reg_vars.usePrior = true;

% Determine which parameters should be estimated
reg_vars.which_vars.beta_0 = true; % intercept
reg_vars.which_vars.beta_1 = true; % PE (fixed learning rate)
reg_vars.which_vars.beta_2 = true; % interaction PE and RU
reg_vars.which_vars.beta_3 = true; % interaction PE and CPP
reg_vars.which_vars.beta_4 = true; % interaction PE and hit
reg_vars.which_vars.beta_5 = true; % interaction PE and noise condition
reg_vars.which_vars.beta_6 = false; % interaction PE and visible
reg_vars.which_vars.beta_7 = false; % interaction EE and visible
reg_vars.which_vars.omikron_0 = true; % motor noise (independent of UP)
reg_vars.which_vars.omikron_1 = true; % learning-rate noise (dependent on UP)
reg_vars.which_vars.uniform = false; % uniform component for outlier predictions
reg_vars.regressionComponents = [reg_vars.which_vars.beta_0, reg_vars.which_vars.beta_1,...
    reg_vars.which_vars.beta_2, reg_vars.which_vars.beta_3, reg_vars.which_vars.beta_4,...
    reg_vars.which_vars.beta_5, reg_vars.which_vars.beta_6, reg_vars.which_vars.beta_7];

% Create regression-object instance
regression = ForRegression(reg_vars);

% Sample random model parameters that we try to recover
df_params = table();

if reg_vars.which_vars.beta_0
    df_params.beta_0 = unifrnd(-0.1,0.1, n_subj, 1);
end

if reg_vars.which_vars.beta_1
    df_params.beta_1 = rand(n_subj, 1);
end

if reg_vars.which_vars.beta_2
    df_params.beta_2 = rand(n_subj,1);
end

if reg_vars.which_vars.beta_3
    df_params.beta_3 = rand(n_subj,1);
end

if reg_vars.which_vars.beta_4
    df_params.beta_4 = rand(n_subj,1);
end

if reg_vars.which_vars.beta_5
    df_params.beta_5 = unifrnd(-0.1,0.1, n_subj, 1);
end

if reg_vars.which_vars.beta_6
    df_params.beta_6 = unifrnd(-0.1,0.1, n_subj, 1);
end

if reg_vars.which_vars.beta_7
    df_params.beta_7 = unifrnd(-0.1,0.1, n_subj, 1);
end

df_params.omikron_0 = unifrnd(3, 10, n_subj,1);

if reg_vars.which_vars.omikron_1
    df_params.omikron_1 = rand(n_subj, 1) * 0.3;
end

df_params.subj_num = (1:n_subj)';

% Simulate updates based on sampled parameters
n_trials = 400;
samples = regression.sample_data(df_params, n_trials);

% ----------------------------
% 2. Estimate regression model
% ----------------------------

% Translate table to structure
samplesStruct = table2struct(samples, 'ToScalar', true);

% Estimate regression model
results = regression.run_estimation(samplesStruct);

% --------------------
% 3. Plot correlations
% --------------------

behavLabels = {reg_vars.beta_0, reg_vars.beta_1, reg_vars.beta_2,...
    reg_vars.beta_3, reg_vars.beta_4, reg_vars.beta_5,...
    reg_vars.beta_6, reg_vars.beta_7, reg_vars.omikron_0,...
    reg_vars.omikron_1};

whichParamsVec = struct2array(reg_vars.which_vars);
behavLabels = behavLabels(whichParamsVec);
gridSize = [3,3];
for_recoverySummary(df_params, results, behavLabels, gridSize)