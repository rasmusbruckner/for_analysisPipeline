% This is a simple illustration of the likelihood function 
%
% This can be extended to other parameters.
% 
% 1. Load data
% 2. Simulate data
% 3. Evaluate likelihood with different hazard rates


% Identify parent directory of this config script
parentDirectory = fileparts(fileparts(mfilename('fullpath')));
cd(parentDirectory)
addpath(genpath(parentDirectory));

% This is the BIDS folder
bidsDir = strcat(parentDirectory, filesep, 'for_data', filesep, 'for_bids_data');

% ------------
% 1. Load data
% ------------

% Run preprocessing to get all behavioral data
allSubBehavData = for_preprocessing(bidsDir);

% Number of subjects
n_subj = length(unique(allSubBehavData.subj_num));

% ----------------
% 2. Simulate data
% ----------------

% Initialize agent variables
agent_vars = ForAgentVars();
agent_vars.max_x = 2 * pi;

% Initialize simulation variables
sim_est_err = NaN(n_subj, 1);
true_params = NaN;

% Independent model parameter initialization
df_model = table();
df_model.omikron_0 = 10;
df_model.omikron_1 = 0; 
df_model.h = 0.5;
df_model.s = 1;
df_model.u = 0;
df_model.sigma_H = 0;
df_model.subj_num = 1;
sub = 1;

% Select current subject
subsetIdx = allSubBehavData.subj_num == sub;

% Apply filtering to each field of the struct
df_subj = structfun(@(x) x(subsetIdx, :), allSubBehavData, 'UniformOutput', false);

% Extract model parameters
sel_coeffs = df_model(df_model.subj_num == sub, :);

% Set agent variables
agent_vars.h = sel_coeffs.h;
agent_vars.s = sel_coeffs.s;
agent_vars.u = exp(sel_coeffs.u);
agent_vars.sigma_H = sel_coeffs.sigma_H;
agent_vars.tau_0 = 0.999;
agent_vars.sigma_0 = 6.1875;
agent_vars.mu_0 = pi;
agent_vars.max_x = 2 * pi;

% Create agent instance
agent = AlAgentRBM(agent_vars);

% Create simple array for RBM
sel_coeffs = table2array(sel_coeffs);

% Simulation
sim = true;
[~, df_sim] = for_task_agent_int(df_subj, agent, agent_vars, sel_coeffs, sim);

% --------------------------------------------------
% 3. Evaluate likelihood with different hazard rates
% --------------------------------------------------

% Turn off simulation
sim = false;

% Initialize variables
h = linspace(0,1,100);
llh_mix = nan(1, 100);

% Cycle over subjects
for i = 1:length(h)

    % Use different hazard rate for each iteration
    df_model.h = h(i);

    % Extract model parameters
    sel_coeffs = df_model(df_model.subj_num == sub, :);
    agent_vars.h = sel_coeffs.h;

    % Create simple array for RBM
    sel_coeffs = table2array(sel_coeffs);

    % Create agent instance
    agent = AlAgentRBM(agent_vars);

    % Simulation
    [llh_rbm, df_data] = for_task_agent_int(df_sim, agent, agent_vars, sel_coeffs, sim);
    
    % Sum negative log-likelihodds
    llh_mix(i) = -1 * sum(llh_rbm);

end

% Plot results
figure()
plot(h, llh_mix)