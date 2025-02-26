classdef ForEstimation
    % FORESTIMATION This class specifies the instance variables and methods
    % for RBM parameter estimation
    
    properties

        % Parameter names
        omikron_0
        omikron_1
        h
        s
        u
        sigma_H
        
        % Fixed starting points
        omikron_0_x0
        omikron_1_x0
        h_x0
        s_x0
        u_x0
        sigma_H_x0
        
        % Range for random starting points
        omikron_0_x0_range
        omikron_1_x0_range
        h_x0_range
        s_x0_range
        u_x0_range
        sigma_H_x0_range

        % Estimation lower boundaries
        lowerBoundaries

        % Estimation upper boundaries
        upperBoundaries 
        
        % Free parameter indexes
        which_vars
        
        % Fixed parameter values
        fixed_mod_coeffs
        
        % Other properties
        n_subj % number of subjects
        rand_sp % random starting points
        n_sp % number of starting points
        use_prior % use prior for regularization

    end
    
    methods

        function obj = ForEstimation(est_vars)
            % Constructor to initialize estimation variables

            % Parameter names
            obj.omikron_0 = est_vars.omikron_0;
            obj.omikron_1 = est_vars.omikron_1;
            obj.h = est_vars.h;
            obj.s = est_vars.s;
            obj.u = est_vars.u;
            obj.sigma_H = est_vars.sigma_H;
            
            % Fixed starting points
            obj.omikron_0_x0 = est_vars.omikron_0_x0;
            obj.omikron_1_x0 = est_vars.omikron_1_x0;
            obj.h_x0 = est_vars.h_x0;
            obj.s_x0 = est_vars.s_x0;
            obj.u_x0 = est_vars.u_x0;
            obj.sigma_H_x0 = est_vars.sigma_H_x0;
            
            % Range for random starting points
            obj.omikron_0_x0_range = est_vars.omikron_0_x0_range;
            obj.omikron_1_x0_range = est_vars.omikron_1_x0_range;
            obj.h_x0_range = est_vars.h_x0_range;
            obj.s_x0_range = est_vars.s_x0_range;
            obj.u_x0_range = est_vars.u_x0_range;
            obj.sigma_H_x0_range = est_vars.sigma_H_x0_range;
            
            % Estimation lower boundaries
            obj.lowerBoundaries = est_vars.lowerBoundaries;

            % Estimation upper boundaries
            obj.upperBoundaries = est_vars.upperBoundaries;
            
            % Parameters that should be estimated
            obj.which_vars = est_vars.which_vars;
            
            % Fixed parameter values
            obj.fixed_mod_coeffs = est_vars.fixed_mod_coeffs;
            
            % Extract other properties
            obj.n_subj = est_vars.n_subj;
            obj.rand_sp = est_vars.rand_sp;
            obj.n_sp = est_vars.n_sp;
            obj.use_prior = est_vars.use_prior;

        end
        
        function results = run_estimation(obj, allSubBehavData, agent_vars)
            % RUN_ESTIMATION This function implements the parallel
            % estimation of the model parameters
            % 
            %   Input
            %       obj: Estimation-object instance
            %       allSubBehavData: Behavioral data
            %       agent_vars: Agent-variables-object instance
            % 
            %   Output
            %       results: Structure containing estimation results
            
            % Free parameters
            which_vars_vec = struct2array(obj.which_vars);

            % Initialize
            results = struct(); % results struct
            results.parameters = nan(obj.n_subj, sum(which_vars_vec)); % estimation coefficients
            results.llh = nan(obj.n_subj,1); % negative log-likelihood
            results.bic = nan(obj.n_subj,1); % BIC
            f = waitbar(0,'Estimating RBM'); % Initialize waitbar

            % Cycle over subjects and fit RBM
            for i = 1:obj.n_subj

                % Update waitbar
                waitbar(i/obj.n_subj,f,'Estimating RBM');

                % Logical index for ID
                indices = (allSubBehavData.ID == i);

                % Extract subset for fields dynamically
                subBehavData = structfun(@(x) x(indices), allSubBehavData, 'UniformOutput', false);

                % Estimate model
                [bestParams, bestNegLogLike, bic] = obj.model_estimation(subBehavData, agent_vars);
                            
                % Store parameters
                results.parameters(i,:) = bestParams;
            
                % Store negative log-likelihood
                results.llh(i) = bestNegLogLike;

                % Store BIC
                results.bic(i) = bic;

            end

            % Close waitbar
            close(f)

            % Define column names
            columnNames = {obj.omikron_0, obj.omikron_1, obj.h, obj.s, obj.u, obj.sigma_H};
            columnNames = columnNames(which_vars_vec);

            % Convert matrix to table and assign column names
            results.parameters = array2table(results.parameters, 'VariableNames', columnNames);

        end
        
        function [bestParams, bestNegLogLike, bic] = model_estimation(obj, df_subj, agent_vars)
            % MODEL_ESTIMATION This function estimates the free
            % parameters of the model for a given subject
            % 
            %   Input
            %       obj: Estimation-object instance
            %       df_subj: Structure containing subject data
            %       agent_vars: Agent variables object instance
            %   
            %   Output
            %       bestParams: Best parameter estimates
            %       bestNegLogLike: Minimum negative log-likelihood
            %       bic: Bayesian information criterion
            
             % Free parameters
            which_vars_vec = struct2array(obj.which_vars);

            % Lower and upper bounds
            lb = struct2array(obj.lowerBoundaries);
            ub = struct2array(obj.upperBoundaries);
            
            % Optimization settings
            options = optimset('Algorithm', 'interior-point', 'Display', 'final'); 

            % Initialize log-likelihood variable
            bestNegLogLike = inf;
            bestParams = NaN;
            
            % Cycle over starting points
            for r = 1:obj.n_sp

                % Get start points
                if obj.rand_sp
                    startPoint = obj.get_starting_point();
                else
                    startPoint = [obj.omikron_0_x0, obj.omikron_1_x0,...
                        obj.h_x0, obj.s_x0, obj.u_x0, obj.sigma_H_x0];   
                end
                
                % Define objective function
                objFun = @(x) obj.llh(x, df_subj, agent_vars);

                % Optimize parameters
                [x_opt, fval] = fmincon(objFun, startPoint(which_vars_vec), [], [], [], [], lb(which_vars_vec), ub(which_vars_vec), [], options);

                % Store best results
                if fval < bestNegLogLike
                    bestNegLogLike = fval;
                    bestParams = x_opt;
                end
            end
                      
            % Compute BIC
            bic = obj.compute_bic(bestNegLogLike, sum(which_vars_vec), length(df_subj.ID));
        
        end
        
        function llh_sum = llh(obj, coeffs, df, agent_vars)
            %LLH_SUM This function computes the cumulative negative
            % log-likelihood of the data under the model
            % 
            %   Input
            %       obj: Estimation-object instance
            %       coeffs: Free parameters
            %       df: Data frame of the current subject
            %       agent_vars: Agent variables object instance
            % 
            %   Output
            %       llh_sum: Cumulative negative log-likelihood
            
            % Get fixed parameters
            fixed_coeffs = obj.fixed_mod_coeffs;
            
            % Initialize parameter list and counters
            sel_coeffs = [];
            idx = 1;
            
            % Assign selected coefficients for model estimation
            keys = fieldnames(obj.which_vars);
            for i = 1:numel(keys)
                key = keys{i};
                if obj.which_vars.(key)
                    sel_coeffs(end+1) = coeffs(idx);
                    idx = idx + 1;
                else
                    sel_coeffs(end+1) = fixed_coeffs.(key);
                end
            end
            
            % Set agent variables for the model
            agent_vars.h = sel_coeffs(3);
            agent_vars.s = sel_coeffs(4);
            agent_vars.u = exp(sel_coeffs(5));
            agent_vars.sigma_H = sel_coeffs(6);

            % Call the AlAgent object
            agent = AlAgentRBM(agent_vars);
            sim = false;

            % Compute the likelihood using task-agent interaction
            llh_mix = for_task_agent_int(df, agent, agent_vars, sel_coeffs, sim);
            
            % Consider prior over uncertainty-underestimation coefficient
            if obj.use_prior
                u_prob = log(normpdf(sel_coeffs(5), 0, 5));
            else
                u_prob = 0;
            end
            
            % Compute cumulative negative log-likelihood
            llh_sum = -1 * (sum(llh_mix) + u_prob);
        end

        function startPoint = get_starting_point(obj)
            %GET_STARTING_POINT This function samples random start points
            % for regression model estimation
            %
            %   Input
            %       obj: Estimation-object instance
            %
            %   Output
            %       startPoint: Sampled start points 
        
            % Sample from uniform distribution within specified range
            startPoint = [unifrnd(obj.omikron_0_x0_range(1), obj.omikron_0_x0_range(2)),...
                unifrnd(obj.omikron_1_x0_range(1), obj.omikron_1_x0_range(2)),...
                unifrnd(obj.h_x0_range(1), obj.h_x0_range(2)),...
                unifrnd(obj.s_x0_range(1), obj.s_x0_range(2)),...
                unifrnd(obj.u_x0_range(1), obj.u_x0_range(2)),...
                unifrnd(obj.sigma_H_x0_range(1), obj.sigma_H_x0_range(2))];
        end

        function bic = compute_bic(~, llh, n_params, n_trials)
            % BIC This function computes the Bayesian Information Criterion (BIC)
            % 
            % See Stephan et al. (2009). Bayesian model selection for
            % group studies. NeuroImage
            %
            %   Input
            %       llh: Negative log-likelihood
            %       n_params: Number of free parameters
            %       n_trials: Number of trials
            % 
            %   Output
            %       bic: Bayesian Information Criterion value

            bic = (-1 * llh) - (n_params / 2) * log(n_trials);
        end
    end
end
