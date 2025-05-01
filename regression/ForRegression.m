classdef ForRegression
    %FORREGRESSION This class specifies the instance variables and
    % methods of the common methods of the circular regression analyses

    properties

        % Coefficient names
        beta_0 % intercept
        beta_1 % PE
        beta_2 % RU*PE
        beta_3 % CPP*PE
        beta_4 % hit*PE
        beta_5 % kappa*PE
        beta_6 % visible*PE
        beta_7 % visible*EE
        omikron_0 % motor noise
        omikron_1 % learning-rate noise
        uniform % uniform component

        % Start point of the parameters
        beta_0_x0
        beta_1_x0
        beta_2_x0
        beta_3_x0
        beta_4_x0
        beta_5_x0
        beta_6_x0
        beta_7_x0
        omikron_0_x0
        omikron_1_x0
        uniform_x0

        % Range of random starting points
        beta_0_x0_range
        beta_1_x0_range
        beta_2_x0_range
        beta_3_x0_range
        beta_4_x0_range
        beta_5_x0_range
        beta_6_x0_range
        beta_7_x0_range
        omikron_0_x0_range
        omikron_1_x0_range
        uniform_x0_range

        % When prior is used: pior mean
        beta_0_prior_mean
        beta_1_prior_mean
        beta_2_prior_mean
        beta_3_prior_mean
        beta_4_prior_mean
        beta_5_prior_mean
        beta_6_prior_mean
        beta_7_prior_mean
        omikron_0_prior_mean
        omikron_1_prior_mean
        uniform_prior_mean
        prior_mean

        % When prior is used: pior width
        beta_0_prior_width
        beta_1_prior_width
        beta_2_prior_width
        beta_3_prior_width
        beta_4_prior_width
        beta_5_prior_width
        beta_6_prior_width
        beta_7_prior_width
        omikron_0_prior_width
        omikron_1_prior_width
        uniform_prior_width
        prior_width

        % Parameters to be estimated
        which_vars

        % Estimation lower boundaries
        lowerBoundaries

        % Estimation upper boundaries
        upperBoundaries

        % Other properties
        n_subj % number of subjects
        rand_sp % random starting points
        n_sp % number of starting points
        regressionComponents % free parameters
        usePrior % use prior for regularization

    end

    methods

        function obj = ForRegression(reg_vars)
            % Constructor to initialize default regression variables

            % Coefficient names
            obj.beta_0 = reg_vars.beta_0;
            obj.beta_1 = reg_vars.beta_1;
            obj.beta_2 = reg_vars.beta_2;
            obj.beta_3 = reg_vars.beta_3;
            obj.beta_4 = reg_vars.beta_4;
            obj.beta_5 = reg_vars.beta_5;
            obj.beta_6 = reg_vars.beta_6;
            obj.beta_7 = reg_vars.beta_7;
            obj.omikron_0 = reg_vars.omikron_0;
            obj.omikron_1 = reg_vars.omikron_1;
            obj.uniform = reg_vars.uniform;

            % Start point of the parameters
            obj.beta_0_x0 = reg_vars.beta_0_x0;
            obj.beta_1_x0 = reg_vars.beta_1_x0;
            obj.beta_2_x0 = reg_vars.beta_2_x0;
            obj.beta_3_x0 = reg_vars.beta_3_x0;
            obj.beta_4_x0 = reg_vars.beta_4_x0;
            obj.beta_5_x0 = reg_vars.beta_5_x0;
            obj.beta_6_x0 = reg_vars.beta_6_x0;
            obj.beta_7_x0 = reg_vars.beta_7_x0;
            obj.omikron_0_x0 = reg_vars.omikron_0_x0;
            obj.omikron_1_x0 = reg_vars.omikron_1_x0;
            obj.uniform_x0 = reg_vars.uniform_x0;

            % Range of random starting points
            obj.beta_0_x0_range = reg_vars.beta_0_x0_range;
            obj.beta_1_x0_range = reg_vars.beta_1_x0_range;
            obj.beta_2_x0_range = reg_vars.beta_2_x0_range;
            obj.beta_3_x0_range = reg_vars.beta_3_x0_range;
            obj.beta_4_x0_range = reg_vars.beta_4_x0_range;
            obj.beta_5_x0_range = reg_vars.beta_5_x0_range;
            obj.beta_6_x0_range = reg_vars.beta_6_x0_range;
            obj.beta_7_x0_range = reg_vars.beta_7_x0_range;
            obj.omikron_0_x0_range = reg_vars.omikron_0_x0_range;
            obj.omikron_1_x0_range = reg_vars.omikron_1_x0_range;
            obj.uniform_x0_range = reg_vars.uniform_x0_range;

            % When prior is used: pior mean
            obj.beta_0_prior_mean = reg_vars.beta_0_prior_mean;
            obj.beta_1_prior_mean = reg_vars.beta_1_prior_mean;
            obj.beta_2_prior_mean = reg_vars.beta_2_prior_mean;
            obj.beta_3_prior_mean = reg_vars.beta_3_prior_mean;
            obj.beta_4_prior_mean = reg_vars.beta_4_prior_mean;
            obj.beta_5_prior_mean = reg_vars.beta_5_prior_mean;
            obj.beta_6_prior_mean = reg_vars.beta_6_prior_mean;
            obj.beta_7_prior_mean = reg_vars.beta_7_prior_mean;
            obj.omikron_0_prior_mean = reg_vars.omikron_0_prior_mean;
            obj.omikron_1_prior_mean = reg_vars.omikron_1_prior_mean;
            obj.uniform_prior_mean = reg_vars.uniform_prior_mean;
            obj.prior_mean = reg_vars.prior_mean;

            % When prior is used: pior width
            obj.beta_0_prior_width = reg_vars.beta_0_prior_width;
            obj.beta_1_prior_width = reg_vars.beta_1_prior_width;
            obj.beta_2_prior_width = reg_vars.beta_2_prior_width;
            obj.beta_3_prior_width = reg_vars.beta_3_prior_width;
            obj.beta_4_prior_width = reg_vars.beta_4_prior_width;
            obj.beta_5_prior_width = reg_vars.beta_5_prior_width;
            obj.beta_6_prior_width = reg_vars.beta_6_prior_width;
            obj.beta_7_prior_width = reg_vars.beta_7_prior_width;
            obj.omikron_0_prior_width = reg_vars.omikron_0_prior_width;
            obj.omikron_1_prior_width = reg_vars.omikron_1_prior_width;
            obj.uniform_prior_width = reg_vars.uniform_prior_width;
            obj.prior_width = reg_vars.prior_width;

            % Parameters to be estimated
            obj.which_vars = reg_vars.which_vars;

            % Estimation lower boundaries
            obj.lowerBoundaries = reg_vars.lowerBoundaries;

            % Estimation upper boundaries
            obj.upperBoundaries = reg_vars.upperBoundaries;

            % Extract other properties
            obj.n_subj = reg_vars.n_subj;
            obj.rand_sp = reg_vars.rand_sp;
            obj.n_sp = reg_vars.n_sp;
            obj.regressionComponents = reg_vars.regressionComponents;
            obj.usePrior = reg_vars.usePrior;

        end

        function results = run_estimation(obj, allSubBehavData)
            %RUN_ESTIMATION This function manages the estimation of the
            % circular regression model
            %
            %   Input
            %       obj: Regression-object instance
            %       allSubBehavData: Behavioral data
            %
            %   Ouptut
            %       results: Structure with coefficients, error terms, and
            %       negative log-likelihood
            
            % Throw error if omikron_0 not free
            if obj.which_vars.omikron_0 == false
                error("Regression model needs omikron_0 parameter. Please add it back in")
            end

            % Free parameters
            which_vars_vec = struct2array(obj.which_vars);

            % Initialize coefficient structure
            results = struct();

            % Regression coefficients
            results.parameters = nan(obj.n_subj, sum(which_vars_vec));

            % Negative log-likelihood
            results.llh = nan(obj.n_subj,1);

            % Initialize waitbar
            f = waitbar(0,'Estimating regression model');

            % Cycle over subjects and fit circular regression model
            for i = 1:obj.n_subj

                % Update waitbar
                waitbar(i/obj.n_subj,f,'Estimating regression model');

                % Logical index for ID
                indices = (allSubBehavData.subj_num == i);

                % Extract subset for fields dynamically
                subBehavData = structfun(@(x) x(indices), allSubBehavData, 'UniformOutput', false);

                % Estimate model
                [bestParams, bestNegLogLike] = obj.estimation(subBehavData);

                % Store regression coefficients
                results.parameters(i,:) = bestParams;

                % Store negative log-likelihood
                results.llh(i) = bestNegLogLike;

            end

            % Close waitbar
            close(f)

            % Define column names
            columnNames = {obj.beta_0, obj.beta_1, obj.beta_2,...
                obj.beta_3, obj.beta_4, obj.beta_5, obj.beta_6,...
                obj.beta_7, obj.omikron_0, obj.omikron_1, obj.uniform};
            columnNames = columnNames(which_vars_vec);

            % Convert matrix to table and assign column names
            results.parameters = array2table(results.parameters, 'VariableNames', columnNames);

        end

        function [bestParams, bestNegLogLike] = estimation(obj, subBehavData)
            %ESTIMATION This function performs the estimation of the
            % circular regression model for a given subject
            %
            %   Input
            %       obj: Regression-object instance
            %       subBehavData: Subject's behavioral data
            %
            %   Output
            %       bestParams: Best parameter estiamtes
            %       bestNegLogLike: Minimum negative log-likelihood

            % Free parameters
            which_vars_vec = struct2array(obj.which_vars);

            % Lower and upper bounds
            lb = struct2array(obj.lowerBoundaries);
            ub = struct2array(obj.upperBoundaries);

            % Data matrix
            datamat = obj.get_datamat(subBehavData);

            % Create structure containing all relevant data for regression
            data = struct();

            % Select trials for regression
            a_t = subBehavData.a_t; % all updates

            % Select finite and updating trials
            sel = isfinite(a_t) & all(isfinite(datamat), 2) & a_t~=0;
            data.Y = a_t(sel);
            data.X = datamat(sel,:);

            % Optimization setting (different options possible):

            % options = optimset('Algorithm','interior-point',...
            % 'MaxPCGIter', 5000, 'MaxProjCGIter', 5000, 'MaxSQPIter',...
            % 5000, 'MaxRLPIter', 5000);
            % options = optimset('Algorithm', 'sqp', 'Display','final');

            % This is the one we're currently using
            options = optimset('Algorithm', 'interior-point', 'Display', 'final');

            % Initialize log-likelihood variable
            bestNegLogLike = inf;
            bestParams = NaN;

            % Cycle over starting points
            for i = 1:obj.n_sp

                % Get start points
                if obj.rand_sp
                    startPoint = obj.get_starting_point();
                else
                    startPoint = [obj.beta_0_x0, obj.beta_1_x0, obj.beta_2_x0,...
                        obj.beta_3_x0, obj.beta_4_x0, obj.beta_5_x0, obj.beta_6_x0,...
                        obj.beta_7_x0, obj.omikron_0_x0, obj.omikron_1_x0,...
                        obj.uniform_x0];
                end

                % Define objective function with data as an extra parameter
                objFun = @(x) obj.llh(x, data);

                % Estimate parameters
                [params, negLogLike] = fmincon(objFun, startPoint(which_vars_vec), [], [], [], [], lb(which_vars_vec), ub(which_vars_vec), [], options);

                % Store best results
                if negLogLike < bestNegLogLike
                    bestParams = params;
                    bestNegLogLike = negLogLike;
                end
            end
        end

        function negLogLike = llh(obj, params, data)
            % NEGLOGLIKE This function compute the negative log-likelihood
            % conditional on the regression parameters
            %
            %   Input
            %       obj: Regression-object instance
            %       params: Regression parameters
            %       data: Subject data
            %
            %   Output
            %       negLogLike: Negative log-likelihood

            % Extract free regression parameters
            coeffs = params(1:sum(obj.regressionComponents));
            coeffs_prior_mean = obj.prior_mean(obj.regressionComponents);
            coeffs_prior_width = obj.prior_width(obj.regressionComponents);

            % Predicted updates
            yHat = data.X*coeffs';

            % Residuals
            if obj.which_vars.omikron_1

                % Combining motor and learning-rate noise
                abs_dist = abs(yHat);
                motor_noise = params(sum(obj.regressionComponents)+1);
                lr_noise = params(sum(obj.regressionComponents)+2);
                concentration = residual_fun(abs_dist, motor_noise, lr_noise);

            else
                % Motor noise only
                concentration = repmat(params(sum(obj.regressionComponents)+1), length(data.X),1);
            end

            % Initialize likelihood vector
            allVM_likelihoods = nan(length(data.Y),1);

            % Compute probability of update conditional on regression model
            for i = 1:length(data.Y)
                allVM_likelihoods(i) = circ_vmpdf(data.Y(i), yHat(i), concentration(i));
            end

            % Add uniform mixture component:
            % Accounts for trials where subjects predict randomly so
            % points that violate the gaussian don't break the model
            if obj.which_vars.uniform == 1
                % If there is a uniform mixture, adjust the likelihoods
                % accordingly
                allVM_likelihoods = allVM_likelihoods.*(1-params(end))+ (1./(2.*pi)).*(params(end));
            end

            % Negative log-likelihood
            negLogLike = -1.*sum(log(allVM_likelihoods));

            % Favour estimates closer to prior mean
            % Currently implemented for key parameters; consider adding
            % priors for error terms and uniform component
            if obj.usePrior == 1
                priorProb = sum(log(normpdf(coeffs, coeffs_prior_mean, coeffs_prior_width)));
                if priorProb < 1e-300
                    priorProb = 1e-300; % set some minimum prior probability (flat outside some range)
                end
                negLogLike = negLogLike-priorProb;
            end

            if ~isfinite(negLogLike) || isnan(negLogLike)
                error('log-likelihood incorrect');
            end
        end

        function startPoint = get_starting_point(obj)
            %GET_STARTING_POINT This function samples random start points
            % for regression model estimation
            %
            %   Input
            %       obj: Regression-object instance
            %
            %   Output
            %       startPoint: Sampled start points

            % Sample from uniform distribution within specified range
            startPoint = [unifrnd(obj.beta_0_x0_range(1), obj.beta_0_x0_range(2)),...
                unifrnd(obj.beta_1_x0_range(1), obj.beta_1_x0_range(2)),...
                unifrnd(obj.beta_2_x0_range(1), obj.beta_2_x0_range(2)),...
                unifrnd(obj.beta_3_x0_range(1), obj.beta_3_x0_range(2)),...
                unifrnd(obj.beta_4_x0_range(1), obj.beta_4_x0_range(2)),...
                unifrnd(obj.beta_5_x0_range(1), obj.beta_5_x0_range(2)),...
                unifrnd(obj.beta_6_x0_range(1), obj.beta_6_x0_range(2)),...
                unifrnd(obj.beta_7_x0_range(1), obj.beta_7_x0_range(2)),...
                unifrnd(obj.omikron_0_x0_range(1), obj.omikron_0_x0_range(2)),...
                unifrnd(obj.omikron_1_x0_range(1), obj.omikron_1_x0_range(2))...
                unifrnd(obj.uniform_x0_range(1), obj.uniform_x0_range(2))];
        end

        function datamat = get_datamat(obj, subBehavData)
            %GET_DATAMAT This function creates the design matrix for the
            % circular regression model
            %
            %   Input
            %       obj: Regression-object instance
            %       subBehavData: Subject's behavioral data
            %
            %   Output
            %       datamat: Design matrix

            % Initialize matrix
            datamat = [];

            % Fill depending on specificed free parameters
            if obj.which_vars.beta_0 == true
                datamat = [datamat, ones(size(subBehavData.delta_t))];
            end

            if obj.which_vars.beta_1 == true
                datamat = [datamat, subBehavData.delta_t];
            end

            % Note: This is just one option and we could implement other
            % combinations. This one ensures that the LR is not > 1
            if obj.which_vars.beta_2 == true
                datamat = [datamat, (subBehavData.tau_t .* (1-subBehavData.omega_t)) .* subBehavData.delta_t];
            end

            if obj.which_vars.beta_3 == true
                datamat = [datamat, subBehavData.omega_t .* subBehavData.delta_t];
            end

            if obj.which_vars.beta_4 == true
                datamat = [datamat, subBehavData.delta_t .* subBehavData.hit_dummy];
            end

            if obj.which_vars.beta_5 == true
                datamat = [datamat, subBehavData.delta_t .* subBehavData.sigma_dummy];
            end

            if obj.which_vars.beta_6 == true
                datamat = [datamat, subBehavData.delta_t .* subBehavData.visible_dummy];
            end

            if obj.which_vars.beta_7 == true
                datamat = [datamat, subBehavData.e_t .* subBehavData.visible_dummy];
            end
        end

        function df_sim = sample_data(obj, df_params, n_trials, allSubBehavData)
            %SAMPLE_DATA This function randomly samples regression data for
            % model validation purposes
            %
            %   Input
            %       obj: Regression-object instance
            %       df_params: Regression paramters for simulation
            %       n_trials: Number of trials
            %       allSubBehavData: Optional subject behavioral data
            %
            %   Output
            %       df_sim: Sampled regression updates

            % Check if file name suffix is provided
            if ~exist('allSubBehavData', 'var') || isempty(allSubBehavData)
                sim_data = true;
            else
                sim_data = false;
            end

            % Number of simulations
            n_sim = length(df_params.omikron_0);

            % Initialize
            df_data = table(); % regression predictors
            df_sim = table(); % simulated data

            % Cycle over simulations
            for i = 1:n_sim

                % Extract regression coefficients
                coeffs = table2array(df_params(i,:));

                % Regression variables
                if sim_data

                    % Randomly generate data
                    df_data.delta_t = unifrnd(-pi, pi, n_trials, 1);
                    df_data.e_t = unifrnd(-pi, pi, n_trials, 1);
                    df_data.visible_dummy = binornd(1, 0.1, n_trials, 1);
                    df_data.hit_dummy = binornd(1, 0.5, n_trials, 1);
                    df_data.sigma_dummy = [zeros(n_trials/2, 1); ones(n_trials/2, 1)];
                    df_data.tau_t = rand(n_trials, 1);
                    df_data.omega_t = rand(n_trials, 1);
                else

                    % Optionally based on subject data:

                    % Logical index for ID
                    indices = (allSubBehavData.subj_num == i);

                    % Extract subset for fields dynamically
                    subBehavData = structfun(@(x) x(indices), allSubBehavData, 'UniformOutput', false);

                    % Extract regression data
                    df_data.delta_t = subBehavData.delta_t;
                    df_data.e_t = subBehavData.e_t;
                    df_data.visible_dummy = subBehavData.visible_dummy;
                    df_data.hit_dummy = subBehavData.hit_dummy;
                    df_data.sigma_dummy = subBehavData.sigma_dummy;
                    df_data.tau_t = subBehavData.tau_t;
                    df_data.omega_t = subBehavData.omega_t;
                end

                % Create design matrix
                datamat = obj.get_datamat(df_data);

                % Initialize parameter list and counters
                sel_coeffs = [];
                idx = 1;

                % Assign selected coefficients for model estimation
                keys = fieldnames(obj.which_vars);
                for c = 1:numel(keys)
                    key = keys{c};
                    if obj.which_vars.(key)
                        sel_coeffs(end+1) = coeffs(idx);
                        idx = idx + 1;
                    end
                end

                % Extract free regression parameters
                coeffs = sel_coeffs(1:sum(obj.regressionComponents));

                % Predicted updates
                yHat = datamat * coeffs';

                % Residuals
                if obj.which_vars.omikron_1

                    % Compute updating noise based on common function
                    abs_dist = abs(yHat);
                    motor_noise = sel_coeffs(sum(obj.regressionComponents)+1);
                    lr_noise = sel_coeffs(sum(obj.regressionComponents)+2);
                    concentration = residual_fun(abs_dist, motor_noise, lr_noise);

                else
                    % Motor noise only
                    concentration = repmat(sel_coeffs(sum(obj.regressionComponents)+1), length(datamat),1);
                end

                % Sample updates from Gaussian using standard deviation
                % yHat = normrnd(yHat, sqrt(1./concentration));
                for j = 1:length(yHat)
                    yHat(j) = circ_vmrnd(yHat(j), concentration(j), 1);
                end
                % todo: van mises!

                % Store update and ID
                df_data.a_t = yHat;
                df_data.subj_num = repmat(i, height(df_data), 1);

                % Combine all data
                df_sim = [df_sim; df_data];

            end
        end
    end
end