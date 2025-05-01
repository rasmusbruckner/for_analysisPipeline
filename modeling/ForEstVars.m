classdef ForEstVars
    % This class defines the instance variables for estimation

    properties

        % Parameter names for data frame
        omikron_0 = 'omikron_0'; % Motor noise
        omikron_1 = 'omikron_1'; % Learning-rate noise
        h = 'h'; % Hazard rate
        s = 's'; % Surprise sensitivity
        u = 'u'; % Uncertainty underestimation
        sigma_H = 'sigma_H'; % Catch trial

        % Fixed starting points (used if random starting points are disabled)
        omikron_0_x0 = 10;
        omikron_1_x0 = 0.01;
        h_x0 = 0.1;
        s_x0 = 0.999;
        u_x0 = 0.0;
        sigma_H_x0 = 1;

        % Range for random starting-point values (used if random starting points are enabled)
        omikron_0_x0_range = [3, 20];
        omikron_1_x0_range = [0.001, 1];
        h_x0_range = [0.001, 0.99];
        s_x0_range = [0.001, 0.99];
        u_x0_range = [1, 10];
        sigma_H_x0_range = [0, 0.5];

        % Free parameter indexes
        which_vars = struct();

        % Fixed values for fixed parameters
        fixed_mod_coeffs;

        % Estimation lower boundaries
        lowerBoundaries = struct();

        % Estimation lower boundaries
        upperBoundaries = struct();

        % Other attributes
        n_subj = NaN; % Number of participants
        rand_sp = true; % Use of random starting points during estimation
        n_sp = 10; % Number of random starting points
        use_prior = true; % Prior of uncertainty-underestimation parameter
    end

    methods
        function obj = ForEstVars()
            % Constructor: Initializes default estimation variables

            % Determine which parameters should be estimated
            obj.which_vars.omikron_0 = true;
            obj.which_vars.omikron_1 = true;
            obj.which_vars.h = true;
            obj.which_vars.s = true;
            obj.which_vars.u = true;
            obj.which_vars.sigma_H = true;

            % Estimation lower boundaries
            obj.lowerBoundaries.omikron_0_lb = 3;
            obj.lowerBoundaries.omikron_1_lb = 0.001;
            obj.lowerBoundaries.h_lb = 0.001;
            obj.lowerBoundaries.s_lb = 0.001;
            obj.lowerBoundaries.u_lb = -2;
            obj.lowerBoundaries.sigma_H_lb = 0;

            % Estimation upper boundaries
            obj.upperBoundaries.omikron_0_ub = 20; %100;
            obj.upperBoundaries.omikron_1_ub = 1;
            obj.upperBoundaries.h_ub = 0.999;
            obj.upperBoundaries.s_ub = 0.999;
            obj.upperBoundaries.u_ub = 15;
            obj.upperBoundaries.sigma_H_ub = 0.5;

            % Fixed and normative model coefficients
            obj.fixed_mod_coeffs = struct(...
                obj.omikron_0, 5, ...
                obj.omikron_1, 0.0, ...
                obj.h, 0.1, ...
                obj.s, 1.0, ...
                obj.u, 0.0, ...
                obj.sigma_H, 0.0);
        end
    end
end
