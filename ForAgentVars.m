classdef ForAgentVars
    %FORAGENTVARS This class specifies the instance variables
    % of the reduced Bayesian model

    properties

        % Default agent variables
        s = 1; % surprise sensitivity
        h = 0.1; % hazard rate
        u = 0.0; % uncertainty underestimation
        q = 0; % reward bias
        sigma = 10; % noise in the environment (standard deviation)
        sigma_0 = 6.1875; % initial variance of predictive distribution
        sigma_H = 1; % catch-trial standard deviation (uncertainty) of helicopter cue
        tau_0 = 0.5; % initial relative uncertainty
        omega_0 = 1; % initial change-point probability
        mu_0 = 3.1416; % initial belief about mean
        max_x = 2*pi; % maximum outcome

    end

    methods(Static)
        function obj = AgentVars()
            % Constructor to initialize default agent variables
        end
    end
end
