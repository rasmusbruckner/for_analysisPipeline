function concentration = residual_fun(abs_dist, motor_noise, lr_noise)
%RESIDUAL_FUN This function computes updating noise (residuals)
% as a combination of two noise components
%
% The noise variable is returned in terms of van-Mises concentration,
% which is some kind of precision, where var = 1/concentration.
%
%   Input
%       abs_dist: Absolute distance -- predicted update or prediction error
%       motor_noise: Motor-noise parameter (imprecise motor control)
%       lr_noise: Learning-rate-noise parameter (more noise for larger update)
%
%   Output
%       concentration: Updating noise expressed as von-Mises concentration

% Compute updating noise expressed as variance, where
% (1) 1/motor noise is updating variance due to imprecise motor control and
% (2) learning-rate noise models more noise for larger updates
var = (1./motor_noise) + abs_dist * lr_noise;

% Translate back to concentration
concentration = 1./var;

end