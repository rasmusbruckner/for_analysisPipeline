function kappa_up = residual_fun(abs_dist, motor_noise, lr_noise)
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
%       kappa_up: Updating noise expressed as von-Mises concentration

% Compute updating noise expressed as variance
% (1) motor noise is updating variance due to imprecise motor control and
% (2) learning-rate noise models more noise for larger updates
up_noise = motor_noise + lr_noise * (rad2deg(abs_dist));

% Convert std of update distribution to radians and kappa
up_noise_radians = deg2rad(up_noise);
up_var_radians = up_noise_radians.^2;
kappa_up = 1 ./ up_var_radians;

end