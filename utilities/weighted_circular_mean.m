function mean_angle = weighted_circular_mean(angles, weights)
%WEIGHTED_CIRCULAR_MEAN Computes the weighted mean of circular data
%
% Input:
%   angles - vector of angles in radians
%   weights - vector of weights (same size as angles)
%
% Output:
%   mean_angle - weighted mean angle in [0, 2*pi)

if nargin < 2
    error('Two input arguments required: angles and weights.');
end

angles = angles(:);  % Ensure column vector
weights = weights(:);

if length(angles) ~= length(weights)
    error('Angles and weights must be the same length.');
end

z = weights .* exp(1i * angles);
mean_angle = angle(sum(z) / sum(weights));

% Wrap to [0, 2*pi)
mean_angle = mod(mean_angle, 2*pi);
end
