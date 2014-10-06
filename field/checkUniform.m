function [isUniform] = checkUniform(measurementPoints)
% function [isUniform] = checkUniform(measurementPoints)
%
% If mesh is non-uniform, then need to use calculated element volumes to
% scale point loads (otherwise simple element volume calculation can be done
% for all uniform elements).
%

x = unique(measurementPoints(:, 1));
y = unique(measurementPoints(:, 2));
z = unique(measurementPoints(:, 3));

isUniform = all(abs(diff(x)/(x(2)-x(1))-1 < 10^-9)) &&...
            all(abs(diff(y)/(y(2)-y(1))-1 < 10^-9)) &&...
            all(abs(diff(z)/(z(2)-z(1))-1 < 10^-9));
