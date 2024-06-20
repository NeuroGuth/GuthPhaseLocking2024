function out = TG_RandomPointsInCircle(myR)
%
% Creates n random points drawn from a circular area.
%
% Input
%   maxR    --> maximum radius
%   minR    --> minimum radius
%   N       --> number of points to create
%   centerX --> center of circle
%   centerY --> center of circle

% get limits in range [0, 1]
limits  = ([myR.minR, myR.maxR] ./ myR.maxR) .^ 2; % cave: square!

% calculate random points within circle in polar coordinates
r       = sqrt(rand(myR.N, 1) .* range(limits) + min(limits));
theta   = rand(myR.N, 1) .* 2 .* pi;

% transform into desired range
rTransf = r .* myR.maxR;

% transform into Cartesian coordinates
x       = myR.centerX + rTransf .* cos(theta);
y       = myR.centerY + rTransf .* sin(theta);

% output
out     = [x, y];