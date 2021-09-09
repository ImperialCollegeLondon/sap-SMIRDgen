function [ r, theta, phi ] = mycart2sph( x , y, z )
%This function convert the cartesian coordinate system into spherical
%coordinate
% Implemented by Sina Hafezi. Nov 2014. EEE Dept. Imperial College London.

%   INPUTS
%------------
%   x:  x (either vector or value)
%   y:  y (either vector or value)
%   z:  z (either vector or value)

%   OUTPUTS
%------------
%   r:  radius (range)
%   theta:  azimuth angle (radian) [0 2pi)
%   phi:  elevation angle (radian) [-pi/2 pi/2]


hypotxy = hypot(x,y);
r = hypot(hypotxy,z);
phi = asin(z./r);
theta = wrapTo2Pi(atan2(y,x));
end

