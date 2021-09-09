function [ x y z ] = mysph2cart( r, theta, phi )
%This function convert the spherical coordinate system into cartesian
%coordinate
% Implemented by Sina Hafezi. Nov 2014. EEE Dept. Imperial College London.

%   INPUTS
%------------
%   r:  radius or range (either vector or value)
%   theta:  azimuth angle (radian) (either vector or value)
%   phi:  elevation angle (radian) (either vector or value)


%   OUTPUTS
%------------
%   x:  x 
%   y:  y 
%   z:  z

% In matlab Pi is not preciesly Pi. E.g. sin(pi) is not absolutely zero
if (theta==pi)
    theta=sym(pi);
end
if (phi==pi)
    phi=sym(phi);
end

x=r.*cos(phi).*cos(theta);
y=r.*cos(phi).*sin(theta);
z=r.*sin(phi);


end

