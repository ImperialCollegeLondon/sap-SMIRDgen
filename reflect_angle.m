function reflected_angle = reflect_angle( angle, p )
%This function finds the emission angle source-receiver based on the
%position of image and the emission angle of image-reciever
%   Implemented by Sina Hafezi. Nov 2014. EEE Dept. Imperial College London.
%   Input(s)
%--------------
%   angle:  spherical angle [azimuth elevation] (in radian) 
%   p:      image position vector [px py pz] (ps are binary)
%   NOTE THE RANGE OF ANGLES AS FOLLOW:
%   azimuth: [0 2pi)
%   elevation: [-pi/2 +pi/2]


%   Outpus(s)
%--------------
%   reflected_angle:    spherical angle after all reflections [azimuth
%   elevation] (in radian)

theta=angle(1);
phi=angle(2);

% check if input angles are within the range
[theta phi]=angle_sanity(theta, phi);

if (p(1)==0)
    if (p(2)==1)
        theta=2*pi-theta;
    end
else
    if (p(2)==0)
        theta=wrapTo2Pi(pi-theta);
    else
        theta=wrapTo2Pi(pi+theta);
    end
end

if (p(3)==1)
    phi=-phi;
end
reflected_angle=[theta phi];

end

