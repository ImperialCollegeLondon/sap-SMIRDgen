function [azimuth elevation] = angle_sanity( azimuth, elevation )
%This function checks the range of the input spheircal angles and wraps to
%the correct range if needed.
% Implemented by Sina Hafezi. Nov 2014. EEE Dept. Imperial College London.

%%       Input(s)
%-------------------------------------
%   azimuth:    azimuth angle in radian
%   elevation:  elevation angle in radian

%%       Output(s)
%-------------------------------------
%   azimuth:    azimuth angle in radian ranged within [0 2Pi)
%   elevation:  elevation angle in radian ranged within [-Pi/2 Pi/2]

if ~( (elevation>=-pi/2) && (elevation<=pi/2) )
    elevation=wrapToPi(elevation);
    if (elevation>pi/2)
        elevation=pi-elevation;
        azimuth=azimuth+pi;
    elseif (elevation<-pi/2)
        elevation=-pi-elevation;
        azimuth=azimuth+pi;
    end
end
if ~( (azimuth>=0) && (azimuth<2*pi) )
    azimuth=wrapTo2Pi(azimuth);
    azimuth(azimuth==2*pi)=0;
end

end

