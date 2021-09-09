function gain = gain_function( emission, frequency, type_or_table, interp_method )
%This function retuns the freq-dependant gain associated with a particular
%direction for a standard or a cusotmised directivity pattern
% Implemented by Sina Hafezi. Nov 2014. EEE Dept. Imperial College London.

%       Input(s)
%-------------------------------------
%   emission:      emission angles [azimuth elevation] w/ resp. to the source orientaton (radian)
%   frequency:      vector containing the frequency samples
%   type_or_table:  either type of the directivity or the table containing
%   the customised directivity pattern. deflaut='omni'
%   %interp_method: method for interpolation (linear, nearest, spline,
%   cubic) default='linear'

%         STANDARD DIRECTIVITY PATTERNS (Char)
%               Type    |   Symbol
% -----------------------------------
%    omnidirectional    |   omni
%      bidirectional    |   bi
% 	   hypercardioid    |   hyper
%           cardioid    |   card
%        subcardioid    |   sub
%      hemispherical    |   hemi
%              delta    |   delta
%       invers-delta    |   idelta

%        CUSTOMISED DIRECTIVITY PATTERN (1x4 Cell)
%  { [azi] (1xm) , [ele] (1xn) , [freq] (1xk) , [gain] (nxmxk) }


%       Output(s)
%-------------------------------------
%	gain:	freq-dependant directivity gain vector with the length same as
%	the input frequency vector


% Sanity checks
narginchk(1,4);
if (nargin<3)
    type_or_table='omni'; %default if not specified
end
if (nargin<4)
    interp_method='linear';
end
% check if input angles are within the range
theta=emission(1);
phi=emission(2);
[theta phi]=angle_sanity(theta, phi);

if ((iscell(type_or_table))||(isequal(type_or_table,'head')))
    	% Head or Arbitrary (Cell look up)
	% either head or customised (look up table)
                if(isequal(type_or_table,'head'))
                    % standrd head manikin
                    % table= initialise it w/ custom table [mxnxk]
                else
                    % arbitrary values
                    table=type_or_table;
                end
                
                % look up or interpolate
                % table cell: { theta/azi [1xm] , phi/ele [1xn] , freq [1xk] , gain [nxmxk] }
                if (~iscell(table))
                    error('Table must be a cell {1x4}');
                end
                if (~isequal(size(table),[1,4]))
                    error('Table must have a size of {1x4} { [azimuths 1xm] , [elevations 1xn] , [freqs 1xk] , [gains nxmxk] }');
                end
                if (min(size(table{4}))<2)
                    error('gain must have at least two samples in each dimension');
                end
                if ((min(size(table{1}))>1)||(min(size(table{2}))>1)||(min(size(table{3}))>1))
                    error('first, second and the third table matrix must be 1D vector');
                end
                if (~isequal(size(table{4}),[length(table{2}) length(table{1}) length(table{3})]))
                    error('size of gain table must be n x m x k where m, n and k are respectively the length of the first, second and third elements in table');
                end
                
                % angle range sanity check for table given by the user
                if ((max(table{1})>=2*pi) || (min(table{1})<0))
                    error('Azimuth angle in table must be within [0 2Pi) radian');
                end
                if ((max(table{2})>pi/2) || (min(table{1})<-pi/2))
                    error('Elevation angle in table must be within [-Pi/2 Pi/2] radian');
                end
                
                
                
                switch interp_method
                        case 'linear'
                        case 'nearest'
                        case 'spline'
                        case 'cubic'
                        otherwise
                            error('interpolation method not found');
                end
                gain=interp3(table{1},table{2},table{3},table{4},theta,phi,frequency,interp_method);
                gain=squeeze(gain);
                
else
alpha=-1; %assuming it's not a standard source directivity type
% standard source types

switch type_or_table
	case 'omni'
        %omnidirectional
        alpha=1;
	case 'bi'
        % bidirectional
        alpha=0;
	case 'hyper'
        % hypercardiod
        alpha=0.25;
	case 'card'
        % cardioid
        alpha=0.5;
	case 'sub'
        % subcardioid
        alpha=0.75;
end
if (alpha~=-1)
	% standard source directivity type
    cos_theta=cos(theta);
    cos_phi=cos(phi);
    % In matlab pi is not really pi. E.g. sin(pi) is not absolutely zero
    if ((theta==pi/2)||(theta==3*pi/2))
        cos_theta=0;
    end
    if ((phi==pi/2)||(phi==3*pi/2))
        cos_phi=0;
    end
    
	strength = cos_phi * cos_theta;
	gain = abs(alpha + (1-alpha) * strength);
    
else
	% non-standard source directivity type (other special types or customised)
        
	% other types
	switch type_or_table
        case 'hemi'
            % hemispherical
            if ((theta>pi/2)&&(theta<3*pi/2)) 
                gain=0;
            else
                gain=1; 
            end
        case 'delta'
            %delta
            if ((theta==0)&&(phi==0))
                gain=1;
            else
                gain=0; 
            end
        case 'idelta'
                %inverse-delta
                if ((theta==0)&&(phi==0))
                    gain=0;
                else
                    gain=1; 
                end
        otherwise
                error('No Standard Source Directivity Type found!');

    end 
    

end
gain=repmat(gain,size(frequency));
end

end

