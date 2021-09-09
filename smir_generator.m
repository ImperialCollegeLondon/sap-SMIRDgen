function [h, H, beta_hat] = smir_generator(c, fs, r, s, orient, L, b, sphType, sphRadius, mic, N_harm, n_sample, oversampling_factor, order, type_or_table, interp_method)
%This function generate the room impulse response for a spherical microhpne
%array and a directional source in a rectangular room.
% Implemented by Sina Hafezi based on a Pseudocode by Jarret, Habets, and
% Naylor ("Rigid sphere room impulse response simulation: Algorithm and
% applications" 2012 paper)

%       Input(s)
%-------------------------------------
%     c:           speed of sound in m/s
%     fs:          processing sampling frequency in Hz
%     r:           1 x 3 vector specifying the (x,y,z) coordinates of the 
%                 centre of the array in m
%     s:           1 x 3 vector specifying the (x,y,z) coordinates of the 
%                 source in m

%     orient:      sound source orientation [yaw pitch roll] (in radian)
%     L:           1 x 3 vector specifying the room dimensions (x,y,z) in m
%     b:           2 x 3 walls reflection coefficient matrix [betha_x1 betha_y1 betha_z1;
%   betha_x2 betha_y2 betha_z2] or Reverberation Time T_60 in s
%     sphType:     type of spherical microphone array ['open', 'rigid']
%     sphRadius:   radius of the spherical microphone array in m
%     mic:         M x 2 matrix specifying the angles of the microphones
%                 (azimuth,elevation) in radians
%     N_harm:      maximum spherical harmonic order to use in spherical
%                 harmonic decomposition
%     oversampling_factor:           oversampling factor
%     n_sample:     number of samples of the RIR to calculate 
%                 
%     order:       maximum reflection order (default=-1, all possible reflections)
%     type_or_table:     either type of the source directivity pattern or the table containing
%                        the customised source directivity pattern (more info see
%                       gain_function.m) default='omni'
%     interp_method:    method for interpolation ('linear', 'nearest', 'spline',
%                       'cubic') default= 'linear'
%
%       Output(s):
%------------------
%     h           M x nsample matrix containing the calculated RIR(s)
%     H           M x K*nsample/2+1 matrix containing the calculated RTF(s)
%     beta_hat    If beta is the reverberation time, the calculated
%                 reflection coefficient is returned.

% Compute wall reflection coefficients from reverberation time
if (length(b) == 1)
    V = L(1)*L(2)*L(3);
    S = 2*(L(1)*L(3)+L(2)*L(3)+L(1)*L(2));
    TR = b;
    % Sabin-Franklin's formula
    alfa = 24*V*log(10)/(c*S*TR);
    if (alfa > 1)
        error('Error: The reflection coefficients cannot be calculated using the room parameters (room dimensions and reverberation time) supplied. Please supply the reflection coefficients or change the reverberation time/room dimensions.');
    end
    beta_hat = sqrt(1-alfa);
    b = repmat(beta_hat, 2, 3);
else
    beta_hat = b;
end


% rotating the source-directivity cartesina coordiante
ax_c0=eye(3); % default xyz axes (per row)
q_quart=quaternion(orient,'euler','zyx','frame');
ax_cr=rotatepoint(q_quart,ax_c0); % rotated coordiante system for directional source


% maximum distance of image to receiver based on the number of RIR samples
% specified by the user
max_dist=n_sample*c/fs;

Nm=[ceil(max_dist/(2*L(1))) ceil(max_dist/(2*L(2))) ceil(max_dist/(2*L(3)))];

total_mics = size(mic,1);            % Number of microphones
N_FFT = oversampling_factor * n_sample;                 % Oversampling

% spherical to cartesian for mic positions
mic_cart=zeros(total_mics,3);
[mic_cart(:,1) mic_cart(:,2) mic_cart(:,3)]=mysph2cart(sphRadius,mic(:,1),mic(:,2));

% Cell storing all images info ( order, p, m, Rpm, reflected_orientation, overall reflection co )
Images = get_images_info( s, ax_cr, L, b, Nm, order );
total_images=size(Images,1);

%SHD (Spherical Harmonic Decomposition)
max_k=N_FFT/2+1;   % maximum number of wavenumbers
wave_n=(([0:(max_k-1)]/N_FFT)*fs)*2*pi/c; %wavenumbers
%Delta=zeros(max_k,N_harm+1);
%Gamma=zeros(max_k,N_harm+1);

no_overflow = 0;
overflow_warning = 0;
while (no_overflow == 0)
    l = 0 : N_harm;  
        if isequal(sphType,'rigid')
            %Delta(k,l+1)= i / ( ( 0.5 * (besselh(l-1,wave_n(k)*sphRadius) - besselh(l+1,wave_n(k)*sphRadius)) )*(wave_n(k)*sphRadius)^2);
            Hankel_derivative=repmat(sqrt(pi./(2*wave_n*sphRadius)),N_harm+1,1) .* ( (bsxfun(@besselh,l+0.5-1,wave_n'*sphRadius).' - bsxfun(@besselh,l+0.5+1,wave_n'*sphRadius).')/2 - bsxfun(@besselh,l+0.5,wave_n'*sphRadius).' ./ repmat((2*wave_n*sphRadius),N_harm+1,1) );
            mode_strength= i./(Hankel_derivative .* repmat((wave_n*sphRadius).^2,N_harm+1,1)); % Wronskian relation
            if (any(any(isinf(bsxfun(@besselh,l+0.5-1,wave_n'*sphRadius)))))
                N_harm = N_harm - 1;
                if (overflow_warning == 0)
                    warning('The chosen value of N_harm is too high; trying a lower number...');
                    overflow_warning = 1;
                end
            else
                no_overflow = 1;
            end
        else
            %Delta(k,l+1)= besselj(l,wave_n(k)*sphRadius);
            mode_strength=repmat(sqrt(pi./(2*wave_n*sphRadius)),N_harm+1,1) .* bsxfun(@besselj,l+0.5,wave_n'*sphRadius).';
            no_overflow = 1;
        end
        %Gamma(k,l+1)= i*wave_n(k)*Delta(k,l+1)/(4*pi);
end

shd_k_l_dependent_all_sources = 1i * mode_strength.' .* repmat(wave_n.',1,N_harm+1);
shd_angle_l_dependent_all_sources = (2*l+1);    


H=zeros(total_mics,max_k);
h=zeros(total_mics,n_sample);

for image_id=1:total_images
    image_order=Images{image_id,1};
	p=Images{image_id,2};
	m=Images{image_id,3};
	R_pm=Images{image_id,4}; %origin-image vector
    image_orient=Images{image_id,5}; % [3x3] reflected_oriented coordinate system (each row is an axis), xyz
	overal_betha=Images{image_id,6};
    
    % arrayCentre-image vector
    R_pm=R_pm-r;
    
    
    if (norm(R_pm) + sphRadius <max_dist)
        Upsilon=zeros(total_mics,(N_harm+1));
        for n_mic=1:total_mics
            Cos_Theta=dot(R_pm, mic_cart(n_mic,:)) / (norm(R_pm) * norm(mic_cart(n_mic,:))); % cosine of angle between mic and image w/ resp. to array centre
            %Sometimes the dot product is slightly smaller/larger than one, which causes problems for legendre()
			if (Cos_Theta < -1)
                Cos_Theta = -1;
            elseif (Cos_Theta > 1)
				Cos_Theta = 1;
            end
            Psi=sphLegendre(N_harm,Cos_Theta);
            Upsilon(n_mic,:)=Psi.*shd_angle_l_dependent_all_sources ;
        end
        
        shd_k_l_dependent=zeros(size(shd_k_l_dependent_all_sources));
        for k=1:max_k
            sphbesselh_out=sphBesselh(N_harm,wave_n(k)*norm(R_pm));
            shd_k_l_dependent(k,:)=shd_k_l_dependent_all_sources(k,:) .* sphbesselh_out;
        end
        
        for n_mic=1:total_mics
            
            % emission angle calculation
            image_receiver=mic_cart(n_mic,:)-R_pm;
            
            % getting the cartesian coordiante of image-receiver vector in
            % the imgage oreinted coordiante system
            new_x=dot(image_receiver,image_orient(1,:));
            new_y=dot(image_receiver,image_orient(2,:));
            new_z=dot(image_receiver,image_orient(3,:));
            [rrad,az elev]=mycart2sph(new_x,new_y,new_z); % spherical coordinate of emission angle
            
            % emission angle w/ resp. to orientation
            emission_angle=[az elev];
            directivity_gain=gain_function(emission_angle,wave_n*c/(2*pi), type_or_table, interp_method );
            
            for k=1:max_k
                for ll=1:(N_harm+1)
                    H(n_mic,k)=H(n_mic,k)+overal_betha*directivity_gain(k)*Upsilon(n_mic,ll)*shd_k_l_dependent(k,ll);
                end
            end
            
        end
        
    end
        

end


% NOTE: The room transfer function is undefined for k = 0.
% This results in a DC offset when taking the inverse
% transform to obtain the RIR.
H(find(isnan(H))) = 0;

% To get the RIR the right way around (with MATLAB's Fourier  transform
% definition), we need an exponential in the form exp(-ikR)/kR not 
% exp(ikR)/kR, so we take the conjugate here. Alternatively Williams 
% defines the Fourier transform differently in Fourier Acoustics and 
% keeps exp(ikR).
H = conj(H);

h = ifft([H conj(H(:,N_FFT/2:-1:2))], N_FFT, 2, 'symmetric');
% Truncate oversampled RIRs
h = h(:,1:n_sample);


% High pass filtering (as proposed by Allen and Berkley) at 0.01 of
% sampling freq
%{
W = 2*pi*100/fs;
R1 = exp(-W);
B1 = 2*R1*cos(W);
B2 = -R1 * R1;
A1 = -(1+R1);

Y=zeros(1,3);

for p=1:n_sample
    X0=h(p);
    Y(2:end)=Y(1:end-1);
    Y(1)=B1*Y(2)+B2*Y(3)+X0;
    h(p)=Y(1)+A1*Y(2)+R1*Y(3);
end

%}
end

