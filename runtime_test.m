% Test Run script for SMIRDgen 
%   Sina Hafezi, Imperial College London 2021.
close all;
%% Setup
procFs = 8000;                      % Sampling frequency (Hz)
c = 343;                            % Sound velocity (m/s)
nsample = 300;                      % Length of desired RIR
N_harm = 20;                        % Maximum order of harmonics to use in SHD
K = 2;                              % Oversampling factor

L = [8 4 8];                        % Room dimensions (x,y,z) in m
sphLocation = [5+0.042 2 4];            % Receiver location (x,y,z) in m
%s = [2.37 4.05 4.4];                % Source location(s) (x,y,z) in m
s = [2 2 4];
%beta = repmat(0.5,2,3);       % Room reflection coefficients [\beta_x_1 \beta_y_1 \beta_z_1 ; \beta_x_2 \beta_y_2  \beta_z_2] 
beta=repmat([.8 .4 .2],2,1);
order = 1;                         % Reflection order (-1 is maximum reflection order)

sphRadius = 0.042;                  % Radius of the sphere (m)
sphType = 'rigid';                  % Type of sphere (open/rigid)

mic = [pi 0];			% Microphone positions (azimuth, elevation) in radian

orient=deg2rad([0 0 0]);           % (radian) Source Orientation [yaw pitch roll] yaw(azimuth) pitch(-elevatoin) roll

% Directivity pattern

% Standard pattern
source_directivity_pattern='omni';

% Customized pattern
%{
%Note for Customized pattern: The sample azimuth, elevation and frequency must fully cover their range due to interpolation)
azimuth_samples=deg2rad(0:5:(360-5));
elevation_samples=deg2rad(-90:5:90);
frequency_samples=0:100:(procFs*K/2);
gain=ones(length(elevation_samples),length(azimuth_samples),length(frequency_samples)); % Sampled Omnidirectinal
source_directivity_pattern={azimuth_samples,elevation_samples,frequency_samples,gain}; % Constructing a cell for pattern
%}

interp_method='linear';             % interpolation method used in directivity pattern in case of customised pattern


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

%% Run spherical array simulation
tic;
[h, H, beta_hat] = smir_generator(c, procFs, sphLocation, s,orient, L, beta, sphType, sphRadius, mic, N_harm, nsample, K, order,source_directivity_pattern,interp_method);
time=toc
% Plotting
%%{
mic_to_plot = 1;

figure;
subplot(211);
    plot([0:nsample-1]/procFs,h(mic_to_plot,1:nsample), 'r')
    xlim([0 (nsample-1)/procFs]);
    title(['Room impulse response at microphone ', num2str(mic_to_plot)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
subplot(212);
    plot((0:1/(K*nsample):1/2)*procFs,mag2db(abs(H(mic_to_plot,1:K*nsample/2+1))), 'r');
    title(['Room transfer function magnitude at microphone ', num2str(mic_to_plot)]);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude (dB)');
%}


