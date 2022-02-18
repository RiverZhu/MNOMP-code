clc; clear; close all;
set(0,'DefaultLineMarkerSize',12);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesFontSize',14);
% rng(1)
%% Define Scenario
N = 256; % Length of Sinusoid
% number of sinusoids in the mixture of sinusoids
K = 5;    
T = 40;
% SNR OF EACH MEASUREMENT
% effective SNR for each sinusoid in dB scale 
% (if all N measurements put together)
SNR = [44; 42; 40; 38; 36]+6; 
% SNR = [44; 42; 40; 38; 36; 38*ones(3,1)]+6; 

% (add 6 dB to compensate power loss due to 
% compressive measurements M = N/4 - we take one-fourth of 
% measurements and we loose one fourth of the power) 

% noise level
sigma = 1;

% sample gains and frequencies
omega_true = 2*pi*rand(K,1);
gain_true = bsxfun(@times,sqrt(sigma) * (10.^(SNR/20)), exp(1j*2*pi*rand(K,T)));  % K*T

% original signal
y_full = exp(1j* (0:(N-1)).' * omega_true.')/sqrt(N) * gain_true;

%% Choose measurement type
% Compressive measurements
% Number of compressive measurements
M = round(N/4); 
% type of measurement matrix
measure_type = 'cmplx_bernoulli'; 
% windowing weights
% options 'all_ones' (default), 'hamming' and 'hann'
window_type = 'all_ones'; 
S = generateMeasMat(N, M, measure_type, window_type);


%% Measurements
% add measurement noise
noise = sigma * (randn(N,T) + 1j*randn(N,T))/sqrt(2);
y_noisy = y_full + noise;

% take **compressive** measurements
y = S * y_noisy;

%% NOISE STATS
% *** noise is not white ***
% noise covariance is sigma^2*S*S';
covar = sigma * (S*S');

%% USING **extractSpectrum** WHEN NOISE IS COLORED
% whiten measurements using noise covariance matrix
whitenMAT = (covar)^(-1/2);
y_white = whitenMAT * y; % noise is white in y_white
S_white = whitenMAT * S; % new measurement matrix
sigma_white = 1; % after whitening - noise power is 0 dB

%% Algorithm 
% stopping criterion parameter tau
% when we know the noise level sigma^2,
% we choose tau to ensure that the false alarm rate is fixed:
% ie., the probability that we over estimate the size of the
% support set is fixed to some pre-determined number p_fa
p_fa = 1e-2;
% tau1 = sigma_white^2 * ( log(N) - log( log(1/(1-p_fa)) ) );
sigma = sigma_white^2;
tau = sigma*chi2inv((1-p_fa)^(1/N), 2*T)/2;   % threshold is not accurate in this setting
% [omega_est, gain_est, residue] = MNOMP(y_white, S_white, tau, overSamplingRate, R_s , R_c);
[omega_est, gain_est, residue] = MNOMP(y_white, S_white, tau);

% abs(omega_est1-omega_est)
%% Plot results
figure; 
subplot(1,2,1); polar(omega_true, abs(gain_true(:,1)), 'bo');
hold on; polar(omega_est, abs(gain_est(:,1)),'rx');
title(sprintf('Magnitude and frequency\nof each sinusoid'));
legend({'Truth', 'Estimate'}, 'Location', 'SouthOutside');

subplot(1,2,2); polar(angle(gain_true(:,1)), abs(gain_true(:,1)), 'bo');
hold on; polar(angle(gain_est(:,1)), abs(gain_est(:,1)), 'rx');
title(sprintf('Magnitude and phase\nof each sinusoid'));
legend({'Truth', 'Estimate'}, 'Location', 'SouthOutside');

figure; stem(0:(length(residue)-1), 10*log10(residue/sigma_white^2/M/T));
xlabel('Number of detected sinusoids');
ylabel('Residual power (relative to \sigma^2) in dB');
hold on; plot(K, 0, 'rx');
legend({'Estimate', 'Truth'});
