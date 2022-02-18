clc; clear; close all;
set(0,'DefaultLineMarkerSize',12);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesFontSize',14);
% rng(1)

%% Define Scenario
N = 256; % Length of Sinusoid
% number of sinusoids in the mixture of sinusoids
K = 5;    

% SNR OF EACH MEASUREMENT
% SNR for each sinusoid in dB scale 
% (all N measurements put together)
SNR = [44; 42; 40; 38; 36];

sigma = 1;              % noise variance sigma^2, instead of sigma
T = 60;
% sample gains and frequencies
omega_true = 2*pi*rand(K,1);
gain_true = bsxfun(@times,sqrt(sigma) * (10.^(SNR/20)), exp(1j*2*pi*rand(K,T)));  % K*T

% original signal
y_full = exp(1j* (0:(N-1)).' * omega_true.')/sqrt(N) * gain_true;

%% Choose measurement type
% normal measurements
M = N; % number of measurements = N
S = eye(N);

%% Measurements
% add measurement noise
noise = sigma*(randn(N,T) + 1j*randn(N,T))/sqrt(2);
y_noisy = y_full + noise;

% take measurements
y = S * y_noisy;
    
%% Algorithm 
% stopping criterion parameter tau
% when we know the noise level sigma^2,
% we choose tau to ensure that the false alarm rate is fixed:
% ie., the probability that we over estimate the size of the
% support set is fixed to some pre-determined number p_fa
p_fa = 1e-2;
% tau = sigma^2 * ( log(N) - log( log(1/(1-p_fa)) ) );
% T = 1;
tau= sigma*chi2inv((1-p_fa)^(1/N), 2*T)/2;

[omega_est, gain_est, residue] = MNOMP(y, S, tau);

% [omega_est, gain_est, residue] = extractSpectrum(y, S, tau);

%% Plot results
figure; 
subplot(1,2,1); polar(omega_true, sum(abs(gain_true).^2,2), 'bo');
hold on; polar(omega_est, sum(abs(gain_est).^2,2),'rx');
title(sprintf('Magnitude and frequency\nof each sinusoid'));
legend({'Truth', 'Estimate'}, 'Location', 'SouthOutside');

subplot(1,2,2); polar(angle(gain_true(:,1)), abs(gain_true(:,1)), 'bo');
hold on; polar(angle(gain_est(:,1)), abs(gain_est(:,1)), 'rx');
title(sprintf('Magnitude and phase\nof each sinusoid'));
legend({'Truth', 'Estimate'}, 'Location', 'SouthOutside');

figure; stem(0:(length(residue)-1), 10*log10(residue/sigma^2/M/T));
xlabel('Number of detected sinusoids');
ylabel('Residual power (relative to \sigma^2) in dB');
hold on; plot(K, 0, 'rx');
legend({'Estimate', 'Truth'});
