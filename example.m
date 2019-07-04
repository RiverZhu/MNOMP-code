% MNOMP for solving line spectrum estimation problem in the MMV setting
% Code is written by Lin Han and Jiang Zhu. If you have any
% problems, please contact jiangzhu16@zju.edu.cn
% Date: July 04 , 2019

clc; 
clear; 
close all;
rng(1)
%% Parameter initialization
N = 50;                                                        % Length of Sinusoid
K = 5;                                                        % number of sinusoids in the mixture of sinusoids 
T = 6;                                                          % number of snapshots  
gamma = 4;                                              % the oversampling ratio for MNOMP
SNR = 30;
sigma = 1;                                                  % noise variance
grid_interval = 2*pi/N;                               % interval of grids
ratio = 1;                                                 % frequency_interval_min = 2.5*grid_interval
opt = 0;                                                   % Refine step, opt=0, the method adopted by Madhow, opt = 1, another method which calculates the deriative directly.

%% Model generation
w_true = inner_space(grid_interval, K, ratio);      % the true value
w_true_sort = sort(w_true);
Y_noiseless = zeros(N,T);
A = zeros(N,K);
X = zeros(K,T);
for i = 1:K
     A(:,i) = exp(1j* (0:(N-1)).' * w_true(i))/sqrt(N);
     x_hat = randn(1, T) + 1j*randn(1, T);
     k = sqrt( 10.^(SNR/10)*T*sigma/(norm(A(:,i)*x_hat, 'fro'))^2 );   
     X(i, :) = k*x_hat;
     Y_noiseless = Y_noiseless + A(:,i) * X(i,:);    
end       
Noise = sqrt(sigma)*(randn(size(Y_noiseless)) + 1j*randn(size(Y_noiseless)))/sqrt(2);   % noise
Y = Y_noiseless + Noise;

%% MNOMP
M = N; 
R = gamma*N;
Phi = eye(N);
R_c = 1;
R_s = 3;
P_false_nominal = 0.01;
tau_mnomp= sigma*chi2inv((1-P_false_nominal)^(1/N), 2*T)/2;
[omegaList, xList, res_inf_normSq_rot]  = extract_mnomp(Y, R, Phi, M, R_s, R_c, tau_mnomp,opt);


%% plot
figure(1)
lw = 2;
msz = 8;
fsz = 14;
amp = zeros(K, 1);
amp_est = zeros(K, 1);

for k = 1:K
    amp(k) = norm(X(k, :));
    amp_est(k) = norm(xList(:, k));
end

polar(w_true, amp, 'bo');
hold on; 
polar(omegaList, amp_est,'rx');
title(sprintf('Magnitude and frequency\nof each sinusoid'));
legend({'Truth', 'Estimate'}, 'Location', 'SouthOutside');
set(gca, 'FontSize', fsz,'FontName','Times New Roman', 'LineWidth',lw);


