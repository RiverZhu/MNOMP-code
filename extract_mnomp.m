function [omegaList, xList, res_inf_normSq_rot] = extract_mnomp(Y, R, A, M, R_s, R_c, tau, opt)
%%%%% Multiple snapshot newtonied orthogonal matching pursuit algorithm
% Code is written by Lin Han and Jiang Zhu. If you have any
% problems, please contact jiangzhu16@zju.edu.cn
% Date: July 04 , 2019
% Input:
% Y: the uncompressed measurements
% R: the number of grid frequencies = gamma*N_sensor
% A: the compression matrix
% M: the number of compressed measurements
% R_s: the number of single refinements
% R_c: the number of cyclic refinements
% tau: the thresholding value

% Output:
% omegaList: the estimated frequencies of sinusoids
% xList: the estimated x vectors
% res_inf_normSq_rot: the residue energy

[N_sensors, N_snapshots] = size(Y);
omegaList = [];
xList  = [];
Y_r = Y;

coarseOmega = 2*pi*(0:(R-1))/R;       % the grid frequency
a_ori = zeros(M,R);
ant_idx = 0:N_sensors-1;
ant_idx = ant_idx - (N_sensors-1)/2;
ant_idx = ant_idx(:);
for i = 1:R
    a_ori(:,i) = A*exp(1j* ant_idx * coarseOmega(i))/sqrt(N_sensors);
end
res_inf_normSq_rot = inf;
while (res_inf_normSq_rot >= tau)
    
    [omega, x, Y_r] = detectnew(Y_r, R, a_ori, coarseOmega, ant_idx, A);            % detecting a new frequency
    for i = 1:R_s
        [omega, x, Y_r] = refineone(ant_idx, omega, N_sensors, Y_r, x, A, opt);           % refine (omega, x) using single frequency Newton update algorithm
    end
    omegaList = [omegaList;omega];
    xList  = [xList; x];
    [omegaList, xList, Y_r] = refineAll(Y_r, omegaList, xList, N_sensors, ant_idx, R_s, R_c, A,opt);           % refine parameters of already detected (omega, x) one at a time
    [omegaList, xList, Y_r]= solveleastsquares(A, ant_idx, omegaList, Y);            % update all x vector estimate by least squares
    res_inf_normSq_rot = residue(Y_r);
end
omegaList = wrap_2pi(omegaList);         % restricting frequencies to [0, 2*pi)
for i=1:N_snapshots
    xList(:,i) = xList(:,i).*exp(1j*ant_idx(1)*omegaList);           % return to the state before symmetry
end
xList = xList.';
end

function [omega, x, Y_r] = detectnew(Y, R, a_ori, coarseOmega, ant_idx, A)
%%%%% detect the frequency and the corresponding x vector estimate
% Input:
% Y: the uncompressed measurements
% R: the number of grid frequencies = gamma*N_sensor
% a_ori
% coarseOmega: the grid frequency
% ant_idx: 0:N_sensors-1
% A: the compression matrix

% Output:
% omega: the detected frequency
% x: the corresponding x vector estimate
% Y_r: the residual energy
[N_sensors, ~] = size(Y);
S = zeros(1,R);
for k = 1:R
    S(:,k) = trace(real(Y'*a_ori(:,k)*a_ori(:,k)'*Y/(a_ori(:,k)'*a_ori(:,k))));
end
[~,IDX] = max(S);
omega = coarseOmega(IDX);                               % find the frequency
a_omega = A*exp(1j* ant_idx * omega)/sqrt(N_sensors);
x = a_omega'*Y*(a_omega'*a_omega)^(-1);          % update x vector estimate
Y_r = Y - a_omega*x;                                            % generate the residual measurements
end

function [omega, x, Y_r] = refineone(ant_idx, omega, N_sensors, Y_r, x, A,opt)
%%%%%% single refinement
% Iuput:
% ant_idx: 0:N_sensors-1
% omega: the detected frequency in the detection step
% N_sensors: the number of sensors
% Y_r: the residual measurements after the detection step
% x: the x vector estimate in the detection step
% A: the compression matrix

% Output:
% omega: the refined frequency
% x: the refined x vector estimate
% Y_r: the refined residual measurements
[~, T] = size(Y_r);
a_omega = A*exp(1j*ant_idx*omega)/sqrt(N_sensors);
da_omega = A*1j*ant_idx.*a_omega;
d2a_omega = -A*ant_idx.^2.*a_omega;
Y = Y_r + a_omega*x;

if opt ==1
    der1 = trace(real(Y'*(da_omega*a_omega' + a_omega*da_omega')*Y/(a_omega'*a_omega) - ...
        Y'*(a_omega*a_omega')*Y*(da_omega'*a_omega + a_omega'*da_omega)/(a_omega'*a_omega)^2));
    d1 = Y'*(d2a_omega*a_omega' + 2*(da_omega*da_omega') + a_omega*d2a_omega')*Y/(a_omega'*a_omega);
    d2 = -2*Y'*(da_omega*a_omega' + a_omega*da_omega')*Y*(da_omega'*a_omega + a_omega'*da_omega)/(a_omega'*a_omega)^2;
    d3 = Y'*(a_omega*a_omega')*Y*(d2a_omega'*a_omega + 2*(da_omega'*da_omega) + a_omega'*d2a_omega)/(a_omega'*a_omega)^2;
    d4 = -2*Y'*(a_omega*a_omega')*Y*(da_omega'*a_omega + a_omega'*da_omega)^2/(a_omega'*a_omega)^3;
    der2 = trace(real((d1+d2-d3-d4)));
    
    if der2< 0
        omega1_next = omega - der1/der2;             % the estimated frequency by using the Newton method
    else
        omega1_next = omega;
    end
else
    der1 = zeros(1, T);
    der2 = zeros(1, T);
    for ii = 1:T
        der1(ii) = -2*real(x(ii) * Y_r(:, ii)'*da_omega);
        der2(ii) = -2*real(x(ii) * Y_r(:, ii)'*d2a_omega) +...
            2*abs(x(ii))^2*(da_omega'*da_omega);
    end
    
    % UPDATE OMEGA
    der1_sum = sum(der1);
    der2_sum = sum(der2);
    if der2_sum > 0
        omega1_next = omega - der1_sum/der2_sum;
    else
        omega1_next = omega - sign(der1)*(1/4)*(2*pi/N)*rand(1);
    end
end



a_theta_next = A*exp(1j*ant_idx*omega1_next)/sqrt(N_sensors);
x_next = (a_theta_next'*a_theta_next)^(-1)*a_theta_next'*Y;         % update the corresponding x vector estimate
Y_r_next = Y - a_theta_next*x_next;
if trace(Y_r_next'*Y_r_next)<=trace(Y_r'*Y_r)         % update the residual measurements
    omega = omega1_next;
    x = x_next;
    Y_r = Y_r_next;
else
    x = a_omega'*Y*(a_omega'*a_omega)^(-1);       % update the corresponding x vector estimate
    Y_r = Y - a_omega*x;        % the residual measurements
end
end

function [omegaList, xList, Y_r] = refineAll(Y_r, omegaList, xList, N_sensors, ant_idx, R_s, R_c, A, opt)
%%%%% cyclic refinements
% Input:
% Y_r: the residual measurements after the single refinements
% omegaList: all the detected frequencies
% xList: all the detected x vector estimate
% N_sensors: the number of sensors
% ant_idx: 0:N_sensors-1
% R_s: the number of single refinements
% R_c: the number of cyclic refinements
% A: the compression matrix

% Output:
% omegaList: the refined frequencies after cyclic refinements
% xList: the refined x vector estimate after cyclic refinements
% Y_r: the refined residual measurements

K = length(omegaList);          % the number of detected frequencies
order = 1:K;
for i = 1:R_c
    for j = 1:K
        l = order(j);
        omega = omegaList(l);
        x = xList(l,:);
        for kk = 1:R_s
            [omega, x, Y_r] = refineone(ant_idx, omega, N_sensors, Y_r, x, A, opt);
        end
        omegaList(l) = omega;
        xList(l,:) = x;
    end
end
end

function omega_prime= wrap_2pi(omega)
%%%%% restricting frequencies to [0, 2*pi)
% Input:
% omega: to be restricted frequency
% Output:
% omega_prime: restricted frequency
omega_prime = angle(exp(1j*omega));
omega_prime(omega_prime < 0) = omega_prime(omega_prime < 0) + 2*pi;
end

function [omegaList, xList, Y_r] = solveleastsquares(A, ant_idx, omegaList, Y)
%%%%% update all x vector estimate by least squares
% Input:
% A: the compression matrix
% ant_idx: 0:N_sensors-1
% omegaList: the refined frequencies
% Y: the uncompressed measurements

% Output:
% omegaList: the refined frequencies
% xList: the refined x vector estimate
% Y_r: the residual measurements
N_sensor = length(ant_idx);
B = A*exp(1j*ant_idx*omegaList.')/sqrt(N_sensor);
xList = (B'*B)\(B'*Y);
Y_r = Y - B*xList;
end

function [res_inf_normSq_rot] = residue(Y_r)
%%%%% the computation of residual energy
[N_sensor, N_snapshots] = size(Y_r);
A = eye(N_sensor);
a_ori = A*dftmtx(N_sensor)/sqrt(N_sensor);
g = zeros(N_snapshots, 1);
g_all = zeros(N_sensor, 1);
for k = 1:N_sensor
    for i = 1:N_snapshots
        g(i) = (abs(a_ori(:,k)'*Y_r(:,i))).^2;
    end
    g_all(k) = sum(g);
end
res_inf_normSq_rot = max(g_all);
end