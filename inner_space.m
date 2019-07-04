function [omega_true] = inner_space(delta_dft, N_signal, ratio)
% compute the interval between any two frequency of omega_true
% Code is written by Lin Han and Jiang Zhu. If you have any
% problems, please contact jiangzhu16@zju.edu.cn
% Date: July 04 , 2019
    while true
        omega_true = 2*pi*rand(N_signal,1);
        omega_true = sort(omega_true);
        interval = ones(N_signal, N_signal);
        for j =1:N_signal -1
             for k = j+1:N_signal
                 a = abs(omega_true(k) - omega_true(j));
                 b = abs(2*pi - a);
                 interval(j,k) = min(a,b);
             end
        end     
        if min(min(interval)) > ratio*delta_dft
            break;
        end
        
    end
end

