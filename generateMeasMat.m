function S = generateMeasMat(N, M, proj_style, window_style)
% INPUTS
% N is the length of sinusoid
% M is the number of measurements (optional) (default is N)
% proj_style (optional) type of measurement matrix
%       SUPPORTED TYPES: subsample randomly : subsample (M < N)
%                        complex gaussian   : cmplx_gauss 
%                        real gaussian      : gauss
%                        bernoulli {pm 1}   : bernoulli
%                        complex bernoulli {pm 1,pm j} : cmplx_bernoulli
%                        regular            : full (default)
% window_style (optional) : all_ones (default), hamming, hann
% 
% https://bitbucket.org/wcslspectralestimation/continuous-frequency-estimation

if ~exist('M','var'), M = N;
elseif isempty(M), M = N; end

if ~exist('proj_style','var'), proj_style = 'full';
elseif isempty(proj_style), proj_style = 'full'; end

if ~exist('window_style','var'), window_style = 'all_ones';
elseif isempty(window_style), window_style = 'all_ones'; end



switch proj_style
    case 'subsample'
        samp = randperm(N);
        samp = samp(1:M);
        S = sparse(1:M, samp, ones(1, M), M, N);
    case 'cmplx_gauss'
        S = (randn(M,N) + 1j*randn(M,N))/sqrt(2*N);
    case 'gauss'
        S = randn(M,N)/sqrt(N);
    case 'bernoulli'
        S = (2 * (rand(M,N) > 0.5) - 1)/sqrt(N);
    case 'cmplx_bernoulli'
        S = randn(M,N) + 1j*randn(M,N);
        S = exp(1j*floor(angle(S)/(pi/2))*pi/2)/sqrt(N);
    case 'full'
        if M ~= N
            error('Number of measurements M must be equal to the length of the sinusoid N');
        end
        S = eye(N);
    otherwise
        error('Measurement matrix not supported');
end


switch window_style
    case 'all_ones'
        window = ones(N,1);
    case 'hamming'
        window = hamming(N);
    case 'hann'
        window = hann(N + 2);
        window = window(2:end-1);
    otherwise
        error('Windowing function not supported');
end
S = S * diag(window);
