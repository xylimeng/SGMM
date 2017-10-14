function [obs, CompactLabel] = pars2obs(template, mu, sigma, labels)
% PARS2OBS generates observations from parameters including the template maps,
%     mu's and sigma's. When the latent labels are provided, the templates
%     will not be used.
%
%     'template' - matrix $n$ by $K$, where $n$ is the sample size (one
%         dimensional), and $K$ is the number of clusters. Each row of
%         $template$ must sum to one;
%     'mu' - the mean vector of size 1 by $K$;
%     'sigma' - the standard deviation vector of size 1 by $K$;
%     'labels' - optional; the underlying labels. When 'labels' is provided,
%         'template' will not used.
%
% % Examples for pars2obs.m
%     n = 100; % sample size 
%     K = 2; % number of components
%     mu = [0.1, 0.2]'; 
%     sigma = [0.1, 0.1]'; 
% 
%     location = (1:n)/n; % scale to be in [0,1]
%     b1 = cdf('norm', location, 0.4, 0.1); % sigmoid function 
%     b2 = 1 - b1; 
%     template = [b1', b2']; 
% 
%     [obs, labels] = pars2obs(template, mu, sigma); 
% 
%     h = figure; 
%     subplot(1,2,1)
%     plot(location,b1, 'red')
%     hold on
%     plot(location, b2, 'blue')
%     title('b function'); 
%     hleg1 = legend('b1', 'b2');
%     set(hleg1, 'Location', 'NorthEast')
%     hold off
% 
%     subplot(1,2,2) 
%     scatter(location(labels == 1), obs(labels == 1), 'red')
%     title(sprintf('mu = [%.2f, %.2f]\nsigma = [%.2f, %.2f]', mu, sigma))
%     hold on 
%     scatter(location(labels == 2), obs(labels == 2), 'blue')
%     hold off
%
%   Version 0.2; 12-Mar-2014 12:34pm
%       add reshape for (mu, sigma) in case the input is row vector; 
%   Version 0.1; 23-Nov-2013 14:05:49
%       fix the bug regarding the order of CompactLabel (sort by rowNumber)
%   Version 0.0; Date: 22-Nov-2013 13:30:30

% To do: 1. delete par2obs.m, P2z_obs.m
%        2. incorporate convolution transformation in the furture but not
%        here; the convolution will be applied when generationg template or
%        after obs, which means neither has nothing to do with this
%        function PAR2OBS
%        3. do simulation to check the histgram of obs 

% reshape first 
mu = reshape(mu, [length(mu), 1]); 
sigma = reshape(sigma, [length(sigma), 1]); 

if nargin < 4
    n = size(template, 1);
    SparseLabel = mnrnd(1, template); % size same as template : sparse
    vec_sigma = SparseLabel * sigma;
    vec_mean = SparseLabel * mu;
    obs = randn([n, 1]) .* vec_sigma + vec_mean;
    [rowNumber, CompactLabel] = find(SparseLabel);
    tmp = sortrows([rowNumber, CompactLabel]); % sort by rowNumbers
    CompactLabel = tmp(:,2); 
end

if nargin == 4
    vec_sigma = sigma(labels);
    vec_mean =  mu(labels);
    obs = randn([n, 1]) .* vec_sigma + vec_mean;
    CompactLabel = labels;
end

end


