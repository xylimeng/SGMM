% function to apply ITM to observation
% The input includes:
%     obs - observations (vectorized);
%     B_mat - prior probability maps (n by K) - b_{ik};
%     par_ini - initial estimates of (mu, sigma, gamma);
%     maxIter - max number of iterations;
%     TolFun - tolerance for convergence;
% Output:
%     par - estimated parameters values;
%     labels - estimated labels (hard segmentation);
%
% Versions:
%     0.3 - edited @ 12/3/14 - change 'maxIter' to 1000 (default). 
%     0.2 - edited @ 5/5/14 - standize gamma's to ensure their sum == 1
%     0.1 - edited at 12-Mar-2014
%     0.0 - created at 23-Nov-2013

% To-do:
%     1. set up initial values

% debug utility: 
%   obs = obs_slice;
%   B_mat = B_mat_slice;
%   par_ini = Par_true;
%   TolFun = 1e-5;
%   maxIter = 100;

function [par, labels, posterior] = obs2ITM(obs, B_mat, par_ini, maxIter, TolFun)

if nargin < 4
    TolFun = 1e-5;
    maxIter = 1000;
end

if nargin < 5
    TolFun = 1e-5;
end

% 'obs': n by 1
% 'B_mat' : n by K
n = numel(obs);

if (size(obs,2) > 1)
    error('input observation must have the dimension n by 1')
end


K = size(B_mat, 2);

mu = zeros([K, maxIter]);
sigma = zeros([K, maxIter]);
pi_k = zeros([K, maxIter]);
gamma_k = zeros([K, maxIter]);

weight = zeros([n, K]);
L = zeros([n, K, maxIter]);
log_L = zeros([1, maxIter]);

%% initialize pars using the input par_ini
i = 1;
mu(:,1) = par_ini.mu';
sigma(:,1) = par_ini.sigma';
pi_k(:,1) = par_ini.p';
gamma_k(:,1) = par_ini.gamma';
convergence = false;

%% calculate summary for ith iter: PI_mat, weight, log_L (E-step)
PI_mat = P2Pt(B_mat, gamma_k(:,i));
for j = 1:K
    L(:,j,i) = pdf('norm', obs, mu(j,i), sigma(j,i));
    s = PI_mat(:,j); % eq.(19) in the United Seg paper
    weight(:,j) = s .* L(:,j,i);
end
den = sum(weight,2);
weight = weight ./ repmat(den, [1, K]);  % q_{ik} in eq. (18) in the United Seg paper
log_L(i) = sum(log(den)); % surprising: it's just den!

%% update and judge convergence
while ((i <= maxIter) && (~convergence))
    %%% prepare to update gamma_k
    update_gamma_tmp = B_mat * gamma_k(:,i);
    %%% update pi_k mu and sigma for (i+1)th iter (M-step)
    for j = 1:K
        pi_k(j,i+1) = mean(weight(:,j));
        sum_weight = sum(weight(:,j));
        mu(j, i+1) = weight(:,j)' * obs / sum_weight;
        sigma(j, i+1) = sqrt(weight(:,j)' * ((obs(:) - mu(j,i+1)).^2) / sum_weight);
        gamma_k(j, i+1) = sum(sum_weight) ./ sum(B_mat(:,j) ./ update_gamma_tmp);
    end        
    
    % standize gamma's 
    gamma_k(:,i+1) = gamma_k(:,i+1) ./ sum(gamma_k(:,i+1)); 
    
    %%% calculate summary for (i+1)th iter: weight, log_L (E-step)    
    PI_mat = P2Pt(B_mat, gamma_k(:,i+1));
    for j = 1:K
        L(:,j,i+1) = pdf('norm', obs, mu(j,i+1), sigma(j,i+1));
        s = PI_mat(:,j); % eq.(19) in the United Seg paper
        weight(:,j) = s .* L(:,j,i);
    end
    den = sum(weight,2);
    weight = weight ./ repmat(den, [1, K]);  % q_{ik} in eq. (18) in the United Seg paper    
    log_L(i+1) = sum(log(den)); % surprising: it's just den!       
    
    %%% converge or not
    ratio = abs((log_L(i + 1) - log_L(i))/log_L(i)); %relative change
    convergence = (ratio < TolFun);
    i = i + 1;
end

%% save estimates and display messages
iter = i - 1;

if (iter == maxIter) && (~convergence)
    fprintf('maxIter reached but not converge')
end
if (iter < maxIter)
    fprintf(['Converged! Iter Num = ', num2str(iter), '\n'])
end

par.mu = mu(:,i)';
par.sigma = sigma(:,i)';
par.p = pi_k(:,i)';
par.gamma = gamma_k(:,i)';
posterior = weight; 
[~, labels] = max(weight,[], 2);

end
