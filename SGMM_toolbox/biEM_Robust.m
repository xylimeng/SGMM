% function to implement robust bivariate EM algorithm
% create the function obs2ITM (similar)
% input: obs, B_mat, par_ini, maxIter, TolFun
function [par, weight, labels] = biEM_Robust(y, B_mat, par_ini, tunning, Display)
% Versions: 
%     v0.1 - 5/26/2015
%          - update the robust function using Section 2.8 in the book by Geoffrey Mclachlan and Kaye Basford
%          - disregard Bai's approach completely because we have multivariate obs.
%          - play with 'tunning' for your problem 
%
% Last modified: 5/26/2015
TolFun = 1e-5;
maxIter = 1000;

if nargin < 4
    tunning = sqrt(2 .* (1 - 1/9  + 1/3 * 1.645)^3);
    Display = 'off'; 
end 

if nargin < 5
    Display = 'off';
end

[n, K] = size(B_mat);
d = size(y, 2); % dimension of y's : d = 2 here

% pre-create matrics
mu = zeros([d, K, maxIter]);
sigma = zeros([d, d, K, maxIter]);
pi_k = zeros([K, maxIter]);
gamma_k = zeros([K, maxIter]);

weight = zeros([n, K]);
r1 = weight; % re-weight for robust EM - mean
r2 = weight; % re-weight for robust EM - covariance matrix
L = zeros([n, K, maxIter]);
log_L = zeros([1, maxIter]);
% %% weighting functions for robust EM used in Bai's paper - disregarded 
% if strcmp(robust_mean, 'huber') == 1
%     RobustMean = @(x) max(-1.345, min(1.345, x));
% end
% 
% if strcmp(robust_mean, 'tukey') == 1
%     RobustMean = @(x) x .* max((1 - (x./4.685).^2), 0).^2;
% end
% 
% RobustScale = @(x) min(1 - (1 - (x./1.56).^2).^3, 1) ./ (x.^2);

%% Weight functions for robust EM 
% tunning = sqrt(chi2inv(0.90, 2)); % Devlin, Gnanadesikan and Kettenring(1981)
% tunning = sqrt(2 .* (1 - 1/9  + 1/3 * 1.645)^3); % Campbell (1984, 1985)
% tunning = 3; % test
u1 = @(s) min(1, tunning./s); 
%% initialize pars using the input par_ini
i = 1;
mu(:,:,1) = par_ini.mu;
sigma(:,:,:, 1) = par_ini.sigma;
pi_k(:,1) = par_ini.p';
gamma_k(:,1) = par_ini.gamma';
convergence = false;

%% calculate summary for ith iter: PI_mat, weight, log_L (E-step)
PI_mat = P2Pt(B_mat, gamma_k(:,i));
for j = 1:K
    L(:,j,i) = mvnpdf(y, mu(:,j,i)', sigma(:, :, j, i));
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
        
        %%% calculate re-weights for robust EM - mean
        r1(:, j) = mahal_dist(y, mu(:, j, i), sigma(:, :, j, i));
        weight_mean = weight(:,j) .* u1(r1(:, j));
        sum_weight_mean = sum(weight_mean);
        % differs from d = 1 to d = 2
        mu(:, j, i+1) = y' * weight_mean / sum_weight_mean;
        
        res = y - repmat(mu(:, j, i+1)', [n, 1]);        
        %%% calculate re-weights for robust EM - covariance
        r2(:, j) = mahal_dist(y, mu(:, j, i + 1), sigma(:, :, j, i));
        weight_scale = weight(:, j) .* (u1(r2(:, j))).^2;
        sum_weight_scale = sum(weight_scale);
        sigma(:, :, j, i+1) = res' * (res .* repmat(weight_scale, [1, d])) / sum_weight_scale;
        sigma(:, :, j, i+1) = (sigma(:, :, j, i+1) + sigma(:, :, j, i+1)')./2;
        % two purposes: make it symmetric and multiplied by 2 - robust EM
        gamma_k(j, i+1) = sum(sum_weight) ./ sum(B_mat(:,j) ./ update_gamma_tmp);
    end
    
    % standize gamma's
    gamma_k(:,i+1) = gamma_k(:,i+1) ./ sum(gamma_k(:,i+1));
    
    %%% calculate summary for (i+1)th iter: weight, log_L (E-step)
    PI_mat = P2Pt(B_mat, gamma_k(:,i+1));
    for j = 1:K
        L(:,j,i+1) = mvnpdf(y, mu(:,j,i+1)', sigma(:, :, j, i+1));
        s = PI_mat(:,j); % eq.(19) in the United Seg paper
        weight(:,j) = s .* L(:,j,i+1); %NOTE: double-check use 'i + 1' or 'i'?
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

if strcmp(Display, 'on')
    if (iter == maxIter) && (~convergence)
        fprintf('maxIter reached but not converge')
    end
    if (iter < maxIter)
        fprintf(['Converged! Iter Num = ', num2str(iter), '\n'])
    end
end

i = iter;
par.mu = mu(:,:,i);
par.sigma = sigma(:,:,:,i);
par.p = pi_k(:,i)';
par.gamma = gamma_k(:,i)';
% posterior = weight;
[~, labels] = max(weight,[], 2);
end


