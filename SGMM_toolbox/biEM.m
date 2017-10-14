% function to implement bivariate EM algorithm
% create the function obs2ITM (similar)
% input: obs, B_mat, par_ini, maxIter, TolFun

function [par, weight, labels] = biEM(y, B_mat, par_ini, Display)
TolFun = 1e-5;
maxIter = 1000;
if nargin < 4
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
L = zeros([n, K, maxIter]);
log_L = zeros([1, maxIter]);

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
        % differs from d = 1 to d = 2
        mu(:, j, i+1) = y' * weight(:,j) / sum_weight;
        res = y - repmat(mu(:, j, i+1)', [n, 1]);
        sigma(:, :, j, i+1) = res' * (res .* repmat(weight(:,j), [1, d])) / sum_weight;
        sigma(:, :, j, i+1) = (sigma(:, :, j, i+1) + sigma(:, :, j, i+1)')./2;
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


