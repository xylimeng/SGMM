% %% function to obtain scores based on the estimated parameters
% par - struct with fields mu, sigma;
% obs - noisy observations;
% Pi_mat - prior probability maps
% posterior - posterior probability maps
% 
% OUTPUT: Z (scores) with the following types
% type - 'QT' for Quantile Transformation
%        'soft1' for soft & moment demoninator
%        'soft2' for soft & intuitive demoninator
%        'hard' for hard segmentaion

% Version:
%     0.2 edited at 05/05/14
%         - previous 'QT' transformation is wrong - it used the posterior
%         probability;
%         - preivously, there is no 'soft2'.
%     0.1 documentated at 03/12/14
%     0.0 created at 11/19/13

function Z = pars2score(par, obs, Pi_mat, posterior)

% reshape first: work with column vector
par.mu = reshape(par.mu, [length(par.mu), 1]); 
par.sigma = reshape(par.sigma, [length(par.sigma), 1]); 
par.gamma = reshape(par.gamma, [length(par.gamma), 1]); 

% QT-pi: quantile transformation using Pi_mat. 
[n,K] = size(Pi_mat);
u = zeros([n, K]);
for k = 1:K
    u(:,k) = cdf('normal', obs, par.mu(k), par.sigma(k));
end

copula = sum(u .* Pi_mat, 2); 
Z.QTvPi = icdf('normal', copula, 0, 1);

% % QT2: quantile transformation using w_mat. 
% w_mat = posterior; 
% [n,K] = size(w_mat);
% u = zeros([n, K]);
% for k = 1:K
%     u(:,k) = cdf('normal', obs, par.mu(k), par.sigma(k));
% end
% 
% copula = sum(u .* w_mat, 2); 
% Z.QTvW = icdf('normal', copula, 0, 1);

% % soft1: soft segmentation with marginal sd of Y as the denominator
% obs_mean = posterior * par.mu; % applicable for both 'homo' or 'spatical' cases
% obs_2moment = posterior* (par.mu.^2 + par.sigma.^2);
% obs_sigma = sqrt(obs_2moment - obs_mean.^2);
% Z.soft1 = (obs - obs_mean)./obs_sigma;

% soft2: soft segmentation with wighted sigma's as the denominator
obs_mean = posterior * par.mu; % applicable for both 'homo' or 'spatical' cases
obs_sigma = sqrt(posterior * ((par.sigma.^2)));
Z.soft2 = (obs - obs_mean)./obs_sigma;

% hard segmentation
[~, labels] = max(posterior,[], 2);
obs_mean = par.mu(labels);
obs_sigma = par.sigma(labels);
Z.hard = (obs - obs_mean)./obs_sigma;



end