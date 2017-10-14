function [obs, idx_label] = parsobs_bi(template, par) 
% generate obs (n by 2) from template and parameters 
% template: n by K
% par: has fields mu and sigma 

% generate bivariate observation
SparseLabel = mnrnd(1, template); % size same as template : sparse
idx_label = logical(SparseLabel);
n_Normal = sum(SparseLabel);
[n, K] = size(template); 
obs = zeros([n, 2]);
for i = 1:3
    MU = par.mu(:,i)';
    SIGMA = par.sigma(:, :, i);
    test = mvnrnd(MU, SIGMA, n_Normal(i));
    obs(idx_label(:, i), :) = test;
end
end 
