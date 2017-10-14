%% function to evaluate the posterior probability based on prior + parameters + observations

% Version: 
%   0.0 - created at 3/12/14

function [posterior] = pars2posterior(B_mat, obs, par)

[n,K] = size(B_mat);
L = zeros([n,K]);
weight = zeros([n,K]);

PI_mat = P2Pt(B_mat, par.gamma);
for j = 1:K
    L(:,j) = pdf('norm', obs, par.mu(j), par.sigma(j));
    s = PI_mat(:,j); % eq.(19) in the United Seg paper
    weight(:,j) = s .* L(:,j);
end
den = sum(weight,2);
weight = weight ./ repmat(den, [1, K]);  % q_{ik} in eq. (18) in the United Seg paper
posterior = weight;
end
