function [score] = par2score_2D_marginal(obs, mu, sigma, posterior,type)

K = size(posterior, 2); 
d = size(obs, 2); % univariate or bivariate or higher dimension 

if strcmp(type, 'hard') 
    [~, labels] = max(posterior,[], 2);
    posterior = zeros(size(posterior)); 
    for j = 1:K
        posterior((labels == j), j) = 1; 
    end 
end 


bk.mean = posterior * mu'; 
bk.center = obs - bk.mean; 
marginal_var = zeros(size(mu)); % d by K
for j = 1:K
    marginal_var(:,j) = diag(sigma(:,:,j));
end
bk.sigma = sqrt(posterior * marginal_var');
score = bk.center ./ bk.sigma; 


end 