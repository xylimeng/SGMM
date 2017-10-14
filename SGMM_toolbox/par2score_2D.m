function [contrast_score] = par2score_2D(obs, mu, sigma, posterior, contrast_coef, type)

K = size(posterior, 2); 
d = size(obs, 2); % univariate or bivariate or higher dimension 

if strcmp(type, 'hard') 
    [~, labels] = max(posterior,[], 2);
    posterior = zeros(size(posterior)); 
    for j = 1:K
        posterior((labels == j), j) = 1; 
    end 
end 

% normalize the contrast coef 
contrast_coef = reshape(contrast_coef, [d, 1])./sqrt(sum(contrast_coef(:).^2));

bk.mean = posterior * mu'; 
bk.center = obs - bk.mean; 

% re-scaled using sigma 
temp_scores = zeros([size(obs), K]); 
for j = 1:K 
    this_sigma = sigma(:,:,j); 
    root_sigma = sqrtm(inv(this_sigma)); 
    temp_scores(:,:,j) = bk.center * root_sigma; 
end 

scores = zeros(size(obs)); 
for j = 1:K
    scores = temp_scores(:, :, j) .* repmat(posterior(:,j), [1, d]) + scores;
end
 
contrast_score = scores * contrast_coef; 
end 