% a function from par estimation, obs to marginal scores
function score = par_struct2score(par_struct, obs, B_mat, contrast_coef)

contrast_coef = reshape(contrast_coef, [length(contrast_coef), 1]); 
[n, K] = size(B_mat);
weight = zeros([n, K]);
L = zeros([n, K]);

mu = par_struct.mu;
sigma = par_struct.sigma;


if length(fieldnames(par_struct)) == 4
    gamma = par_struct.gamma;
    PI_mat = P2Pt(B_mat, gamma);
    for j = 1:K
        L(:,j) = mvnpdf(obs, mu(:,j)', sigma(:, :, j));
        s = PI_mat(:,j); % eq.(19) in the United Seg paper
        weight(:,j) = s .* L(:,j);
    end
    den = sum(weight,2);
    weight = weight ./ repmat(den, [1, K]);  % q_{ik} in eq. (18) in the United Seg paper
    
end

if length(fieldnames(par_struct)) == 3
    p = par_struct.p;
    obj_GM = gmdistribution(mu',sigma,p);
    weight = posterior(obj_GM, obs);
end

score.soft = par2score_2D_marginal(obs, mu, sigma, weight, 'soft') * contrast_coef; 
score.soft = score.soft ./ sqrt(var(score.soft)); 
score.hard = par2score_2D_marginal(obs, mu, sigma, weight, 'hard') * contrast_coef; 
score.hard = score.hard ./ sqrt(var(score.hard)); 
end 