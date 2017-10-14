% unit patch for all simulations with argument (tumor, convolution) or not
% Last modified: 5/29/2015

function [all, two_tail, left_tail, right_tail, count] = unitPatch(tumor_indicator, convolution_indicator, output_mat, patch_start, patch_end)

addpath(genpath('../../SGMM_toolbox/'))
load('../setting.mat')
all = cell([patch_end - patch_start + 1, 1]);
% nrep_patch = patch_end - patch_start + 1;
count = 0;
cutoff = abs(norminv(0.001./2));
cutoff_one_side = abs(norminv(0.001));
[n, K] = size(B_mat);
two_tail = zeros([n, 6]); % 6 cols
left_tail = zeros([n, 6]);
right_tail = zeros([n, 6]);
ret_all = zeros([3, 6, patch_end - patch_start + 1]);
for ith = patch_start:patch_end
    rng(ith);
    count = count + 1;
    obs = parsobs_bi(PI_mat, par);
    % add tumor or not
    if tumor_indicator
        tumor_center = [221, 127];
        tumor_radius = 10;
        tumor_mean = 15; % increase it to have better contrast
        tumor_sd = 1;
        % @ 5/12 we observe that: original EM is robust to the mean
        % the tumor increases the variability of the original EM est.
        % to look at this: use tumor_mean = 25;
        
        [n_row, n_col] = size(mask_brain);
        row_idx = repmat(1:n_row, [n_col, 1])';
        col_idx = repmat(1:n_col, [n_row, 1]);
        tumor_idx = (((row_idx - tumor_center(1)).^2 + (col_idx - tumor_center(2)).^2) <= tumor_radius.^2);
        n_tumor = sum(tumor_idx(:));
        tumor_idx_vec = find(tumor_idx(mask_brain));
        
        tumor_obs = obs;
        tumor_obs(tumor_idx_vec, 1) = randn([n_tumor, 1]) .* tumor_sd + tumor_mean;
        tumor_obs(tumor_idx_vec, 2) = randn([n_tumor, 1]) .* tumor_sd + tumor_mean;
        clear obs;
        obs = tumor_obs;
    end
    
    
    % add convolution or NOT
    if convolution_indicator
        sigma_kernel = 1;
        hsize = 6 .* sigma_kernel;
        gaussian = fspecial('gaussian', hsize, sigma_kernel);
        % h = fspecial('gaussian', hsize, sigma) returns a rotationally symmetric
        % Gaussian lowpass filter of size hsize with standard deviation sigma
        % (positive). hsize can be a vector specifying the number of rows and
        % columns in h, or it can be a scalar, in which case h is a square matrix.
        % Not recommended. Use imgaussfilt or imguassfilt3 instead.
        % We keep using fspecial here since we need 'nanconv' to address NaN's.
        obs1_conv = nanconv(vec2im(obs(:, 1), mask_brain), gaussian, 'nanout');
        obs2_conv = nanconv(vec2im(obs(:, 2), mask_brain), gaussian, 'nanout');
        % has NaN's at the same place of input with the option 'nanout'
        clear obs; 
        obs = [obs1_conv(mask_brain(:)), obs2_conv(mask_brain(:))];
    end
    %% Estimation: No-Robust EM
    par_ini = obs2ini(obs, B_mat);
    [par_hat, weight_hat, ~] = biEM(obs, B_mat, par_ini) ;
    score_soft = par2score_2D(obs, par_hat.mu, par_hat.sigma, weight_hat, [1, -1], 'soft');
    score_hard = par2score_2D(obs, par_hat.mu, par_hat.sigma, weight_hat, [1, -1], 'hard');
    
    %% Robust EM
    tunning = sqrt(chi2inv(0.99, 2)); % Devlin, Gnanadesikan and Kettenring(1981)
    % tunning = sqrt(2 .* (1 - 1/9  + 1/3 * icdf('norm', 0.95))^3); % Campbell (1984, 1985)
    % tunning = 3; % test
    [par_rb, weight_rb, ~] = biEM_Robust(obs, B_mat, par_ini, tunning) ;
    score_rb_soft = par2score_2D(obs, par_rb.mu, par_rb.sigma, weight_rb, [1, -1], 'soft');
    score_rb_hard = par2score_2D(obs, par_rb.mu, par_rb.sigma, weight_rb, [1, -1], 'hard');
    
    %% Tradition GM - GM
    options = statset('MaxIter', 1000);
    obj_GM = gmdistribution.fit(obs, K, 'Replicates', 2, 'Options', options);
    % extract estimates
    par_GM.mu = obj_GM.mu';
    par_GM.sigma = obj_GM.Sigma;
    par_GM.p = obj_GM.PComponents;
    clear posterior %bad naming
    posterior_GM = posterior(obj_GM, obs); % posterior prob
    score_GM_soft = par2score_2D(obs, par_GM.mu, par_GM.sigma, posterior_GM, [1, -1], 'soft');
    score_GM_hard = par2score_2D(obs, par_GM.mu, par_GM.sigma, posterior_GM, [1, -1], 'hard');
    
    % GMM.hard, GMM.soft, SGMM.hard, SGMM.soft
    score = cat(2, score_GM_hard, score_GM_soft, score_hard, score_soft, score_rb_hard, score_rb_soft);
    two_tail = (abs(score) > cutoff) + two_tail;
    left_tail = (score < -cutoff_one_side) + left_tail;
    right_tail = (score > cutoff_one_side) + right_tail;
    
    all{count} = {par_GM, par_hat, par_rb};
    
    ret = zeros([3, 6]);
    ret(1, 1:3) = sqrt(sum((par_hat.mu - par.mu).^2))./sqrt(sum(par.mu.^2));
    ret(2, 1:3) = sqrt(sum((par_rb.mu - par.mu).^2))./sqrt(sum(par.mu.^2));
    ret(3, 1:3) = sqrt(sum((par_GM.mu - par.mu).^2))./sqrt(sum(par.mu.^2));
    
    for k = 1:3
        norm1 = norm(par.sigma(:, :, k));
        ret(1,k+3) = norm(par_hat.sigma(:,:,k) - par.sigma(:, :, k))./norm1;
        ret(2,k+3) = norm(par_rb.sigma(:,:,k) - par.sigma(:, :, k))./norm1;
        ret(3,k+3) = norm(par_GM.sigma(:,:,k) - par.sigma(:, :, k))./norm1;
    end
    
    ret_all(:, :, count) = ret;
    
    save(output_mat, 'all', 'ret_all', 'two_tail', 'right_tail', 'left_tail', 'count', '-v7.3')
end







