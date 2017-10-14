% marginal score and the corresponding contrast
tumor_indicator = 0;
convolution_indicator = 0;
save_mat_file = sprintf('%d%d_all_par', tumor_indicator, convolution_indicator);
load(save_mat_file)
addpath(genpath('../../SGMM_toolbox/'))
load('../setting.mat')

[n, K] = size(B_mat); 
two_tail = zeros([n, 6]); % 6 cols
left_tail = zeros([n, 6]);
right_tail = zeros([n, 6]);

ith = 1; 
rng(ith);
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

% all{count} = {par_GM, par_hat, par_rb};

score_GM = par_struct2score(par_all{ith}{2}, obs, B_mat); % marginal scores 

