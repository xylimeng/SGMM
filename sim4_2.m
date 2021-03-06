all_Nrep = 10^3;
idx_start = 1:(10):all_Nrep; 
n_patch = length(idx_start); 
idx_start(1) = 1;
idx_end = [idx_start(2:end) - 1, all_Nrep]; 
patch_size = idx_end - idx_start + 1; 
fprintf('%d patches: size from %d to %d \n', n_patch, min(patch_size), max(patch_size));
fprintf('idx from %d to %d \n', min(idx_start), max(idx_end)); 

tumor_indicator = 1; 
convolution_indicator = 0; 
output_mat_file = sprintf('output_mat_%d%d', tumor_indicator, convolution_indicator); 
mkdir(output_mat_file)
% run this in hpc cluster 
matlabpool 
parfor ith_patch = 1:n_patch
    fprintf('%dth patch start:', ith_patch)
    output_mat = sprintf('%s/patch_%d.mat', output_mat_file, ith_patch); 
    obs_dir = output_mat_file; 
    unitPatch_mvt(tumor_indicator, convolution_indicator, output_mat, obs_dir, idx_start(ith_patch), idx_end(ith_patch));
    fprintf('%dth patch end. \n', ith_patch)
end    
matlabpool



