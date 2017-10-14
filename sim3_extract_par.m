all_Nrep = 10^5;
idx_start = 1:(10^3):all_Nrep;
n_patch = length(idx_start);
idx_start(1) = 1;
idx_end = [idx_start(2:end) - 1, all_Nrep];
patch_size = idx_end - idx_start + 1;
fprintf('%d patches: size from %d to %d \n', n_patch, min(patch_size), max(patch_size));
fprintf('idx from %d to %d \n', min(idx_start), max(idx_end));
% tumor_indicator = 0;
% convolution_indicator = 0;
for tumor_indicator = 0:1;
    for convolution_indicator = 0:1;
        par_all = cell([all_Nrep, 1]);
        output_mat_file = sprintf('output_mat_%d%d', tumor_indicator, convolution_indicator);
        save_mat_file = sprintf('%d%d_all_par', tumor_indicator, convolution_indicator);
        %% compare two par's
        addpath('../../SGMM_toolbox')
        load('../setting.mat')
        % load('patch_26.mat', 'all'); % testing
        for ith_patch = 1:n_patch
            file_mat = sprintf('%s/patch_%d.mat', output_mat_file, ith_patch);
            if exist(file_mat) == 2
                fprintf('%s loaded \n', file_mat)
                load(file_mat)
                % transform to a summary
                par_all(idx_start(ith_patch):idx_end(ith_patch)) = all;
                % end of summary
            end
        end
        save(save_mat_file, 'par_all');
        clear par_all;
    end
end

