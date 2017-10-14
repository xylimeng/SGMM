all_Nrep = 10^5;
idx_start = 1:(10^3):all_Nrep;
n_patch = length(idx_start);
idx_start(1) = 1;
idx_end = [idx_start(2:end) - 1, all_Nrep];
patch_size = idx_end - idx_start + 1;
fprintf('%d patches: size from %d to %d \n', n_patch, min(patch_size), max(patch_size));
fprintf('idx from %d to %d \n', min(idx_start), max(idx_end));

% % run this block in cluster to obtain summary .mat
% for tumor_indicator = 0:1;
%     for convolution_indicator = 0:1
%         output_mat_file = sprintf('output_mat_%d%d', tumor_indicator, convolution_indicator);
%         load('../setting.mat')
%         [n, K] = size(B_mat);
%         final_left_tail = zeros([n, 6]);
%         final_right_tail = zeros([n, 6]);
%         final_two_tail = zeros([n, 6]);
%         final_count = 0;
%         for ith_patch = 1:n_patch
%             file_mat = sprintf('%s/patch_%d.mat', output_mat_file, ith_patch);
%             if exist(file_mat) == 2
%                 fprintf('%s loaded \n', file_mat)
%                 load(file_mat)
%                 final_left_tail = final_left_tail + left_tail;
%                 final_right_tail = final_right_tail + right_tail;
%                 final_two_tail = final_two_tail + two_tail;
%                 final_count = final_count + count;
%             end
%         end
%         final_left_tail  = final_left_tail ./final_count;
%         final_right_tail  = final_right_tail ./final_count;
%         final_two_tail  = final_two_tail ./final_count;
%
%         save_mat_file = sprintf('%d%d_summary_tail.mat', tumor_indicator, convolution_indicator);
%         save(save_mat_file, 'final_two_tail', 'final_left_tail', 'final_right_tail', 'final_count');
%     end
% end

%% run this block to plot
% GMM.hard, GMM.soft, SGMM.hard, SGMM.soft, Robust.hard, Robust.soft
addpath('../../SGMM_toolbox')
load('../setting.mat')
for tumor_indicator = 0:1;
    for convolution_indicator = 0:1;
        save_mat_file = sprintf('%d%d_summary_tail.mat', tumor_indicator, convolution_indicator);
        load(save_mat_file)
        name = {'GMM Hard', 'GMM Soft', 'SGMM Hard', 'SGMM Soft', 'RB-SGMM Hard', 'RB-SGMM Soft'};
        
        plot_names = {'left_tail', 'right_tail', 'two_tail'};
        plot_date = {final_left_tail, final_right_tail, final_two_tail};
        plot_id = sprintf('%d%d', tumor_indicator, convolution_indicator);
        
        for ith_plot = 1:3
            tod = plot_date{ith_plot} ./ 0.001;
            % tod(tod > 1) = 1.2;
            tod = (tod > 1);
            m=32; a = hot(m); b = 1-a; diffmap = [b(1:2:m,:); a(1:2:m,:)];
            
            plot_idx = [1,3,5,2,4,6];
            for jth = 1:6
                subplot(2, 3, jth)
                j = plot_idx(jth);
                imagesc(vec2im(tod(:, j), mask_brain), [0, 2]);
                colormap hot; axis image xy off; colorbar
                title(name{j})
                % colormap(diffmap);
                % set(gcf, 'Colormap', redbluecmap)
                colormap hot(8)
            end
            saveas(gcf, sprintf('output_png/%s_%s_test.png', plot_id, plot_names{ith_plot}))
        end
    end
end



