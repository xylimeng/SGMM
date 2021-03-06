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
        output_mat_file = sprintf('output_mat_%d%d', tumor_indicator, convolution_indicator);
        
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
                output = zeros([3, 7, length(all)]);
                for rep = 1:length(all)
                    for method = 1:3
                        par2 = all{rep}{method};
                        output(method,1:3, rep) = sqrt(sum((par.mu - par2.mu).^2));
                        for j = 1:3
                            output(method,3 + j, rep) = norm(par2.sigma(:,:,j) - par.sigma(:,:,j));
                        end
                        
                        if method == 1
                            output(method, 7, rep) = norm(repmat(par2.p, [size(B_mat, 1), 1]) - PI_mat);
                        else
                            output(method, 7, rep) = norm(P2Pt(B_mat, par2.gamma) - PI_mat);
                        end
                    end
                end
                % end of summary
                save(sprintf('%s/summary_patch_%d.mat', output_mat_file, ith_patch), 'output');
            end
        end
        
    end
end

