%% Generate tables for bivariate simulation
% First run 'sim3_summary.m' to obtain all summary.mat
%% Table 1. parameter estimation of three methods in various scenarios
all_data = {'output_mat_00/summary_patch_1.mat', ...
    'output_mat_10/summary_patch_18.mat', ...
    'output_mat_01/summary_patch_1.mat'};
all_file = {'output_table/table_neither.txt', ...
    'output_table/table_tumor.txt', ...
    'output_table/table_conv.txt'};

for ith = 1:3
    clearvars -except ith all_data all_file
    load(all_data{ith})
    
    report.mean = mean(output, 3);
    report.N = size(output, 3);
    report.sd = sqrt(var(output, 1, 3)./report.N);
    
    [report.row, report.col] = size(report.mean);
    % use fprintf to generate latex tables
    name.row = {'GMM', 'SGMM', 'RB-SGMM'};
    name.col = {'Method', '$\mu_1$', '$\mu_2$', '$\mu_3$', ...
        '$\Sigma_1$', '$\Sigma_2$','$\Sigma_3$', '$\pi$'};
    
    tab = cell([1, report.row + 1]); % including title row
    % title row
    title = name.col{1};
    for col = 2:length(name.col) % name.col is 1 longer than report.col
        title = sprintf('%s & %s', title, name.col{col});
    end
    title = sprintf('%s \\\\ \\hline', title);
    tab{1} = title;
    
    for row = 1:report.row
        tab_row = name.row{row};
        for col = 1:report.col
            tab_add = sprintf('& %.2f (%.2f)', report.mean(row, col), report.sd(row, col));
            tab_row = sprintf('%s %s', tab_row, tab_add);
            if col == report.col
                tab_row = sprintf('%s \\\\', tab_row);
            end
        end
        tab{row + 1} = tab_row;
    end
    
    latex = sprintf('\\begin{tabular}{*{%d}{c}} \\hline \n', length(name.col));
    % latex = '';
    for row = 1:report.row + 1
        latex = sprintf('%s %s \n', latex, tab{row});
    end
    latex = sprintf('%s \\hline \\end{tabular} \n', latex);
    fileID = fopen(all_file{ith}, 'w');
    fprintf(fileID, '%s', latex);
    fclose(fileID);
end

%% Update: Round 2 @ 10/8/17
all_data = {'output_mat_00/summary_patch_1.mat', ...
    'output_mat_10/summary_patch_18.mat'};
all_file = {'R2_neither.txt', 'R2_tumor.txt'};

ret = cell([1, 2]);
for ith = 1:2
    clearvars -except ith all_data all_file ret
    load(all_data{ith})
    
    report.mean = mean(output, 3);
    report.N = size(output, 3);
    report.sd = sqrt(var(output, 1, 3)./report.N);
    
    ret{ith} = report;
end

table = [ret{1}.mean, ret{2}.mean; ret{1}.sd, ret{2}.sd];
table = table([1,4,2,5,3,6], :); 

name.row = {'GMM', 'SE', 'SGMM', 'SE', 'RB-SGMM', 'SE'};
name.col = {'Method', '$\mu_1$', '$\mu_2$', '$\mu_3$', ...
    '$\Sigma_1$', '$\Sigma_2$','$\Sigma_3$', '$\pi$'};
tab = cell([1, 6]); % including title row
for row = 1:6
    tab_row = name.row{row};
    for col = 1:14
        if mod(row, 2) == 1
            tab_add = sprintf('& %.2f', table(row, col));
        else
            tab_add = sprintf('& $(%.2f)$', table(row, col) * 100);
            % tab_add = sprintf('& $(%.2f)$', table(row, col) * 100);
        end
        tab_row = sprintf('%s %s', tab_row, tab_add);
        if col == 14
            tab_row = sprintf('%s \\\\', tab_row);
        end
    end
    tab{row} = tab_row;
end

latex = '';
for row = 1:6
    latex = sprintf('%s %s \n', latex, tab{row});
end
