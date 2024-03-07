file_path = ['ecg/resource/MIT_100_combined_f.mat'];

% ---加载.mat文件中的数据，获取存储在其中的初级融合特征---
load(file_path);

[num_rows, num_cols] = size(combined_f); % 获取数组大小
sliding_window = floor(num_rows/10000) + 1; % 滑动窗口尺寸
n_sliding_window = ceil(num_rows/sliding_window);% 滑动窗口总数

% column_1 = combined_f(:, 2);
% full_column_1 = full(column_1);

% column_combined_f = column_1(1: 14, 1);
% max_val = max(column_combined_f);

mf = NaN(n_sliding_window, num_cols);
% 进行MaxFusion融合降维
for i = 1:num_cols
    column_vector = NaN(n_sliding_window, 1);
    for j = 1:n_sliding_window
        if j * sliding_window > num_rows
            column_combined_f = combined_f(((j - 1) * sliding_window + 1): num_rows, i);
            max_val = max(column_combined_f);
            column_vector(j, 1) = max_val;
        else
            column_combined_f = combined_f((j - 1) * sliding_window + 1: j * sliding_window, i);
            max_val = max(column_combined_f);
            column_vector(j, 1) = max_val;
        end
    end
    mf(:, i) = column_vector;
end
% 将融合后的数据保存到.mat文件中，方面后续操作
folder_path = 'ecg/resource';
file_name = ['MIT_100_mf.mat'];

% 保存 segments 数组到 .mat 文件
save(fullfile(folder_path, file_name), 'mf');
