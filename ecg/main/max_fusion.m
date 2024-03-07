file_path = ['ecg/resource/MIT_100_combined_f.mat'];

% ---����.mat�ļ��е����ݣ���ȡ�洢�����еĳ����ں�����---
load(file_path);

[num_rows, num_cols] = size(combined_f); % ��ȡ�����С
sliding_window = floor(num_rows/10000) + 1; % �������ڳߴ�
n_sliding_window = ceil(num_rows/sliding_window);% ������������

% column_1 = combined_f(:, 2);
% full_column_1 = full(column_1);

% column_combined_f = column_1(1: 14, 1);
% max_val = max(column_combined_f);

mf = NaN(n_sliding_window, num_cols);
% ����MaxFusion�ںϽ�ά
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
% ���ںϺ�����ݱ��浽.mat�ļ��У������������
folder_path = 'ecg/resource';
file_name = ['MIT_100_mf.mat'];

% ���� segments ���鵽 .mat �ļ�
save(fullfile(folder_path, file_name), 'mf');
