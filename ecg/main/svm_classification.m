%% 
%--- 读取训练集以及测试集 ---
files = {'100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '111', '112', '113', '114', '115', '116', '117', '118', '119', '121', '122', '123', '124', '200', '201', '202', '203', '205', '207', '208', '209', '210', '212', '213', '214', '215', '217', '219', '220', '221', '222', '223', '228', '230', '231', '232', '233', '234'};

data = [];
label = [];

for file_index = 1:length(files)
    file_path = ['ecg/resource/MIT_' files{file_index} '_mf.mat'];

    % ---加载.mat文件中的数据，获取存储在其中的初级融合特征---
    load(file_path);
    
    [n_rows, n_cols] = size(mf);
    
    % 生成标签
    person_labels = ones(1, n_cols) * str2double(files{file_index});
    label = [label person_labels];
    
    data = [data mf];
end

% 将分割的数据保存到.mat文件中，方面后续操作
folder_path = 'ecg/resource';
data_file_name = 'MIT_train_data.mat';
label_file_name = 'MIT_label_data.mat';

% 保存 segments 数组到 .mat 文件
save(fullfile(folder_path, data_file_name), 'data');
save(fullfile(folder_path, label_file_name), 'label');
%% 

% --- SVM训练 ---
% 加载数据
data_file_path = 'ecg/resource/MIT_train_data.mat';
label_file_path = 'ecg/resource/MIT_label_data.mat';

% ---加载.mat文件中的数据，获取训练集测试集---
load(data_file_path);
load(label_file_path);

label = label';

% 随机产生训练集和测试集
data_size = length(label);
n = randperm(data_size);
% 训练集——1/5个样本
train_data = data(:, n(1:data_size / 5));
train_data = sparse(train_data);
train_label = label(n(1:data_size / 5), :);
% 测试集——4/5个样本
test_data = data(:, n(data_size / 5 + 1: end));
test_data = sparse(test_data);
test_data = test_data';
test_label = label(n(data_size / 5 + 1: end), :);

% ---训练---
tic;
models = train(train_label, train_data', '-s 1'); % we use linear SVM classifier (C = 1), calling libsvm library
LinearSVM_TrnTime = toc;
%%

% ---预测---
[predicted_label, accuracy, decision_values] = predict(test_label, test_data, models, ''); % label predictoin by libsvm

% 将融合后的数据保存到.mat文件中，方面后续操作
folder_path = 'ecg/resource';
file_name = 'model.mat';

% 保存训练好的模型到 .mat 文件
save(fullfile(folder_path, file_name), 'models');