%% 
%--- ��ȡѵ�����Լ����Լ� ---
files = {'100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '111', '112', '113', '114', '115', '116', '117', '118', '119', '121', '122', '123', '124', '200', '201', '202', '203', '205', '207', '208', '209', '210', '212', '213', '214', '215', '217', '219', '220', '221', '222', '223', '228', '230', '231', '232', '233', '234'};

data = [];
label = [];

for file_index = 1:length(files)
    file_path = ['ecg/resource/MIT_' files{file_index} '_mf.mat'];

    % ---����.mat�ļ��е����ݣ���ȡ�洢�����еĳ����ں�����---
    load(file_path);
    
    [n_rows, n_cols] = size(mf);
    
    % ���ɱ�ǩ
    person_labels = ones(1, n_cols) * str2double(files{file_index});
    label = [label person_labels];
    
    data = [data mf];
end

% ���ָ�����ݱ��浽.mat�ļ��У������������
folder_path = 'ecg/resource';
data_file_name = 'MIT_train_data.mat';
label_file_name = 'MIT_label_data.mat';

% ���� segments ���鵽 .mat �ļ�
save(fullfile(folder_path, data_file_name), 'data');
save(fullfile(folder_path, label_file_name), 'label');
%% 

% --- SVMѵ�� ---
% ��������
data_file_path = 'ecg/resource/MIT_train_data.mat';
label_file_path = 'ecg/resource/MIT_label_data.mat';

% ---����.mat�ļ��е����ݣ���ȡѵ�������Լ�---
load(data_file_path);
load(label_file_path);

label = label';

% �������ѵ�����Ͳ��Լ�
data_size = length(label);
n = randperm(data_size);
% ѵ��������1/5������
train_data = data(:, n(1:data_size / 5));
train_data = sparse(train_data);
train_label = label(n(1:data_size / 5), :);
% ���Լ�����4/5������
test_data = data(:, n(data_size / 5 + 1: end));
test_data = sparse(test_data);
test_data = test_data';
test_label = label(n(data_size / 5 + 1: end), :);

% ---ѵ��---
tic;
models = train(train_label, train_data', '-s 1'); % we use linear SVM classifier (C = 1), calling libsvm library
LinearSVM_TrnTime = toc;
%%

% ---Ԥ��---
[predicted_label, accuracy, decision_values] = predict(test_label, test_data, models, ''); % label predictoin by libsvm

% ���ںϺ�����ݱ��浽.mat�ļ��У������������
folder_path = 'ecg/resource';
file_name = 'model.mat';

% ����ѵ���õ�ģ�͵� .mat �ļ�
save(fullfile(folder_path, file_name), 'models');