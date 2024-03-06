file_path = 'ecg/resource/MIT_100_primary_feature.mat';

% ---加载.mat文件中的数据，获取存储在其中的初级融合特征---
load(file_path);

% PCANet参数
PCANet.NumStages = 2;         % 层数
PCANet.PatchSize = [7 7];     % 
PCANet.NumFilters = [9 9];
PCANet.HistBlockSize = [7 7]; 
PCANet.BlkOverLapRatio = 0.5;
PCANet.Pyramid = [];

% [ftrain, V, BlkIdx] = PCANet_train(primary_features,PCANet,1); % BlkIdx serves the purpose of learning block-wise DR projection matrix; e.g., WPCA

V = cell(PCANet.NumStages,1); 

% ---第一卷积层---
out_primary_features = primary_features;
Num_primary_features = length(primary_features);
primary_features_Idx = (1:Num_primary_features)';

V{1} = PCA_FilterBank(out_primary_features, PCANet.PatchSize(1), PCANet.NumFilters(1)); % compute PCA filter banks

[out_primary_features, primary_features_Idx] = PCA_output(out_primary_features, primary_features_Idx, PCANet.PatchSize(1), PCANet.NumFilters(1), V{1});

[f, BlkIdx] = HashingHist(PCANet,primary_features_Idx,out_primary_features); % compute the feature of image "idx"

% ---第二卷积层---
second_out_primary_features = out_primary_features;
second_primary_features_Idx = primary_features_Idx;

V{2} = PCA_FilterBank(second_out_primary_features, PCANet.PatchSize(2), PCANet.NumFilters(2)); % compute PCA filter banks

[second_out_primary_features, second_primary_features_Idx] = PCA_output(second_out_primary_features, second_primary_features_Idx, PCANet.PatchSize(2), PCANet.NumFilters(2), V{2});

[second_f, second_BlkIdx] = HashingHist(PCANet,second_primary_features_Idx,second_out_primary_features); % compute the feature of image "idx"

% ---对第一层 第二层 提取特征进行级联---
combined_f = [f; second_f];

% 将分割的数据保存到.mat文件中，方面后续操作
folder_path = 'ecg/resource';
file_name = 'MIT_100_combined_f.mat';

% 保存 segments 数组到 .mat 文件
save(fullfile(folder_path, file_name), 'combined_f');