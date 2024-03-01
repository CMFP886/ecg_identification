clc; clear all;
% [signal,Fs,tm]=rdsamp('data_set/ptb-diagnostic-ecg-database-1.0.0//ptb-diagnostic-ecg-database-1.0.0/patient002/s0015lre',1,20000);
[signal,Fs,tm]=rdsamp('data_set/mit-bih-arrhythmia-database-1.0.0//mit-bih-arrhythmia-database-1.0.0/100',1,15000);
[ann]=rdann('data_set/ptb-diagnostic-ecg-database-1.0.0//ptb-diagnostic-ecg-database-1.0.0/patient002/s0015lre', 'hea');
% plot(tm,signal)

scales = 6; % 尺度范围，这里选择了 6 层尺度分解
wavelet = 'mexh'; % 选择 Mexican hat wavelet

% 执行连续小波变换
coeffs = cwt(signal, 1:scales, wavelet);

% 可视化小波变换结果
figure;
imagesc((1:length(signal))/Fs, scales, abs(coeffs));
set(gca,'YDir','normal');
xlabel('Time (s)');
ylabel('Scale');
title('Continuous Wavelet Transform with Mexican Hat Wavelet (6 Levels)');
colorbar;

figure;
for i = 1:scales
    subplot(scales, 1, i);
    plot(abs(coeffs(i,:))); % 绘制小波系数的绝对值
    title(['Scale ', num2str(i)]);
end

% 求取绝对值
abs_coeffs = abs(coeffs(6,:));

% 初始化存储最大值和最小值的向量
max_values = [];
min_values = [];

% 遍历绝对值系数中相邻两个点
for i = 1:length(abs_coeffs)-1
    % 比较相邻两个点的大小
    if abs_coeffs(i) > abs_coeffs(i+1)
        max_values = [max_values abs_coeffs(i)];
        min_values = [min_values abs_coeffs(i+1)];
    else
        max_values = [max_values abs_coeffs(i+1)];
        min_values = [min_values abs_coeffs(i)];
    end
end

% % 从第二个点开始遍历到倒数第二个点
% for i = 2:length(abs_coeffs)-1
%     % 取出当前点及其前后两个点的值
%     previous_value = abs_coeffs(i-1);
%     current_value = abs_coeffs(i);
%     next_value = abs_coeffs(i+1);
%     
%     % 比较并找出最大值和最小值
%     max_val = max([previous_value, current_value, next_value]);
%     min_val = min([previous_value, current_value, next_value]);
%     
%     % 将最大值和最小值存储到相应的数组中
%     max_values = [max_values, max_val];
%     min_values = [min_values, min_val];
% end

% 取最大阈值和最小阈值
% 对最大值和最小值进行排序
sorted_max_values = sort(max_values, 'descend');
sorted_min_values = sort(min_values, 'ascend');

% 计算最大阈值
if length(sorted_max_values) >= 8
    max_threshold = mean(sorted_max_values(1:8));
else
    % 如果最大值数量不足8个，可以选择其他处理方式，比如使用全部的最大值
    max_threshold = mean(sorted_max_values);
end

% 计算最小阈值
if length(sorted_min_values) >= 100
    min_threshold = mean(sorted_min_values(1:100));
else
    % 如果最小值数量不足100个，可以选择其他处理方式，比如使用全部的最小值
    min_threshold = mean(sorted_min_values);
end

% 计算阈值
threshold = (max_threshold - min_threshold) * 0.25;

% 初步标注 R 波波峰
[r_peaks_row, r_peaks_col] = find(coeffs(6,:) > threshold);

%计算相邻 R 波的间距
R_peak_intervals = diff(r_peaks_col);

%计算 R 波的平均间期
mean_R_interval = mean(R_peak_intervals);

%提高定位准确度
min_index_arr = [];
for i = 1:length(r_peaks_col)-1
    %当标注的两个相邻 R 波间距小于 0.3 倍平均间期时
    if R_peak_intervals(i) < 0.3 * mean_R_interval
        %去除幅度较小的一个 R 波
        [~, min_index] = min(abs_coeffs(r_peaks_col(i:i+1)));
        min_index_arr = [min_index_arr min_index+i-1];
%     %当标注的两个相邻 R 波间距大于 1.5 倍平均间期时
%     elseif R_peak_intervals(i) > 1.5 * mean_R_interval
%         %降低初始阈值
%         threshold = threshold * 0.99;
%         [r_peaks_row, r_peaks_col] = find(abs_coeffs > threshold);
% %         %计算相邻 R 波的间距
% %         R_peak_intervals = diff(r_peaks_col);
% % 
% %         %计算 R 波的平均间期
% %         mean_R_interval = mean(R_peak_intervals);
%         
%         min_index_arr = [];
%         i = 1;
    end
end

r_peaks_col(min_index_arr) = [];

% 在原始信号上标记 R 波波峰
figure;
plot(signal); hold on;
scatter(r_peaks_col, signal(r_peaks_col), 'ro', 'LineWidth', 1);
title('原始 ECG 信号与 R 波波峰标记');
xlabel('样本点');
ylabel('幅值');
legend('ECG 信号', 'R 波波峰');