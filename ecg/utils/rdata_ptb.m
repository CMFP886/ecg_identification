clc; clear all;
% [signal,Fs,tm]=rdsamp('data_set/ptb-diagnostic-ecg-database-1.0.0//ptb-diagnostic-ecg-database-1.0.0/patient002/s0015lre',1,20000);
[signal,Fs,tm]=rdsamp('data_set/mit-bih-arrhythmia-database-1.0.0//mit-bih-arrhythmia-database-1.0.0/100',1,15000);
[ann]=rdann('data_set/ptb-diagnostic-ecg-database-1.0.0//ptb-diagnostic-ecg-database-1.0.0/patient002/s0015lre', 'hea');
% plot(tm,signal)

scales = 6; % �߶ȷ�Χ������ѡ���� 6 ��߶ȷֽ�
wavelet = 'mexh'; % ѡ�� Mexican hat wavelet

% ִ������С���任
coeffs = cwt(signal, 1:scales, wavelet);

% ���ӻ�С���任���
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
    plot(abs(coeffs(i,:))); % ����С��ϵ���ľ���ֵ
    title(['Scale ', num2str(i)]);
end

% ��ȡ����ֵ
abs_coeffs = abs(coeffs(6,:));

% ��ʼ���洢���ֵ����Сֵ������
max_values = [];
min_values = [];

% ��������ֵϵ��������������
for i = 1:length(abs_coeffs)-1
    % �Ƚ�����������Ĵ�С
    if abs_coeffs(i) > abs_coeffs(i+1)
        max_values = [max_values abs_coeffs(i)];
        min_values = [min_values abs_coeffs(i+1)];
    else
        max_values = [max_values abs_coeffs(i+1)];
        min_values = [min_values abs_coeffs(i)];
    end
end

% % �ӵڶ����㿪ʼ�����������ڶ�����
% for i = 2:length(abs_coeffs)-1
%     % ȡ����ǰ�㼰��ǰ���������ֵ
%     previous_value = abs_coeffs(i-1);
%     current_value = abs_coeffs(i);
%     next_value = abs_coeffs(i+1);
%     
%     % �Ƚϲ��ҳ����ֵ����Сֵ
%     max_val = max([previous_value, current_value, next_value]);
%     min_val = min([previous_value, current_value, next_value]);
%     
%     % �����ֵ����Сֵ�洢����Ӧ��������
%     max_values = [max_values, max_val];
%     min_values = [min_values, min_val];
% end

% ȡ�����ֵ����С��ֵ
% �����ֵ����Сֵ��������
sorted_max_values = sort(max_values, 'descend');
sorted_min_values = sort(min_values, 'ascend');

% ���������ֵ
if length(sorted_max_values) >= 8
    max_threshold = mean(sorted_max_values(1:8));
else
    % ������ֵ��������8��������ѡ����������ʽ������ʹ��ȫ�������ֵ
    max_threshold = mean(sorted_max_values);
end

% ������С��ֵ
if length(sorted_min_values) >= 100
    min_threshold = mean(sorted_min_values(1:100));
else
    % �����Сֵ��������100��������ѡ����������ʽ������ʹ��ȫ������Сֵ
    min_threshold = mean(sorted_min_values);
end

% ������ֵ
threshold = (max_threshold - min_threshold) * 0.25;

% ������ע R ������
[r_peaks_row, r_peaks_col] = find(coeffs(6,:) > threshold);

%�������� R ���ļ��
R_peak_intervals = diff(r_peaks_col);

%���� R ����ƽ������
mean_R_interval = mean(R_peak_intervals);

%��߶�λ׼ȷ��
min_index_arr = [];
for i = 1:length(r_peaks_col)-1
    %����ע���������� R �����С�� 0.3 ��ƽ������ʱ
    if R_peak_intervals(i) < 0.3 * mean_R_interval
        %ȥ�����Ƚ�С��һ�� R ��
        [~, min_index] = min(abs_coeffs(r_peaks_col(i:i+1)));
        min_index_arr = [min_index_arr min_index+i-1];
%     %����ע���������� R �������� 1.5 ��ƽ������ʱ
%     elseif R_peak_intervals(i) > 1.5 * mean_R_interval
%         %���ͳ�ʼ��ֵ
%         threshold = threshold * 0.99;
%         [r_peaks_row, r_peaks_col] = find(abs_coeffs > threshold);
% %         %�������� R ���ļ��
% %         R_peak_intervals = diff(r_peaks_col);
% % 
% %         %���� R ����ƽ������
% %         mean_R_interval = mean(R_peak_intervals);
%         
%         min_index_arr = [];
%         i = 1;
    end
end

r_peaks_col(min_index_arr) = [];

% ��ԭʼ�ź��ϱ�� R ������
figure;
plot(signal); hold on;
scatter(r_peaks_col, signal(r_peaks_col), 'ro', 'LineWidth', 1);
title('ԭʼ ECG �ź��� R ��������');
xlabel('������');
ylabel('��ֵ');
legend('ECG �ź�', 'R ������');