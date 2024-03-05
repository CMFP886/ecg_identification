file_path = 'ecg/resource/MIT_100.mat';

% ---加载.mat文件中的数据，获取存储在其中的已划分的ecg信号---
load(file_path);
signal = segments{445}; % 获取某段信号

% ---对信号进行Hilbert变换---
hilbert_signal = hilbert(signal); % 进行Hilbert变换
hilbert_imag_signal = imag(hilbert_signal); % 取虚部

% 绘制图像查看
time_vector = (0:length(signal)-1) / 360; % 将采样点索引转换为时间（秒）
figure;
subplot (211);
plot (time_vector, signal);
title ( '原始ecg信号' );
subplot (212);
plot (time_vector, hilbert_imag_signal);
title ( 'Hilbert变换信号' );
ylim ([-1,1]);

% ---对Hilbert变换后的信号进行重采样，取500个采样点---
% ---论文为600个，但是MIT-BIH数据库一段信号并没有那么长，现在对其缩短，看效果咋样

% 新的采样点数量
hilbert_samples = 600;
estimated_samples = 500;

% 进行重采样
resampled_hilbert_signal = resample(hilbert_imag_signal, hilbert_samples, length(hilbert_imag_signal));

% 绘制图像查看
time_vector = (0:length(signal)-1) / 360; % 将采样点索引转换为时间（秒）
resampled_time_vector = (0:length(resampled_hilbert_signal)-1) / 360; % 将采样点索引转换为时间（秒）
figure;
subplot (211);
plot (time_vector, hilbert_imag_signal);
title ( 'Hilbert变换信号' );
subplot (212);
plot (resampled_time_vector, resampled_hilbert_signal);
title ( '重采样后信号' );
ylim ([-1,1]);

% ---对信号进行经典谱估计（间接法）---
Fs = length(signal); % 信号的采样频率
nfft = 512;
cxn = xcorr(signal, 'unbiased');
estimated_signal = fft(cxn, nfft);
Pxx=abs(estimated_signal);% 计算功率谱密度

% 可视化
index=0:round(nfft/2-1);
k=index*Fs/nfft;
plot_Pxx=10*log10(Pxx(index+1));
figure;
plot(k,plot_Pxx);
title('进行频率谱估计后的图像');

% ---对功率谱估计结果进行重采样---
plot_Pxx = normalize(plot_Pxx, 'range', [0, 1]);
resampled_estimated_signal = resample(plot_Pxx, estimated_samples, length(plot_Pxx));

% ---将两个向量进行拼接---
combined_vector = [resampled_hilbert_signal; resampled_estimated_signal];

% ---进行归一化操作---
resampled_combined_vector = normalize(combined_vector, 'range', [0, 1]);

figure;
plot(resampled_combined_vector);
title('归一化后的图像');