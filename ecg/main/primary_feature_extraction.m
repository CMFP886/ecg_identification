file_path = 'ecg/resource/MIT_100.mat';

% ---����.mat�ļ��е����ݣ���ȡ�洢�����е��ѻ��ֵ�ecg�ź�---
load(file_path);
signal = segments{445}; % ��ȡĳ���ź�

% ---���źŽ���Hilbert�任---
hilbert_signal = hilbert(signal); % ����Hilbert�任
hilbert_imag_signal = imag(hilbert_signal); % ȡ�鲿

% ����ͼ��鿴
time_vector = (0:length(signal)-1) / 360; % ������������ת��Ϊʱ�䣨�룩
figure;
subplot (211);
plot (time_vector, signal);
title ( 'ԭʼecg�ź�' );
subplot (212);
plot (time_vector, hilbert_imag_signal);
title ( 'Hilbert�任�ź�' );
ylim ([-1,1]);

% ---��Hilbert�任����źŽ����ز�����ȡ500��������---
% ---����Ϊ600��������MIT-BIH���ݿ�һ���źŲ�û����ô�������ڶ������̣���Ч��զ��

% �µĲ���������
hilbert_samples = 600;
estimated_samples = 500;

% �����ز���
resampled_hilbert_signal = resample(hilbert_imag_signal, hilbert_samples, length(hilbert_imag_signal));

% ����ͼ��鿴
time_vector = (0:length(signal)-1) / 360; % ������������ת��Ϊʱ�䣨�룩
resampled_time_vector = (0:length(resampled_hilbert_signal)-1) / 360; % ������������ת��Ϊʱ�䣨�룩
figure;
subplot (211);
plot (time_vector, hilbert_imag_signal);
title ( 'Hilbert�任�ź�' );
subplot (212);
plot (resampled_time_vector, resampled_hilbert_signal);
title ( '�ز������ź�' );
ylim ([-1,1]);

% ---���źŽ��о����׹��ƣ���ӷ���---
Fs = length(signal); % �źŵĲ���Ƶ��
nfft = 512;
cxn = xcorr(signal, 'unbiased');
estimated_signal = fft(cxn, nfft);
Pxx=abs(estimated_signal);% ���㹦�����ܶ�

% ���ӻ�
index=0:round(nfft/2-1);
k=index*Fs/nfft;
plot_Pxx=10*log10(Pxx(index+1));
figure;
plot(k,plot_Pxx);
title('����Ƶ���׹��ƺ��ͼ��');

% ---�Թ����׹��ƽ�������ز���---
plot_Pxx = normalize(plot_Pxx, 'range', [0, 1]);
resampled_estimated_signal = resample(plot_Pxx, estimated_samples, length(plot_Pxx));

% ---��������������ƴ��---
combined_vector = [resampled_hilbert_signal; resampled_estimated_signal];

% ---���й�һ������---
resampled_combined_vector = normalize(combined_vector, 'range', [0, 1]);

figure;
plot(resampled_combined_vector);
title('��һ�����ͼ��');