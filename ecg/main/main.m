%% 
clear;
[signal,Fs,tm]=rdsamp('data_set/mit-bih-arrhythmia-database-1.0.0/mit-bih-arrhythmia-database-1.0.0/109',1,10000);
points = length(signal);

level=8; wavename='bior2.6';
ECGsignalM1 = signal(:, 1);
ecgdata=ECGsignalM1;
% figure(2);
% plot(ecgdata(1:points));grid on ;axis tight;axis([1,points,-2,5]);
% title('原始ECG信号');
%%%%%%%%%%进行小波变换8层
[C,L]=wavedec(ecgdata,level,wavename);
%%%%%%%提取尺度系数，
A1=appcoef(C,L,wavename,1);
A2=appcoef(C,L,wavename,2);
A3=appcoef(C,L,wavename,3);
A4=appcoef(C,L,wavename,4);
A5=appcoef(C,L,wavename,5);
A6=appcoef(C,L,wavename,6);
A7=appcoef(C,L,wavename,7);
A8=appcoef(C,L,wavename,8);
%%%%%%%提取细节系数
D1=detcoef(C,L,1);
D2=detcoef(C,L,2);
D3=detcoef(C,L,3);
D4=detcoef(C,L,4);
D5=detcoef(C,L,5);
D6=detcoef(C,L,6);
D7=detcoef(C,L,7);
D8=detcoef(C,L,8);
%%%%%%%%%%%%重构
A8=zeros(length(A8),1); %去除基线漂移,8层低频信息
RA7=idwt(A8,D8,wavename);
RA6=idwt(RA7(1:length(D7)),D7,wavename);
RA5=idwt(RA6(1:length(D6)),D6,wavename);
RA4=idwt(RA5(1:length(D5)),D5,wavename);
RA3=idwt(RA4(1:length(D4)),D4,wavename);
RA2=idwt(RA3(1:length(D3)),D3,wavename);
D2=zeros(length(D2),1); %去除高频噪声，2层高频噪声
RA1=idwt(RA2(1:length(D2)),D2,wavename);
D1=zeros(length(D1),1);%去除高频噪声，1层高频噪声
DenoisingSignal=idwt(RA1,D1,wavename);
% figure(3);
% plot(DenoisingSignal);
% title('去除噪声的ECG信号'); grid on; axis tight;axis([1,points,-2,5]);
clear ecgdata;

level=4;    
sr=360; 
%读入ECG信号
%load ecgdata.mat;
%load ECGsignalM1.mat;
%load Rsignal.mat
mydata = DenoisingSignal;
ecgdata=mydata';
swa=zeros(4,points);%存储概貌信息
swd=zeros(4,points);%存储细节信息
signal=ecgdata(0*points+1:1*points); %取点信号

%算小波系数和尺度系数
%低通滤波器 1/4 3/4 3/4 1/4
%高通滤波器 -1/4 -3/4 3/4 1/4
%二进样条小波

for i=1:points-3
   swa(1,i+3)=1/4*signal(i+3-2^0*0)+3/4*signal(i+3-2^0*1)+3/4*signal(i+3-2^0*2)+1/4*signal(i+3-2^0*3);
   swd(1,i+3)=-1/4*signal(i+3-2^0*0)-3/4*signal(i+3-2^0*1)+3/4*signal(i+3-2^0*2)+1/4*signal(i+3-2^0*3);
end
j=2;
while j<=level
   for i=1:points-24
     swa(j,i+24)=1/4*swa(j-1,i+24-2^(j-1)*0)+3/4*swa(j-1,i+24-2^(j-1)*1)+3/4*swa(j-1,i+24-2^(j-1)*2)+1/4*swa(j-1,i+24-2^(j-1)*3);
     swd(j,i+24)=-1/4*swa(j-1,i+24-2^(j-1)*0)-3/4*swa(j-1,i+24-2^(j-1)*1)+3/4*swa(j-1,i+24-2^(j-1)*2)+1/4*swa(j-1,i+24-2^(j-1)*3);
   end
   j=j+1;
end
%画出原信号和尺度系数。小波系数
%figure(10);
%subplot(level+1,1,1);plot(ecgdata(1:points));grid on ;axis tight;
%title('ECG信号在j=1,2,3,4尺度下的尺度系数对照');
%for i=1:level
%    subplot(level+1,1,i+1);
%    plot(swa(i,:));axis tight;grid on; xlabel('time');ylabel(strcat('a  ',num2str(i)));
%end
%figure(11);
%subplot(level,1,1); plot(ecgdata(1:points)); grid on;axis tight;
%title('ECG信号及其在j=1,2,3,4尺度下的尺度系数及小波系数');
%for i=1:level
%    subplot(level+1,2,2*(i)+1);
%    plot(swa(i,:)); axis tight;grid on;xlabel('time');
%    ylabel(strcat('a   ',num2str(i)));
%    subplot(level+1,2,2*(i)+2);
%    plot(swd(i,:)); axis tight;grid on;
%    ylabel(strcat('d   ',num2str(i)));
%end

%画出原图及小波系数
%figure(12);
%subplot(level,1,1); plot(real(ecgdata(1:points)),'b'); grid on;axis tight;
%title('ECG信号及其在j=1,2,3,4尺度下的小波系数');
%for i=1:level
%    subplot(level+1,1,i+1);
%    plot(swd(i,:),'b'); axis tight;grid on;
%    ylabel(strcat('d   ',num2str(i)));
%end

%**************************************求正负极大值对**********************%
ddw=zeros(size(swd));
pddw=ddw;
nddw=ddw;
%小波系数的大于0的点
posw=swd.*(swd>0);
%斜率大于0
pdw=((posw(:,1:points-1)-posw(:,2:points))<0);
%正极大值点
pddw(:,2:points-1)=((pdw(:,1:points-2)-pdw(:,2:points-1))>0);
%小波系数小于0的点
negw=swd.*(swd<0);
ndw=((negw(:,1:points-1)-negw(:,2:points))>0);
%负极大值点
nddw(:,2:points-1)=((ndw(:,1:points-2)-ndw(:,2:points-1))>0);
%或运算
ddw=pddw|nddw;
ddw(:,1)=1;
ddw(:,points)=1;
%求出极值点的值,其它点置0
wpeak=ddw.*swd;
wpeak(:,1)=wpeak(:,1)+1e-10;
wpeak(:,points)=wpeak(:,points)+1e-10;

%画出各尺度下极值点
%figure(13);
%for i=1:level
%    subplot(level,1,i);
%    plot(wpeak(i,:)); axis tight;grid on;
%ylabel(strcat('j=   ',num2str(i)));
%end
%subplot(4,1,1);
%title('ECG信号在j=1,2,3,4尺度下的小波系数的模极大值点');

interva2=zeros(1,points);
intervaqs=zeros(1,points);
Mj1=wpeak(1,:);
Mj3=wpeak(3,:);
Mj4=wpeak(4,:);
%画出尺度3极值点
% figure(14);
% plot (Mj3);
%title('尺度3下小波系数的模极大值点');

posi=Mj3.*(Mj3>0);
%求正极大值的平均
thposi=(max(posi(1:round(points/4)))+max(posi(round(points/4):2*round(points/4)))+max(posi(2*round(points/4):3*round(points/4)))+max(posi(3*round(points/4):4*round(points/4))))/4;
posi=(posi>thposi/3);
nega=Mj3.*(Mj3<0);
%求负极大值的平均
thnega=(min(nega(1:round(points/4)))+min(nega(round(points/4):2*round(points/4)))+min(nega(2*round(points/4):3*round(points/4)))+min(nega(3*round(points/4):4*round(points/4))))/4;
nega=-1*(nega<thnega/4);
%找出非0点
interva=posi+nega;
loca=find(interva);
for i=1:length(loca)-1
    if abs(loca(i)-loca(i+1))<80
       diff(i)=interva(loca(i))-interva(loca(i+1));
    else
       diff(i)=0;
    end
end
%找出极值对
loca2=find(diff==-2);
%负极大值点
interva2(loca(loca2(1:length(loca2))))=interva(loca(loca2(1:length(loca2))));
%正极大值点
interva2(loca(loca2(1:length(loca2))+1))=interva(loca(loca2(1:length(loca2))+1));
intervaqs(1:points-10)=interva2(11:points);
countR=zeros(1,1);
countQ=zeros(1,1);
countS=zeros(1,1);
mark1=0;
mark2=0;
mark3=0;
i=1;
j=1;
Rnum=0;
%*************************求正负极值对过零。即R波峰值，并检測出QRS波起点及终点*******************%
while i<points
    if interva2(i)==-1
       mark1=i;
       i=i+1;
       while(i<points&interva2(i)==0)
          i=i+1;
       end
       mark2=i;
%求极大值对的过零点
       mark3= round((abs(Mj3(mark2))*mark1+mark2*abs(Mj3(mark1)))/(abs(Mj3(mark2))+abs(Mj3(mark1))));
%R波极大值点
       R_result(j)=mark3-10;%为何-10？经验值吧
       countR(mark3-10)=1;
%求出QRS波起点
       kqs=mark3-10;
       markq=0;
     while (kqs>1)&&( markq< 3)
         if Mj1(kqs)~=0
            markq=markq+1;
         end
         kqs= kqs -1;
     end
  countQ(kqs)=-1;
  
%求出QRS波终点  
  kqs=mark3-10;
  marks=0;
  while (kqs<points)&&( marks<3)
      if Mj1(kqs)~=0
         marks=marks+1;
      end
      kqs= kqs+1;
  end
  countS(kqs)=-1;
  i=i+60;
  j=j+1;
  Rnum=Rnum+1;
 end
i=i+1;
end

num2=1;
while(num2~=0)
   num2=0;
%j=3,过零点
   R=find(countR);
%过零点间隔
   R_R=R(2:length(R))-R(1:length(R)-1);
   RRmean=mean(R_R);
%当两个R波间隔小于0.4RRmean时,去掉值小的R波
for i=2:length(R)
    if (R(i)-R(i-1))<=0.4*RRmean
        num2=num2+1;
        if signal(R(i))>signal(R(i-1))
            countR(R(i-1))=0;
        else
            countR(R(i))=0;
        end
    end
end
end

num1=2;
while(num1>0)
   num1=num1-1;
   R=find(countR);
   R_R=R(2:length(R))-R(1:length(R)-1);
   RRmean=mean(R_R);
%当发现R波间隔大于1.6RRmean时,减小阈值,在这一段检測R波
for i=2:length(R)
    if (R(i)-R(i-1))>1.6*RRmean
        Mjadjust=wpeak(4,R(i-1)+80:R(i)-80);
        points2=(R(i)-80)-(R(i-1)+80)+1;
%求正极大值点
        adjustposi=Mjadjust.*(Mjadjust>0);
        adjustposi=(adjustposi>thposi/4);
%求负极大值点
        adjustnega=Mjadjust.*(Mjadjust<0);
        adjustnega=-1*(adjustnega<thnega/5);
%或运算
        interva4=adjustposi+adjustnega;
%找出非0点
        loca3=find(interva4);
        diff2=interva4(loca3(1:length(loca3)-1))-interva4(loca3(2:length(loca3)));
%假设有极大值对,找出极大值对
        loca4=find(diff2==-2);
        interva3=zeros(points2,1)';
        for j=1:length(loca4)
           interva3(loca3(loca4(j)))=interva4(loca3(loca4(j)));
           interva3(loca3(loca4(j)+1))=interva4(loca3(loca4(j)+1));
        end
        mark4=0;
        mark5=0;
        mark6=0;
    while j<points2
         if interva3(j)==-1;
            mark4=j;
            j=j+1;
            while(j<points2&interva3(j)==0)
                 j=j+1;
            end
            mark5=j;
%求过零点
            mark6= round((abs(Mjadjust(mark5))*mark4+mark5*abs(Mjadjust(mark4)))/(abs(Mjadjust(mark5))+abs(Mjadjust(mark4))));
            countR(R(i-1)+80+mark6-10)=1;
            j=j+60;
         end
         j=j+1;
     end
    end
 end
end
%% 
% --- 根据R波位置分割信号 ---
% 确定每个时间点对应的信号采样点索引
sample_indices = R_result'; % 四舍五入取最接近的采样点索引
if sample_indices(1) == 0
    sample_indices(1) = 1;
end
% % 获取信号的第一通道
% signal_1 = signal(:,1);
% 划分段落
segments = cell(floor((length(sample_indices) - 1) / 2), 1); % 创建一个单元数组以存储不同段
for i = 1:2:length(sample_indices) - 2
    start_index = sample_indices(i);
    end_index = sample_indices(i+2) - 1; % 结束索引要减1，以保证不重叠

    % 检查索引是否超出信号范围
    if end_index > length(signal)
        end_index = length(signal);
    end

    % 从信号中提取段
    segments{(i + 1) / 2} = signal(start_index:end_index);
end

% % 绘制第一段数据
% % 提取第一段数据
% first_segment = segments{1};
% 
% % 生成对应的时间向量
% time_vector = (0:length(first_segment)-1) / 360; % 将采样点索引转换为时间（秒）
% 
% figure;
% plot(time_vector, first_segment);
% xlabel('时间（秒）');
% ylabel('信号幅值');
% title('第一段信号');
%%
% --- 对分割的信号进行初级特征提取 ---
primary_features = cell(length(segments), 1);
for i = 1:length(segments)
    segment_signal = segments{i}; % 获取某段信号

    % ---对信号进行Hilbert变换---
    hilbert_signal = hilbert(segment_signal); % 进行Hilbert变换
    hilbert_imag_signal = imag(hilbert_signal); % 取虚部

%     % 绘制图像查看
%     time_vector = (0:length(signal)-1) / 360; % 将采样点索引转换为时间（秒）
%     figure;
%     subplot (211);
%     plot (time_vector, signal);
%     title ( '原始ecg信号' );
%     subplot (212);
%     plot (time_vector, hilbert_imag_signal);
%     title ( 'Hilbert变换信号' );
%     ylim ([-1,1]);

    % ---对Hilbert变换后的信号进行重采样，取500个采样点---
    % ---论文为600个，但是MIT-BIH数据库一段信号并没有那么长，现在对其缩短，看效果咋样

    % 新的采样点数量
    hilbert_samples = 600;
    estimated_samples = 500;

    % 进行重采样
    resampled_hilbert_signal = resample(hilbert_imag_signal, hilbert_samples, length(hilbert_imag_signal));

%     % 绘制图像查看
%     time_vector = (0:length(signal)-1) / 360; % 将采样点索引转换为时间（秒）
%     resampled_time_vector = (0:length(resampled_hilbert_signal)-1) / 360; % 将采样点索引转换为时间（秒）
%     figure;
%     subplot (211);
%     plot (time_vector, hilbert_imag_signal);
%     title ( 'Hilbert变换信号' );
%     subplot (212);
%     plot (resampled_time_vector, resampled_hilbert_signal);
%     title ( '重采样后信号' );
%     ylim ([-1,1]);

    % ---对信号进行经典谱估计（间接法）---
    Fs = length(segment_signal); % 信号的采样频率
    nfft = 512;
    cxn = xcorr(segment_signal, 'unbiased');
    estimated_signal = fft(cxn, nfft);
    Pxx=abs(estimated_signal);% 计算功率谱密度

    % 可视化
    index=0:round(nfft/2-1);
    k=index*Fs/nfft;
    plot_Pxx=10*log10(Pxx(index+1));
%     figure;
%     plot(k,plot_Pxx);
%     title('进行频率谱估计后的图像');

    % ---对功率谱估计结果进行重采样---
    plot_Pxx = normalize(plot_Pxx, 'range', [0, 1]);
    resampled_estimated_signal = resample(plot_Pxx, estimated_samples, length(plot_Pxx));

    % ---将两个向量进行拼接---
    combined_vector = [resampled_hilbert_signal, resampled_estimated_signal];

    % ---进行归一化操作---
    resampled_combined_vector = normalize(combined_vector, 'range', [0, 1]);

%     figure;
%     plot(resampled_combined_vector);
%     title('归一化后的图像');

    % 将向量折叠成10*110的矩阵
    combined_matrix = reshape(resampled_combined_vector, 110, 10)';

    primary_features{i} = combined_matrix;
end
%%
% --- 对初级特征进行深层特征提取 ---
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
%%
% --- 进行MaxFusion算法对特征进行降维 ---
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
%%
% ---对降维后的信号进行预测
rand_label = rand(1, length(mf(1, :)))';
mf = mf';
mf = sparse(mf);

% ---加载.mat文件中的数据，获取训练好的模型---
file_path = 'ecg/resource/model.mat';

load(file_path);

[predicted_label, accuracy, decision_values] = predict(rand_label, mf, models, '-q');

% 投票表决
most_voted_label = mode(predicted_label);

disp(['预测身份为person' num2str(most_voted_label)]);