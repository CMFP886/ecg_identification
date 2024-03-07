%% 
clear;
[signal,Fs,tm]=rdsamp('data_set/mit-bih-arrhythmia-database-1.0.0/mit-bih-arrhythmia-database-1.0.0/109',1,10000);
points = length(signal);

level=8; wavename='bior2.6';
ECGsignalM1 = signal(:, 1);
ecgdata=ECGsignalM1;
% figure(2);
% plot(ecgdata(1:points));grid on ;axis tight;axis([1,points,-2,5]);
% title('ԭʼECG�ź�');
%%%%%%%%%%����С���任8��
[C,L]=wavedec(ecgdata,level,wavename);
%%%%%%%��ȡ�߶�ϵ����
A1=appcoef(C,L,wavename,1);
A2=appcoef(C,L,wavename,2);
A3=appcoef(C,L,wavename,3);
A4=appcoef(C,L,wavename,4);
A5=appcoef(C,L,wavename,5);
A6=appcoef(C,L,wavename,6);
A7=appcoef(C,L,wavename,7);
A8=appcoef(C,L,wavename,8);
%%%%%%%��ȡϸ��ϵ��
D1=detcoef(C,L,1);
D2=detcoef(C,L,2);
D3=detcoef(C,L,3);
D4=detcoef(C,L,4);
D5=detcoef(C,L,5);
D6=detcoef(C,L,6);
D7=detcoef(C,L,7);
D8=detcoef(C,L,8);
%%%%%%%%%%%%�ع�
A8=zeros(length(A8),1); %ȥ������Ư��,8���Ƶ��Ϣ
RA7=idwt(A8,D8,wavename);
RA6=idwt(RA7(1:length(D7)),D7,wavename);
RA5=idwt(RA6(1:length(D6)),D6,wavename);
RA4=idwt(RA5(1:length(D5)),D5,wavename);
RA3=idwt(RA4(1:length(D4)),D4,wavename);
RA2=idwt(RA3(1:length(D3)),D3,wavename);
D2=zeros(length(D2),1); %ȥ����Ƶ������2���Ƶ����
RA1=idwt(RA2(1:length(D2)),D2,wavename);
D1=zeros(length(D1),1);%ȥ����Ƶ������1���Ƶ����
DenoisingSignal=idwt(RA1,D1,wavename);
% figure(3);
% plot(DenoisingSignal);
% title('ȥ��������ECG�ź�'); grid on; axis tight;axis([1,points,-2,5]);
clear ecgdata;

level=4;    
sr=360; 
%����ECG�ź�
%load ecgdata.mat;
%load ECGsignalM1.mat;
%load Rsignal.mat
mydata = DenoisingSignal;
ecgdata=mydata';
swa=zeros(4,points);%�洢��ò��Ϣ
swd=zeros(4,points);%�洢ϸ����Ϣ
signal=ecgdata(0*points+1:1*points); %ȡ���ź�

%��С��ϵ���ͳ߶�ϵ��
%��ͨ�˲��� 1/4 3/4 3/4 1/4
%��ͨ�˲��� -1/4 -3/4 3/4 1/4
%��������С��

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
%����ԭ�źźͳ߶�ϵ����С��ϵ��
%figure(10);
%subplot(level+1,1,1);plot(ecgdata(1:points));grid on ;axis tight;
%title('ECG�ź���j=1,2,3,4�߶��µĳ߶�ϵ������');
%for i=1:level
%    subplot(level+1,1,i+1);
%    plot(swa(i,:));axis tight;grid on; xlabel('time');ylabel(strcat('a  ',num2str(i)));
%end
%figure(11);
%subplot(level,1,1); plot(ecgdata(1:points)); grid on;axis tight;
%title('ECG�źż�����j=1,2,3,4�߶��µĳ߶�ϵ����С��ϵ��');
%for i=1:level
%    subplot(level+1,2,2*(i)+1);
%    plot(swa(i,:)); axis tight;grid on;xlabel('time');
%    ylabel(strcat('a   ',num2str(i)));
%    subplot(level+1,2,2*(i)+2);
%    plot(swd(i,:)); axis tight;grid on;
%    ylabel(strcat('d   ',num2str(i)));
%end

%����ԭͼ��С��ϵ��
%figure(12);
%subplot(level,1,1); plot(real(ecgdata(1:points)),'b'); grid on;axis tight;
%title('ECG�źż�����j=1,2,3,4�߶��µ�С��ϵ��');
%for i=1:level
%    subplot(level+1,1,i+1);
%    plot(swd(i,:),'b'); axis tight;grid on;
%    ylabel(strcat('d   ',num2str(i)));
%end

%**************************************����������ֵ��**********************%
ddw=zeros(size(swd));
pddw=ddw;
nddw=ddw;
%С��ϵ���Ĵ���0�ĵ�
posw=swd.*(swd>0);
%б�ʴ���0
pdw=((posw(:,1:points-1)-posw(:,2:points))<0);
%������ֵ��
pddw(:,2:points-1)=((pdw(:,1:points-2)-pdw(:,2:points-1))>0);
%С��ϵ��С��0�ĵ�
negw=swd.*(swd<0);
ndw=((negw(:,1:points-1)-negw(:,2:points))>0);
%������ֵ��
nddw(:,2:points-1)=((ndw(:,1:points-2)-ndw(:,2:points-1))>0);
%������
ddw=pddw|nddw;
ddw(:,1)=1;
ddw(:,points)=1;
%�����ֵ���ֵ,��������0
wpeak=ddw.*swd;
wpeak(:,1)=wpeak(:,1)+1e-10;
wpeak(:,points)=wpeak(:,points)+1e-10;

%�������߶��¼�ֵ��
%figure(13);
%for i=1:level
%    subplot(level,1,i);
%    plot(wpeak(i,:)); axis tight;grid on;
%ylabel(strcat('j=   ',num2str(i)));
%end
%subplot(4,1,1);
%title('ECG�ź���j=1,2,3,4�߶��µ�С��ϵ����ģ����ֵ��');

interva2=zeros(1,points);
intervaqs=zeros(1,points);
Mj1=wpeak(1,:);
Mj3=wpeak(3,:);
Mj4=wpeak(4,:);
%�����߶�3��ֵ��
% figure(14);
% plot (Mj3);
%title('�߶�3��С��ϵ����ģ����ֵ��');

posi=Mj3.*(Mj3>0);
%��������ֵ��ƽ��
thposi=(max(posi(1:round(points/4)))+max(posi(round(points/4):2*round(points/4)))+max(posi(2*round(points/4):3*round(points/4)))+max(posi(3*round(points/4):4*round(points/4))))/4;
posi=(posi>thposi/3);
nega=Mj3.*(Mj3<0);
%�󸺼���ֵ��ƽ��
thnega=(min(nega(1:round(points/4)))+min(nega(round(points/4):2*round(points/4)))+min(nega(2*round(points/4):3*round(points/4)))+min(nega(3*round(points/4):4*round(points/4))))/4;
nega=-1*(nega<thnega/4);
%�ҳ���0��
interva=posi+nega;
loca=find(interva);
for i=1:length(loca)-1
    if abs(loca(i)-loca(i+1))<80
       diff(i)=interva(loca(i))-interva(loca(i+1));
    else
       diff(i)=0;
    end
end
%�ҳ���ֵ��
loca2=find(diff==-2);
%������ֵ��
interva2(loca(loca2(1:length(loca2))))=interva(loca(loca2(1:length(loca2))));
%������ֵ��
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
%*************************��������ֵ�Թ��㡣��R����ֵ������y��QRS����㼰�յ�*******************%
while i<points
    if interva2(i)==-1
       mark1=i;
       i=i+1;
       while(i<points&interva2(i)==0)
          i=i+1;
       end
       mark2=i;
%�󼫴�ֵ�ԵĹ����
       mark3= round((abs(Mj3(mark2))*mark1+mark2*abs(Mj3(mark1)))/(abs(Mj3(mark2))+abs(Mj3(mark1))));
%R������ֵ��
       R_result(j)=mark3-10;%Ϊ��-10������ֵ��
       countR(mark3-10)=1;
%���QRS�����
       kqs=mark3-10;
       markq=0;
     while (kqs>1)&&( markq< 3)
         if Mj1(kqs)~=0
            markq=markq+1;
         end
         kqs= kqs -1;
     end
  countQ(kqs)=-1;
  
%���QRS���յ�  
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
%j=3,�����
   R=find(countR);
%�������
   R_R=R(2:length(R))-R(1:length(R)-1);
   RRmean=mean(R_R);
%������R�����С��0.4RRmeanʱ,ȥ��ֵС��R��
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
%������R���������1.6RRmeanʱ,��С��ֵ,����һ�μ�yR��
for i=2:length(R)
    if (R(i)-R(i-1))>1.6*RRmean
        Mjadjust=wpeak(4,R(i-1)+80:R(i)-80);
        points2=(R(i)-80)-(R(i-1)+80)+1;
%��������ֵ��
        adjustposi=Mjadjust.*(Mjadjust>0);
        adjustposi=(adjustposi>thposi/4);
%�󸺼���ֵ��
        adjustnega=Mjadjust.*(Mjadjust<0);
        adjustnega=-1*(adjustnega<thnega/5);
%������
        interva4=adjustposi+adjustnega;
%�ҳ���0��
        loca3=find(interva4);
        diff2=interva4(loca3(1:length(loca3)-1))-interva4(loca3(2:length(loca3)));
%�����м���ֵ��,�ҳ�����ֵ��
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
%������
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
% --- ����R��λ�÷ָ��ź� ---
% ȷ��ÿ��ʱ����Ӧ���źŲ���������
sample_indices = R_result'; % ��������ȡ��ӽ��Ĳ���������
if sample_indices(1) == 0
    sample_indices(1) = 1;
end
% % ��ȡ�źŵĵ�һͨ��
% signal_1 = signal(:,1);
% ���ֶ���
segments = cell(floor((length(sample_indices) - 1) / 2), 1); % ����һ����Ԫ�����Դ洢��ͬ��
for i = 1:2:length(sample_indices) - 2
    start_index = sample_indices(i);
    end_index = sample_indices(i+2) - 1; % ��������Ҫ��1���Ա�֤���ص�

    % ��������Ƿ񳬳��źŷ�Χ
    if end_index > length(signal)
        end_index = length(signal);
    end

    % ���ź�����ȡ��
    segments{(i + 1) / 2} = signal(start_index:end_index);
end

% % ���Ƶ�һ������
% % ��ȡ��һ������
% first_segment = segments{1};
% 
% % ���ɶ�Ӧ��ʱ������
% time_vector = (0:length(first_segment)-1) / 360; % ������������ת��Ϊʱ�䣨�룩
% 
% figure;
% plot(time_vector, first_segment);
% xlabel('ʱ�䣨�룩');
% ylabel('�źŷ�ֵ');
% title('��һ���ź�');
%%
% --- �Էָ���źŽ��г���������ȡ ---
primary_features = cell(length(segments), 1);
for i = 1:length(segments)
    segment_signal = segments{i}; % ��ȡĳ���ź�

    % ---���źŽ���Hilbert�任---
    hilbert_signal = hilbert(segment_signal); % ����Hilbert�任
    hilbert_imag_signal = imag(hilbert_signal); % ȡ�鲿

%     % ����ͼ��鿴
%     time_vector = (0:length(signal)-1) / 360; % ������������ת��Ϊʱ�䣨�룩
%     figure;
%     subplot (211);
%     plot (time_vector, signal);
%     title ( 'ԭʼecg�ź�' );
%     subplot (212);
%     plot (time_vector, hilbert_imag_signal);
%     title ( 'Hilbert�任�ź�' );
%     ylim ([-1,1]);

    % ---��Hilbert�任����źŽ����ز�����ȡ500��������---
    % ---����Ϊ600��������MIT-BIH���ݿ�һ���źŲ�û����ô�������ڶ������̣���Ч��զ��

    % �µĲ���������
    hilbert_samples = 600;
    estimated_samples = 500;

    % �����ز���
    resampled_hilbert_signal = resample(hilbert_imag_signal, hilbert_samples, length(hilbert_imag_signal));

%     % ����ͼ��鿴
%     time_vector = (0:length(signal)-1) / 360; % ������������ת��Ϊʱ�䣨�룩
%     resampled_time_vector = (0:length(resampled_hilbert_signal)-1) / 360; % ������������ת��Ϊʱ�䣨�룩
%     figure;
%     subplot (211);
%     plot (time_vector, hilbert_imag_signal);
%     title ( 'Hilbert�任�ź�' );
%     subplot (212);
%     plot (resampled_time_vector, resampled_hilbert_signal);
%     title ( '�ز������ź�' );
%     ylim ([-1,1]);

    % ---���źŽ��о����׹��ƣ���ӷ���---
    Fs = length(segment_signal); % �źŵĲ���Ƶ��
    nfft = 512;
    cxn = xcorr(segment_signal, 'unbiased');
    estimated_signal = fft(cxn, nfft);
    Pxx=abs(estimated_signal);% ���㹦�����ܶ�

    % ���ӻ�
    index=0:round(nfft/2-1);
    k=index*Fs/nfft;
    plot_Pxx=10*log10(Pxx(index+1));
%     figure;
%     plot(k,plot_Pxx);
%     title('����Ƶ���׹��ƺ��ͼ��');

    % ---�Թ����׹��ƽ�������ز���---
    plot_Pxx = normalize(plot_Pxx, 'range', [0, 1]);
    resampled_estimated_signal = resample(plot_Pxx, estimated_samples, length(plot_Pxx));

    % ---��������������ƴ��---
    combined_vector = [resampled_hilbert_signal, resampled_estimated_signal];

    % ---���й�һ������---
    resampled_combined_vector = normalize(combined_vector, 'range', [0, 1]);

%     figure;
%     plot(resampled_combined_vector);
%     title('��һ�����ͼ��');

    % �������۵���10*110�ľ���
    combined_matrix = reshape(resampled_combined_vector, 110, 10)';

    primary_features{i} = combined_matrix;
end
%%
% --- �Գ��������������������ȡ ---
% PCANet����
PCANet.NumStages = 2;         % ����
PCANet.PatchSize = [7 7];     % 
PCANet.NumFilters = [9 9];
PCANet.HistBlockSize = [7 7]; 
PCANet.BlkOverLapRatio = 0.5;
PCANet.Pyramid = [];

% [ftrain, V, BlkIdx] = PCANet_train(primary_features,PCANet,1); % BlkIdx serves the purpose of learning block-wise DR projection matrix; e.g., WPCA

V = cell(PCANet.NumStages,1); 

% ---��һ�����---
out_primary_features = primary_features;
Num_primary_features = length(primary_features);
primary_features_Idx = (1:Num_primary_features)';

V{1} = PCA_FilterBank(out_primary_features, PCANet.PatchSize(1), PCANet.NumFilters(1)); % compute PCA filter banks

[out_primary_features, primary_features_Idx] = PCA_output(out_primary_features, primary_features_Idx, PCANet.PatchSize(1), PCANet.NumFilters(1), V{1});

[f, BlkIdx] = HashingHist(PCANet,primary_features_Idx,out_primary_features); % compute the feature of image "idx"

% ---�ڶ������---
second_out_primary_features = out_primary_features;
second_primary_features_Idx = primary_features_Idx;

V{2} = PCA_FilterBank(second_out_primary_features, PCANet.PatchSize(2), PCANet.NumFilters(2)); % compute PCA filter banks

[second_out_primary_features, second_primary_features_Idx] = PCA_output(second_out_primary_features, second_primary_features_Idx, PCANet.PatchSize(2), PCANet.NumFilters(2), V{2});

[second_f, second_BlkIdx] = HashingHist(PCANet,second_primary_features_Idx,second_out_primary_features); % compute the feature of image "idx"

% ---�Ե�һ�� �ڶ��� ��ȡ�������м���---
combined_f = [f; second_f];
%%
% --- ����MaxFusion�㷨���������н�ά ---
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
%%
% ---�Խ�ά����źŽ���Ԥ��
rand_label = rand(1, length(mf(1, :)))';
mf = mf';
mf = sparse(mf);

% ---����.mat�ļ��е����ݣ���ȡѵ���õ�ģ��---
file_path = 'ecg/resource/model.mat';

load(file_path);

[predicted_label, accuracy, decision_values] = predict(rand_label, mf, models, '-q');

% ͶƱ���
most_voted_label = mode(predicted_label);

disp(['Ԥ�����Ϊperson' num2str(most_voted_label)]);