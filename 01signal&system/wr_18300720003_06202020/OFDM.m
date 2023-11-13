function [FFToutput,output,error] = OFDM( N,input,n_channel,SN )


%输入信号的串行/并行变换
count=0;
concurrent=zeros(N,256);
for i=1:1:N
    for j=1:1:256
        concurrent(i,j)=input(1,count+j);
    end
    count=count+256;
end

%OFDM系统
%对并行输入信号进行IFFT
OFDMinput=zeros(N,256);
for i=1:1:N
    OFDMinput(i,:)=sqrt(256)*ifft(concurrent(i,:));
end
%对并行输入信号插入保护间隔
OFDMprotect=zeros(N,272);
for i=1:1:N
    for j=1:1:16
        OFDMprotect(i,j)=OFDMinput(i,240+j);
    end
    for j=1:1:256
        OFDMprotect(i,j+16)=OFDMinput(i,j);
    end
end

%对OFDM的输出信号进行多径传播：其中ak^2满足均匀分布
h=zeros(N,n_channel);
for i=1:1:N
   h(i,:)=0.5.^0.5.*complex((randn(1,n_channel)),randn(1,n_channel));
end

y0=zeros(N,272);
temp=zeros(N,272+n_channel-1);
for i=1:1:N
    temp(i,:)=conv(h(i,:),OFDMprotect(i,:));
    y0(i,:)=temp(i,1:272);
end

%Pn=10^((0-SN)/10);%外加噪声的平均功率

y0=awgn(y0,SN,'measured');%将噪声加入到多径传播后产生的y0中

%反OFDM系统
%去除保护间隔，并进行串/并行转换
FOFDM=zeros(N,256);
for i=1:1:N
    for j=1:1:256
        FOFDM(i,j)=y0(i,16+j);
    end
end

%对FOFDM进行DFT或FFT
FFTout=zeros(N,256);
for i=1:1:N
    FFTout(i,:)=1/sqrt(256)*fft(FOFDM(i,:));
end
%均衡
H=zeros(N,256);
for i=1:1:N
    H(i,:)=fft(h(i,:),256);
end
for i=1:1:N
    FFTout(i,:)=FFTout(i,:)./H(i,:);
end


%并/串转换
FFToutput=zeros(1,(256*N));
count=0;
for i=1:1:N
    for j=1:1:256
        FFToutput(1,count+j)=FFTout(i,j);
    end
    count=count+256;
end

%判决
output=zeros(1,(256*N));
for i=1:(256*N)
    if real(FFToutput(1,i))>0
        output(1,i)=1;
    else
        output(1,i)=0;
    end
end

%计算误码率
error=sum(abs(output-input))/(256*N);
end