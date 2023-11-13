%% F大调各音展示
dt=1/(8*1000);
t=0:dt:0.5;
f=[349.23 392 440 466.16 523.25 587.33 659.25];

T=[];
for i =1:7
    tmp=sin(2*pi*f(i).*t);
    T=[T tmp];
end


sound(T);
%% 演奏乐曲
tmp6b=sin(2*pi*f(6)/2.*t)
tmp1=sin(2*pi*f(1).*t)
tmp2=sin(2*pi*f(2).*t)
tmp5=sin(2*pi*f(5).*t);
tmp6=sin(2*pi*f(6).*t);
T2=[];
T2=[T2 sin(2*pi*f(5).*[t t t]) tmp6 sin(2*pi*f(2).*[t t t t])  sin(2*pi*f(1).*[t t t]) tmp6b  sin(2*pi*f(2).*[t t t t])];

sound(T2);
%% 包络
bl=t;
len=size(bl);
bl=t;
for i =1:len(2)/3
  bl(i)=i/(len(2)/3);
end
for i =floor(len(2)/3):floor(len(2)/3*2)
  bl(i)=1;
end
for i =floor(len(2)*2/3):len(2)
  bl(i)=3/len(2)*(len(2)-i);
end

a=ones(16,1);
BL=kron(a,bl');
%包络后的声音
T3=BL'.*T2;
figure
subplot(2,1,1);
plot(T3);
ylim([-1.5 1.5]);
title('包络后的时域波形');
Y1=fft(T3);
L=size(T3);
L=L(2);
P2 = abs(Y1/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
[max1,loc]=max(P1);
Fs=1/dt;
f = Fs*(0:(L/2))/L;
subplot(2,1,2);
plot(f,P1);
title('包络后的频率波形');
xlim([250 500]);
sound(T3);
%% 降低八度
% tmp26b=sin(2*pi*f(6).*t)
% tmp21=sin(2*pi*f(1)*2.*t)
% tmp22=sin(2*pi*f(2)*2.*t)
% tmp25=sin(2*pi*f(5)*2.*t);
% tmp26=sin(2*pi*f(6)*2.*t);
% T4=[];
% 
% T4=[T4 tmp25 tmp25 tmp25 tmp26 tmp22 tmp22 tmp22 tmp22 tmp21 tmp21 tmp21 tmp26b tmp22 tmp22 tmp22 tmp22 ];
% sound(T4);
f=[349.23 392 440 466.16 523.25 587.33 659.25];
f=f./2;
tmp6b=sin(2*pi*f(6)/2.*t);
tmp1=sin(2*pi*f(1).*t);
tmp2=sin(2*pi*f(2).*t);
tmp5=sin(2*pi*f(5).*t);
tmp6=sin(2*pi*f(6).*t);
T2=[];
T2=[T2 tmp5 tmp5 tmp5 tmp6 tmp2 tmp2 tmp2 tmp2 tmp1 tmp1 tmp1 tmp6b tmp2 tmp2 tmp2 tmp2 ];

sound(T2'.*BL);
%% 升高八度
f=[349.23 392 440 466.16 523.25 587.33 659.25];
f=f.*2;
tmp6b=sin(2*pi*f(6)/2.*t);
tmp1=sin(2*pi*f(1).*t);
tmp2=sin(2*pi*f(2).*t);
tmp5=sin(2*pi*f(5).*t);
tmp6=sin(2*pi*f(6).*t);
T2=[];
T2=[T2 tmp5 tmp5 tmp5 tmp6 tmp2 tmp2 tmp2 tmp2 tmp1 tmp1 tmp1 tmp6b tmp2 tmp2 tmp2 tmp2 ];
sound(T2'.*BL);
%% 自选音乐
f=[349.23 392 440 466.16 523.25 587.33 659.25];
t4=0:dt:2;
tmp3=sin(2*pi*f(3).*t)+0.2*sin(2*pi*f(3)*2.*t)+0.1*sin(2*pi*f(3)*4.*t);
tmp5=sin(2*pi*f(5).*t)+0.2*sin(2*pi*f(5)*2.*t)+0.1*sin(2*pi*f(5)*4.*t);
tmp5_4=sin(2*pi*f(5).*t4)+0.2*sin(2*pi*f(5)*2.*t4)+0.1*sin(2*pi*f(5)*4.*t4);
tmp6=sin(2*pi*f(6).*t)+0.2*sin(2*pi*f(6)*2.*t)+0.1*sin(2*pi*f(6)*4.*t);
tmp1u=sin(2*pi*f(1)*2.*t)+0.2*sin(2*pi*f(1)*2*2.*t)+0.1*sin(2*pi*f(1)*2*4.*t);
T2=[];
T2=[T2 tmp3 tmp3 tmp3 tmp5 tmp6 tmp1u tmp1u tmp6 tmp5 tmp5 tmp5 tmp6 tmp5_4];
sound(T2'.*BL(4:end));
%% 加入谐波
f=[349.23 392 440 466.16 523.25 587.33 659.25];
tmp6b=sin(2*pi*f(6)/2.*t)+0.2*sin(2*pi*f(6).*t)+0.1*sin(2*pi*f(6)*1.5.*t);
tmp1=sin(2*pi*f(1).*t)+0.2*sin(2*pi*f(1)*2.*t)+0.1*sin(2*pi*f(1)*3.*t);
tmp2=sin(2*pi*f(2).*t)+0.2*sin(2*pi*f(2)*2.*t)+0.1*sin(2*pi*f(2)*3.*t);
tmp5=sin(2*pi*f(5).*t)+0.2*sin(2*pi*f(5)*2.*t)+0.1*sin(2*pi*f(5)*3.*t);
tmp6=sin(2*pi*f(6).*t)+0.2*sin(2*pi*f(6)*2.*t)+0.1*sin(2*pi*f(6)*3.*t);
T2=[];
T2=[T2 tmp5 tmp5 tmp5 tmp6 tmp2 tmp2 tmp2 tmp2 tmp1 tmp1 tmp1 tmp6b tmp2 tmp2 tmp2 tmp2 ];

sound(T2'.*BL);

Y=fft(T2);
L=size(T2);
L=L(2);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs=1/dt;
f = Fs*(0:(L/2))/L;
figure
plot(f,P1) 
grid on
[maxv,maxl]= findpeaks(P1,'minpeakheight',0.03,'minpeakdistance',20);
len=size(maxl);
len=len(2);
for i =1:len
text(f(maxl(i)),P1(maxl(i)),['(',num2str(floor(f(maxl(i)))),')'],'color','b');
end
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
xlim([290 1500]);
ylabel('|P1(f)|')

%% 输入信号

clear all
fall=[];
dt=1/(8*1000);
filename = 'fmt.wav';
%找极值分割
y=audioread(filename);
[maxv,maxl]= findpeaks(y,'minpeakheight',0.2,'minpeakdistance',1200);
figure
hold on
plot(y) ;
plot(maxl,maxv,'*','color','R');   
maxnum=size(maxl);
maxnum=maxnum(1);
y1=y(1:maxl(1));
Y1=fft(y1);
L=size(y1);
L=L(1);
P2 = abs(Y1/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
[max1,loc]=max(P1);
Fs=1/dt;
f = Fs*(0:(L/2))/L;
figure
plot(f,P1);
fall=[fall f(loc)];
for i =1:(maxnum-1)
    len=size(y);
    len=len(1);
    yi=y(maxl(i):maxl(i+1));
    Yi=fft(yi);
    L=size(yi);
    L=L(1);
    P2 = abs(Yi/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    [max1,loc]=max(P1);
    Fs=1/dt;
    f = Fs*(0:(L/2))/L;
    
%     figure
%    
%    plot(f,P1);
%    f1=174
%    text(f1,0.03,['(',num2str(f1),')'],'color','b');
%    text(2*f1,0.01,['(',num2str(2*f1),')'],'color','b');
%    text(3*f1,0.013,['(',num2str(3*f1),')'],'color','b');
%    text(4*f1,0.0,['(',num2str(4*f1),')'],'color','b');
%   
%     title(i+1);
    fall=[fall f(loc)];
end
yfin=y(maxl(maxnum):len);
Yfin=fft(yfin);
L=size(yfin);
L=L(1);
P2 = abs(Yfin/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
[max1,loc]=max(P1);
Fs=1/dt;
f = Fs*(0:(L/2))/L;
fall=[fall f(loc)];


sound(y);
guitar=cell2mat(struct2cell(load('Guitar.MAT')));
figure
plot(guitar);
title('吉他');
y1=y(100:124);
Y=fft(y1);
L=size(y1);
L=L(1);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs=1/dt;
f = Fs*(0:(L/2))/L;
figure
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')

ylabel('|P1(f)|')
%%
dt=1/(8*1000);
t=0:dt:0.5;
f=[349.23 392 440 466.16 523.25 587.33 659.25];
f=f./2;
tmp6b=sin(2*pi*f(6)/2.*t)+0.28*sin(2*pi*f(6)/2*2.*t)+0.11*sin(2*pi*f(6)/2*3.*t)+0.1*sin(2*pi*f(6)/2*4.*t);
tmp1=sin(2*pi*f(1).*t)+0.28*sin(2*pi*f(1)*2.*t)+0.11*sin(2*pi*f(1)*3.*t)+0.1*sin(2*pi*f(1)*4.*t);
tmp2=sin(2*pi*f(2).*t)+0.28*sin(2*pi*f(2)*2.*t)+0.11*sin(2*pi*f(2)*3.*t)+0.1*sin(2*pi*f(2)*4.*t);
tmp5=sin(2*pi*f(5).*t)+0.28*sin(2*pi*f(5)*2.*t)+0.11*sin(2*pi*f(5)*3.*t)+0.1*sin(2*pi*f(5)*4.*t);
tmp6=sin(2*pi*f(6).*t)+0.28*sin(2*pi*f(6)*2.*t)+0.11*sin(2*pi*f(6)*3.*t)+0.1*sin(2*pi*f(6)*4.*t);
T2=[];
T2=[T2 tmp5 tmp5 tmp5 tmp6 tmp2 tmp2 tmp2 tmp2 tmp1 tmp1 tmp1 tmp6b tmp2 tmp2 tmp2 tmp2 ];
sound(T2);
%加入包络
x1=1:4001;
x2=1:8002;
x4=1:4*4001;
BL2=100*exp(-x2./300);
BL1=100*exp(-x1./300);
BL4=100*exp(-x4./300);
BLg=[];
BLg=[BL2 BL1 BL1 BL4 BL2 BL1 BL1 BL4];
sound(BLg.*T2);
plot(BLg.*T2);
%% 各个音符的调整
dt=1/(8*1000);
t=0:dt:0.5;
f=[349.23 392 440 466.16 523.25 587.33 659.25];
f=f./2;
tmp6b=sin(2*pi*f(6)/2.*t)+sin(2*pi*f(6)/2*2.*t)+0.43*sin(2*pi*f(6)/2*3.*t)+0.43*sin(2*pi*f(6)/2*4.*t);
tmp1=sin(2*pi*f(1).*t)+0.22*sin(2*pi*f(1)*2.*t)+0.44*sin(2*pi*f(1)*3.*t)+0.15*sin(2*pi*f(1)*4.*t);
tmp2=sin(2*pi*f(2).*t)+0.33*sin(2*pi*f(2)*2.*t)+0.43*sin(2*pi*f(2)*3.*t);
tmp5=sin(2*pi*f(5).*t)+0.22*sin(2*pi*f(5)*2.*t)+0.44*sin(2*pi*f(5)*3.*t)+0.15*sin(2*pi*f(5)*4.*t);
tmp6=sin(2*pi*f(6).*t)+sin(2*pi*f(6)*2.*t)+0.43*sin(2*pi*f(6)*3.*t)+0.43*sin(2*pi*f(6)*4.*t);
T2=[];
T2=[T2 tmp5 tmp5 tmp5 tmp6 tmp2 tmp2 tmp2 tmp2 tmp1 tmp1 tmp1 tmp6b tmp2 tmp2 tmp2 tmp2 ];
sound(T2);
%加入包络
x1=1:4001;
x2=1:8002;
x4=1:4*4001;
BL2=100*exp(-x2./300);
BL1=100*exp(-x1./300);
BL4=100*exp(-x4./300);
BLg=[];
BLg=[BL2 BL1 BL1 BL4 BL2 BL1 BL1 BL4];
sound(BLg.*T2);
plot(BLg.*T2);
