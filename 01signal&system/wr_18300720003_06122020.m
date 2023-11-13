%% 匹配滤波器的输出信号就是输入信号的自相关函数
dt=0.01;
t=-1:dt:3;
x1=heaviside(t)-heaviside(t-2);
x2=heaviside(t)-2*heaviside(t-1)+heaviside(t-2);
x3=sin(2*pi*t.*t);
[y1 lags1]=xcorr(x1);
[y2 lags2]=xcorr(x2);
[y3 lags3]=xcorr(x3);
figure
subplot(3,1,1);
plot(lags1,y1);

subplot(3,1,2);
plot(lags2,y2);

subplot(3,1,3);
plot(lags3,y3);
%%
clear all
dt=0.01;
t=0:dt:10;
x1=sin(2*pi*1.6*t);
x2=square(2*pi*0.8*t);
x3=sawtooth(2*pi*0.4*t);
p=zeros([1 50]);
p(1:6)=[1 1 1 1 1 1];
x4li=zeros([1 1001]);
x1li=zeros([1 1001]);
x2li=zeros([1 1001]);
x3li=zeros([1 1001]);

for i= 1:1000
    
    x1li(i)=x1(i)*(mod(i,50)==1);
    x2li(i)=x2(i)*(mod(i,50)==13);
    x3li(i)=x3(i)*(mod(i,100)==52);
    x4li(i)=x1(i)*(mod(i,50)==38);

end
x1pam=conv(x1li,p);
x2pam=conv(x2li,p);
x3pam=conv(x3li,p);
x4pam=conv(x4li,p);
xpam=x1pam+x2pam+x3pam+x4pam;
t2=0:dt:10.49;
figure
subplot(5,1,1);
plot(t2,x1pam);
xlim([0 10]);
title('正弦信号')
subplot(5,1,2);
plot(t2,x2pam)
xlim([0 10]);
title('矩形波')
subplot(5,1,3);
plot(t2,x3pam)
title('锯齿波')
xlim([0 10]);
subplot(5,1,4);
plot(t2,x4pam)
title('正弦信号再采样');
xlim([0 10]);
subplot(5,1,5);
plot(t2,xpam)
title('复用');
xlim([0 10]);