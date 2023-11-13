dt=0.01;
dw=0.01;
t=-1:dt:3;
w=-50:dw:50;
x1=heaviside(t)+heaviside(t-1)-2*heaviside(t-2);
x2=heaviside(t)-heaviside(t-1);
R=xcorr(x1,x2);
subplot(5,1,1);
t2=-3:dt:5;
plot(t2-1,R);
title('corr');
xlim([-1 3]);
Fw1=dt*exp(w'*t*(-j))*x1';
Fw2=dt*exp(w'*t*(-j))*x2';
subplot(5,1,2);
plot(w,Fw1);
title('x1频谱')
Fr=Fw1.*Fw2;
subplot(5,1,3);
plot(w,Fw2);
title('x2频谱');
subplot(5,1,4);
plot(w,Fr);
title('互相关频谱')
ht=dw*exp(j*t'*w)*Fr/(2*pi);
subplot(5,1,5);
plot(t-1,ht);
title('频率反傅立叶的互相关')
xlim([-1 3]);
ylim([0 2]);
%%
clear all
dt=0.01;
dw=0.01;
t=-3:dt:3;
w=-10:dw:10;
x=sin(t)+randn(1,601);
subplot(6,1,1);
plot(t,x);
title('信号时域谱');
subplot(6,1,2);
[pxx,w2] = periodogram(x);
plot(w2,10*log10(pxx));
title('信号功率谱');
subplot(6,1,3);
H=10./(10+j.*w);
plot(w,H);
title('单位冲激响应频谱');
ht=dw*exp(j*t'*w)*H'/(2*pi);
subplot(6,1,4);
plot(t,ht);
title('单位冲击响应时域谱');
Fx=dt*exp(w'*t*(-j))*x';
h=H.*Fx';
ht=dw*exp(j*t'*w)*h'/(2*pi);
subplot(6,1,5);
plot(t,ht);
title('输出信号时域波形');
[pxx,w] = periodogram(ht);
subplot(6,1,6);
plot(w,10*log10(pxx));
title('输出信号功率谱');
xlim([0 3])