dt=0.01;
N1=11;
N2=21;
t=0:dt:2*pi;
y=sin(t);
interval=-1:2/N1:1;
interval2=-1:2/N2:1;
out=quantiz(y,interval);
out2=quantiz(y,interval2);
figure
subplot(2,1,1);
plot(t,2*out./11-1);

hold on
plot(t,y);
plot(t,2*out2./21-1);
title('原信号和量化信号对比');
legend('11采样','原始信号','21采样');
hold off
subplot(2,1,2);
hold on
plot(t,2*out2./21-1-y);
plot(t,2*out./11-1-y);
legend('21采样','11采样');
title('量化误差');
%%
clear all
dt=0.01;
N1=100;
t=-4:dt:4;
y=exp(t);
interval=exp(-4):(exp(4)-exp(-4))/N1:exp(4);
out=quantiz(y,interval);
figure
subplot(2,1,1);
hold on
plot(t,(exp(4)-exp(-4))*out./N1-exp(-4));
plot(t,y);
title('原信号和量化信号对比');
hold off
subplot(2,1,2);
plot(t,(exp(4)-exp(-4))*out./N1-exp(-4)-y);
title('均匀量化误差');

ycom=compand(y,87.56,exp(4),'A/compressor');
interval3=ycom(1):(exp(4)-exp(-4))/N1:ycom(801);
out2=quantiz(ycom,interval3);
out3 = compand(out2,87.56,out2(801),'A/expander')
figure
subplot(2,1,1);
hold on
plot(t,(exp(4)-exp(-4))*out3./N1-exp(-4));
plot(t,y);
title('原信号和非均匀量化信号对比');
hold off
subplot(2,1,2);
plot(t,(exp(4)-exp(-4))*out3./N1-exp(-4)-y);
title('非均匀量化误差');
