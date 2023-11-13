%%
t = (-10:0.05:10)';
x=(t>-0.5).*(t<1);
h=t./2.*(t>0).*(t<2);
y=conv(x,h)*0.05;
subplot(3,1,1);
t2=(-20:0.05:20)';
plot(t,x);title("x");axis([ -1 4 0 2]);

subplot(3,1,2);
plot(t2,y);title("y");axis([ -1 4 0 2]);

subplot(3,1,3);
plot(t,h);title("h");axis([ -1 4 0 2]);
%%
a=[1 -0.9];
b=[0.05];
t=0:20;
x=(t>=0);
x2=impz(b,a,t);
y1=filtic(b,a,1);
y2=filter(b,a,x);
y3=filter(b,a,x,y1);
subplot(4,1,1)
stem(t,x);
title("输入序列");
xlabel('n');
subplot(4,1,2)
stem(t,y2);
xlabel('n')
title('响应序列1');
subplot(4,1,3)
stem(t,y3);
xlabel('n')
title('响应序列2');
subplot(4,1,4)
stem(t,x2);
xlabel('n')
ylabel('h(n)');
title('响应序列');