%%
N1=2;
N=10;
t=(-N/2):1:N/2;
x=(t>=-N1).*(t<=N1);
a=exp((t'*t).*(-2*pi*j/N))*x'./N;
a=kron([1,1,1]',a)';
w=kron([1,1,1],t);

subplot(3,1,1);
stem(a);
grid on;
N1=1;
N=10;
t=(-N/2):1:N/2;
x=(t>=-N1).*(t<=N1);
a2=exp((t'*t).*(-2*pi*j/N))*x'./N;
a2=kron([1,1,1]',a2)';
w=kron([1,1,1],t);
subplot(3,1,2);
stem(a2);
grid on;

N1=2;
N=20;
t=(-N/2):1:N/2;
x=(t>=-N1).*(t<=N1);
a3=exp((t'*t).*(-2*pi*j/N))*x'./N;
a3=kron([1,1,1]',a3)';
w=kron([1,1,1],t);
subplot(3,1,3);
stem(a3);


%%
%init
dt=0.01;
T=1;T1=1/4;
t=-T/2:dt:T/2;
N=size(t);
N=N(2);
M=80;
%ak
A=exp((-M/2:1:M/2)'*(ceil(-N/2):1:floor(N/2)).*(-dt*j*2*pi/T));
x=(t>=-T1).*(t<=T1);
a=(A*x').*1/N;
subplot(3,1,1);
stem(a);
title('ak');

%recover
N2=20;
B1=exp((-N2/2:1:N2/2)'*(ceil(-N/2):1:floor(N/2)).*(dt*j*2*pi/T));
a1= a(31:51);
x1=B1'*a1;
B2=exp((-M/2:1:M/2)'*(ceil(-N/2):1:floor(N/2)).*(dt*j*2*pi/T));
x2=B2'*a;
subplot(3,1,2);
plot(t,x1);
title('N=20的复原');
subplot(3,1,3);
plot(t,x2);
title('N=80的复原');
