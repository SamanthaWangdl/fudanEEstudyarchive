%%
%实验一
t = (-1:0.01:1)';
impulse = t==0;
unitstep = t>=0;
plot(t,[impulse unitstep])
%%
%实验二
t = (-1:0.01:2)';
z1=2.*func1(t);
z2=func1(t-0.5);
z3=func1(2.*t);
z4=func1(t)+func2(t);
z5=func1(t).*func2(t);
plot(t,[z1,z2,z3,z4,z5]);
function x=func1(t)
x=sin(2*pi*t).*(t>=0);
end
function y=func2(t)
y=exp(-t).*(t>=0);
end
%%
t = (-0.5:0.01:1)';