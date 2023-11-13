function plotwuma( NO,input,FFToutput,SER )
%绘制星座图
figure(2*NO-1);
subplot(2,3,1);
scatter(real(input),imag(input));
title('星座图：输入信号');
subplot(2,3,2);
scatter(real(FFToutput(1,1:5*256)),imag(FFToutput(1,1:5*256)));
title('星座图：输出信号，信噪比0dB');
subplot(2,3,3);
scatter(real(FFToutput(2,1:5*256)),imag(FFToutput(2,1:5*256)));
title('星座图：输出信号，信噪比5dB');
subplot(2,3,4);
scatter(real(FFToutput(3,1:5*256)),imag(FFToutput(3,1:5*256)));
title('星座图：输出信号，信噪比10dB');
subplot(2,3,5);
scatter(real(FFToutput(4,1:5*256)),imag(FFToutput(4,1:5*256)));
title('星座图：输出信号，信噪比15dB');
subplot(2,3,6);
scatter(real(FFToutput(5,1:5*256)),imag(FFToutput(5,1:5*256)));
title('星座图：输出信号，信噪比20dB');

%绘制误码率曲线
figure(2*NO);
x=[0,5,10,15,20];
semilogy(x,SER,'-o');grid on;
end


