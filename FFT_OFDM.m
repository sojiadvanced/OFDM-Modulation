%%IFFT & FFT using FFT OFDM for the sub-caarier signals

fs = 1000; %Sampling frequency 
L = 1500; %Signal duration
t = (0:1500)./fs; % Intervals of sampling period


S = 1*sin(2*pi*100*t) + 1*sin(2*pi*200*t) + 1*sin(2*pi*300*t);

Y = fft(S);
P2 = abs(Y/L);
%P1 = P2(1:L/2+1);
%P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:L)/L;
figure (1)
plot(1000*t(1:50),S(1:50))
title('Original Signal')
xlabel('t (milliseconds)')
ylabel('S(t)')

figure(2)
plot(f,P2)
title('Transformed Signal')
xlabel('f (Hz)')
ylabel('Y(w)')

%figure(3)
%orig = ifft(Y);
%orig = orig(1:L/2+1);
%plot(1000*t(1:50),orig(1:50))
%title('Original Signal')
%xlabel('t (milliseconds)')
%ylabel('S(t)')


