
%Adesoji Bello
%IMPLEMENTATION OF FFT OFDM WITH 16 BLOCK SIZE OF 64 BITS/BLOCK

%BPSK Modulation Parameters------------------------
no_Subcarriers = 64;            %No of subcarriers
no_Databits = 48;               %Data subcarriers
no_Pilotbits = 4;
bw = 20*10^6;                   %20MHz Bandwidth

%------------Derived Parameters-----------------

deltaF = bw/no_Subcarriers;                 %Sub-carrier bandwidth
symbol_Period = 1/deltaF;                           %Total symbol period for 64 sub carriers
cp_Delay = symbol_Period/4;                              %Cyclic Prefix delay period
total_Delay = cp_Delay + symbol_Period;                      %Total symbol delay
no_cyclicPrefix = no_Subcarriers * (cp_Delay/symbol_Period);         %No of symbols for CP
data_Subcarriers = no_Databits + no_Pilotbits;                                          %Total no of used sub-carriers


%%
%%-----------Random Bits at the Transmitter-------
original_data = [];
for i =1:1:16

    s = 2 * randi([0 1], 1, no_Subcarriers) - 1;      %Generate N random bits with BPSK Modulation
    original_data = [original_data s];


end
time_signal =[];
for i=1:1:16
    time_signal = [time_signal (no_Subcarriers/sqrt(data_Subcarriers)) * ifft(original_data(((i-1)*no_Subcarriers)+1:no_Subcarriers*i))];
end
% time_signal = (1/sqrt(no_Subcarriers)) * ifft(original_data);
cyclicPref_signal = [];
for i=1:1:16
    xblock = time_signal(((i-1)*no_Subcarriers)+1:no_Subcarriers*i);
    cyclicPref_signal =[cyclicPref_signal [xblock(no_Subcarriers-no_cyclicPrefix+1:no_Subcarriers), xblock]];        %Time modulated baseband signal after adding CP

end

%%
%--------------Baseband to Passband---------------
subcarrier_delay = total_Delay/(no_Subcarriers+no_cyclicPrefix);
t = subcarrier_delay/50: subcarrier_delay/50: subcarrier_delay;      %Iteration of time reduced to 10
fc = 2.4 * 10^9;                    %Center frequency of WLAN 2.4GHz
car_real =[];
car_imag =[];
carrier = [];
[M,P] = size(cyclicPref_signal);                  %Size of the modulated signal with the addition of Cyclic Prefix

%P = Length of the OFDM Signal after adding the Cyclic Prefix   #1280

for i =1:P
   car_real = real(cyclicPref_signal(i))* cos(2*pi*fc*t);   %Modulate real part (In Phase Part)
    car_imag = imag(cyclicPref_signal(i))* sin(2*pi*fc*t);   %Modulate imaginary part (Quadrature Part)
  carrier = [carrier car_real+car_imag];        %Concatenation of real & imaginary part for carrier transmission
end



%%
%------------Channel Modelling & Receiver System  ----------------------------
%An AWGN channel is modelled with a varying SNR
n = 11;
 num_errors = [];
 ber = [];
 snr = zeros(1,n);
 rcv_wocp = [];
 rcv_data_true = [];
 k =0;
 rcvsignal_freq = [];
 for m =1:1:n
    snr = k;
    rcv_carrier = awgn(carrier,snr,'measured');

    %Inphase & Quadrature down-conversion to obtain the baseband signal
    r = [];
    r_real = [];
    r_imag = [];

    for j = 1:1:P
        %In Phase Part
        r_in = rcv_carrier((j-1)*length(t)+1: j*length(t)) .* cos(2*pi*fc*t);
        r_in_int = (trapz(t,r_in)) * (2/subcarrier_delay);  %Integration over half a symbol period for the real part
        r_real = r_in_int;
    
         %Quadrature Part
         r_qd = rcv_carrier((j-1)*length(t)+1: j*length(t)) .* sin(2*pi*fc*t);
         r_qd_int = (trapz(t,r_qd)) * (2/subcarrier_delay);  %Integration over half a symbol period for the real part
         r_imag = r_qd_int;
    
         r = [r r_real+ i*r_imag];               %Received Signal Vector with CP
    end
    v = no_Subcarriers+16;
    
  %Removal of CYCLIC PREFIX
    for i=1:1:16
         rblock = r(((i-1)*v)+1:v*i);    
         rclip = rblock(no_cyclicPrefix+1:no_cyclicPrefix+no_Subcarriers);  
         rcv_wocp = [rcv_wocp rclip];
    end

    rcv_wo_cyclicprefix = rcv_wocp;             %Receiced signal without cyclic prefix
            

    %FFT Block
    for i=1:1:16
        rcvsignal_freq = [rcvsignal_freq  sqrt(data_Subcarriers)/no_Subcarriers*fft(rcv_wo_cyclicprefix(((i-1)*no_Subcarriers)+1:no_Subcarriers*i))];
    end

    rcv_data = rcvsignal_freq;
    %Demodulation of the recovered signal
    rcv_data(rcv_data<0) = 0;
    rcv_data(rcv_data>0) = 1;
                                     %Signal after Demodulation
    original_data(original_data<0) = 0;
%     received_signal(received_signal<0) = 0;

    [num_errors(m), ber(m)] = biterr(original_data,rcv_data);
    
  
    rcv_data(1:end)=[];
    rcv_wocp(1:end)=[];
    rcvsignal_freq = [];
 k = k+1;
end    
SNR = 0:1:10;
%Plot of BER vs SNR
plot(SNR,ber, 'linewidth',2.0);
xlabel('Eb/N0')
ylabel('BER')
title('Sensitivity of BER vs CFO = 0 for FFT OFDM')
grid on

T = table(SNR',ber');
T.Properties.VariableNames = {'SNR','BER'};
T

