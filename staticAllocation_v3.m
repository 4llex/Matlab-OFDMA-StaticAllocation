%%% Simulação de alocação estatica de usuarios em simbolo OFDM
%%% OFDMA with static allocation

TargetSer = 1e-3;                   %% SER Alvo
SNR = 10:2:24;                         %% XXX
N = 1584;                       %% Numero de Subportadoras
b = zeros(1,N);                     %% Vetor de Bits das portadoras / Numerologia 3
Total_bits = zeros(1,length(SNR));  %% Total de bits em um simbolo
Total_bits_medio = zeros(1,length(SNR));
quantizar = 'no';                  %%

%% SNR gap para constelação M-QAM:
Gamma=(1/3)*qfuncinv(TargetSer/4)^2; % Gap to channel capacity M-QAM


%% 
%subPower = 20/1854; % 20 seria a potencia max do sistema de transmissao
% LTE EVA CHANNEL
freq_sample = N*15e3; %30.72e6; sample rate do LTE
EVA_SR3072_Delay           =[0 30 150 310 370 710 1090 1730 2510].*1e-9;
EVA_SR3072_PowerdB_Gain    = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]; %  20*log10(0.39)= -8.1787 => -8.1787 dB -> Voltage-ratio = 0.398107

chan_EVA = rayleighchan((1/(freq_sample)),0,EVA_SR3072_Delay,EVA_SR3072_PowerdB_Gain);        
impulse= [1; zeros(N - 1,1)];  


num_itr = 10000;
for i=1:length(SNR)
     j=0;
    while j<num_itr 
        %H = ones(1,N);%% H ideal
        h = filter(chan_EVA, impulse)';
        H = fft(h,N);
        
        

        % Ideal Channel
        b = log2(1 + ((abs(H).^2)* 10^(SNR(i)/10) )/Gamma);

        % Quantização
        if strcmp(quantizar,'yes')
            b(b<2) = 0;
            b((b>2)&(b<4)) = 2;
            b((b>4)&(b<6)) = 4;
            b((b>6)&(b<8)) = 6;
            b(b>8) = 8;
        end

        Total_bits(i) = Total_bits(i) + sum(b);
        %bm(j) = b;
        
        j = j+1;
    end  
    Total_bits(i) = Total_bits(i)/num_itr;
end


%% Gera graficos de Bits/SNR
figure;
plot(SNR, Total_bits/N, '-o');
title('quantidade bits por SNR');
xlabel('SNR [dB]'); 
ylabel('Bits'); 
grid on;
grid minor;
