%%% Simulação de alocação estatica de usuarios em simbolo OFDM
%%% OFDMA with static allocation

TargetSer = 1e-3;                           %% SER Alvo
SNR = 5:2:30;                               %% XXX
N = 1584;                                   %% Numero de Subportadoras
b = zeros(1,N);                             %% Vetor de Bits das portadoras / Numerologia 3
Total_bits = zeros(1,length(SNR));          %% Total de bits em um simbolo
bits_per_rb = zeros(1,length(SNR));         %% qtd media de Bits por RB 
quantizar = 'yes';                          %%
RB = 132;                                   %% qtd de RB
sc_per_rb = 12;                             %% SubCarriers per RB, depends numerology    
%% SNR gap para constelação M-QAM:
Gamma=(1/3)*qfuncinv(TargetSer/4)^2; % Gap to channel capacity M-QAM


%% 
%subPower = 20/1854; % 20 seria a potencia max do sistema de transmissao
% LTE EVA CHANNEL
freq_sample = 23.76e6;     %N*15e3; %30.72e6; sample rate do LTE
EVA_SR3072_Delay           =[0 30 150 310 370 710 1090 1730 2510].*1e-9;
EVA_SR3072_PowerdB_Gain    = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]; %  20*log10(0.39)= -8.1787 => -8.1787 dB -> Voltage-ratio = 0.398107

chan_EVA = rayleighchan((1/(freq_sample)),0,EVA_SR3072_Delay,EVA_SR3072_PowerdB_Gain);        
impulse= [1; zeros(N - 1,1)];  


num_itr = 10000;
for i=1:length(SNR)
    j=0;
    while j<num_itr 
        %H = ones(1,N);%% H ideal
        h1 = filter(chan_EVA, impulse)';
        %H1 = fft(h1,N/3);
        H1 = fft(h1,N);
        %reset(chan_EVA);
        
        h2 = filter(chan_EVA, impulse)';
        %H2 = fft(h2,N/3);
        H2 = fft(h2,N);
        %reset(chan_EVA);
        
        h3 = filter(chan_EVA, impulse)';
        %H3 = fft(h3,N/3);
        H3 = fft(h3,N);
        %reset(chan_EVA);
        
        H = horzcat(H1(1:528),H2(529:1056),H3(1057:1584));
        
        % Calcula Resposta em frequencia média para os 132 RB's
        H_rb = rb_h_media(H, sc_per_rb);

        % Ideal Channel
        b = log2(1 + ((abs(H_rb).^2)* 10^(SNR(i)/10) )/Gamma);

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
    
    bits_per_rb(i) = (Total_bits(i)/RB)*12; 
end



%% Gera graficos de Bits/SNR
figure;
plot(SNR, bits_per_rb, '-o');
title('Alocação Estática em sistema de multiplo acesso Ortogonal');
xlabel('SNR [dB]'); 
ylabel('Bits/RB'); 
grid on;
grid minor;
