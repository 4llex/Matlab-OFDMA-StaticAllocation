%%% Simulação de alocação estatica de usuarios em simbolo OFDM
%%% OFDMA with static allocation

%% Define Numerology
Numerology = 3;

if (Numerology == 1)
     N = 6336;
     sc_per_rb = 48;
     RE = 720;
else
     N = 1584;
     sc_per_rb = 12;
     RE = 720;
end

%%
TargetSer = 1e-3;                           %% SER Alvo
SNR = 6:2:30;                               %% XXX
%N = 6336;                                   %% Numero de Subportadoras
b = zeros(1,N);                             %% Vetor de Bits das portadoras / Numerologia 3
Total_bits = zeros(1,length(SNR));          %% Total de bits em um simbolo
bits_per_rb = zeros(1,length(SNR));         %% qtd media de Bits por RB 
quantizar = 'yes';                          %%
RB = 132;                                   %% qtd de RB
%sc_per_rb = 48;                            %% SubCarriers per RB, depends numerology    
%% SNR gap para constelação M-QAM:
Gamma=(1/3)*qfuncinv(TargetSer/4)^2; % Gap to channel capacity M-QAM


%% 
%subPower = 20/1854; % 20 seria a potencia max do sistema de transmissao
% LTE EVA CHANNEL
freq_sample = 23.76e6;     %N*15e3; %30.72e6; sample rate do LTE
EVA_SR3072_Delay           = [0 30 150 310 370 710 1090 1730 2510].*1e-9;
EVA_SR3072_PowerdB_Gain    = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]; %  20*log10(0.39)= -8.1787 => -8.1787 dB -> Voltage-ratio = 0.398107

chan_EVA = rayleighchan((1/(freq_sample)),0,EVA_SR3072_Delay,EVA_SR3072_PowerdB_Gain);        
impulse= [1; zeros(N - 1,1)];  

mask1 = [ones(1, 44), zeros(1, 132-44)];
mask2 = circshift(mask1, 44);
mask3 = circshift(mask2, 44);


num_itr = 1000;
for i=1:length(SNR)
    i
    j=0;
    while j<num_itr 
        h1 = filter(chan_EVA, impulse)';
        H1 = fft(h1,N);
        %reset(chan_EVA);
        
        h2 = filter(chan_EVA, impulse)';
        H2 = fft(h2,N);
        %reset(chan_EVA);
        
        h3 = filter(chan_EVA, impulse)';
        H3 = fft(h3,N);
        %reset(chan_EVA);
        
        % Calcula Resposta em frequencia média para os 132 RB's
        H1 = rb_h_media(H1, sc_per_rb);
        H2 = rb_h_media(H2, sc_per_rb);
        H3 = rb_h_media(H3, sc_per_rb);
        SNRLIN = 10^(SNR(i)/10);
        P  = 20;
        Pu = P/3;
        
        [subPower1,~, subCapacity1 ] = fcn_waterfilling(Pu, P/(SNRLIN*RB), Gamma, H1, mask1);
        [subPower2,~, subCapacity2 ] = fcn_waterfilling(Pu, P/(SNRLIN*RB), Gamma, H2, mask2);
        [subPower3,~, subCapacity3 ] = fcn_waterfilling(Pu, P/(SNRLIN*RB), Gamma, H3, mask3);
    
        b = subCapacity1 + subCapacity2 +  subCapacity3;

        % Quantização
        if strcmp(quantizar,'yes')
            b(b<2) = 0;
            b((b>2)&(b<4)) = 2;
            b((b>4)&(b<6)) = 4;
            b((b>6)&(b<8)) = 6;
            b(b>8) = 8;
        end

%         figure(2);
%         cla;
%         plot(subPower1); hold on; plot(subPower2); plot(subPower3);
%         drawnow();

%         figure(3);
%         cla;
%         plot(abs(H1)); hold on; plot(abs(H2)); plot(abs(H3));
%         drawnow();
%         pause(1);
        
        Total_bits(i) = Total_bits(i) + sum(b);
        %bm(j) = b;
        
        j = j+1;
    end  
    
    Total_bits(i) = Total_bits(i)/num_itr;
    
    bits_per_rb(i) = (Total_bits(i)/RB)*RE; 
end

%% Saving Vector in a File
if (Numerology == 1)
     Sim.DataSNR = SNR;   
     Sim.DataBPRB = bits_per_rb;
     FileName = strcat('C:\Users\alexrosa\Documents\MATLAB\DynamicAllocation\staticAllocation_num1.mat'); 
     save(FileName,'Sim');  
else     
     Sim.DataSNR = SNR;   
     Sim.DataBPRB = bits_per_rb;
     FileName = strcat('C:\Users\alexrosa\Documents\MATLAB\DynamicAllocation\staticAllocation_num3.mat'); 
     save(FileName,'Sim');
end    

%% Gera graficos de Bits/SNR
figure;
plot(SNR, bits_per_rb, '-xr');
title('Alocação Estática em sistema de multiplo acesso Ortogonal');
xlabel('SNR [dB]'); 
ylabel('Bits/RB'); 
grid on;
grid minor;
