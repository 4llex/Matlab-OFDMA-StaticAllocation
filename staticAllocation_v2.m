%%% Simulação de alocação estatica de usuarios em simbolo OFDM
%%% OFDMA with static allocation

TargetSer = 1e-3;                   %% SER Alvo
SNR = 5:30;                         %% XXX
N = 16;%1584;                       %% Numero de Subportadoras
b = zeros(1,N);                     %% Vetor de Bits das portadoras / Numerologia 3
Total_bits = zeros(1,length(SNR));  %% Total de bits em um simbolo

quantizar = 'yes';                  %%

%% SNR gap para constelação M-QAM:
Gamma=(1/3)*qfuncinv(TargetSer/4)^2; % Gap to channel capacity M-QAM


%% 
%subPower = 20/1854; % 20 seria a potencia max do sistema de transmissao



for i=1:length(SNR)
        
    H = ones(1,N);%% H ideal
    

    % Ideal Channel
    b = log2(1 + ((abs(H).^2)* 10^(SNR(i)/10) )/Gamma);
    
    % EPA Channel
    %b(i) = log2(1 + ((abs(H).^2).*subPower)/(Gamma*SigmaSqr));
    
    if strcmp(quantizar,'yes')
        b(b<2) = 0;
        b((b>2)&(b<4)) = 2;
        b((b>4)&(b<6)) = 4;
        b((b>6)&(b<8)) = 6;
        b(b>8) = 8;
    end
    
    Total_bits(i) = sum(b);
        
end


%% Gera graficos de Bits/SNR
figure;
plot(SNR, Total_bits/N, '-o');
title('quantidade bits por SNR');
xlabel('SNR [dB]'); 
ylabel('Bits'); 
grid on;
grid minor;
