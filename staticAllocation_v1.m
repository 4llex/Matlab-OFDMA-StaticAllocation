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
        
    H = ones(1,1584);%% H ideal
    
    
    for j=1:N
        
        % Ideal Channel
        b(j) = log2(1 + ((abs(H(j))^2)* 10^(SNR(i)/10) )/Gamma);
        if strcmp(quantizar,'yes')
            if(b(j) < 2)
                b(j) = 0;
            elseif (b(j) < 4)
                b(j) = 2;
            elseif (b(j) < 6)
                b(j) = 4;
            elseif (b(j) < 8)
                b(j) = 6;
            else
                b(j) = 8;
            end
        
        end
        Total_bits(i) = Total_bits(i) + b(j);
        
        % EPA Channel
        %b(i) = log2(1 + ((abs(H).^2).*subPower)/(Gamma*SigmaSqr));
    end
    
    
end


%% Gera graficos de Bits/SNR
figure;
plot(SNR, Total_bits/N, '-o');
title('quantidade bits por SNR');
xlabel('SNR [dB]'); 
ylabel('Bits'); 
grid on;
grid minor;
