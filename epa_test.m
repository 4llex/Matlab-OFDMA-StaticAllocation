%%%%%%%% Rayleigh channel(t) with multi tap %%%%%%%%
%h1 = ((randn(1,1)+1i*randn(1,1)))/ (sqrt(2) * sqrt(1));
           
   epa = [0.8 0.3 0.1];
   epa_pos = [1 5 8];
   N = 8;
   
   h1 = zeros(1,N);
   h1(epa_pos) = epa .* (randn(1,3)+randn(1,3)*1i)
   