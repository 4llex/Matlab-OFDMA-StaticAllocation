simbolo = rand(1, 1584);

n_RB = 132;
RB = 1;
Hrb = zeros(1,n_RB);

for i=1:12:n_RB*12
   
    Hrb(RB) = sum(simbolo(i:RB*12))/12;
    
    RB=RB+1;
    
end
