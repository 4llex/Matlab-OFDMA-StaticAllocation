function [ Hrb ] = rb_h_media( H )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    qtd_RB = 132;
    RB = 1;
    Hrb = zeros(1,qtd_RB);

    for i=1:12:qtd_RB*12

        Hrb(RB) = sum(H(i:RB*12))/12; %abs??

        RB=RB+1;

    end


end

