function alpha=alpha_func(d,K,x,t)

a=0.2;

C_k=((tanh(d)+1)/2).^a.*log((tanh(d)+1)/2);

C1=min(abs(C_k));

C1_bar=(1/C1)^(2/3);

alpha=C1_bar*K*(log(2*K/x)/(2*t))^(3/4);

end

