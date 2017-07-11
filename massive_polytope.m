clear all
M = 50;
S = 10;                                                                                                                                                                                                                                                   ;
k = 1:S;
k = fliplr(k);
coeff = (factorial(S)./(factorial(k).*factorial(S-k))).*(((M-S+1)/S).^k) -... 
(factorial(S+1)./(factorial(k).*factorial(S-k+1))).*(((M-S)/(S+1)).^k);
coeff = [-((M-S)/(S+1))^(S+1) coeff];
roots = roots(coeff);   
real_roots = real(roots(imag(roots)==0))