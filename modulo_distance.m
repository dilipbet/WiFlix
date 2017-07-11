function [Dab] = modulo_distance(a,b,X,Y)
% function [Dab] = modulo_distance(a,b,X,Y)
% Given two vectors of locations, a, and b, and a wrap-around tile X by Y
% This function produces a matrix of the pairwise modulo-distances between
% elements of a and b, 
j=sqrt(-1);

a=a(:); b=b(:); b=b.';
am = modulo_location(a,X,Y);
bm = modulo_location(b,X,Y);

A = am * ones(1,length(bm));
B = ones(length(am),1)  * bm;

Dxdir = abs(real(A-B));  Dxwrap = X-abs(real(A-B));  Dx = min(Dxdir, Dxwrap);

Dydir = abs(imag(A-B));  Dywrap = Y-abs(imag(A-B));  Dy = min(Dydir, Dywrap);

Dab = abs(Dx+j*Dy);

end

