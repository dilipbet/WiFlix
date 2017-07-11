function am = modulo_location(a,X,Y)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


a_real = real(a);

ar_min = min(a_real);

if ar_min<0
    da = ceil(-ar_min/X);
    a_real = a_real+ (da+1)*X;
end

am_x = rem(a_real,X);



a_imag = imag(a);

ai_min = min(a_imag);

if ai_min<0
    da = ceil(-ai_min/Y);
    a_imag = a_imag+ (da+1)*Y;
end

am_y = rem(a_imag,Y);

am = am_x +sqrt(-1)* am_y;


end

