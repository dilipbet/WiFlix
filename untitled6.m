arriv = [3 4 5 13 14 15 16 18 19];
N = length(arriv)
blue = zeros(1,N);
a = 1;
for i = 2: 25
    num_arriv = (arriv > i-1)&(arriv <= i);
    num = sum(num_arriv,2);
    blue(i) = blue(i-1)+num;
end
stairs(blue)
hold all
green = zeros(1,15);
green(13) = 1;
green(14) = 2;
green(15) = 3;
stairs(green)