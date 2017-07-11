function[pathloss] = winner(helper_user_distance)
d = helper_user_distance;
pathloss = zeros(size(d));
f = 5;
% LOS
A = 18.7;
B = 46.8;
C = 20;
shad_var = 9;
dd = d;
dd(dd<3) = 3;
PL_los = A*log10(dd) + B + C*log10(f/5) + sqrt(shad_var)*randn(1);

%NLOS
A_nlos = 36.8;
B_nlos = 43.8;
C_nlos = 20;
shad_var_nlos = 16;
PL_nlos = A_nlos*log10(dd) + B_nlos + C_nlos*log10(f/5) + sqrt(shad_var_nlos)*randn(1);

los_prob = zeros(size(d));
lessthan = (d <= 2.5);
morethan = (d > 2.5);
los_prob(lessthan) = 1;
los_prob(morethan) = 1 - 0.9*((1 - (1.24 - 0.6*log10(d(morethan))).^3).^(1/3));
max(max(los_prob))

N = ones(size(d));
los_components = binornd(N,los_prob);
pathloss((los_components==1)) = PL_los((los_components==1));
pathloss((los_components==0)) = PL_nlos((los_components==0));

