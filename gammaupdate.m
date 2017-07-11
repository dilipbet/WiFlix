function[gamma] = gammaupdate(virtq,V)
gamma = V./virtq;
gamma(gamma>1) = 1;
gamma(gamma<0) = 0;
