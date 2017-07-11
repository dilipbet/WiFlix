function [Lxy_t,Lxy,X,Y] = arch_layout(D,N,M)
%function [Lxy_t,Anchorsxy, Non_Anchors,X,Y] = arch_layout(D,N,M)
%  This program generates 
% -- a substrate of size X=N by Y=M*sqrt(3)/2 (with wrap-around)
% --  a (hypothetical) hexagonal lattice X by M elements (Lxy_t)
% -- D nodes arethe uniformly dropped in the X by Y tile
% -- A representative node (out of the D nodes) is chose for each of the
%    hypothetical (MN) nodes --> these representatives correspond to the
%    anchor nodes(Anchorsxy)
% -- the rese of the D nodes are non_anchors (Non_Anchors)
%
X=N;
Y=M*sqrt(3)/2;

Lxy_t =zeros(M,N);
j=sqrt(-1);
for m=0:M-1
    if rem(m,2)==0
        xo = 1/4+j*sqrt(3)/2 * (1/2+(m));
    else
        xo = 3/4+j*sqrt(3)/2 * (1/2+(m));
    end
    Lxy_t(m+1,:) = xo+[0:N-1];  
end
FS=16;
h1=figure(1); clf;
plot([0, X; X X+j*Y; X+j*Y, j*Y; j*Y, 0], 'r'); hold on;
plot(Lxy_t, 'ro');
xlabel('X coordinate','FontSize',FS); ylabel('Y coordinate','FontSize',FS);
title('AP Locations (+)','FontSize',FS);  


%Lxy = (0.9 + 0.3*rand(D,1)) + sqrt(-1)*(0.9 + 0.2*rand(D,1));
Lxy(1:D) = 0.7 + 0.2*rand(D,1) + sqrt(-1)*(0.35 + 0.1*rand(D,1));
%Lxy(D-9:D) = X*rand(10,1)+ sqrt(-1)*Y*rand(10,1);
%Lxy = X*rand(D,1) + sqrt(-1)*Y*rand(D,1);

plot(Lxy, 'b+');
hold off;
axis([0 X 0 Y]);
axis('square')





end
