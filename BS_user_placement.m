function[BS_locations, user_locations,num_cells,Dab] = BS_user_placement(macro_side,cell_side,num_users,M,Nx,Ny)
%%********This function places base stations in whichever way you want*****
%num_cells = macro_side/cell_side; % for simple 1-D model
num_cells = Nx*Ny;
[Lxy_t,Lxy,X,Y] = arch_layout(M,Nx,Ny);
% small cell base stations are located at the centre of every voronoi cell
%****positions of femto BS's *****
[Dab] = modulo_distance(Lxy_t,Lxy,X,Y);
BS_locations = reshape(Lxy_t.',1,Nx*Ny);      

 %BS_locations = [0 80]; %for simple 1-D example
%****** USER POSITIONS*************
user_locations = Lxy.';
%***** Users are generated from either a uniform distribution or Gaussian
%distribution*****

%user_locations(1:150) = linspace(1,10,150);
%user_locations(151:200) = linspace(40,79,50);

%user_locations = 20*ones(1,num_users); % simple 1-D model
%user_locations(1:num_users - 10) = 5;
%user_locations(num_users-9:num_users) = 35;







