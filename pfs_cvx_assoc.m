function[activity_fractions] = pfs_cvx_assoc(schedulable_rates_matrix,num_cells,num_users,max_users_served)
rates = schedulable_rates_matrix;  
factor = 1./cumsum(ones(num_cells,num_users,max_users_served),3);
cvx_begin

variable alph(num_cells,num_users,max_users_served)
maximize sum(log(sum(sum(alph.*rates,3),1)),2)
0 <= alph <= 1;
sum(sum(alph.*factor,2),3) <= 1;
sum(sum(alph,3),1) <= 1;

cvx_end
activity_fractions = alph; 

