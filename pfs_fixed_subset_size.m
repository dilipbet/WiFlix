function[activity_fixedS] = pfs_fixed_subset_size(schedulable_rates_matrix,num_cells,num_users,max_users_served)
rates = schedulable_rates_matrix(:,:,max_users_served);  
%factor = 1./cumsum(ones(num_cells,num_users,max_users_served),3);
cvx_begin

variable alph(num_cells,num_users)
maximize sum(log(sum(alph.*rates,1)),2)
0 <= alph <= 1;
sum(alph,2) <= max_users_served;
sum(alph,1) <= 1;

cvx_end
activity_fixedS = alph;