function[maxmin_fixedS] = mms_fixed_subset_size(schedulable_rates_matrix,num_cells,num_users,max_users_served)
rates = schedulable_rates_matrix(:,:,max_users_served);
param = 1 - 1000*ones(1,num_users);
%factor = 1./cumsum(ones(num_cells,num_users,max_users_served),3);
cvx_begin

variable alph(num_cells,num_users)
maximize sum_smallest(sum(alph.*rates,1),1)
0 <= alph <= 1;
sum(alph,2) <= max_users_served;
sum(alph,1) <= 1;

cvx_end
maxmin_fixedS = alph;