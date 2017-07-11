function[activity_fractions] = pfs_cvx(schedulable_rates_matrix,num_users,max_users_served)
rates = permute(schedulable_rates_matrix(1,:,:),[3 2 1]);
factor = 1./cumsum(ones(max_users_served, num_users),1);
cvx_begin

variable alph(max_users_served, num_users)
maximize sum(log(sum(alph.*rates,1)),2)
0 <= alph <= 1;
sum(sum(alph.*factor)) <= 1;

cvx_end
activity_fractions = alph;   