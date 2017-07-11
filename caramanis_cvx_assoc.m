function[caramanis_fractions] = caramanis_cvx_assoc(schedulable_rates_matrix,num_cells,num_users,max_users_served)
rates = schedulable_rates_matrix(:,:,1);
size(rates)
num_cells;
num_users;
cvx_begin

variable alph(num_cells,num_users)
maximize sum(sum(alph.*log(rates)))+sum(entr(sum(alph,2)),1)
0 <= alph <= 1;
sum(alph,1) == 1;

cvx_end
caramanis_fractions = alph;