function[served_user_rates,served_user_index] = WSRmax(schedulable_rates_matrix,weights,...
    max_users_served,num_users,num_cells,max_rate,zz,i)
%*******3D VISUALIZATION OF OPERATIONS*******
%1st row of 1st page lists the weighted rates that all users can get
%from BS1 if scheduled *alone*.
%1st row of 2nd page lists the weighted rates that all users can get from
%BS1 if scheduled in subsets of size 2. 
% and so on.... 'i'th row of kth page lists weighted rates of users from BS'i'
%when scheduled in subsets of size k.
if(max_rate==0)
weights_rep = repmat(weights,[num_cells 1 max_users_served]);
%size(weights_rep)
%size(schedulable_rates_matrix)
[sorted_weighted_rates,indices] = sort(weights_rep.*schedulable_rates_matrix,2,'descend');
%%
% temp is a matrix which is the indicator of the best *i* users in the *i*
% th page for every BS
temp = true(num_cells,num_users,max_users_served);
for i = 1:max_users_served
  temp(:,i+1:num_users,i)= false;  
end
%%
sorted_subset = sorted_weighted_rates.*temp;% only the best *i* users in the *i*th page
indices_subset = indices.*temp; % indices indicate the user positions of the *i* sorted users (given by sort function in matlab)in the *i* th page
weighted_sums = sum(sorted_subset.*temp,2);%sum the weighted rates of the *i* sorted users in the *i* th page
max_sum = max(weighted_sums,[],3);%calculate the maximum along the pages. This gives the max weighted sum rate
%%
% the goal here is to find the best *cardinality* to serve and then choose
% the *best subset of users* with that cardinality.
cardinality_choice_ind = (weighted_sums == repmat(max_sum,[1 1 max_users_served]));%best *cardinality* chosen by the BS
cardinality_choice_ind (cumsum(cardinality_choice_ind,3) > 1) = false;% in case, more than one cardinality is chosen, choose the smallest among the best
%cardinality_choice = subset_cardinality.*cardinality_choice_ind;
%best_rate_subsets =  sum(sorted_subset.*repmat(cardinality_choice_ind,[1 num_users 1]),3);
best_subset_indices = sum(indices_subset.*repmat(cardinality_choice_ind, [1 num_users 1]),3);%a row vector corresponding to every BS which indicates the indices of the best set of users of best cardinality
cardinality = repmat(sum((best_subset_indices>0),2),1,num_users);% a row vector (corresponding to every BS) with elements equal to the best cardinality 
actual_matrix_indices = (cardinality-1)*num_cells*num_users + (best_subset_indices - 1)*num_cells+cumsum(ones(num_cells,num_users),1);
%actual matrix index gets you the actual index in the 3-D matrix. first
%summand takes you to the correct page. 2nd summand takes you to the
%correct column and the 3rd summand takes you to the correct row.
actual_matrix_indices(best_subset_indices == 0)= 0;%temporarily equate to zero all the users which are not served
tempora = max(max(actual_matrix_indices));
actual_matrix_indices(best_subset_indices == 0)= tempora; %equate all zero entries in the matrix to the max index chosen
served_user_index = false(num_cells,num_users,max_users_served);%an indicator matrix which has a *1/true* at the location where a user is served.
served_user_index(actual_matrix_indices)= true;
%activity(actual_matrix_indices) = activity(actual_matrix_indices) + 1;
%activity1 = activity;
served_user_rates = sum(schedulable_rates_matrix.*served_user_index,3);
else
    served_user_index = false(num_cells,num_users);
    scheduled_user_ind = zeros(num_cells,1);
    for j = 1 : num_cells
        user_ind_list = find(zz(j,:));
        jj = mod(i,length(user_ind_list));
        if(jj==0)
           scheduled_user_ind(j) = user_ind_list(length(user_ind_list));
        else
           scheduled_user_ind(j) = user_ind_list(jj); 
           served_user_index(j,scheduled_user_ind(j)) = true;
        end
    end
    served_user_rates = schedulable_rates_matrix.*served_user_index;
    
end


