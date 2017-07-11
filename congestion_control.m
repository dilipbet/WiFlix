function[d_choice, rate_choice] = congestion_control(xx,rand_start,total_gop,i,virtq1,virtq2,num_users,num_slots_chunk)
if(mod(i,num_slots_chunk)==1) 
rand_start_rep = repmat(rand_start,8,1);
rand_start_rep = mod((rand_start_rep+i-1),total_gop*4);
rand_start_rep(rand_start_rep==0) = total_gop*4;
oppa = cumsum(repmat(ones(8,1),1,num_users),1);
lin_indices1 = 2*total_gop*4*8 + 8*(rand_start_rep-1)+oppa;
lin_indices2 = 3*total_gop*4*8 + 8*(rand_start_rep-1)+oppa;
dvec_t_rep = xx(lin_indices1);
chunk_size_rep = xx(lin_indices2);

objec = repmat(virtq1,8,1).*chunk_size_rep-...
    repmat(virtq2,8,1).*dvec_t_rep;
opt_id = (repmat(min(objec),8,1)==objec);
temp_id = cumsum(opt_id);
opt_id(temp_id>1) = 0;
d_choice = sum(dvec_t_rep.*opt_id);
rate_choice = sum(chunk_size_rep.*opt_id);
else 
    d_choice = zeros(1,num_users);
    rate_choice = zeros(1,num_users);
end