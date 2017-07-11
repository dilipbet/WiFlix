function[active_users] = demand(i,active_users,chunk_number_requested)
prob_req = 0.005;
video_req_over = (mod(chunk_number_requested,1000)== 0);
active_users(video_req_over)= 0;
idle_ind = (active_users==0);
n = sum(idle_ind);
active_users(idle_ind)=binornd(1,prob_req,[1,n]);
if(i<1000)
    active_users(1) = 1;
else
    active_users(1) = 0;
end
active_now = (active_users == 1);
newly_active = xor(idle_ind,active_now);
