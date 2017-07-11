function[allchunk_reception_profile,served_bits]...
     = reception_time(rate_evolve,size_evolve,num_users,N,num_symbols_time)
 allchunk_reception_profile = inf*ones(N,num_users);
 size_till_now = cumsum(size_evolve,1);
 
 for i = 1:N-1
   
chunk_reception_time = inf*ones(1,num_users); %sum of the amount of bits of all chunks that have been 
%requested by a user until time $i$
chunk_list = repmat(size_till_now(i,:),N*num_symbols_time,1);

served_bits = cumsum(rate_evolve,1);%the list of total amount of service 
% offered until every slot to the 'active queues' at slot $i$.
served_chunks = (served_bits >= chunk_list);               
service_time_index = cumsum(served_chunks,1);



%service_time_index(service_time_index > 1) = 0;
[row,col] = find(service_time_index==1);
chunk_reception_time(col) = ceil(row./num_symbols_time)+1;
allchunk_reception_profile(i,:) = chunk_reception_time; 
i
 end
%actual_trans_profile(i,:) = sum(actual_transbits.*helper_index,1)