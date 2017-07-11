function[schedulable_rates_matrix, subset_cardinality,received_power,pathloss_SINR,zz]...
    = physical_layer(path_loss,num_cells,num_users,num_antennas,max_users_served,...
    BS_power,noise_power,max_rate)
%***To calculate number of symbols in one chunk duration********




delta = 40;
alpha = 3.5;

%noise_power = 1;
%macro_power = 10^4;
%femto_power = 10^2; % 20 dBm
%BS_user_distance = Dab;


%BS_user_distance = abs(repmat(BS_locations,1,num_users)-repmat(user_locations,num_cells,1)); 
%zz = (BS_user_distance <= 60);
if(max_rate==0)
zz = ones(num_cells,num_users);
else
    zz = (path_loss == repmat(max(path_loss,[],1),num_cells,1));
    if(sum(sum(zz,1)>1)>0)
        stop
    end
end
%path_loss1 = 1./(1+(BS_user_distance/delta).^alpha);
 received_power = repmat(BS_power',1,num_users).*path_loss;
%received_power = 2*ones(1,num_users);
%received_power(1:10) = 12;
%received_power(11:20) = 90;
%received_power = unifrnd(2,20,1,num_users);

%*******RATE CALCULATIONS FROM MASSIVE MIMO***************


signal_plus_int_power = sum(received_power,1);
int_power = repmat(signal_plus_int_power,num_cells,1)- received_power;
pathloss_SINR = received_power./(noise_power + int_power);
pathloss_SINR = pathloss_SINR.*zz;

subset_cardinality = repmat(cumsum(ones(1,1,max_users_served),3),num_cells,num_users); 
dimension_factor = (num_antennas +1 - subset_cardinality)./subset_cardinality;
schedulable_rates_matrix = log2(1 + repmat(pathloss_SINR,[1 1 max_users_served]).*dimension_factor);

%******RATE CALCULATION FROM ERGODIC CAPACITY SINGLE ANTENNA******
% fade = complex(sqrt(1/2)*randn(num_cells,num_users),...
%      sqrt(1/2)*randn(num_cells,num_users));
%  received_power_fade = received_power.*(abs(fade).^2);
%  
%  signal_plus_int_power = sum(received_power_fade,1);
%  int_power = repmat(signal_plus_int_power,num_cells,1)- received_power_fade;
%  pathloss_SINR = received_power_fade./(noise_power + int_power);
%  pathloss_SINR = pathloss_SINR.*zz;
%  
% schedulable_rates_matrix = log2(1+pathloss_SINR);
% subset_cardinality = 1;

% schedulable_rates_matrix = num_symbols*exp(1./pathloss_SINR).*expint(1./pathloss_SINR);
% schedulable_rates_matrix(isnan(schedulable_rates_matrix)) = 0;
% schedulable_rates_matrix(isinf(schedulable_rates_matrix)) = 0;
% subset_cardinality = repmat(cumsum(ones(1,1,max_users_served),3),num_cells,num_users);

