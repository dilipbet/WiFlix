%We assume OFDM parameters identical to LTE.
%hold all
%clear all
%close all
stream = RandStream('mt19937ar','Seed',28);
    RandStream.setGlobalStream(stream);


% This is the code to setup a simple layout for ITA 2014 load balancing
% simulations.

%In this code, we consider a square single macro cell with wrap around
%assumption.

% two-tiers
% single macro cell simulated with the only macro BS located in the center
% of the cell
% wrap around
% square region
% Tier 2 is the small cell BSs, they are located according to PPP. Density
% will be a parameter
% Users will be also located according to PPP, with another denstity
% parameter
% Tier distributions are independent
%
% Small cell ans MAcro cell locations and numbers are chosen once for the
% first frame and kept fixed.

% each frame corresponds to a new sample path for USER PPPs. ie. at each frame
% the number of points is chosen according to the Poisson distribution with
% parameter density*area, and then that many users are uniformly
% distributed in the area.


cell_type = 'simple'
% Define the area
if(strcmp(cell_type,'square'))
    cell_side = 900; %in meters
elseif(strcmp(cell_type,'square9'))
    cell_side = 900; %in meters
    MC_square_number = 1+1*sqrt(-1);
    hotspot_density_multiplier = 10;
elseif(strcmp(cell_type,'simple'))
    cell_side = 80;
    MC_square_number = 1+1*sqrt(-1);
    hotspot_density_multiplier = 10;
else
    missing_cell_type;
end

SC_relative_density = 5;% relative density wrt Macro (which is 1)
UE_relative_density = 500;% relative density wrt Macro (which is 1) % ralted to the toal number of  users
% in Square scenario these many users will be uniformly distributed
% but in square9 scenario, hotspot_density_multiplier = h

% h/(9+h-1) of them will be uniformly distributed in the center small
% square.



alpha = 4; % as in alpha fairness

%---------------------
% subgradient search parameters
a = 1;
b = 0;
imax = 10000;
%----------------------
delta = 0.01; %for reduced LP

%---------------

NumFrames = 1; % A frame corresponds to one layout.

MC_powerdBm = 46; %dBm
SC_powerdBm = 35; %dBm
Noise_powerdBm = -104; %dBM

MC_M = 100;
SC_M = 40;

MC_S = 10;
SC_S = 4;

PathlossModel = 'Simple'

if(strcmp(PathlossModel,'Caramannis'))
    MC_A = 34;
    MC_B = 40;
    SC_A = 37;
    SC_B = 30;
    shad_var = 8;%log normal shadowing deviation
    
    pm1 = MC_A;
    pm2 = MC_B;
    pm3 = SC_A;
    pm4 = SC_B;
    pm5 = shad_var;
elseif(strcmp(PathlossModel,'WINNER'))
    how_to_differentiate_tiers;
    
    pm1 = [];
    pm2 = [];
    pm3 = [];
    pm4 = [];
    pm5 = [];
elseif(strcmp(PathlossModel,'Simple'))
    MC_alpha = 3.5;
    MC_delta = 40;
    SC_alpha = 4;
    SC_delta = 40;
    
    pm1 = MC_alpha;
    pm2 = MC_delta;
    pm3 = SC_alpha;
    pm4 = SC_delta;
    pm5 = [];
else
    missing_pathloss_model;
end

filename = sprintf('alpha%d_celltype_%s_cellside%d_MCs%d+i%d_h%d_SC%d_UE%d_NF%d_MCP%d_SCP%d_NP%d_MCM%d_MCS%d_SCM%d_SCS%d_PLmodel_%s_pm1_%dpm2_%dpm3_%dpm4_%dpm5_%d_subgradient_a%d_b%d_imax%d_delta%d.mat',alpha,cell_type,cell_side,real(MC_square_number),imag(MC_square_number),hotspot_density_multiplier,SC_relative_density,UE_relative_density,NumFrames,MC_powerdBm,SC_powerdBm,Noise_powerdBm,MC_M,MC_S,SC_M,SC_S,PathlossModel,pm1,pm2,pm3,pm4,pm5,a,b,imax,delta)


%%%%%%%%%%%%%% SET SIMULATION PARAMETERS BEFORE THIS LINE  %%%%%%%%%%%%%

MC_power = 10^(0.1*MC_powerdBm);
SC_power = 10^(0.1*SC_powerdBm);
Noise_power = 10^(0.1*Noise_powerdBm);

%[cell_area,MC_location] = MC_related(cell_type,cell_side,MC_square_number);
cell_area = cell_side^2;
SC_density = SC_relative_density/cell_area;
UE_density = UE_relative_density/cell_area;

%Small Cell Locations
%IMPORTANT: At each frame number of SC's might be different.
%SC_NumPoints(f) gives the number of SCs for the frame = f
if(strcmp(cell_type,'square'))
    frame_variant = 1;
    [SC_Locs, SC_NumPoints] =  PointPoissonProcess2D(cell_area,SC_density, frame_variant, cell_side, NumFrames);
elseif(strcmp(cell_type,'square9'))
    frame_variant = 0;
    [SC_Locs, SC_NumPoints] =  PointPoissonProcess2D(cell_area,SC_density, frame_variant, cell_side, NumFrames);
elseif(strcmp(cell_type,'simple'))
    SC_Locs = zeros(1,5);
    SC_Locs(3) = 40+40i;
    SC_Locs(1) = 20+20i;
    SC_Locs(2) = 20+60i;
    SC_Locs(4) = 60+20i;
    SC_Locs(5) = 60+60i;
end

burasi = 1;
%pause
%UE Locations
%IMPORTANT: At each frame number of UE's might be different.
%UE_NumPoints(f) gives the number of UEs for the frame = f

if(strcmp(cell_type,'square'))
    frame_variant = 1;
    [UE_Locs, UE_NumPoints] =  PointPoissonProcess2D(cell_area,UE_density, frame_variant, cell_side, NumFrames);
else
    frame_variant = 1;
    [UE_Locs, UE_NumPoints] =  PointPoissonProcess2D(cell_area,UE_density, frame_variant, cell_side, NumFrames);
    burasi  = 2;
    %pause
    for frame = 1:NumFrames
        % frame = frame
        uni = round(8/(8+hotspot_density_multiplier)*UE_NumPoints(frame));
        hotspot = UE_NumPoints(frame)-uni;
        %update the locs of the initial hotspot many users to the center
        %square
        UE_Locs(frame, 1:hotspot) = cell_side/3+sqrt(-1)*cell_side/3+(rand(1,hotspot)+rand(1,hotspot)*sqrt(-1))*cell_side/3;
    end
end

burasi = 3;



% first BS is the MC
% the rest is the SC's
%BS_Locs = [MC_location*ones(NumFrames,1) SC_Locs];
BS_Locs = SC_Locs;
figure(1)
scatter(real(BS_Locs),imag(BS_Locs),'o','b')
hold on
scatter(real(UE_Locs),imag(UE_Locs),'*','r')
grid
DistanceM = BS_UE_Distance_Calculation(BS_Locs,UE_Locs,NumFrames,cell_side);
path_loss = Pathloss_Calculation(DistanceM,PathlossModel,pm1,pm2,pm3,pm4,pm5);

num_users = length(UE_Locs);
num_cells = length(BS_Locs);
advanced_receiver=0;
max_rate=0;
%**** IMPORT VIDEO DATA*****
total_gop = 200;
v1 = importdata('video1_ssim.csv');
v2 = importdata('video2_ssim.csv');
v3 = importdata('video3_ssim.csv');
v4 = importdata('video4_ssim.csv');
v22 = zeros(total_gop*4,4);
v33 = zeros(total_gop*4,4);
%Since, additional bits of each layer are given in the imported files, to find the total bits
%corresponding to a particular quality, we need to take the cumulative sum.
%That is what we do in the following 'for loop'.
for dd = 1:total_gop
    v1((dd-1)*8+1:(dd-1)*8+8,4) = cumsum(v1((dd-1)*8+1:(dd-1)*8+8,4),1);
    v4((dd-1)*8+1:(dd-1)*8+8,4) = cumsum(v4((dd-1)*8+1:(dd-1)*8+8,4),1);
    v3((dd-1)*4+1:(dd-1)*4+4,4) = cumsum(v3((dd-1)*4+1:(dd-1)*4+4,4),1);
    v2((dd-1)*4+1:(dd-1)*4+4,4) = cumsum(v2((dd-1)*4+1:(dd-1)*4+4,4),1);
end
%v2 and v3 have only 4 layers. To make it compatible with the future code
%in congestion_control.m, we assume all videos have the same fixed number
%of layers = 8. Since videos 2 and 3 have 4 layers, we artificially append
%the 5th 6th 7th and 8th layers which all have the same quality as the 4th
%layer. Thus, the table for v2 and v3 is modified to v22 and v33 in the
%folloing for loop. Now, the new tables v22 and v33 have for each GOP, 4
%more artificial layers of best quality. Last 5 entries of the table will
%be same and will correspond to the best quality.
for tt = 1:total_gop
v22((tt-1)*8+1:(tt-1)*8+4,:) = v2((tt-1)*4+1:(tt-1)*4+4,:);
v22((tt-1)*8+5:(tt-1)*8+8,:) = repmat(v2(tt*4,:),4,1);
v33((tt-1)*8+1:(tt-1)*8+4,:) = v3((tt-1)*4+1:(tt-1)*4+4,:);
v33((tt-1)*8+5:(tt-1)*8+8,:) = repmat(v3(tt*4,:),4,1);
end
%v22 = circshift(v4,100*8+1);
%v33 = circshift(v1,100*8+1);
num_layers = 8*ones(1,num_users); % we make a list of the number of layers in each video
num_layers(1:4) = [8 4 4 8]; 
xx = [v1;v22;v33;v4];% we concatenate the tables of all video trace files to a single table
rand_start = ceil(linspace(16,800,num_users));% This indicates the starting point of
                                                     % video download for
                                                     % every user.
%rand_start(1:4) = [1 201 401 601]; %since we have 4 videos concatenated one after the other
                                   %and there are 4 users each requesting a
                                   %full video with 200 GOP's, we set the
                                   %start points to 1, 201, 401 and 601
                                   %resply


%*****IMPORTING VIDEO DATA ENDS*****************

%*****LTE SYSTEM PARAMETERS*****************

coherence_time = (10)*10^(-3); % in seconds
RB_duration = (0.5)*10^(-3);
chunk_duration = 0.5; %in seconds
schedule_slot = chunk_duration;
system_bandwidth = 18*10^6; %20MHz
subcarrier_bandwidth = 12*(15*10^3);%15 KHz
num_symbols_freq = system_bandwidth/(subcarrier_bandwidth);
num_symbols_time_slot = floor(coherence_time/RB_duration);
num_slots_chunk = floor(chunk_duration/coherence_time);
num_symbols_RB = 7*12;
num_symbols_slot = num_symbols_RB*num_symbols_freq*num_symbols_time_slot;

%************

%macro_side = 40; % length of entire 1D layout
%cell_side = 5; % length of Voronoi cell
num_antennas = 40;
max_users_served = 10;
V = 10*10^14;
N = 201;
a = cumsum(ones(1,num_users),2);
virtq = zeros(1,num_users);
virtq2 = (V/2)*ones(1,num_users);
%[BS_locations, user_locations,num_cells,Dab] = BS_user_placement(macro_side,cell_side,num_users,M,Nx,Ny);
b = repmat(2.^a,num_cells,1);
%***Small topology*****
% load('topology20.mat','z_h','z_u');
% BS_locations = z_h;
% user_locations = z_u;


%*****LArge topology******
BS_locations = BS_Locs;
user_locations = UE_Locs;
BS_power = SC_power*ones(size(BS_locations));
%BS_power(1) = MC_power;

[schedulable_rates_matrix, subset_cardinality,received_power,pathloss_SINR,assoc_mat] =...
    physical_layer(path_loss,num_cells,num_users,num_antennas,max_users_served,...
    BS_power,Noise_power,max_rate);


%R_max = max(max(max(schedulable_rates_matrix)));
%served_users_evolution = zeros(num_cells, num_users,N);     
rate_evolve = zeros(N*num_slots_chunk,num_users);
virtq_evolve = zeros(N*num_slots_chunk,num_users);
virtq2_evolve = zeros(N*num_slots_chunk,num_users);
size_evolve = zeros(N,num_users);
d_evolve = zeros(N,num_users);
subset_evolution = zeros(num_cells,N*num_slots_chunk);
cardinality_evolution = zeros(num_cells,N*num_slots_chunk);
simul_assoc = zeros(N*num_slots_chunk,num_users);
chunk_number_list = zeros(N,num_users);
chunk_number_requested = zeros(1,num_users);
req_time_stamp = zeros(1,num_users);
req_list = zeros(N,num_users);
del_param = 3;
for i = 1:N*num_slots_chunk
    
%     [schedulable_rates_matrix, subset_cardinality,received_power,pathloss_SINR,pathloss,BS_user_distance] =...
%     physical_layer(BS_locations, user_locations,num_cells,...
%     num_users,num_antennas,max_users_served);
%     
     [d_choice,rate_choice] = congestion_control(xx,rand_start,total_gop,i,virtq,virtq2,num_users,num_slots_chunk);
[served_user_rates,served_user_index] = WSRmax(schedulable_rates_matrix,virtq,max_users_served,num_users,...
    num_cells,max_rate,assoc_mat,i); 

served_users = sum(served_user_index,3);
%served_users_evolution(:,:,i) = served_users;
subset_evolution(:,i) = sum(served_users.*b,2);
 simul_assoc(i,:) = sum(served_users,1);
 if(advanced_receiver==1)
 offered = sum(served_user_rates,1);
 else
     offered = max(served_user_rates,[],1);
 end
[gamma] = gammaupdate(virtq2,V);
%actual_served = min(virtq, offered);
virtq = max(virtq-num_symbols_slot*offered + rate_choice,0);
actual_served = min(rate_choice + virtq, num_symbols_slot*offered);
virtq2 = max(virtq2 + gamma - d_choice,0); 

virtq_evolve(i,:) = virtq;
virtq2_evolve(i,:) = virtq2;
if(mod(i,num_slots_chunk)==1)
   chunk_index = floor(i/num_slots_chunk)+1;
   size_evolve(chunk_index,:) = rate_choice;
   chunk_number_requested = chunk_number_requested + 1;
   chunk_number_list(chunk_index,:) = chunk_number_requested;
   req_time_stamp = chunk_index;
   req_list(chunk_index,:)= req_time_stamp;
   d_evolve(chunk_index,:) = d_choice;
end


rate_evolve(i,:) = actual_served;
i

end

[allchunk_reception_profile,served_bits]...
     = reception_time(rate_evolve,size_evolve,num_users,N,num_slots_chunk);
alluser_chunk_delay = allchunk_reception_profile - req_list;
delay_profile = alluser_chunk_delay;

%% 
pb_buffers = zeros(N,num_users);
play_start_time = zeros(N,num_users);
buff_start_time = zeros(N,num_users);
max_delay_window = inf*ones(N,num_users);
start = zeros(1,num_users);
 playable_time = zeros(N,num_users);
 last_playable = zeros(1,num_users);
 arrival_profile = allchunk_reception_profile;
 new_chunks_avail = zeros(N,num_users);
 last_playable_chunk = zeros(N,num_users);
 buffer_status = zeros(N,num_users);
 for i= 2:N
     temporaryy = chunk_number_list.*((allchunk_reception_profile<=i)...
         & (chunk_number_list > repmat(last_playable,N,1))); 
      num_chunk_avail = sum(temporaryy~=0);                                  %number of new chunks which have arrived before ***i***
     new_chunks_avail(i,:) = num_chunk_avail;
     last_playable = last_playable + num_chunk_avail;
     last_playable_chunk(i,:) = last_playable;
      max_delay_window(i,:) = max(delay_profile.*((max(0,i-10) < allchunk_reception_profile)...
         & (allchunk_reception_profile <= i)),[],1);
     max_delay_window(max_delay_window==0)=inf;
     temp1 = pb_buffers(i-1,:);
     temp3 = max_delay_window(i,:);
     new_startup = ((pb_buffers(i-1,:) == 0) & (start==1)& (mod(chunk_number_list(i,:),101)~=0));
     play_start = ((temp1 >= del_param*temp3) &(start == 0));
     temp_start = play_start_time(i,:);
     temp_stop = buff_start_time(i,:);
     temp_start(play_start) = i;
     temp_stop(new_startup) = i;
     play_start_time(i,:) = temp_start;
     buff_start_time(i,:) = temp_stop;
     start(new_startup) = 0;
     start(play_start)=1;
     %start(stop) = 0;
     temp1((start==1)) = max(temp1((start==1))-1,0) + num_chunk_avail((start==1));
     temp1((start==0)) = temp1((start==0)) + num_chunk_avail((start==0));
     pb_buffers(i,:) = temp1;
     buffer_status(i,:) = start;
      %i
 end
 temp = cumsum((play_start_time~=0),1);
 tempor = (temp == 1);
 wanted = (cumsum(tempor) == 1);
 [row, col] = find(wanted);
 row  = row;
 col = col';
 startup_time = zeros(1,num_users);
 startup_time(col) = row';
 startup_time(startup_time == 0) = N;
%% 
figure(2)
hold all
cdfplot((sum((buffer_status==0),1)*100)/N)
xlabel('buffering time (x)')
ylabel('fraction of users with pre-buffering time < x')
%%

% num_rebuff = sum(buff_start_time>0);
% figure(3)
% cdfplot(num_rebuff)
% xlabel('number of re-buffering periods (x)')
% ylabel('fraction of users with no. of re-buffering periods < x')
% hold all
% 
% frac_rebuff_time = 1 - (sum(buffer_status)+startup_time)/N;
% figure(4)
% cdfplot(frac_rebuff_time*100)
% xlabel('percentage of playback time spent in re-buffering mode (x)')
% ylabel('fraction of users with % of re-buffering < x')
%hold all
%%
figure(3)
hold all
timeave_ssim = mean(d_evolve);
cdfplot(timeave_ssim)
xlabel('SSIM averaged over delivered chunks (x)')
ylabel('fraction of users with average SSIM < x')
%% 
figure(4)
hold all
num_received = sum(delay_profile ~= inf);
delay_profile(delay_profile == inf) = 0;
sum_delay = sum(delay_profile);
avg_delay = sum_delay./num_received;
cdfplot(avg_delay)
xlabel('average delay in reception of chunks (x)')
ylabel('fraction of users with average delay < x')
%%
figure(5)
hold all
yy = sum(virtq_evolve+virtq2_evolve,2);
zz = cumsum(yy)./cumsum(ones(size(yy)));
plot(zz)
grid

figure(6)
hold all
plot(sum(log(cumsum(d_evolve,1)./cumsum(ones(size(d_evolve)))),2))
grid



%jj = jj+1;
%hold all
%activity = activity/N;







