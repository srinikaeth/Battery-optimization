
%% Model inputs pre-processing

% S_base = 1 MW

% Getting the distribution network model
mpc = Node123;

% Getting the Ybus from the network model
[Y, YF, YT] = makeYbus(mpc);

% X and R for the voltage model
X = csvread('X_2.csv');
R = csvread('R_2.csv');

% Extracting the buy and sell prices
p_buy = csvread('P_buy_TOU2.csv');
p_sell = csvread('P_sell_TOU2.csv');

% Slack bus
slack = 1;

%%% Adding the residential loads

% Reading the residential loads
unalt_loads = csvread('2300_homes.csv');

% Getting the PV generation data
pv_gen = csvread('PV_energy.csv');
pv_gen = pv_gen';

% Input for each case
T = 72;
n_buses = 123;

home_per_bus = 18;
n_homes = home_per_bus*n_buses;

% Power factor generated randomly
pf = 0.7 + 0.25*rand(n_buses,T);
Q_mult = tan(acos(pf)); 
 
% Theta limits
thetaL = -ones(n_buses,T);
thetaU = ones(n_buses,T);

% Restricting input data to no of homes required
unalt_loads = unalt_loads(1:n_homes,:);

% Restricting load data to no of hours required
unalt_loads = unalt_loads(:, 1:T);

% Extracting only the required PV data
pv_gen = pv_gen(1,1:T);

% Extracting only the required buy and sell prices
p_buy = p_buy(1:T,1);
p_sell = p_sell(1:T,1);

% Converting prices to /MWh
p_buy = 1000*p_buy;
p_sell = 1000*p_sell;

% Getting transformer limits
TP = csvread('Trans_3.csv');

Cost_plot = zeros(89,1);

%% Loop for different penetration limits

for pen = 0:10:110

storage_bus = pen;
n_storage = home_per_bus*storage_bus;

full_loads = unalt_loads;

% Subtracting the PV generation from only homes with storage
for v = 1:n_storage
    full_loads(v,:) = full_loads(v,:) - pv_gen;
end

% Full loads is in kW

% Getting the OPF input in MW
%L_opf = 0.001*full_loads;

% Merging all the loads in the same bus as one column
for m = 1:123
    k = home_per_bus*(m-1) + 1; 
    total_loads(m,:) = sum(full_loads(k:k+home_per_bus-1,:),1);
end

% Getting the OPF input in MW
L_buses = 0.001*total_loads;

% Solving the OPF
[ dcopf_soln ] = dcopf_storage_network_v19(Y,T,p_buy,p_sell,L_buses,slack,thetaL,thetaU,n_storage, ...
                                            pen,n_buses,n_homes,storage_bus,X,R,Q_mult,TP);

% Storing the cost in a separate variable for plotting                                        
Cost_plot(pen+1,1) = dcopf_soln.Cost;
% B_plot = dcopf_soln.b;

end                                        
