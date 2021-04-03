% Apply Data mining for PVST with variable price, 20181202

%% load in data 
clear all;
close all;
load Bornholm.mat;
load DApricethreeyear.mat;
load PVSTwithvariableprice;
% load 2016PC;
% load 2016PD;
% load 2015PC;
% load 2015PD;

%% rearrange input and output
% output: cell(365*24) -> 8760*1
PCharge = reshape(cell2mat(P_C),8760,1);
PDischarge = reshape(cell2mat(P_D),8760,1);
y = abs(PCharge) - abs(PDischarge);  % output(continuos)
% input:
PVinput = PVgeneration(:,1);
Loadinput = tenwithHP(:,1);
%% retail price 
Spot_hourly=DA2017(:,1);
Tf_DSO=DA2017(:,2); % tariff DSO
Tf_TSO=DA2017(:,3);
% Spot_hourly=DA2016(:,1);
% Tf_DSO=DA2016(:,2); % tariff DSO
% Tf_TSO=DA2016(:,3);
% Spot_hourly=DA2015(:,1);
% Tf_DSO=DA2015(:,2); % tariff DSO
% Tf_TSO=DA2015(:,3);
Tf_subscription=0;
RandS=0; % renewable energy support
electricity_tax=0.774;
VAT=0.25;
price=(1+VAT)*(Spot_hourly+Tf_TSO+Tf_DSO+Tf_subscription+RandS+electricity_tax);   % price for buying, 0.2 retail energy price
% price = DA2017(:,2)+DA2017(:,1);  % only DSO tariff
% price = DA2017(:,1);  % spot price

% data per day --> matrix (365*24)
for i_day = 1:365
    PVinput_day(i_day,:) = PVinput(i_day*24-23:i_day*24);
    Loadinput_day(i_day,:) = Loadinput(i_day*24-23:i_day*24);
    price_day(i_day,:) = price(i_day*24-23:i_day*24);    
    y_day(i_day,:) = y(i_day*24-23:i_day*24);
end

PVinput_day_sum = sum(PVinput_day,2);
Loadinput_day_sum = sum(Loadinput_day,2);

%% Kmeans to classify by PV and Load
% number of clusters
N = 4; % PV
M = 4; % load

% matlab内置的Kmeans
% [Idx_PV,C_PV, sumd_PV, D_PV] = kmeans(PVinput_day_sum,N);
% [Idx_load, C_load, sumd_load, D_load] = kmeans(Loadinput_day_sum,M);
[Idx_PV,C_PV, sumd_PV, D_PV] = kmeans(PVinput_day_sum,N,'Start','cluster');
[Idx_load, C_load, sumd_load, D_load] = kmeans(Loadinput_day_sum,M,'Start','cluster');

% 自己写的Kmeans
% [u_PV, Idx_PV] = KMeans_zhao(PVinput_day_sum,N);
% [u_load, Idx_load] = KMeans_zhao(Loadinput_day_sum,M);

%% plot clustered data
col=colormap(jet(N));
for i=1:N
[a{i},b{i}]=find(Idx_PV==i);
PV{i}=PVinput_day_sum(a{i});
plot(PV{i},'o','Color',col(i,:));hold on
end
col1=colormap(jet(M));
figure;
for i=1:M
[a1{i},b1{i}]=find(Idx_load==i);
Load{i}=Loadinput_day_sum(a1{i});
plot(Load{i},'o','Color',col1(i,:));hold on
end
%% classification for all data by labels of PV&load
label = classify_zhao(Idx_PV,Idx_load,N,M);

% A is used to save all types of data (per day) in a cell
A = cell(1,N*M);

for i = 1:N*M
    loc_label = find(label==i);
    A{i} = y_day(loc_label,:);
end

% time = 1:24;
[row_A, col_A] = size(A);

% for i = 1:col_A
    [row_i, col_i] = size(A{i});  
%      figure(i);
    
%     for ii=1 : row_i 
%          plot(time, A{i}(ii, :));
%         hold on;
%     end
    
% end

%% get mean value for classified y
y_mean = zeros(col_A,24); 
for i = 1:col_A
    y_mean(i,:) = mean(A{i},1); 
end

%% apply new control strategy for everyday
for i = 1:365
    y_day_new(i,:) = y_mean(label(i),:);
end
y_new = y_day_new(1,:);
for i = 2:365
    y_new = [y_new y_day_new(i,:)];
end
y_new = y_new';

%% battery setting

N_PV=1;  % numbers of PV
E_PV=[4.5 4.5 4.5 8.5 12.5 4.5 12.5 4.5 8.5 8.5];
battery_size=0.6.*E_PV(1:N_PV);  % unit in kWh
SOC_PV_initial=0.2*ones(1,N_PV).*battery_size;
SOC_PV_min = 0.2.*battery_size;
SOC_PV_max = 0.85.*battery_size;
eta_1=0.9*ones(1,1); % efficiency of charging
eta_2=1./(0.95*ones(1,1)); % efficiency of discharging

% consider efficiency
y_new_eff = zeros(length(y_new),1);
y_eff = zeros(length(y),1);
for i = 1:size(y_new,1)
    if y_new(i)>0
        y_new_eff(i) = y_new(i)*eta_1;
    elseif y_new(i)<0
        y_new_eff(i) = y_new(i)*eta_2;
    end
end

for i = 1:size(y,1)
    if y(i)>0
        y_eff(i) = y(i)*eta_1;
    elseif y(i)<0
        y_eff(i) = y(i)*eta_2;
    end
end

% regulated by battery size
for i = 1:size(y_new_eff,1) %找到y中第一个大于0的即充电的
    if y_new_eff(i)>0
        charge = y_new_eff(i);
        if charge>(SOC_PV_max-SOC_PV_min)
            charge = SOC_PV_max-SOC_PV_min;
            y_new_eff(i) = SOC_PV_max-SOC_PV_min; %保证第一个充电的值不大于SOC_PV_max-SOC_PV_min
        end
        loc = i+1;
        break
    else
        charge = 0;
        y_new_eff(i) = 0;
    end
end
for j = loc:size(y_new_eff,1)
    if charge+y_new_eff(j)>=0
        former = charge;
        charge = charge+y_new_eff(j);
        if charge>(SOC_PV_max-SOC_PV_min)
            y_new_eff(j) = (SOC_PV_max-SOC_PV_min)-former;
            charge = (SOC_PV_max-SOC_PV_min);
        end
    else
        y_new_eff(j) = -charge;
        charge = 0;
    end
end

%% check SOC
SOC = nan(size(y_new_eff,1),1);
sum_SOC = SOC_PV_initial;
for k = 1:size(y_new_eff,1)
    sum_SOC = sum_SOC + y_new_eff(k);
    SOC(k,1) = sum_SOC;
end
figure;
plot(SOC);

%% calculate retail price
T=24; % optimization horizon
DAoneyear=DA2017; % Day ahead price, 1st column DAM, 2rd DSO tariff, 3rd TSO tariff
price_buy = zeros(length(DA2015),1);
price_sell = zeros(length(DA2015),1);
for k=1:length(DA2015)/T
    Spot_hourly=DAoneyear((k-1)*T+1:k*T,1)./1e3;
    Tf_DSO=DAoneyear((k-1)*T+1:k*T,2); % tariff DSO
    Tf_TSO=DAoneyear((k-1)*T+1:k*T,3);
    % used for optimization from prosumer point of view
    alfa_1=(1+VAT)*(Spot_hourly+Tf_TSO+Tf_DSO+Tf_subscription+RandS+electricity_tax); 
    price_buy((k-1)*T+1:k*T) = alfa_1; 
    alfa_2=0.32; % price for selling
    price_sell((k-1)*T+1:k*T) = alfa_2*ones(24,1); 
end

%% calculate the cost reduction
final_est = -PVinput + Loadinput + y_new;
final_opt = -PVinput + Loadinput + y;
cost_est = zeros(length(DA2015),1);
cost_opt = zeros(length(DA2015),1);
cost_load = zeros(length(DA2015),1);

for i = 1:length(DA2015)
    % est cost
    if(final_est(i)<=0)
        cost_est(i) = final_est(i)*price_sell(i);
    else
        cost_est(i) = final_est(i)*price_buy(i);
    end
    % optimized cost
    if (final_opt(i)<=0)
        cost_opt(i) = final_opt(i)*price_sell(i);
    else
        cost_opt(i) = final_opt(i)*price_buy(i);
    end
    % only load cost
    if (Loadinput(i)<=0)
        cost_load(i) = Loadinput(i)*price_sell(i);
    else
        cost_load(i) = Loadinput(i)*price_buy(i);
    end
end

final_cost_est = sum(cost_est);
final_cost_opt = sum(cost_opt);
(final_cost_est-final_cost_opt)/final_cost_opt
final_cost_load = sum(cost_load);