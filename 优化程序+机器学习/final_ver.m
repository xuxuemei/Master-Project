% Apply Data mining for PVST with variable price
%% load in data 
clear all;
close all;
load Bornholm.mat;
load DApricethreeyear.mat;
load PVSTwithvariableprice;

%%% training model %%%
%% rearrange input and output
% output: cell(365*24) -> 8760*1
PCharge = reshape(cell2mat(P_C),8760,1);
PDischarge = reshape(cell2mat(P_D),8760,1);
y = abs(PCharge) - abs(PDischarge);  % output(continuos)
% input: 1st resident as model
PVinput = PVgeneration(:,1); 
Loadinput = tenwithHP(:,1);

%% retail price
Spot_hourly=DA2017(:,1);
Tf_DSO=DA2017(:,2); % tariff DSO
Tf_TSO=DA2017(:,3);
Tf_subscription=0;
RandS=0; % renewable energy support
electricity_tax=0.774;
VAT=0.25;
% price for buying, 0.2 retail energy price
price=(1+VAT)*(Spot_hourly+Tf_TSO+Tf_DSO+Tf_subscription+RandS+electricity_tax);   

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
% er = zeros(10,10);
N = 4; % PV
M = 4; % load

% Kmeans
[Idx_PV,C_PV, sumd_PV, D_PV] = kmeans(PVinput_day_sum,N,'Start','cluster');
[Idx_load, C_load, sumd_load, D_load] = kmeans(Loadinput_day_sum,M,'Start','cluster');

%% plot clustered data
col=colormap(jet(N));
for i=1:N
[a{i},b{i}]=find(Idx_PV==i);
PV{i}=PVinput_day_sum(a{i});
plot(PV{i},'o','Color',col(i,:));hold on
xlabel('Objectives for each category');
ylabel('PV generation (per day)');
end
col1=colormap(jet(M));
figure;
for i=1:M
[a1{i},b1{i}]=find(Idx_load==i);
Load{i}=Loadinput_day_sum(a1{i});
plot(Load{i},'o','Color',col1(i,:));hold on
xlabel('Objectives for each category');
ylabel('Load (per day)');
end

%% classification for all data by labels of PV&load
label = classify_xu(Idx_PV,Idx_load,N,M);

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
%     figure(i);   
%     for ii=1 : row_i 
%          plot(time, A{i}(ii, :));
%         hold on;
%     end    
% end

%% get mean value for classified y
% These are strategy patterns 
y_mean = zeros(col_A,24); 
for i = 1:col_A
    y_mean(i,:) = mean(A{i},1); 
end
for i = 1:col_A
    y_mean(i,:) = mean(A{i},1); 
end
y_all_mean = mean(y_day,1);
nan = NaN(1,24);
for i = 1:col_A
    if isnan(y_mean(i,:))==1
    y_mean(i,:) = y_all_mean(1,:);
    end
end


%%% test %%%
%% use the other residents as test objectives
PVinput_test = PVgeneration(:,10); 
Loadinput_test = tenwithHP(:,10);
% data per day --> matrix (365*24)
for i_day = 1:365
    PVinput_day_test(i_day,:) = PVinput_test(i_day*24-23:i_day*24);
    Loadinput_day_test(i_day,:) = Loadinput_test(i_day*24-23:i_day*24);
end

PVinput_day_sum_test = sum(PVinput_day_test,2);
Loadinput_day_sum_test = sum(Loadinput_day_test,2);

Idx_PV_test = zeros(length(PVinput_day_sum_test),1);
Idx_load_test = zeros(length(Loadinput_day_sum_test),1);

tmp_N = zeros(N,1);
tmp_M = zeros(M,1);
for i_day = 1:365
    for i = 1:N
        tmp_N(i,1) = abs(PVinput_day_sum_test(i_day,1)-C_PV(i,1));
    end
    [C,I] = min(tmp_N);
    Idx_PV_test(i_day,1)=I;
end
for i_day = 1:365
    for i = 1:M
        tmp_M(i,1) = abs(Loadinput_day_sum_test(i_day,1)-C_load(i,1));
    end
    [C,I] = min(tmp_M);
    Idx_load_test(i_day,1)=I;
end

label_test = classify_xu(Idx_PV_test,Idx_load_test,N,M);

for i = 1:365
    y_day_new_test(i,:) = y_mean(label_test(i),:);
end
y_new_test = y_day_new_test(1,:);
for i = 2:365
    y_new_test = [y_new_test y_day_new_test(i,:)];
end
y_new_test = y_new_test';


%% SOC
y_new = y_new_test;
for i = 1:size(y_new,1)%找到y_est中第一个大于0的即充电的
    if y_new(i)>0
        charge = y_new(i);
        if charge>1.755
            charge = 1.755;
            y_new(i) = 1.755;%保证第一个充电的值不大于2.295-0.54
        end
        loc = i+1;
        break
    else
        charge = 0;
        y_new(i) = 0;
    end
end
for j = loc:size(y_new,1)
    if charge+y_new(j)>=0
        former = charge;
        charge = charge+y_new(j);
        if charge>1.755
            y_new(j) = 1.755-former;
            charge = 1.755;
        end
    else
        y_new(j) = -charge;
        charge = 0;
    end
end

%% calculate cost by different methods
T=24; % optimization horizon
price_buy = price;
price_sell = 0.32;

final_est = -PVinput_test + Loadinput_test + y_new;
final_opt = -PVinput_test + Loadinput_test + y;
cost_est = zeros(length(DA2017),1);
cost_opt = zeros(length(DA2017),1);
cost_load = zeros(length(DA2017),1);

for i = 1:length(DA2017)
    % estimation cost
    if(final_est(i)<=0)
        cost_est(i) = final_est(i)*price_sell;
    else
        cost_est(i) = final_est(i)*price_buy(i);
    end
    % optimization cost
    if (final_opt(i)<=0)
        cost_opt(i) = final_opt(i)*price_sell;
    else
        cost_opt(i) = final_opt(i)*price_buy(i);
    end
    % only load cost
    if (Loadinput_test(i)<=0)
        cost_load(i) = Loadinput_test(i)*price_sell;
    else
        cost_load(i) = Loadinput_test(i)*price_buy(i);
    end
end         

final_cost_est = sum(cost_est);
final_cost_opt = sum(cost_opt);
% error rate
er_tmp = (final_cost_est-final_cost_opt)/final_cost_opt;
final_cost_load = sum(cost_load);

