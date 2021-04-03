%% PVST with tariff model, 2015-2016, 2016-2017, 2017-2018
tic
clc
clear all
%% load price and PV generation ,and load
load DApricethreeyear
load Bornholm
load DAprice_2018
load PV_2018
load load_withHP
load load_withoutHP
% load tenwithandwithoutHPs
%% system setting
T=24; % optimization horizon
DAoneyear=DA2015; % Day ahead price, 1st column DAM, 2rd DSO tariff, 3rd TSO tariff
N_PV=1;  % numbers of PV
% E_PV=[4.5 4.5 4.5 8.5 12.5 4.5 12.5 4.5 8.5 8.5]; % the size of PV for each resident
E_PV = 39.44; % case : PV size
battery_size=0.6.*E_PV(1:N_PV);  % unit in kWh, size of the battery or 1.2 1.5
SOC_PV_initial=0.2*ones(1,N_PV).*battery_size; 
SOC_PV_min = 0.2.*battery_size;
SOC_PV_max = 0.85.*battery_size;
%% converter power and efficies
P_converter_ch=0.4*battery_size;  % converter charging power
P_converter_dis=0.4*battery_size;
% predefined variable to increase speed
P_PV_Ch=ones(T,1)*P_converter_ch;
P_PV_Dis=ones(T,1)*P_converter_dis;
% efficiency of charging and discharging
eta_1=0.9*ones(1,N_PV);
eta_2=1./(0.95*ones(1,N_PV));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This part is used for optimal schedule generation
Cost_opt=zeros(length(DA2015)/T,N_PV);
Cost_normal=zeros(length(DA2015)/T,N_PV);
for k=1:length(DA2015)/T
    k
%     Spot_hourly=DAoneyear((k-1)*T+1:k*T,1)./1e3; 
%     P_PV=PVgeneration((k-1)*T+1:k*T,1:N_PV); % PV output
%     P_load=tenwithHP((k-1)*T+1:k*T,1:N_PV); % tenwithHP: load profile with HP;tenwithoutHP: load profile without HP; 
%     Tf_DSO=DAoneyear((k-1)*T+1:k*T,2); % tariff DSO
%     Tf_TSO=DAoneyear((k-1)*T+1:k*T,3);
    Spot_hourly=price_2018((k-1)*T+1:k*T,1)./1e3; 
    P_PV=pv_2018((k-1)*T+1:k*T,1:N_PV); % PV output
    P_load=load_withhp((k-1)*T+1:k*T,1:N_PV);
%     P_load=load_withouthp((k-1)*T+1:k*T,1:N_PV);
    %% used for optimization from prosumer point of view
%     alfa_1=(0.2*ones(24,1)+Tf_TSO+Tf_DSO);   % price for buying, 0.2 retail energy price
%     alfa_2=0.32; % price for selling
    alfa_1=Spot_hourly;
    alfa_2=Spot_hourly;
    %% input SOC of each day
    if k==1
    SOC_int(k,:)=SOC_PV_initial;
    else
    SOC_int(k,:)=SOC_in(k-1,:);
    end
    %% main function, opt for 24 hours
[P_s{k},P_b{k},P_C{k},P_D{k},SOC_in(k,:),SOC{k},obj(k)] = agg_DA_PV_nocurtail(N_PV,P_PV_Ch,....
    P_PV_Dis,SOC_int(k,:),SOC_PV_min,SOC_PV_max,P_PV,eta_1,eta_2,P_load,alfa_1,alfa_2,T);
    %% results
Cost_opt(k,:)=sum(P_s{k}.*alfa_2+P_b{k}.*alfa_1); % revenue of DSO when prosumer has PVST 
Cost_normal(k,:)=sum(P_load.*alfa_1); % without PVST
hour_PLF((k-1)*24+1:k*24,:)=(P_s{k}+P_b{k})./P_load; % hourly peak ratio
daily_PLF(k,:)=sum(P_s{k}+P_b{k})./sum(P_load); % averay daily peak ratio
sum_PV(k,:)=sum(P_s{k}+P_b{k}); % net energy exchange
sum_load(k,:)=sum(P_load);
annal_PLF=sum(sum_PV)./sum(sum_load); % average annual peak load ratio
% figure;bar(SOC{k})
% figure;bar(P_C{k});hold on;bar(-P_D{k});hold on
end
% kk=zeros(365,24);
% for k=1:365
%     Tf_DSO=DAoneyear((k-1)*T+1:k*T,2); % tariff DSO
%     Tf_TSO=DAoneyear((k-1)*T+1:k*T,3);
%     kkk1(k,:)=sum(tenwithHP((k-1)*T+1:k*T,1:N_PV).*(Tf_TSO+Tf_DSO));
%     kkk2(k,:)=sum(P_b{k}.*(Tf_TSO+Tf_DSO));
% end
%% results and plotting
%% cost of prosumer
C_1Y_o=-sum(Cost_opt); % optimized cost
C_1Y_n=sum(Cost_normal); % without opt cost
%% peak load factor
% sum(P_s{})
%% summer case
% kk=166;
% for k=166:172
% %% winter case
% %   kk=291;
% % for k=291:297
%     P_PV_week((k-kk)*24+1:(k-kk+1)*24,1)=PVgeneration((k-1)*24+1:k*24,1:1);
%     P_load_week((k-kk)*24+1:(k-kk+1)*24,1)=tenwithHP((k-1)*24+1:k*24,1:1);
%     P_s_week((k-kk)*24+1:(k-kk+1)*24,1)=P_s{k}(:,1);
%     P_b_week((k-kk)*24+1:(k-kk+1)*24,1)=P_b{k}(:,1);
%     P_Ch_week((k-kk)*24+1:(k-kk+1)*24,1)=P_C{k}(:,1);
%     P_Dis_week((k-kk)*24+1:(k-kk+1)*24,1)=P_D{k}(:,1);
% end
%% plotting
% figure;
% P_p=[P_PV_week P_Dis_week];bar(P_p,'stacked');hold on
% P_n=[-P_load_week -P_Ch_week];bar(P_n,'stacked');hold on
% Balance=P_PV_week+P_Dis_week-P_load_week-P_Ch_week;
% plot(Balance);hold on