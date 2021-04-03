 % PV+ST considering curtailment and tariff model
function [P_s,P_b,P_C,P_D,SOC,Soc,obj] = agg_DA_PV_nocurtail(N_PV,P_PV_Ch,P_PV_Dis,....
    SOC_PV_initial,SOC_PV_min,SOC_PV_max,P_PV,eta_1,eta_2,P_load,alfa_1,alfa_2,T)
%% Define variables
SOC_PV = sdpvar(T+1,N_PV,'full');
P_Ch = sdpvar(T,N_PV,'full');
P_Dis = sdpvar(T,N_PV,'full');
% P_sell = sdpvar(T,N_PV,'full');
% P_buy = sdpvar(T,N_PV,'full');
%% big M variables
% z1 = sdpvar(T,N_PV,'full'); % delta_1*P_market
% z2 = sdpvar(T,N_PV,'full'); % delta_2*P_market
z3 = sdpvar(T,N_PV,'full'); % delta_1*delta_3*P_Ch
z4 = sdpvar(T,N_PV,'full'); % delta_1*delta_4*P_Dis
z5 = sdpvar(T,N_PV,'full'); % delta_2*delta_5*P_Ch
z6 = sdpvar(T,N_PV,'full'); % delta_2*delta_6*P_Dis
%% binary variable
delta_1= binvar(T,N_PV,'full'); % 
delta_2= binvar(T,N_PV,'full'); % 
delta_3= binvar(T,N_PV,'full'); % 
delta_4= binvar(T,N_PV,'full'); % 
delta_5= binvar(T,N_PV,'full'); % 
delta_6= binvar(T,N_PV,'full'); % 
%% Define constraints and objective function
% constraints of SOC
Constraints = [];
%% SOC of PV storage
Constraints = [Constraints, (SOC_PV(1,:) == SOC_PV_initial):['SOC_1']];
% Constraints = [Constraints, SOC_PV(end,:) == SOC_PV(1,:)];
for i=1:N_PV
    for t = 1:T
        Constraints = [Constraints, (SOC_PV(t+1,i) == SOC_PV(t,i)+(z3(t,i)+....
            z5(t,i))*eta_1-(z4(t,i)+z6(t,i))*eta_2):['SOC_PV' num2str(i+t)]];
    end
        Constraints = [Constraints, (SOC_PV_min(i)<=SOC_PV(:,i)<=SOC_PV_max(i)):['SOC_PVL' num2str(i+t)]];
end
Constraints = [Constraints, (0<=P_Ch<=P_PV_Ch):['P_ch']];
Constraints = [Constraints, (0<=P_Dis<=P_PV_Dis):['P_Dis']];
Constraints = [Constraints, (delta_1+delta_2==1):['detla1']]; % surplus energy or not
Constraints = [Constraints, (delta_3+delta_4<=1):['detla3']]; % surplus mode: charging or discharging
Constraints = [Constraints, (delta_5+delta_6<=1):['detla5']]; % shorage mode: charging or discharging
Constraints = [Constraints, (delta_1.*P_PV-delta_1.*P_load-z3+z4>=0):['P_equ1']];
Constraints = [Constraints, (delta_2.*P_PV-delta_2.*P_load-z5+z6<=0):['P_equ2']];
%% big M constraints for charging/discharging
M_3=P_PV_Ch;m_3=zeros(T,N_PV);
M_4=P_PV_Dis;m_4=zeros(T,N_PV);
M3=max(0,-m_3);m3=min(0,-M_3);
M4=max(0,-m_4);m4=min(0,-M_4);n=2;
% P_Ch, sell
Constraints = [Constraints,(z3+M3.*(delta_1+delta_3) <= P_Ch + n*M3):['P_M19']];
Constraints = [Constraints,(z3 <= P_Ch + M3):['P_M20']];
Constraints = [Constraints,(-z3-m3.*(delta_1+delta_3)<= (-P_Ch-n.*m3)):['P_M21']];
Constraints = [Constraints,(-z3 <= -P_Ch-m3):['P_M22']];
Constraints = [Constraints,(z3 +delta_1.*m3 <= 0):['P_M23']];
Constraints = [Constraints,(z3 +delta_3.*m3 <= 0):['P_M24']];
Constraints = [Constraints,(-z3 -delta_1.*M3 <= 0):['P_M23']];
Constraints = [Constraints,(-z3 -delta_3.*M3 <= 0):['P_M24']];
% P_Ch, buy
Constraints = [Constraints,(z5+M3.*(delta_2+delta_5) <= P_Ch + n*M3):['P_M25']];
Constraints = [Constraints,(z5 <= P_Ch + M3):['P_M26']];
Constraints = [Constraints,(-z5-m3.*(delta_2+delta_5) <= -P_Ch -n.*m3):['P_M27']];
Constraints = [Constraints,(-z5 <= -P_Ch - m3):['P_M28']];
Constraints = [Constraints,(z5 +delta_2.*m3 <= 0):['P_M29']];
Constraints = [Constraints,(z5 +delta_5.*m3 <= 0):['P_M30']];
Constraints = [Constraints,(-z5 -delta_2.*M3 <= 0):['P_M29']];
Constraints = [Constraints,(-z5 -delta_5.*M3 <= 0):['P_M30']];
% P_Dis, sell
Constraints = [Constraints,(z4+M4.*(delta_1+delta_4) <= P_Dis + n*M4):['P_M31']];
Constraints = [Constraints,(z4 <= P_Dis + M4):['P_M32']];
Constraints = [Constraints,(-z4-m4.*(delta_1+delta_4) <= -P_Dis -n*m4):['P_M33']];
Constraints = [Constraints,(-z4 <= -P_Dis - m4):['P_M34']];
Constraints = [Constraints,(z4 +delta_1.*m4<= 0):['P_M35']];
Constraints = [Constraints,(z4 +delta_4.*m4<= 0):['P_M36']];
Constraints = [Constraints,(-z4 -delta_1.*M4<= 0):['P_M35']];
Constraints = [Constraints,(-z4 -delta_4.*M4<= 0):['P_M36']];
% P_Dis, buy
Constraints = [Constraints,(z6+M4.*(delta_2+delta_6) <= P_Dis + n*M4):['P_M37']];
Constraints = [Constraints,(z6 <= P_Dis + M4):['P_M38']];
Constraints = [Constraints,(-z6-m4.*(delta_2+delta_6) <= -P_Dis -n*m4):['P_M39']];
Constraints = [Constraints,(-z6 <= -P_Dis - m4):['P_M40']];
Constraints = [Constraints,(z6 +delta_2.*m4 <= 0):['P_M41']];
Constraints = [Constraints,(z6 +delta_6.*m4 <= 0):['P_M42']];
Constraints = [Constraints,(-z6 -delta_2.*M4 <= 0):['P_M41']];
Constraints = [Constraints,(-z6 -delta_6.*M4 <= 0):['P_M42']];
%% objective function
% Objective = -sum(sum(z1-delta_1.*P_load-z3+z4,2).*alfa_2+sum(z2-delta_2.*P_load-z5+z6,2).*alfa_1);
Objective = -sum(sum(delta_1.*P_PV-delta_1.*P_load-z3+z4,2).*alfa_2+....
    sum(delta_2.*P_PV-delta_2.*P_load-z5+z6,2).*alfa_1)+sum(sum(power((delta_1.*P_PV-....
    delta_1.*P_load-z3+z4+delta_2.*P_PV-delta_2.*P_load-z5+z6),2)));
%% Solve the problem
options = sdpsettings('solver','gurobi');
% options = sdpsettings('solver','cplex');
% options = sdpsettings('solver', 'lpsolve');
sol = optimize(Constraints,Objective,options);

% Analyze error flags
if sol.problem == 0
    % Extract and display value
% z1=value(z1);
% z2=value(z2);
z3=value(z3);
z4=value(z4);
z5=value(z5);
z6=value(z6);
Soc=value(SOC_PV); % whole SOC for each hour
SOC=Soc(end,:); % the end of SOC
delta_1=value(delta_1);
delta_2=value(delta_2);
delta_3=value(delta_3);
delta_4=value(delta_4);
delta_5=value(delta_5);
delta_6=value(delta_6);
P_s=delta_1.*P_PV-delta_1.*P_load-z3+z4; % sold energy
P_b=delta_2.*P_PV-delta_2.*P_load-z5+z6; % bought energy
obj=value(Objective);
P_C=z3+z5; % charging power for each hour
P_D=z4+z6; % discharing power for each hour
P_Ch=value(P_Ch);
P_Dis=value(P_Dis);
else
    display('Hmm, something went wrong!');
    sol.info
    yalmiperror(sol.problem)
end
obj_agg=value(Objective);
end