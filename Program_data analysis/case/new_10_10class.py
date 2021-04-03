# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 21:24:57 2019

@author: NIlab
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 22:27:09 2019

@author: xuxue
"""

### special course的方法 应用到case


import scipy.io as sio
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import mpl_toolkits.axisartist as axisartist
from mpl_toolkits.mplot3d import Axes3D #画三维图不可少
from matplotlib import cm  #cm 是colormap的简写

### 导入数据 ###

pv_raw = pd.read_excel('pv_2018_filled.xlsx')
load_raw = pd.read_excel('2018_load.xlsx')
price_raw = pd.read_excel('spot_price_2018.xlsx')
tariff = sio.loadmat('DApricethreeyear.mat')
PC_mat = sio.loadmat('PC_new')
PD_mat = sio.loadmat('PD_new')

PV_df = pd.DataFrame(np.array(pv_raw.Value_perhour).reshape(365,24))

load_df = pd.DataFrame(np.array(load_raw["sum of those w/o hp."]).reshape(365,24))
#load_df = pd.DataFrame(np.array(load_raw["sum of those w hp."]).reshape(365,24))

#price_perhour = pd.DataFrame(np.array(price_raw.iloc[2:,9].dropna()))
#price_perhour = list(price_raw.iloc[2:,9].dropna())
p = price_raw.iloc[2:,9].dropna()

price = []
for i in range(2,len(p)+3):
    if i==1996:
        continue
    num = len(p[i])
    if p[i][0] == '-':
        number = -(int(p[i][1:num-3])+int(p[i][num-2:])*0.01)
    else:
        number = int(p[i][0:num-3])+int(p[i][num-2:])*0.01
    price.append(number)


#price_perhour = pd.DataFrame(np.array(price).reshape(365,24))
spot_price = np.array(price)
price_perhour_spot = spot_price/1000  # euro/mwh to euro/kwh

tariff_raw = tariff["DA2017"] #price包含三列，第一列 spot price 第二列 DSO tariff 第三列 TSO tariff
VAT = 0.25
#Spot_price = tariff_raw[:,0]
DSO_tariff = tariff_raw[:,1]
TSO_tariff = tariff_raw[:,2]
electricity_tax = 0.774
price_perhour_retail = []
for i in range(8760):
    price_perhour_retail.append((1+VAT)*(price_perhour_spot[i] + DSO_tariff[i] + TSO_tariff[i] + electricity_tax))
#    price_perhour_retail.append((1+VAT)*(DSO_tariff[i] + TSO_tariff[i] + electricity_tax))


#prcie_perhour = price_perhour_spot  # 用spot price的情况
price_perhour = np.array(price_perhour_retail) # 用retail price的情况
 
#price_sell = price_perhour_spot
price_sell = np.array([0.32 for i in range(8760)])


# 输出部分处理
PC_raw = PC_mat["P_C"]
PD_raw = PD_mat["P_D"]

# 初始化PC和PD的dataframe 先取第1天数据填入
PC_df = pd.DataFrame(PC_raw[0,0]) #默认是生成一列
PC_df = PC_df.T #转置成行
PD_df = pd.DataFrame(PD_raw[0,0])
PD_df = PD_df.T

# 把剩下的数据依次往dataframe上拼接        
for i in range(1,365):
    PC_temp = pd.DataFrame(PC_raw[0,i])
    PC_temp = PC_temp.T
    PC_df = pd.concat([PC_df,PC_temp])
    PD_temp = pd.DataFrame(PD_raw[0,i])
    PD_temp = PD_temp.T
    PD_df = pd.concat([PD_df,PD_temp])
    
    
    
#PV求和
PVinput_day = []
for i in range(365):
    PVinput_day.append(PV_df.iloc[i,:].sum())
PVinput_day = pd.DataFrame(PVinput_day)
#load求和
Loadinput_day = []
for i in range(365):
    Loadinput_day.append(load_df.iloc[i,:].sum())
Loadinput_day = pd.DataFrame(Loadinput_day)
    
### 应用Kmeans进行分类 ###
#   N：PV分类数    M：load分类数

error_rate_list = [[] for i in range(10)]#根据PV分为10个list,每个list里对应十类load
N_total = 10
M_total = 10

for N in range(1,N_total+1):
    for M in range(1,M_total+1):
        #PV
        estimator_PV = KMeans(n_clusters=N) #构造聚类器
        estimator_PV.fit(PVinput_day) #聚类
        Idx_PV = estimator_PV.labels_  #获取聚类标签
        
        #Load
        estimator_Load = KMeans(n_clusters=M) #构造聚类器
        estimator_Load.fit(Loadinput_day) #聚类
        Idx_load = estimator_Load.labels_  #获取聚类标签
            
        #classification for all data by labels of PV&load
        # classify_xu函数用于将标签分成 N*M 类
        def classify_xu(Idx_PV,Idx_load,N,M):
            label = [None for ii in range(365)]
            for i in range(N):
                for j in range(M):
                    for t in range(365):
                        if Idx_PV[t] == i and Idx_load[t] == j:
                            label[t] = i*M+j
            return label
        
        # label包含所有样本的标签
        label = classify_xu(Idx_PV,Idx_load,N,M);
        
        #A is used to save all types of data (per day) 用于存分过类之后的所有样本
        # y_opt 是365*24的CD值(from optimization)
        y_opt = [[None for ii in range(24)]for jj in range(365)]
        for i in range(365):
            for j in range(24):
                y_opt[i][j] = abs(PC_df.iloc[i,j])-abs(PD_df.iloc[i][j])
                
        A = [[None]for i in range(N*M)]
        for i in range(365):
            A[label[i]].append(y_opt[i])
        for i in range(N*M):
            A[i].pop(0)
        #A = pd.DataFrame(A)
        
        #求每一类的24个小时的每小时平均值
        # CD_classified 分过类之后 N*M 类 profile 
        CD_classified = [[None for i in range(24)]for j in range(N*M)]
        for i in range(N*M):
            if A[i]==[]:
                continue
            for hour in range(24):
                sum_hour = 0
                for j in range(len(A[i])):
                    sum_hour += A[i][j][hour]
                CD_classified[i][hour] = sum_hour/len(A[i])
        CD_classified_df = pd.DataFrame(CD_classified)
        
        # 将分过类的profile按照label再应用到一年当中，即为预测结果y_est
        y_est_raw = [[None for ii in range(24)]for jj in range(365)]
        for i in range(365):
        #    for j in range(24):
                y_est_raw[i] = CD_classified[label[i]]
        
        y_opt_nd = np.array(y_opt)
        y_opt_fl = y_opt_nd.flatten()
        
        y_est_nd = np.array(y_est_raw)
        y_est_test = y_est_nd.flatten()
        y_est = y_est_nd.flatten()
        
        ### SOC limitation ###
        battery_size = 15 #kwh
        SOC_min = 0.2*battery_size
        SOC_max = 0.85*battery_size
        battery_cap = SOC_max - SOC_min
        
        # 找到y_est中第一个大于0的即充电的
        loc = 0
        for i in range(len(y_est)):
            if y_est[i] > 0:
                charge = y_est[i]
                if charge > battery_cap:
                    charge = battery_cap
                    y_est[i] = battery_cap #保证第一个充电的值不大于battery capacity
                loc += 1
                break
            if y_est[i] < 0:
                charge = 0
                y_est[i] = 0
                       
        for j in range(loc,len(y_est)):
            if charge + y_est[j] >= 0:
                former = charge
                charge = charge+y_est[j]
                if charge > battery_cap:
                    y_est[j] = battery_cap-former
                    charge = battery_cap    
            else:
                y_est[j] = -charge
                charge = 0
               
        ### calculate cost by different methods ###
        load_fl = np.array(load_df).flatten()            
        PV_fl = np.array(PV_df).flatten()  
        #price_perhour = np.array(price_perhour)
        final_opt = -PV_fl + load_fl + y_opt_fl    
        final_est = -PV_fl + load_fl + y_est     
        
        cost_opt = []
        cost_est = []   
        cost_load = []
        
        for i in range(len(final_opt)) :  
            if final_opt[i]>=0:  # >0则买电
                cost_opt.append(final_opt[i]*price_perhour[i])
            else:
                cost_opt.append(final_opt[i]*price_sell[i])
                
        for i in range(len(final_est)) :  
            if final_est[i]>=0:  # >0则买电
                cost_est.append(final_est[i]*price_perhour[i])
            else:
                cost_est.append(final_est[i]*price_sell[i])
                
        for i in range(len(load_fl)) :  
            if load_fl[i]>=0:  # >0则买电
                cost_load.append(load_fl[i]*price_perhour[i])
            else:
                cost_load.append(load_fl[i]*price_sell[i])        
                
        cost_opt_sum = sum(cost_opt)
        cost_est_sum = sum(cost_est)
        cost_load_sum = sum(cost_load)
        error_rate = ((cost_est_sum - cost_opt_sum)/cost_opt_sum)*100
        error_rate_list[N-1].append(error_rate)



#plot
#x，y坐标矩阵
X = [[j for i in range(10)] for j in range(1,11)]
Y = [[i for i in range(1,11)] for j in range(10)]
X = np.array(X)
Y = np.array(Y)

error_rate_list_np = np.array(error_rate_list)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X,Y,error_rate_list_np,cmap='rainbow')

ax.invert_yaxis()  #y轴反向
#x轴间距1
x_major_locator=MultipleLocator(1)
#y轴间距1
y_major_locator=MultipleLocator(1)
#设置间距
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
plt.title('Different combinations of N and M',fontsize=20)
ax.set_zlabel('Error Rate (%)',fontsize=20)  # 坐标轴
ax.set_ylabel('M (load)',fontsize=20)
ax.set_xlabel('N (PV)',fontsize=20)

er_df = pd.DataFrame(error_rate_list_np)
er_df.to_excel("er100.xlsx")


