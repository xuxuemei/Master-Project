import scipy.io as sio
import numpy as np
import pandas as pd

# 导入数据
pv_and_load = sio.loadmat('Bornholm.mat')
price = sio.loadmat('DApricethreeyear.mat')
output = sio.loadmat('PVSTwithvariableprice')

res_case = 0 #取第几户家庭的数据，共有十户，可取0-9

PV_raw = pv_and_load["PVgeneration"]
PV_df = pd.DataFrame(PV_raw[:,res_case].reshape(365,24)) #将原本的8760*1重构成365*24，再转换成dataframe

load_raw = pv_and_load["tenwithHP"]
load_df = pd.DataFrame(load_raw[:,res_case].reshape(365,24))

price_raw = price["DA2017"] #price包含三列，第一列 spot price 第二列 DSO tariff 第三列 TSO tariff

PC_raw = output["P_C"]
PD_raw = output["P_D"]

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
        
# 最后得到 365*24 的 PV_df load_df PC_df PD_df
# PC的数据是每小时充电的值，PD是放电值
# 之前是把它们加和之后用正负表示充放电，但是也许这里不用加和，之后可以直接取数值求每小时平均