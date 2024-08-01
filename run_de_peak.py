import os
import sys
import pandas as pd
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt

path = sys.argv[1]

path2 = sys.argv[2]

def cal_deriv(t_x, t_y):
    tx, ty = t_x.values.tolist(), t_y.values.tolist()
    x, y = [], []
    for i in range(len(tx)):
        x.append(tx[i] - 1000*1e-6*ty[i])
        y.append(ty[i])
    #x, y = t_x.values.tolist(), t_y.values.tolist()
    diff_x = []
    for i, j in zip(x[0::], x[1::]):
        diff_x.append(j - i)

    diff_y = []
    for i, j in zip(y[0::], y[1::]):
        diff_y.append(j - i)

    slopes = []
    for i in range(len(diff_y)):
        slopes.append(diff_y[i]/y[i] / (diff_x[i]))
        #print(x[i], diff_y[i], diff_x[i], slopes[i])

    #return x[1::], slopes
    deriv = []
    for i, j in zip(slopes[0::], slopes[1::]):
        deriv.append((0.5 * (i + j))) 
    deriv.insert(0, slopes[0])
    deriv.append(slopes[-1])
    print(len(slopes))
    #for i in slopes:
    #    print(i)
    return x, deriv
files = os.listdir(path)
print(files)
df1 = pd.read_csv(path + '/' + files[0], header=None)

for file in files[1:]:     
  df2 = pd.read_csv(path +'/' +  file,header=None)
  df1 = pd.concat([df1,df2],axis=0,ignore_index=True)

#print (df1)
time = list(set(df1[0].values.tolist()))# df1[0].values.tolist()
print(time)
t_df = df1.loc[(df1[0]==0)&(df1[3]>30)]
df = t_df.sort_values(by = [3], ascending = [True])#sorted(t_df, key = lambda x: x[2],reverse=True)
#print(df.corr())
print (df)

files2 = os.listdir(path2)
print(files2)
df1_2 = pd.read_csv(path2 + '/' + files2[0], header=None)

for file in files2[1:]:
  df2_2 = pd.read_csv(path2 +'/' +  file,header=None)
  df1_2 = pd.concat([df1_2,df2_2],axis=0,ignore_index=True)

t_df2 = df1_2.loc[(df1_2[0]==0)&(df1_2[3]>30)]
df_2 = t_df2.sort_values(by = [3], ascending = [True])

print (df)

peak_V, peak_I = [], []
#for i in range(7,8):
peak_V2, peak_I2 = [], []
predict_peak = []
for i in range(16):
    tem_x, tem_y = cal_deriv(df[3], abs(df[i+4])*1e6)
    #plt.plot(tem_x, tem_y, label='SiPM'+str(i+1))
    peaks, _ = find_peaks(tem_y, prominence=1)
    print(peaks, tem_x[peaks[0]], tem_y[peaks[0]])
    peak_V.append(tem_x[peaks[0]])
    peak_I.append(tem_y[peaks[0]])
    
    tem_x, tem_y = cal_deriv(df_2[3], abs(df_2[i+4])*1e6)
    #plt.plot(tem_x, tem_y, label='SiPM'+str(i+1))
    peaks2, _ = find_peaks(tem_y, prominence=1)
    print(peaks2, tem_x[peaks2[0]], tem_y[peaks2[0]])
    peak_V2.append(tem_x[peaks2[0]])
    peak_I2.append(tem_y[peaks2[0]])
    delta_T  = abs(df[2][peaks[0]] - df_2[2][peaks[0]])
    print('tem',delta_T)
    predict_peak.append(peak_V[i] - 34/1000*delta_T)
    #plt.scatter(tem_x[peaks2], tem_y[peaks2], "ob");
    #plt.scatter(tem_x, tem_y, s = 7, label='SiPM'+str(i+1))
    #plt.scatter(df[2], abs(df[i+3])*1e6, s = 7, label='SiPM'+str(i+1))
tem = str(round(sum(df[2])/len(df[2]), 2))
#print('tem',delta_T)
plt.show()
print(peak_V)
#tem = str(round(sum(df[2])/len(df[2]), 2))
plt.figure(1)
#plt.scatter(peak_V, peak_I);
#sipm = ['SiPm1', 'SiPm2', 'SiPm2', 'SiPm5', 'SiPm5', 'SiPm6', 'SiPm7', 'SiPm8', 'SiPm9', 'SiPm10', 'SiPm11', 'SiPm12', 'SiPm13', 'SiPm14', 'SiPm15', 'SiPm16']
sipm = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
plt.plot(sipm,peak_V,label = '20℃')
plt.plot(sipm,peak_V2, label = '-35℃')
plt.plot(sipm,predict_peak, label = 'prediction -35℃')

plt.title('SiPM breakdown voltage')
plt.xlabel('SiPM')
plt.ylabel('Breakdown Voltage(V)')
d = max(peak_V) - min(peak_V)
#plt.xlim(min(peak_V)-0.01*d, max(peak_V)+0.01*d)
plt.rcParams['font.sans-serif'] = ['sans-serif']
plt.rcParams['axes.unicode_minus'] = False
plt.legend()
plt.grid()
#plt.show()   
plt.savefig("SiPM_V-I_dlnI_peak.png")
plt.show()

