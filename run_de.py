'''
usage: python3 run_de.py file_path 'filename.txt'
'''
import os
import sys
import pandas as pd
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import math

path = sys.argv[1]
filename = sys.argv[2]

def cal_deriv(t_x, t_y):
    x, y = t_x.values.tolist(), t_y.values.tolist()
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

tot = 0
std_dev = 0
peak_V, peak_I, gaus_v = [], [], []
para = {}
#for i in range(7,8):

for i in range(16):
    tem_x, tem_y = cal_deriv(df[3], abs(df[i+4])*1e6)
    plt.plot(tem_x, tem_y, label='SiPM'+str(i+1))
    #peaks2, _ = find_peaks(tem_y, prominence=1)
    #print(i, peaks2, tem_x[peaks2[0]], tem_y[peaks2[0]])
    #peak_V.append(tem_x[peaks2[0]])
    #peak_I.append(tem_y[peaks2[0]])
    peak_I.append(max(tem_y))
    peak_V.append(tem_x[tem_y.index(max(tem_y))])
    #plt.scatter(tem_x[peaks2], tem_y[peaks2], "ob");
    #plt.scatter(tem_x, tem_y, s = 7, label='SiPM'+str(i+1))
    #plt.scatter(df[2], abs(df[i+3])*1e6, s = 7, label='SiPM'+str(i+1))

spot_ave_vol = sum(peak_V)/len(peak_V)#calculate spot average
t_df_gaus = df1.loc[(df1[0]==0)&(df1[3]>spot_ave_vol-0.5)&(df1[3]<spot_ave_vol+0.5)]#saved for gaus fit within 1 voltage
df_gaus = t_df_gaus.sort_values(by = [3], ascending = [True])#saved for gaus fit within 1 voltage

data = open(filename,'w+') 

for j in range(16):
    print('\n SiPM', abs(j-16))
    gaus_x, gaus_y = cal_deriv(df_gaus[3], abs(df_gaus[j+4])*1e6)
    length = len(gaus_x)
    fit_gra = ROOT.TGraph(length, np.array(gaus_x), np.array(gaus_y))
    
    ffit = ROOT.TF1("ffit", "[0]*exp(-((x-[1])/[2])**2)")
    ffit.SetParameter(0, 15)
    ffit.SetParameter(1, 38.5)
    ffit.SetParameter(2, 0.5)

    fit_gra.Fit('ffit')
    print('\n', 'mean: ', ffit.GetParameter(1))
    gaus_v.append(ffit.GetParameter(1))
    tot += ffit.GetParameter(1)

    print('\n SiPM', abs(j-16), ' ', 'peak_I: ', ffit.GetParameter(0),'mean:', ffit.GetParameter(1), 'sigma: ', ffit.GetParameter(2), file = data)

for i in gaus_v:
    std_dev += (tot/16-i)**2/len(gaus_v)

tem = str(round(sum(df[2])/len(df[2]), 2))
print('\n peak:\n', 'I:', peak_I, '\n', '\n', 'V:', peak_V, file = data)
print('\n', 'spot average Vbd:', sum(peak_V)/len(peak_V), file = data)
print('\n', 'temprature', tem, '℃')
print('\n', 'temprature', tem, '℃', file = data)
print('\n', 'gaus average Vbd', tot/16, ' std_dev', math.sqrt(std_dev))
print('\n', 'gaus average Vbd', tot/16, ' std_dev', math.sqrt(std_dev), file = data)
data.close()
plt.figure(1)
plt.title('SiPM V-lnI ('+tem + '℃)')
plt.xlabel('Voltage(V)')
plt.ylabel('dlnI/dV')
plt.xlim(35, 44)
#plt.rcParams['figure.figsize'] = (18,4)
#plt.rcParams['font.sans-serif'] = ['sans-serif']
plt.rcParams['axes.unicode_minus'] = False
#font.family          : sans-serif
#font.sans-serif      : SimHei
plt.legend(loc=2)
plt.grid()
#plt.semilogy()
#plt.show()   
#plt.gcf().set_size_inches(10.5, 10.5)
plt.savefig("SiPM_V-lnI_dlnI.png")
plt.show()

#print(peak_V)
#tem = str(round(sum(df[2])/len(df[2]), 2))
plt.figure(2)
plt.scatter(peak_V, peak_I)
plt.title('SiPM V-lnI ('+tem + '℃)')
plt.xlabel('Voltage(V)')
plt.ylabel('dlnI/dV')
d = max(peak_V) - min(peak_V)
plt.xlim(min(peak_V)-0.01*d, max(peak_V)+0.01*d)
#plt.rcParams['font.sans-serif'] = ['sans-serif']
plt.rcParams['axes.unicode_minus'] = False
#plt.legend(loc=2)
plt.grid()
#plt.show()   
plt.savefig("SiPM_V-I_dlnI_peak.png")
plt.show()
