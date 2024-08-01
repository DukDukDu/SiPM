import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

path = sys.argv[1]

def read_csv(t_path):
    files = os.listdir(t_path)
    print(files)
    t_df1 = pd.read_csv(t_path + '/' + files[0], header=None)

    for file in files[1:]:
        print (file)
        t_df2 = pd.read_csv(t_path +'/' +  file,header=None)
        t_df1 = pd.concat([t_df1,t_df2],axis=0,ignore_index=True)

    print (t_df1[0])
    time = list(set(t_df1[0].values.tolist()))# df1[0].values.tolist()
    print(time)
    #df = df1.loc[(df1[0]==time[6])&(df1[2]>0)]
    t_df = t_df1.loc[(t_df1[0]==0)&(t_df1[3]>0)]
    t_df = t_df.sort_values(by = [3], ascending = [True])
    return t_df


df = read_csv(path)
print (df)
scale = -1
for i in range(16):
    #plt.scatter(df[2], df[1], s = 7, label='SiPM'+str(i+1))
    if sum(df[4]) >0:
        scale = 1
    else:
        scale = -1

    plt.plot(df[3], scale*(df[i+4])*1e6,  label='SiPM'+str(i+1))

    #plt.scatter(df[2], abs(df[i+3])*1e6, s = 7, label='SiPM'+str(i+1))

tem = str(round(sum(df[2])/len(df[2]), 2))
print(tem)
plt.title('SiPM V-I ('+tem + '℃)')
plt.xlabel('Voltage(V)')
plt.ylabel('I($\mu$A)')
d = max(df[3]) - min(df[3])
if scale >0:
    plt.xlim(min(df[3] - 0.01*d), 0.45)
else:
    plt.xlim(min(df[3] - 0.01*d), max(df[3]) +0.01*d)
#plt.xlim(min(df[3] - 0.01*d), max(df[3]) +0.01*d)
#plt.rcParams['figure.figsize'] = (18,4)
#plt.rcParams['font.sans-serif'] = ['sans-serif']
plt.rcParams['axes.unicode_minus'] = False
#font.family          : sans-serif
#font.sans-serif      : SimHei
plt.legend(loc=2)
plt.grid()
plt.semilogy()
#plt.show()   
#plt.gcf().set_size_inches(10.5, 10.5)
plt.savefig("SiPM_V-I.png")
plt.show()
