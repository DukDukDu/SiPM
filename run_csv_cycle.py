import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

path = sys.argv[1]
files = os.listdir(path)
#print(files)
df1 = pd.read_csv(path + '/' + files[0], header=None)

for file in files[1:]:     
    df2 = pd.read_csv(path +'/' +  file,header=None)
    #print(file)
    df1 = pd.concat([df1,df2],axis=0,ignore_index=True)
    #print(file)
#print (df1[0])
time = list(set(df1[0].values.tolist()))# df1[0].values.tolist()
#print(time)
#df = df1.loc[(df1[0]==time[6])&(df1[2]>0)]
t_df = df1.loc[df1[3]>0]
df = t_df.sort_values(by = [1], ascending = [True])
#print (df)

plt.subplot(211)
for i in range(16):
    tmp_df = df.loc[df[5]==i+1].sort_values(by = [1], ascending = [True])
    print(tmp_df)
    tx, ttx = [],tmp_df[1].values.tolist()
    min_t = min(ttx)
    for k in range(len(ttx)):
        tx.append(ttx[k] - min_t)
    work_v = str(round(sum(tmp_df[3])/len(tmp_df[3]), 2))
    #plt.scatter(df[2], df[1], s = 7, label='SiPM'+str(i+1))
    plt.plot(tx, -1*(tmp_df[4])*1e9,  label='SiPM'+str(i+1)+' '+str(work_v)+'V')
    #plt.plot(tx, tmp_df[2],  label='SiPM'+str(i+1)+' '+str(work_v)+'V')

    #plt.scatter(df[2], abs(df[i+3])*1e6, s = 7, label='SiPM'+str(i+1))

#tem = str(round(sum(df[2])/len(df[2]), 2))
#print(tem)
#plt.title('SiPM V-I ('+tem + '℃)')
plt.xlabel('time(S)')
plt.ylabel('I(nA)')
#d = max(df[3]) - min(df[3])
#plt.xlim(min(df[3] - 0.01*d), max(df[3]) +0.01*d)
#plt.rcParams['figure.figsize'] = (18,4)
plt.rcParams['font.sans-serif'] = ['sans-serif']
plt.rcParams['axes.unicode_minus'] = False
#font.family          : sans-serif
#font.sans-serif      : SimHei
#plt.legend(bbox_to_anchor=(1., 0), loc=2)
plt.legend( loc=2)
plt.grid()
plt.semilogy()

#temperature
plt.subplot(212)
for i in range(16):
    tmp_df = df.loc[df[5]==i+1].sort_values(by = [1], ascending = [True])
    print(tmp_df)
    tx, ttx = [],tmp_df[1].values.tolist()
    min_t = min(ttx)
    for k in range(len(ttx)):
        tx.append(ttx[k] - min_t)
    work_v = str(round(sum(tmp_df[3])/len(tmp_df[3]), 2))
    #plt.scatter(df[2], df[1], s = 7, label='SiPM'+str(i+1))
    #plt.plot(tx, -1*(tmp_df[4])*1e9,  label='SiPM'+str(i+1)+' '+str(work_v)+'V')
    plt.plot(tx, tmp_df[2],  label='SiPM'+str(i+1)+' '+str(work_v)+'V')

    #plt.scatter(df[2], abs(df[i+3])*1e6, s = 7, label='SiPM'+str(i+1))

#tem = str(round(sum(df[2])/len(df[2]), 2))
#print(tem)
#plt.title('SiPM V-I ('+tem + '℃)')
plt.xlabel('time(S)')
plt.ylabel('T(℃)')
#d = max(df[3]) - min(df[3])
#plt.xlim(min(df[3] - 0.01*d), max(df[3]) +0.01*d)
#plt.rcParams['figure.figsize'] = (18,4)
plt.rcParams['font.sans-serif'] = ['sans-serif']
plt.rcParams['axes.unicode_minus'] = False
#font.family          : sans-serif
#font.sans-serif      : SimHei
#plt.legend(bbox_to_anchor=(1., 0), loc=2)
#plt.legend( loc=2)
plt.grid()
#plt.semilogy()
#plt.show()   
#plt.gcf().set_size_inches(10.5, 10.5)
plt.savefig("SiPM_V-I.png")
plt.show()
