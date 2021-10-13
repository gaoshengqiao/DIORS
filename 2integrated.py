import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import matplotlib as mpl
import matplotlib.pylab as pylab
def plot(metric,alpha,c):
    ax.scatter(metric.iloc[:,0],metric.iloc[:,1],s=50,alpha=alpha,norm=0.1,c = c)
    
measurement = "Jaccard" #"Signed Jaccard" or "Jaccard"

result_MD = pd.read_csv('docking score.txt',sep='\t')
result_SM = pd.read_csv('Jaccard score.csv',sep=',')
result_MD1 = -result_MD['delta_E']
result_MD1.index = result_MD['pert_iname']
result_SM1 = result_SM[measurement]
result_SM1.index = result_SM['pert_iname']

tmp = pd.concat([result_MD1,result_SM1],axis=1)
tmp.columns = ['-delta_E','Jaccard_score']

tmp = tmp[tmp['-delta_E']>=0]

tmp = tmp.sort_values(by='-delta_E',ascending=False)
tmp.loc[:,'Rd']= range(1,len(tmp)+1)
tmp = tmp.sort_values(by='Jaccard_score',ascending=False)
tmp.loc[:,'Rs']= range(1,len(tmp)+1)
tmp.loc[:,'Rs+Rd']= tmp.loc[:,'Rd']+tmp.loc[:,'Rs']
tmp = tmp.sort_values(by='Rs+Rd',ascending=True)
tmp.loc[:,'Ri']= range(1,len(tmp)+1)

font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 12,
         }

fig = plt.figure()
num=10 #top10标记
tmp = tmp.sort_values(by='-delta_E',ascending=False)
tmp0 = tmp.iloc[:100,:]
tmp1 = tmp.iloc[:num,:]
tmp2 = tmp.sort_values(by='Ri',ascending=True)
tmp3 = tmp2.iloc[:num,:]

ax = fig.add_subplot(111)
tmp0 = tmp.iloc[:100,:]

plot(tmp.loc[:,['-delta_E','Jaccard_score']],0.6,'grey')
plot(tmp0.loc[:,['-delta_E','Jaccard_score']],0.6,'g')
plot(tmp1.loc[:,['-delta_E','Jaccard_score']],0.8,'orangered')
plot(tmp3.loc[:,['-delta_E','Jaccard_score']],0.8,'royalblue')

plt.ylabel('Jaccard score', font2)
plt.xlabel('-ΔE', font2)
ax.spines['top'].set_visible(False) 
ax.spines['right'].set_visible(False) 
plt.legend(['All drugs','Top100','Top10 (molecular docking)','Top10 (intergrated approach)'],loc='upper left')

plt.show()
