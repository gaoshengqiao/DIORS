import pandas as pd
import numpy as np
import os
from time import time

#填入计算类型
way = "reverse" #"mimic" or "reverse"

#填入计算方法
measurement = "Jaccard" #"Signed Jaccard" or "Jaccard"

#填入细胞系
cell='' 

#填入基因形式
genetype ='gene symbol ID' #'entrez ID' or "gene symbol ID"

def Jaccard(S_up,S_up2):
    INT =S_up.intersection(S_up2)
    UNI= S_up.union(S_up2)
    jaccard =len(INT)/len(UNI)
    return jaccard

data1 = pd.read_csv('./signature matching/911 drugs_part0.csv')
data2 = pd.read_csv('./signature matching/911 drugs_part1.csv')
data3 = pd.read_csv('./signature matching/911 drugs_part2.csv')
data4 = pd.read_csv('./signature matching/911 drugs_part3.csv')
data = pd.concat([data1,data2,data3,data4])
data.index=range(len(data))
query =pd.read_csv('./signature matching/input signature.csv')
gene = pd.read_csv('./signature matching/Transcribed_Gene_Metadata.txt',sep='\t')

if len(cell)>0:
    data = data[data['cell']==cell]
    
if genetype =='gene symbol ID':
    if str(query.loc[0,'up_GENE'])!='nan':    
        query_DE_up = pd.Index(set(query.loc[0,'up_GENE'].split(",")))
        query_DE_up = pd.DataFrame(query_DE_up)
        query_DE_up2 = np.unique(query_DE_up.merge(gene,left_on=0,right_on='pr_gene_symbol')['pr_gene_id'])
        query_DE_up = pd.Index(query_DE_up2)
    else:
        query_DE_up = pd.Index([])

    if str(query.loc[0,'down_GENE'])!='nan':
        query_DE_down =  pd.Index(set(query.loc[0,'down_GENE'].split(",")))
        query_DE_down = pd.DataFrame(query_DE_down)
        query_DE_down2 = np.unique(query_DE_down.merge(gene,left_on=0,right_on='pr_gene_symbol')['pr_gene_id'])
        query_DE_down = pd.Index(query_DE_down2)
    else:
        query_DE_down = pd.Index([])

if genetype =='entrez ID':
    if str(query.loc[0,'up_GENE'])!='nan':    
        query_DE_up = pd.Index(set(query.loc[0,'up_GENE'].split(",")))
    else:
        query_DE_up = pd.Index([])

    if str(query.loc[0,'down_GENE'])!='nan':
        query_DE_down =  pd.Index(set(query.loc[0,'down_GENE'].split(",")))
    else:
        query_DE_down = pd.Index([])

n_up = len(query_DE_up)
n_down = len(query_DE_down)
out =np.zeros((len(data),1))
overlapped_query_DE_up=list()
overlapped_query_DE_down=list()

for J in range(len(data)):
    tmp=data.loc[J,'up_GENE']
    tmp2=data.loc[J,'down_GENE']
    DEG_up = pd.Index(set((tmp).split(",")))
    DEG_down = pd.Index((set((tmp2).split(","))))
        
    if (len(query_DE_up)!=0,len(query_DE_down)!=0) == (True,True):  #输入的上下调基因均不为空
        if way == "reverse":
            if measurement == "Signed Jaccard":
                J_score = (Jaccard(query_DE_up,DEG_down)+Jaccard(query_DE_down,DEG_up)-Jaccard(query_DE_up,DEG_up)-Jaccard(query_DE_down,DEG_down))/2
            elif measurement == "Jaccard":
                J_score = (Jaccard(query_DE_up,DEG_down)+Jaccard(query_DE_down,DEG_up))/2
            overlapped_query_DE_up.append(str(list(query_DE_up.intersection(DEG_down))).replace("'","").replace(" ","").replace("]","").replace("[",""))
            overlapped_query_DE_up.append(str(list(query_DE_down.intersection(DEG_up))).replace("'","").replace(" ","").replace("]","").replace("[",""))

        elif way == "mimic":
            if measurement == "Signed Jaccard":
                J_score = (Jaccard(query_DE_up,DEG_up)+Jaccard(query_DE_down,DEG_down)-Jaccard(query_DE_up,DEG_down)-Jaccard(query_DE_down,DEG_up))/2
            elif measurement == "Jaccard":
                J_score = (Jaccard(query_DE_up,DEG_up)+Jaccard(query_DE_down,DEG_down))/2
            overlapped_query_DE_up.append(str(list(query_DE_up.intersection(DEG_up))).replace("'","").replace(" ","").replace("]","").replace("[",""))
            overlapped_query_DE_up.append(str(list(query_DE_down.intersection(DEG_down))).replace("'","").replace(" ","").replace("]","").replace("[",""))
            
    elif (len(query_DE_up)!=0,len(query_DE_down)!=0) == (True,False):  #输入的上调基因不为空，下调基因为空
        if way == "reverse":
            if measurement == "Signed Jaccard":
                J_score = (Jaccard(query_DE_up,DEG_down)-Jaccard(query_DE_up,DEG_up))
            elif measurement == "Jaccard":
                J_score = (Jaccard(query_DE_up,DEG_down))
            overlapped_query_DE_up.append(str(list(query_DE_up.intersection(DEG_down))).replace("'","").replace(" ","").replace("]","").replace("[",""))
            
        elif way == "mimic":
            if measurement == "Signed Jaccard":
                J_score = (Jaccard(query_DE_up,DEG_up)-Jaccard(query_DE_up,DEG_down))
            elif measurement == "Jaccard":
                J_score = (Jaccard(query_DE_up,DEG_up))                  
            overlapped_query_DE_up.append(str(list(query_DE_up.intersection(DEG_up))).replace("'","").replace(" ","").replace("]","").replace("[",""))
      
    elif (len(query_DE_up)!=0,len(query_DE_down)!=0) == (False,True):  #输入的下调基因不为空，上调基因为空
        if way == "reverse":
            if measurement == "Signed Jaccard":
                J_score = (Jaccard(query_DE_down,DEG_up)-Jaccard(query_DE_down,DEG_down))
            elif measurement == "Jaccard":
                J_score = (Jaccard(query_DE_down,DEG_up))
            overlapped_query_DE_up.append(str(list(query_DE_down.intersection(DEG_up))).replace("'","").replace(" ","").replace("]","").replace("[",""))
            
        elif way == "mimic":
            if measurement == "Signed Jaccard":
                J_score = (Jaccard(query_DE_down,DEG_down)-Jaccard(query_DE_down,DEG_up))
            elif measurement == "Jaccard":
                J_score = (Jaccard(query_DE_down,DEG_down))
            overlapped_query_DE_up.append(str(list(query_DE_down.intersection(DEG_down))).replace("'","").replace(" ","").replace("]","").replace("[",""))
            
    out[J,0]=J_score
        
out = pd.DataFrame(out,columns = [measurement])
out=pd.concat([out,pd.DataFrame(np.array(overlapped_query_DE_up),columns=["instersection with queried up"])],axis=1)
out=pd.concat([out,pd.DataFrame(np.array(overlapped_query_DE_down),columns=["instersection with queried down"])],axis=1)
out=pd.concat([data[['pert_iname','sig_id','cell','pert_idose','pert_itime(h)']],out],axis=1)
out = out.sort_values(by=measurement,ascending=False)

#每个药只保留最高分数
drugs = np.unique(out['pert_iname'])
SM_dock_df_filter = pd.DataFrame()
for drug in drugs:
    tmp = out[out['pert_iname']==drug]
    tmp = tmp.sort_values(by = measurement,ascending=False).iloc[:1,:]
    SM_dock_df_filter = pd.concat([SM_dock_df_filter,tmp])
SM_dock_df_filter = SM_dock_df_filter.sort_values(by =measurement,ascending=False)   
SM_dock_df_filter.to_csv("Jaccard score.csv",index=False)

