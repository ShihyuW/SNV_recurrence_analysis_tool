# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:03:13 2020

@author: Jacob-Lab
"""

import pandas as pd
import numpy as np
import collections as col

# import data
#union_variant = pd.read_table("union_of_variants_test.txt")

# select the column
#union_variant = union_variant[["Chr", "Start", "End", "Func.refGene", "Gene.refGene", "AAChange.refGene" ,"AF_eas.2", "avsnp150"]]

# MDF016 & MDF048 family member ANNOVAR file
Sample_list = ["MDD0026_hg38_annotation","MDD0027_hg38_annotation","MDD0131_hg38_annotation","MDD0132_hg38_annotation","MDD0188_hg38_annotation"]
ID_list = ["MDD0026", "MDD0027","MDD0131","MDD0132","MDD0188"]
# generate info tags for N/T
# Otherinfo12 (there are three format...)
# GT:AD:AF:DP:F1R2:F2R1:GQ:PL:GP:PRI:SB:MB:PS
# GT:AD:AF:DP:F1R2:F2R1:GQ:PL:GP:PRI:SB:MB
# GT:AD:AF:F1R2:F2R1:DP:SB:MB
N_info_tags=['GT','AD','AF','DP','F1R2','F2R1','GQ','PL','GP','PRI','SB','MB','PS']
for i in range(13):
    N_info_tags[i]+="_N"
    
# processing columns
for i in Sample_list:
    file=pd.read_csv(f"{i}.txt",sep="\t")
    normal_info=file["Otherinfo13"]

    
    n_list=[]
    for j in range(len(normal_info)):
        n=normal_info[j].split(":")
        n_list.append(n)   
         
    df1=pd.DataFrame(n_list, columns=N_info_tags)
    final_Annotation_file=pd.concat([file,df1],axis=1)
    final_Annotation_file.to_excel(f"{i}.xlsx", index=False)

#Calculate VAF
def calculate_AF_T_N(excel_file):
    val_list=[]
    df=pd.read_excel(excel_file)
    for i in range(len(df)):
        val=df.loc[i]["GQ_N"]-df.loc[i]["AF_N"]
        val_list.append(val)
    valdf=pd.DataFrame(val_list,columns=["GQ_N-AF_N"])
    finaldf=pd.concat([df,valdf],axis=1)
    return finaldf

# create index
def create_idx(df):
    idx_list=[]
    L=df
    L["Start"]=L["Start"].apply(str)
    L["End"]=L["End"].apply(str)
  
    
    for i in range(len(L)):
        idx=L.loc[i]["Chr"]+"/"+L.loc[i]["Start"]+"/"+L.loc[i]["End"]+"/"+L.loc[i]["Ref"]+"/"+L.loc[i]["Alt"]
        idx_list.append(idx)
    idx_df=pd.DataFrame(idx_list,columns=["Idx"])
    finaldf=pd.concat([L,idx_df], axis=1)
    return finaldf

for i in Sample_list:
    create_idx(pd.read_excel(f"{i}.xlsx")).to_excel(f"Indexed_{i}.xlsx", index=False)

# Process union file
pd.read_csv("union_of_variants.txt",sep="\t").to_excel("union_of_variants.xlsx", index=False)
create_idx(pd.read_excel("union_of_variants.xlsx")).to_excel("Indexed_union_of_variants.xlsx", index=False)

# Generate PASSed and 3 filtered array by Idx
U=pd.read_excel("Indexed_union_of_variants.xlsx", index_col="Idx")[["Chr","Start","End","Func.refGene", "ExonicFunc.refGene" ,"Gene.refGene", "AAChange.refGene" ,"AF_eas.1", "avsnp150"]]

for i in Sample_list:
    S=pd.read_excel(f"Indexed_{i}.xlsx",index_col="Idx")
    
    #Filters (Add if needed)
    #filter1=(S['Func.refGene'].isin(["exonic","splicing","exonic;splicing"]))&(S['ExonicFunc.refGene'] != "synonymous SNV")
    #filter2=(S['AF_eas.1'].isin(["0","."]))&(S['TaiwanBiobank-official_Illumina1000-AF'].isin(["0","."]))
    #filter3=(S.AF_T >= 0.05)&(S.AF_T-S.AF_N >= 0.05)
    #GATK_filtered_S=S[(S.Otherinfo10 == "PASS")&(filter1&filter2&filter3)]
    
    ary=U.join(S["GT_N"]).rename(columns={'GT_N':f'{i}'})
    U=ary

ary["Counts"]=ary.iloc[:,8:13].count(axis=1)
#print(ary.iloc[1,8:13])
ary.to_excel("VaraintBasedArray.xlsx")
