#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:07:24 2020

@author: shihyu
"""


import pandas as pd

Sample list = ["Test01", "Test02", "Test03"]

# =============================================================================
# #Generate info tags for N/T
# N_info_tags=['GT','SQ','AD','AF','F1R2','F2R1','DP','SB','MB','PS']
# T_info_tags=['GT','SQ','AD','AF','F1R2','F2R1','DP','SB','MB','PS']
# for i in range(10):
#     N_info_tags[i]+="_N"
#     T_info_tags[i]+="_T"
# 
# #processing columns
# 
# for i in CLL_WGS_list:
#     file=pd.read_csv(f"ANNOVAR_{i}_DRAGEN.hg19_multianno.txt",sep="\t")
#     normal_info=file["Otherinfo13"]
#     tumor_info=file["Otherinfo14"]
#     n_list=[]
#     for j in range(len(normal_info)):
#         n=normal_info[j].split(":")
#         n_list.append(n)   
#      
#     t_list=[]
#     for k in range(len(tumor_info)):
#         t=tumor_info[k].split(":")
#         t_list.append(t)
#         
#     df1=pd.DataFrame(n_list, columns=N_info_tags)
#     df2=pd.DataFrame(t_list, columns=T_info_tags)
#     final_Annotation_file=pd.concat([file,df1,df2],axis=1)
#     final_Annotation_file.to_excel(f"{i}_DRAGEN.nofilter_annotation_vcfinfo.xlsx", index=False)
# 
# =============================================================================

#Calculate VAF
def calculate_AF_T_N(excel_file):
    val_list=[]
    df=pd.read_excel(excel_file)
    for i in range(len(df)):
        val=df.loc[i]["AF_T"]-df.loc[i]["AF_N"]
        val_list.append(val)
    valdf=pd.DataFrame(val_list,columns=["AF_T-AF_N"])
    finaldf=pd.concat([df,valdf],axis=1)
    return finaldf


#create index
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

# =============================================================================
# for i in CLL_WGS_list:
#     create_idx(calculate_AF_T_N(f"{i}_DRAGEN.nofilter_annotation_vcfinfo.xlsx")).to_excel(f"Indexed_{i}_DRAGEN_nofilter_annotation.xlsx", index=False)
# 
# =============================================================================



#Process union file
# =============================================================================
# pd.read_csv("nofilter_union_of_CLL_WGS_DRAGEN_PASSvariants.txt","\t").to_excel("nofilter_union_of_CLL_WGS_DRAGEN_PASSvariants.xlsx", index=False)
# create_idx(pd.read_excel("nofilter_union_of_CLL_WGS_DRAGEN_PASSvariants.xlsx")).to_excel("Indexed_nofilter_union_of_CLL_WGS_DRAGEN_PASSvariants.xlsx", index=False)
# =============================================================================

#Generate PASSed and 3 filtered array by Idx
U=pd.read_excel("Indexed_nofilter_union_of_CLL_WGS_DRAGEN_PASSvariants.xlsx", index_col="Idx")[["Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","avsnp150","cosmic_coding_GRCh37_v91","cosmic_noncoding_GRCh37_v91","AF_eas.1","TaiwanBiobank-official_Illumina1000-AF", "TaiwanBiobank993WGS_AF"]]

for i in CLL_WGS_list:
    S=pd.read_excel(f"Indexed_{i}_DRAGEN_nofilter_annotation.xlsx",index_col="Idx")
    filter1=(S['Func.refGene'].isin(["exonic","splicing","exonic;splicing"]))&(S['ExonicFunc.refGene'] != "synonymous SNV")
    filter2=(S['AF_eas.1'].isin(["0","."]))&(S['TaiwanBiobank-official_Illumina1000-AF'].isin(["0","."]))
    filter3=(S.AF_T >= 0.05)&(S.AF_T-S.AF_N >= 0.05)
    GATK_filtered_S=S[(S.Otherinfo10 == "PASS")&(filter1&filter2&filter3)]
    ary=U.join(GATK_filtered_S[["Otherinfo14","Otherinfo13"]]).rename(columns={'Otherinfo14':f'{i}_N','Otherinfo13':f'{i}_T'})
    U=ary
ary["Counts"]=ary.iloc[:,11:50:2].count(axis=1)
ary[(ary.Counts>0)].to_excel("VaraintBasedArray_3Filter_PASS_0901.xlsx")





#Make list of unique Genes
ary=pd.read_excel("VaraintBasedArray_3Filter_PASS_0901.xlsx").set_index("Gene.refGene")
idx_list=ary.index.unique().to_list()

#Make Gene-based recurrence array
df=pd.concat([pd.DataFrame(idx_list, columns=['Gene.refGene']),pd.DataFrame(columns=CLL_WGS_list)],axis=1)
df.set_index("Gene.refGene", inplace=True)
for i in idx_list:
    for j in CLL_WGS_list:
        n=0
        if isinstance(ary.loc[i][f"{j}_T"],float): #no data in ary[i][j]
            df.at[i,j]=0
        elif isinstance(ary.loc[i][f"{j}_T"],str): #only one position
            df.at[i,j]=1
        else:
            for k,v in ary.loc[i][f"{j}_T"].items(): #more than 2 positions
                if isinstance(v,str):
                    n+=1
            df.at[i,j]=n
    #calculate the occurrence of position for each gene
    if isinstance(ary.loc[i], pd.DataFrame):
        df.at[i,'POS_SUM']=len(ary.loc[i])
    else:
        df.at[i,'POS_SUM']=1
#calculate the sum of occurrence among samples for each gene
df["Occurrence(variants)"]=df.loc[:,CLL_WGS_list].sum(axis=1)
df["Counts"]=(df.loc[:,CLL_WGS_list] != 0).astype(int).sum(axis=1)
df.to_excel("GeneBasedArray_3Filter_PASS_0901.xlsx")


    
