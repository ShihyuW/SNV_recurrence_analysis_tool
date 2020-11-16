#!/usr/bin/env python3
# coding=utf-8

"""
Created on Wed Oct 28 16:03:13 2020

@author: YenlingPeng
"""

import pandas as pd
import os as os

# set the input/output path 
input_path='/work2/lynn88065/script/github/SNV_recurrence_analysis_tool/inputs/'
output_path='/work2/lynn88065/script/github/SNV_recurrence_analysis_tool/outputs/'

# change the directory
os.chdir('{}'.format(input_path))

# family member ID from ANNOVAR file
Sample_list = ["germline_test01","germline_test02"]
union_file_name = 'germline_union_of_variants_test.txt'
union_file = pd.read_csv('{}'.format(union_file_name),sep="\t")

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

# generate info tags for Otherinfo12
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
# create index on each annotation file
    create_idx(final_Annotation_file)

# create index on union file

union_file = create_idx(union_file)

# Generate PASSed and 3 filtered array by Idx
U=pd.DataFrame(union_file)[["Chr","Start","End","Func.refGene", "ExonicFunc.refGene" ,"Gene.refGene", "AAChange.refGene" ,"AF_eas.1", "AF_popmax.1" , "avsnp150"]]

# Process union file
for i in Sample_list:
    S=pd.DataFrame(final_Annotation_file)
    
    #Filters (Add if needed)
    #filter1=(S['Func.refGene'].isin(["exonic","splicing","exonic;splicing"]))&(S['ExonicFunc.refGene'] != "synonymous SNV")
    #filter2=(S['AF_eas.1'].isin(["0","."]))&(S['TaiwanBiobank-official_Illumina1000-AF'].isin(["0","."]))
    #filter3=(S.AF_T >= 0.05)&(S.AF_T-S.AF_N >= 0.05)
    #GATK_filtered_S=S[(S.Otherinfo10 == "PASS")&(filter1&filter2&filter3)]
    
    ary=U.join(S["GT_N"]).rename(columns={'GT_N':"{i}"})
    U=ary

ary["Counts"]=ary.iloc[:,-3:-1].count(axis=1)
ary = ary.fillna(".")

# check the count column
#print(ary.iloc[1,-3:-1])

# export the result to output_path
ary.to_excel('{}'.format(output_path) + "/germline_VaraintBasedArray.xlsx")


