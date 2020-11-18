#!/bin/bash
#PBS -q ntu192G
#PBS -l select=1:ncpus=10
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N find_union_SNV
#PBS -j oe
#PBS -M s0890003@gmail.com
#PBS -m e

INPUT_PATH=/work2/lynn88065/script/github/SNV_recurrence_analysis_tool/inputs
cd $INPUT_PATH

# make the sample list for finding union variant 
ls -l ${INPUT_PATH}/germline_test* |awk '{print $9}' > germline_variants_files_list.txt

# find the union variant by position
# Write every varints with info (for DRAGEN, $1-$66) to a nofilter_combine_variants_info.txt
while read variants_file; do
awk -F "\t" 'BEGIN {OFS = "\t"} NR == 1 {next} {NF=66}1' ${variants_file} >> germline_combine_variants_info.txt
done<${INPUT_PATH}/germline_variants_files_list.txt

# Get header from one of the target file
head -n 1 ${INPUT_PATH}/germline_test01.txt \
         |awk -F "\t" 'BEGIN {OFS = "\t"} {NF=66}1' > germline_union_of_variants_test.txt

# Get the union of all variants with information
sort -n -k1.4 germline_combine_variants_info.txt |uniq >>  germline_union_of_variants_test.txt

