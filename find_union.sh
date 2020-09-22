#!/bin/bash

#Get target file names and make name list(optional)
cd /work2/u1067478/Test/
ls -l /work2/u1067478/Test/Test*.txt|awk '{print $9}'>variants_files_list.txt

#Write every varints with info (for DRAGEN, $1-$80) to a nofilter_combine_variants_info.txt

while read variants_file;do
awk -F "\t" 'BEGIN {OFS = "\t"} NR == 1 {next} {NF=80}1' ${variants_file}>>nofilter_combine_variants_info.txt
done<variants_files_list.txt

#Get header from one of the target file
head -n 1 /work2/u1067478/Test/Test01.txt|awk -F "\t" 'BEGIN {OFS = "\t"} {NF=80}1' > nofilter_union_of_CLL_WES_variants.txt

#Get the union of all variants with information
sort nofilter_combine_variants_info.txt |uniq >> nofilter_union_of_CLL_WES_variants.txt
