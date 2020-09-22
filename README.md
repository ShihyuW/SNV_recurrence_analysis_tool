# SNV_recurrence_analysis_tool
### Environment: python 3.0+ with "openpyxl" module installed ###
This is a recurrence analysis tool for SNP/Indel analysis.
Due to the unique vcf formats among different platforms, this analysis is TOOL SPECIFIC.

Before running the pipeline, the following files are required:
(1) ANNOVAR annotation files from DRAGEN somatic pipeline(.txt) :Test01 ,Test2, Test3.txt were extract from ANNOVAR annotation files of DRAGEN WGS somatic pipeline results.
(2) A Union txt file of variants form (1) (remove duplicates), which also include every necessary variant-specific information(for DRAGEN vcf files is column 1 to 80) : union_of_test_variants.txt were made from those ANNOVAR.txt files with "find_union.sh" on linux.


Currently, the object-oriented function is still under development.
In this version, users still need to modify the recurrence_analysis.py script for each analyis, including a list of sample file names, union variants file name, and each xlsx file name below.
Also, adding filter (i.e. funtion, VAF, MAF....)to the sample file is possible and can be adjust through line 90-93.
