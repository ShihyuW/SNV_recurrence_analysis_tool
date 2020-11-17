#!/usr/bin/bash
#PBS -q ntu192G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N python_test
#PBS -M s0890003@gmail.com
#PBS -j oe
#PBS -m e

export PATH=/pkg/biology/Python/Python3_default/bin:$PATH
PYTHON3=/pkg/biology/Python/Python3_default/bin/python3
SCRIPT_PATH="/work2/lynn88065/script/github/SNV_recurrence_analysis_tool/scripts"
$PYHTON3 $SCRIPT_PATH/germline_snv_recurrence.py

