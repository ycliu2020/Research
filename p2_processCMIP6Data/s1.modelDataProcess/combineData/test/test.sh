###
# @Author       : LYC
# @Date         : 2020-06-29 10:50:00
 # @LastEditTime : 2020-06-29 11:17:20
 # @LastEditors  : LYC
# @Description  :
 # @FilePath     : /code/p2_processCMIP6Data/s1.modelDataProcess/combineData/test/test.sh
# @
###

#!/bin/bash      # shell type
shopt -s extglob # enable extended globbing
#===========================================================================
# Some of the models cut one ensemble member into several files,
#  which include data of different time periods.
# We'd better concatenate them into one at the beginning so that
#  we won't have to think about which files we need if we want
#  to retrieve a specific time period later.
#
# Method:
#	- Make sure 'time' is the record dimension (i.e., left-most)
#	- ncrcat
#
# Input files like:
# /data/cmip5/snc_LImon_bcc-csm1-1_historical_r1i1p1_185001-190012.nc
# /data/cmip5/snc_LImon_bcc-csm1-1_historical_r1i1p1_190101-200512.nc
#
# Output files like:
# /data/cmip5/snc_LImon_bcc-csm1-1_historical_r1i1p1_185001-200512.nc
#
# Online: http://nco.sourceforge.net/nco.html#Combine-Files
#
# Execute this script: bash cmb_fl.sh
#===========================================================================
# Attention
# 1.dont change this script during runing
# 2.dont break the run during runing, if break, please rm the temp file first
# 3.amip: first move /amip/CMIP/*.nc to /amip, then run the script
#
# Name statement
# exp: Experiment
# mdl: model
# esm: ensemble member
# grd: grid_label

locPath='/data1/liuyincheng/CMIP6-mirror/' # Directory of input files
exp_set='ssp245 ssp370'              # Experiment ( could be more )
echo ${#exp_set[@]}
tab='Amon'                               # table_ID
var_set=("hus" "ta" "ts" "ps" "hfls" "hfss"\
    "rlus" "rsus" "rsds" "rlds" "rsuscs" "rsdscs" "rldscs" \
"rlut" "rsut" "rsdt" "rlutcs" "rsutcs" "clisccp" "ua" "va" "wap") # 22 variables


ii=$((${#exp_set[@]}))
echo$ii
for exp in ${exp_set}; do
    echo $exp
    loc=$locPath$exp"/"
    jj=$((${#var_set[@]} - 1))
    echo $jj
done


# for exp in ${exp_set}; do
#     loc=$locPath$exp"/"
#     jj=$((${#var_set[@]} - 1))
#     echo $jj
#     for var_id in `eval echo {0..$jj}`; do            # Loop over two variables
#         # Names of all the models (ls [get file names];
#         #  cut [get model names];
#         #  sort; uniq [remove duplicates]; awk [print])
#         mdl_set=$( ls ${loc}${var_set[var_id]}_${tab}_*_${exp}_*.nc | \
#         cut -d '_' -f 3 | sort | uniq -c | awk '{print $2}' )
#         # Number of models (echo [print contents]; wc [count])
#         mdl_nbr=$( echo ${mdl_set} | wc -w )
#         echo "=============================="
#         echo "There are" ${mdl_nbr} "models for" ${var_set[var_id]}.
        
#         for mdl in ${mdl_set}; do	        # Loop over models
#             # Names of all the ensemble members
#             esm_set=$( ls ${loc}${var_set[var_id]}_${tab}_${mdl}_${exp}_*.nc | \
#             cut -d '_' -f 5 | sort | uniq -c | awk '{print $2}' )
#             # Number of ensemble members in each model
#             esm_nbr=$( echo ${esm_set} | wc -w )
#             echo "------------------------------"
#             echo "Model" ${mdl} "includes" ${esm_nbr} "ensemble member(s):"
#             echo ${esm_set}"."
            
#             for esm in ${esm_set}; do	      # Loop over ensemble members
#                 # Names of all the grid_labels
#                 grd_set=$( ls ${loc}${var_set[var_id]}_${tab}_${mdl}_${exp}_${esm}_*.nc | \
#                 cut -d '_' -f 6 | sort | uniq -c | awk '{print $2}' )
#                 # Number of ensemble members in each model
#                 grd_nbr=$( echo ${grd_set} | wc -w )
#                 echo "------------------------------"
#                 echo "ensemble member" ${esm} "includes" ${grd_nbr} "grid label(s):"
#                 echo ${grd_set}"."
                
#             done
#         done
#     done
# done