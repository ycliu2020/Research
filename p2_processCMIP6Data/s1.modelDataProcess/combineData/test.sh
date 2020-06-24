#!/bin/bash      # shell type
shopt -s extglob # enable extended globbing
###
 # @Author       : LYC
 # @Date         : 2020-06-23 10:38:14
 # @LastEditTime : 2020-06-23 10:51:57
 # @LastEditors  : LYC
 # @Description  : 
 # @FilePath     : /code/p2_processCMIP6Data/s1.modelDataProcess/combineData/test.sh
 # @ 
### 
loc='/data1/liuyincheng/CMIP6-mirror/amip/CMIP/' # Directory of input files
xpt=( 'amip' )                # Experiment ( could be more )
echo "Model Experiment:" ${xpt}.

var=("hus" "ta" "ts" "ps" "hfls" "hfss"\
    "rlus" "rsus" "rsds" "rlds" "rsuscs" "rsdscs" "rldscs" \
"rlut" "rsut" "rsdt" "rlutcs" "rsutcs" "clisccp" "ua" "va" "wap") # 22 variables
rlm='Amon'                         # Realm
echo "Model realm:" ${rlm}.

for var_id in {0..0}; do            # Loop over N variables
    # Names of all the models (ls [get file names];
    #  cut [get model names];
    #  sort; uniq [remove duplicates]; awk [print])
    mdl_set=$( ls ${loc}${var[var_id]}_${rlm}_*_${xpt[0]}_*.nc | \
    cut -d '_' -f 3 | sort | uniq -c | awk '{print $2}' )
    # Number of models (echo [print contents]; wc [count])
    mdl_nbr=$( echo ${mdl_set} | wc -w )
    echo "=============================="
    echo "There are" ${mdl_nbr} "models for" ${var[var_id]}.
    echo "Model list:" ${mdl_set}.
    mdl_nbr_id[var_id]=${mdl_nbr}
    for mdl in ${mdl_set}; do	        # Loop over models
        # Names of all the ensemble members
        nsm_set=$( ls ${loc}${var[var_id]}_${rlm}_${mdl}_${xpt[0]}_*.nc | \
        cut -d '_' -f 5 | sort | uniq -c | awk '{print $2}' )
        # Number of ensemble members in each model
        nsm_nbr=$( echo ${nsm_set} | wc -w )
        echo "------------------------------"
        echo "Model" ${mdl} "includes" ${nsm_nbr} "ensemble member(s):"
        echo ${nsm_set}"."
        
    done
done
