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
exp_set='ssp245 ssp370'               # Experiment ( could be more )
tab='Amon'                               # table_ID
var_set=("hus" "ta" "ts" "ps" "hfls" "hfss"\
    "rlus" "rsus" "rsds" "rlds" "rsuscs" "rsdscs" "rldscs" \
"rlut" "rsut" "rsdt" "rlutcs" "rsutcs" "clisccp" "ua" "va" "wap") # 22 variables

for exp in ${exp_set}; do
    loc=$locPath$exp"/"
    jj=$((${#var_set[@]} - 1))
    echo $exp "start: there are " $jj "vars"
    for var_id in `eval echo {0..$jj}`; do            # Loop over two variables
        # Names of all the models (ls [get file names];
        #  cut [get model names];
        #  sort; uniq [remove duplicates]; awk [print])
        mdl_set=$( ls ${loc}${var_set[var_id]}_${tab}_*_${exp}_*.nc | \
        cut -d '_' -f 3 | sort | uniq -c | awk '{print $2}' )
        # Number of models (echo [print contents]; wc [count])
        mdl_nbr=$( echo ${mdl_set} | wc -w )
        echo "=============================="
        echo "There are" ${mdl_nbr} "models for" ${var_set[var_id]}.
        
        for mdl in ${mdl_set}; do	        # Loop over models
            # Names of all the ensemble members
            esm_set=$( ls ${loc}${var_set[var_id]}_${tab}_${mdl}_${exp}_*.nc | \
            cut -d '_' -f 5 | sort | uniq -c | awk '{print $2}' )
            # Number of ensemble members in each model
            esm_nbr=$( echo ${esm_set} | wc -w )
            echo "------------------------------"
            echo "Model" ${mdl} "includes" ${esm_nbr} "ensemble member(s):"
            echo ${esm_set}"."
            
            for esm in ${esm_set}; do	      # Loop over ensemble members
                # Names of all the grid_labels
                grd_set=$( ls ${loc}${var_set[var_id]}_${tab}_${mdl}_${exp}_${esm}_*.nc | \
                cut -d '_' -f 6 | sort | uniq -c | awk '{print $2}' )
                # Number of ensemble members in each model
                grd_nbr=$( echo ${grd_set} | wc -w )
                echo "------------------------------"
                echo "ensemble member" ${esm} "includes" ${grd_nbr} "grid label(s):"
                echo ${grd_set}"."

                for grd in ${grd_set}; do	      # grid label
                    # Number of files in this ensemble member
                    fl_nbr=$( ls ${loc}${var_set[var_id]}_${tab}_${mdl}_${exp}_${esm}_${grd}_*.nc \
                    | wc -w )
                    
                    # If there is only 1 file, continue to next loop
                    if [ ${fl_nbr} -le 1 ]
                    then
                        echo "There is only 1 file in" ${esm} ${grd}.
                        continue
                    fi
                    
                    # If there is only 0 file, mean dont exist var_set continue to next loop
                    if [ ${fl_nbr} -le 0 ]
                    then
                        echo "There is only 0 file in" ${esm} ${grd}.
                        continue
                    fi
                    
                    echo "There are" ${fl_nbr} "files in" ${esm}  ${grd}.
                    
                    # Starting date of data
                    #   (sed [the name of the first file includes the starting date])
                    yyyymm_str=$( ls ${loc}${var_set[var_id]}_${tab}_${mdl}_${exp}_${esm}_${grd}_*.nc\
                    | sed -n '1p' | cut -d '_' -f 7 | cut -d '-' -f 1 )
                    # Ending date of data
                    #   (sed [the name of the last file includes the ending date])
                    yyyymm_end=$( ls ${loc}${var_set[var_id]}_${tab}_${mdl}_${exp}_${esm}_${grd}_*.nc\
                    | sed -n "${fl_nbr}p" | cut -d '_' -f 7 | cut -d '-' -f 2 )
                    
                    # Concatenate one ensemble member files
                    #   into one along the record dimension (now is time)
                    ncrcat -O ${loc}${var_set[var_id]}_${tab}_${mdl}_${exp}_${esm}_${grd}_*.nc \
                    ${loc}${var_set[var_id]}_${tab}_${mdl}_${exp}_${esm}_${grd}_${yyyymm_str}-${yyyymm_end}
                    
                    # Remove useless files
                    rm ${loc}${var_set[var_id]}_${tab}_${mdl}_${exp}_${esm}_${grd}_!(${yyyymm_str}-${yyyymm_end})
                    
                done
                
            done
        done
    done
done
# old scrip back up
# 1) amip-hist (not the  Primary)
# locPath='/data1/liuyincheng/CMIP6-mirror/amip-hist/' # Directory of input files
# exp_set=( 'amip-hist' )                # Experiment ( could be more )
#
# 2) amip(CMIP)
# locPath='/data1/liuyincheng/CMIP6-mirror/amip/CMIP/' # Directory of input files
# exp_set=( 'amip' )                # Experiment ( could be more )
#
# 3) amip(CFMIP)
# modlist=("BCC-CSM2-MR")
# locPath='/data1/liuyincheng/CMIP6-mirror/amip/CFMIP/' # Directory of input files
# exp_set=( 'amip' )                # Experiment ( could be more )
#
# 4) ssp245
# locPath='/data1/liuyincheng/CMIP6-mirror/ssp245/' # Directory of input files
# exp_set=( 'ssp245' )                # Experiment ( could be more )
#
# 5) ssp370
# locPath='/data1/liuyincheng/CMIP6-mirror/ssp370/' # Directory of input files
# exp_set=( 'ssp370' )                # Experiment ( could be more )
#
# 6) abrupt-4xCO2
# locPath='/data1/liuyincheng/CMIP6-mirror/abrupt-4xCO2/' # Directory of input files
# exp_set=( 'abrupt-4xCO2' )                # Experiment ( could be more )
#
# 7) piControl
# locPath='/data1/liuyincheng/CMIP6-mirror/piControl/' # Directory of input files
# exp_set=( 'piControl' )
#
# 8) cliccp
# locPath='/data1/liuyincheng/CMIP6-mirror/clisccp/'
# exp_set=('amip' 'ssp245' 'ssp370' 'abrupt-4xCO2' 'piControl')