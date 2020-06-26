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
# # Determine whether the directory already exists. If it doesn't, create it. If it does, output "dir exist".
# loc2=/data1/liuyincheng/CMIP6-mirror/amip-hist-combine
# echo "the dir name is $loc2"
# if [ ! -d $loc2  ];then
#     mkdir $loc2
# else
#     echo dir exist
# fi
# 2) amip(CMIP)
locAll='/data1/liuyincheng/CMIP6-mirror/' # Directory of input files
xptAll=( 'ssp245' 'ssp370')                # Experiment ( could be more )
rlm='Amon'                         # Realm
var=("hus" "ta" "ts" "ps" "hfls" "hfss"\
    "rlus" "rsus" "rsds" "rlds" "rsuscs" "rsdscs" "rldscs" \
"rlut" "rsut" "rsdt" "rlutcs" "rsutcs" "clisccp" "ua" "va" "wap") # 22 variables

ii=${#xptAll}
for modelNum in `eval echo {0..$ii}`; do
    xpt=${xptAll[modelNum]}
    loc=$locAll$xpt"/"
    jj=$((${#var[@]} - 1))
    echo $jj
    for var_id in `eval echo {0..$jj}`; do            # Loop over two variables
        # Names of all the models (ls [get file names];
        #  cut [get model names];
        #  sort; uniq [remove duplicates]; awk [print])
        mdl_set=$( ls ${loc}${var[var_id]}_${rlm}_*_${xpt}_*.nc | \
        cut -d '_' -f 3 | sort | uniq -c | awk '{print $2}' )
        # Number of models (echo [print contents]; wc [count])
        mdl_nbr=$( echo ${mdl_set} | wc -w )
        echo "=============================="
        echo "There are" ${mdl_nbr} "models for" ${var[var_id]}.
        
        for mdl in ${mdl_set}; do	        # Loop over models
            # Names of all the ensemble members
            nsm_set=$( ls ${loc}${var[var_id]}_${rlm}_${mdl}_${xpt}_*.nc | \
            cut -d '_' -f 5 | sort | uniq -c | awk '{print $2}' )
            # Number of ensemble members in each model
            nsm_nbr=$( echo ${nsm_set} | wc -w )
            echo "------------------------------"
            echo "Model" ${mdl} "includes" ${nsm_nbr} "ensemble member(s):"
            echo ${nsm_set}"."
            
            for nsm in ${nsm_set}; do	      # Loop over ensemble members
                # Number of files in this ensemble member
                fl_nbr=$( ls ${loc}${var[var_id]}_${rlm}_${mdl}_${xpt}_${nsm}_*.nc \
                | wc -w )
                
                # If there is only 1 file, continue to next loop
                if [ ${fl_nbr} -le 1 ]
                then
                    echo "There is only 1 file in" ${nsm}.
                    continue
                fi
                                
                # If there is only 0 file, mean dont exist var continue to next loop
                if [ ${fl_nbr} -le 0 ]
                then
                    echo "There is only 0 file in" ${nsm}.
                    continue
                fi
                
                echo "There are" ${fl_nbr} "files in" ${nsm}.
                
                # Starting date of data
                #   (sed [the name of the first file includes the starting date])
                yyyymm_str=$( ls ${loc}${var[var_id]}_${rlm}_${mdl}_${xpt}_${nsm}_*.nc\
                | sed -n '1p' | cut -d '_' -f 7 | cut -d '-' -f 1 )
                # Ending date of data
                #   (sed [the name of the last file includes the ending date])
                yyyymm_end=$( ls ${loc}${var[var_id]}_${rlm}_${mdl}_${xpt}_${nsm}_*.nc\
                | sed -n "${fl_nbr}p" | cut -d '_' -f 7 | cut -d '-' -f 2 )
                
                # Concatenate one ensemble member files
                #   into one along the record dimension (now is time)
                ncrcat -O ${loc}${var[var_id]}_${rlm}_${mdl}_${xpt}_${nsm}_*.nc \
                ${loc}${var[var_id]}_${rlm}_${mdl}_${xpt}_${nsm}_gn_${yyyymm_str}-${yyyymm_end}
                
                # Remove useless files
                rm ${loc}${var[var_id]}_${rlm}_${mdl}_${xpt}_${nsm}_gn_!(${yyyymm_str}-${yyyymm_end})
            done
        done
    done
done
# old scrip back up
# 1) amip-hist (not the  Primary)
# locAll='/data1/liuyincheng/CMIP6-mirror/amip-hist/' # Directory of input files
# xptAll=( 'amip-hist' )                # Experiment ( could be more )
#
# 2) amip(CMIP)
# locAll='/data1/liuyincheng/CMIP6-mirror/amip/CMIP/' # Directory of input files
# xptAll=( 'amip' )                # Experiment ( could be more )
#
# 3) amip(CFMIP)
# modlist=("BCC-CSM2-MR")
# locAll='/data1/liuyincheng/CMIP6-mirror/amip/CFMIP/' # Directory of input files
# xptAll=( 'amip' )                # Experiment ( could be more )
#
# 4) ssp245
# locAll='/data1/liuyincheng/CMIP6-mirror/ssp245/' # Directory of input files
# xptAll=( 'ssp245' )                # Experiment ( could be more )
#
# 5) ssp370
# locAll='/data1/liuyincheng/CMIP6-mirror/ssp370/' # Directory of input files
# xptAll=( 'ssp370' )                # Experiment ( could be more )
#
# 6) abrupt-4xCO2
# locAll='/data1/liuyincheng/CMIP6-mirror/abrupt-4xCO2/' # Directory of input files
# xptAll=( 'abrupt-4xCO2' )                # Experiment ( could be more )
#
# 7) piControl
# locAll='/data1/liuyincheng/CMIP6-mirror/piControl/' # Directory of input files
# xptAll=( 'piControl' )
# 
# 8) cliccp
# locAll='/data1/liuyincheng/CMIP6-mirror/clisccp/'
# xptAll=('amip' 'ssp245' 'ssp370' 'abrupt-4xCO2' 'piControl')