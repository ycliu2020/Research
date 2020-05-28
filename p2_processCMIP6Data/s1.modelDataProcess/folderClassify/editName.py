#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File  : editName.py
# Author: liuyincheng
# Date  : 2020/2/27
#  used to auto edit model name 
#
# picontrol
# name='ACCESS-CM2 AWI-CM-1-1-MR BCC-CSM2-MR BCC-ESM1 CanESM5 CESM2 CESM2-FV2 CESM2-WACCM CESM2-WACCM-FV2 FGOALS-g3 FIO-ESM-2-0 GISS-E2-1-G GISS-E2-1-G-CC GISS-E2-1-H GISS-E2-2-G HadGEM3-GC31-LL HadGEM3-GC31-MM MIROC6 MPI-ESM-1-2-HAM MPI-ESM1-2-HR MRI-ESM2-0 NESM3 NorCPM1 NorESM1-F SAM0-UNICON'

# amip
# name='ACCESS-CM2 ACCESS-ESM1-5 BCC-CSM2-MR BCC-ESM1 CanESM5 CESM2 CESM2-WACCM FGOALS-g3 GISS-E2-1-G MIROC6 MPI-ESM1-2-HR MRI-ESM2-0 NESM3 NorCPM1 SAM0-UNICON'

# ssp245
# name='ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CESM2 CESM2-WACCM FGOALS-g3 MIROC6 MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3'

# ssp370
# name='ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR BCC-ESM1 CanESM5 CESM2 CESM2-WACCM FGOALS-g3 MIROC6 MPI-ESM-1-2-HAM MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0'

# adrupt-4co2
name='AWI-CM-1-1-MR BCC-CSM2-MR BCC-ESM1 CanESM5 CESM2 CESM2-WACCM GISS-E2-1-G GISS-E2-1-H GISS-E2-2-G MIROC6 MPI-ESM1-2-HR MRI-ESM2-0 NESM3 NorCPM1 SAM0-UNICON'

add_term=''
temp = name.split(' ') 
print (temp)
print ('There are '+str(len(temp))+' models')

for i,item in enumerate(temp):
    temp[i] = temp[i]+add_term
    newtemp = temp
print (temp)

