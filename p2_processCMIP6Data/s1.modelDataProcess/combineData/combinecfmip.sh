
## combine ncfiles and save to the directory
# file name
modlist=("Amon_BCC-CSM2-MR" \
 "Amon_CAMS-CSM1-0" \
 "Amon_MIROC6" \
 "Amon_MRI-ESM2-0" \
 "Amon_CESM2")

vars=("hus" "ta" "ts" \
"rlus" "rsus" "rsds" "rlds" "rsuscs" "rsdscs" "rldscs" \
"rlut" "rsut" "rsdt" "rlutcs" "rsutcs" ) # 15 variables


loc='/data1/liuyincheng/CMIP6_mirror/amip-hist'#'/ipcc/cfmip/2xco2/atm/mo/'
loc2=''#'/ipcc/cfmip/slabcntl/atm/mo'

for i in {0..4}
do
	echo $i
	modname=${modlist[i]}
	for j in {0..14}
	do
		cd $loc
		cd ${vars[j]}
		cd $modname
		cd run1
		ncrcat *.nc ~/$modname/${vars[j]}_2co2.nc
		cd $loc2
		cd ${vars[j]}
		cd $modname
		cd run1
		ncrcat *.nc ~/$modname/${vars[j]}_control.nc
	done
done
	




