'''
@Author       : LYC
@Date         : 2020-06-23 11:17:58
@LastEditTime : 2020-06-23 20:52:50
@LastEditors  : LYC
@Description  : find all model mdlNames in a specical path
@FilePath     : /code/p2_processCMIP6Data/s1.modelDataProcess/folderClassify/modelList.py
@ 
'''
import os
import shutil
import os.path as op
mdl_vars = ['hus', 'ta', 'ts', 'ps',\
            'rlus', 'rsus', 'rsds', 'rlds', 'rsuscs', 'rsdscs', 'rldscs',\
            'rlut', 'rsut', 'rsdt', 'rlutcs', 'rsutcs']
mdl_rlm = 'Amon'  # model realm
mdl_grid = 'gn'  # model realm

locPath = '/data1/liuyincheng/CMIP6-mirror/amip/CMIP/'
os.chdir('{0}'.format(str(locPath)))  # 将解释器的工作路径切换到要处理的文件夹的路径
names = os.listdir('{0}'.format(locPath))   # 获取当前目录下所有要批量处理的文件名namesAll
# 获取所有文件夹名字(模式名)
mdlNames = []
for name in names:
    if os.path.splitext(name)[1] == '':
        mdlNames.append(name)

# 进入路径读取所有ensemble member
for mdlName in mdlNames:
    print('Model: '+'{0}'.format(mdlName))

    mdlPth = locPath + mdlName+'/'
    os.chdir('{0}'.format(str(mdlPth)))  # 将解释器的工作路径切换到要处理的文件夹的路径
    names = os.listdir('{0}'.format(mdlPth))   # 获取当前目录下所有要批量处理的文件名namesAll
    # 获取所有ensemble member名
    esmNames = []
    for name in names:
        if os.path.splitext(name)[1] == '':
            esmNames.append(name)
    print('ensemble member: '+'{0}'.format(esmNames))

    # 统计变量个数
    for esmName in esmNames:
        esmPth = mdlPth+esmName+'/'
        os.chdir('{0}'.format(str(esmPth)))  # 将解释器的工作路径切换到要处理的文件夹的路径
        # 获取当前目录下所有要批量处理的文件名names
        names = os.listdir('{0}'.format(esmPth))
        print('{0}'.format(esmName)+' have ' + '{0}'.format(len(names))+' vars')
        varsNames=[]
        for name in names:
            temp = name.split('_')
            varsNames.append(temp[0])
        for mdl_var in mdl_vars:
            if '{0}'.format(mdl_var) not in varsNames:
                print('!!!!!!!!!!!!!! '+'{0}'.format(esmName)+ ' lack '+'{0}'.format(mdl_var))
                

            # # 截取所有的模型名 eg:'va_Amon_CESM2-WACCM-FV2_amip_r2i1p1f1_gn_200001-201412.nc'
            # mdlNames = []
            # for i in mdlNames:
            #     temp = i.split('_')
            #     mdlNames.append(temp[2])
            # mdlNames=list(set(mdlNames))
            # print(mdlNames)
