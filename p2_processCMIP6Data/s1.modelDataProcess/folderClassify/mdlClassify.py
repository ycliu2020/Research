#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File  : modelClassify.py
# Author: liuyincheng
# Date  : 2020/2/27
# '''
# 文档整理脚本
#   - 获取批量文档所在的文件夹的地址
#         - 手动查阅
#         -是全部文档的主目录，主路径path
#   - 对每一个文件名进行检索并分类
#         - 具有相同关键字的一批文�?
#   - 对有关键字文件名进行获取地址
#         - 'path{0}'.format(filename)
#         - join(path,filename)
#   - 判断文件名中是否是自己需要处理的文件名，是的话进行移动操作，即整�?
#         - 使用if语句进行判断
#         - 文件的删除、移动、复制、获取路径使用Python的os和shutil模块
# '''

# 批量将模式数据按不同的模式分类
import os
import shutil
import os.path as op

# 编写子函数

# path function
def MainPath(xx):
    main_path = data_path+xx+'/'
    return main_path

# classify function
def file_process(main_path, allIn1):
    os.chdir('{0}'.format(str(main_path)))  # 将解释器的工作路径切换到要处理的文件夹的路径
    names = os.listdir('{0}'.format(main_path))   # 获取当前目录下所有要批量处理的文件名names
    for myIn in allIn1:  # 你所要进行归类的关键字：

        myDst = main_path + myIn  # 目标路径：

        # 如果不存在目录则创建该目录
        if not op.exists(myDst):
            os.mkdir(myDst)
        for name in names:  # 遍历所有的文件夹
            if '{0}'.format(myIn)+'_' in name and '{0}'.format(mdl_rlm) in name and '{0}'.format(mdl_grid) in name:
                # 将上一级路径与文件名组合，得到文件的绝对路径，os.path.join(path,path)
                myScr = op.join(main_path, name)
                shutil.move(myScr, myDst)  # 进行文件移动 原来的路径--> 目标路径
                print('Done...')
                for root, dirs, names1 in os.walk(myScr):
                    if not os.listdir(root):
                         print('Warning! there is no file in ' + '{0}'.format(root))

# anti-classify function
def antifile_process(main_path):
    os.chdir('{0}'.format(str(main_path)))  # 将解释器的工作路径切换到要处理的文件夹的路径
    fileList = os.listdir('{0}'.format(main_path))   # 获取当前目录下所有要批量处理的文件名names
    pathList = [i for i in fileList if os.path.isdir(os.path.join(main_path, i))]  # 绝对目录
    for pp in pathList:
        for root, dirs, names in os.walk(pp):
            for name in names:
                ext = os.path.splitext(name)[1]  # 获取后缀名
                if ext == '.nc':
                    fromdir = os.path.join(root, name)  # ncl文件原始地址
                    moveto = os.path.join(os.path.dirname(root), name)  # dirname 上一层目录
                    os.rename(fromdir, moveto)  # 移动文件
            if not os.listdir(root):
                os.rmdir(root)
 
# auto get model names function
def autoGet_mdlName(main_path):
    os.chdir('{0}'.format(str(main_path)))  # 将解释器的工作路径切换到要处理的文件夹的路径
    namesAll = os.listdir('{0}'.format(main_path))   # 获取当前目录下所有要批量处理的文件名namesAll
    # 获取所有nc文件
    names = []
    for dataname in namesAll:
        if  '{0}'.format(mdl_rlm) in dataname and '{0}'.format(mdl_grid) in dataname:
            names.append(dataname)
    # 截取所有的模型名 eg:'va_Amon_CESM2-WACCM-FV2_amip_r2i1p1f1_gn_200001-201412.nc'
    mdlNames = []
    for i in names:
        temp = i.split('_')
        mdlNames.append(temp[2])
    mdlNames=list(set(mdlNames))
    return mdlNames
# auto get classifed model/ensemble member names function
def autoGet_clsfiedName(main_path):
    os.chdir('{0}'.format(str(main_path)))  # 将解释器的工作路径切换到要处理的文件夹的路径
    namesAll = os.listdir('{0}'.format(main_path))   # 获取当前目录下所有要批量处理的文件名namesAll
    # 获取所有文件夹的名称
    names = []
    for dataname in namesAll:
        if os.path.splitext(dataname)[1] == '':
            names.append(dataname)
    return names

# auto get ensemble member names function
def autoGet_esmName(main_path):
    os.chdir('{0}'.format(str(main_path)))  # 将解释器的工作路径切换到要处理的文件夹的路径
    namesAll = os.listdir('{0}'.format(main_path))   # 获取当前目录下所有要批量处理的文件名namesAll
    # 获取所有nc文件
    names = []
    for dataname in namesAll:
        if os.path.splitext(dataname)[1] == '.nc':
            names.append(dataname)
    # 截取所有的模型名 eg:'va_Amon_CESM2-WACCM-FV2_amip_r2i1p1f1_gn_200001-201412.nc'
    mdlNames = []
    for i in names:
        temp = i.split('_')
        mdlNames.append(temp[4])
    mdlNames=list(set(mdlNames))
    return mdlNames


mdl_xpt = ['amip/CMIP']
# ['abrupt-4xCO2', 'ssp370', 'ssp245',
#            'amip/CMIP', 'amip-hist', 'piControl']  # experiment

# mdl_nsm = 'r1i1p1f1'  # ensemble member
mdl_rlm = 'Amon'  # model realm
mdl_grid = 'gn'  # model realm
data_path = '/data1/liuyincheng/CMIP6-mirror/'



if __name__ == '__main__':
    
    p_1 = input("开始分类还是反分类?(Y/N)：")
    if p_1 == 'y' or p_1 == 'Y':
        for idx, my_xpt in enumerate(mdl_xpt):
            print('====================================================')
            # classific models
            file_process(main_path=MainPath(my_xpt), allIn1=autoGet_mdlName(MainPath(my_xpt)))
            print('{0}'.format(my_xpt)+' Done')
            # model List
            mdlNames=autoGet_clsfiedName(main_path=MainPath(my_xpt))
            print('there are '+'{0}'.format(mdlNames)+' models')
            # classific ensemble member
            for mdlName in mdlNames:
                mdlPth=data_path+my_xpt+'/'+mdlName+'/'
                file_process(main_path=mdlPth, allIn1=autoGet_esmName(mdlPth))
                print('{0}'.format(mdlName)+' Done')
                # ensemble member List
                esmNames=autoGet_clsfiedName(main_path=mdlPth)
                print('there are '+'{0}'.format(esmNames)+' ensemble members')


    elif p_1 == 'n' or p_1 == 'N':
        for idx, my_xpt in enumerate(mdl_xpt):
            # model List
            mdlNames=autoGet_clsfiedName(main_path=MainPath(my_xpt))
            # anti-ensemble member
            for mdlName in mdlNames:
                mdlPth=data_path+my_xpt+'/'+mdlName+'/'
                antifile_process(main_path=mdlPth)
            # anti-model list
            antifile_process(main_path=MainPath(my_xpt))
            print('anti-classify Done')


# # ################### modify this first!##################################################################
# allIn = [[] for _ in range(len(mdl_xpt))]
# # 4co2
# allIn[0] = ['AWI-CM-1-1-MR', 'BCC-CSM2-MR', 'BCC-ESM1', 'CanESM5', 'CESM2', 'CESM2-WACCM', 'GISS-E2-1-G', 'GISS-E2-1-H', 'GISS-E2-2-G', 'MIROC6', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NESM3', 'NorCPM1', 'SAM0-UNICON']

# # ssp370
# allIn[1] = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'AWI-CM-1-1-MR', 'BCC-CSM2-MR', 'BCC-ESM1', 'CanESM5', 'CESM2', 'CESM2-WACCM', 'FGOALS-g3', 'MIROC6', 'MPI-ESM-1-2-HAM', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0']

# # ssp245
# allIn[2] = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'AWI-CM-1-1-MR', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CESM2-WACCM', 'FGOALS-g3', 'MIROC6', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NESM3']

# # amip-cmip
# allIn[3] = ['CAS-ESM2-0', 'MIROC6', 'ACCESS-CM2', 'MIROC-ES2L', 'FGOALS-g3', 'FIO-ESM-2-0', 'EC-Earth3', 'IPSL-CM6A-LR', 'EC-Earth3-Veg', 'NorESM2-LM', 'CESM2-WACCM-FV2', 'BCC-CSM2-MR', 'BCC-ESM1', 'GISS-E2-2-G', 'CIESM', 'GISS-E2-1-G', 'HadGEM3-GC31-LL', 'CESM2', 'UKESM1-0-LL', 'INM-CM5-0', 'ACCESS-ESM1-5', 'HadGEM3-GC31-MM', 'CAMS-CSM1-0', 'NESM3', 'CESM2-FV2', 'MRI-ESM2-0', 'KACE-1-0-G', 'NorCPM1', 'INM-CM4-8', 'MPI-ESM1-2-HR', 'CESM2-WACCM', 'CanESM5', 'SAM0-UNICON', 'E3SM-1-0', 'TaiESM1', 'FGOALS-f3-L']
# # ['ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'BCC-ESM1', 'CanESM5', 'CESM2', 'CESM2-WACCM', 'FGOALS-g3', 'GISS-E2-1-G', 'MIROC6', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NESM3', 'NorCPM1', 'SAM0-UNICON'] 

# # amip-hist
# allIn[4] = ['BCC-CSM2-MR', 'CESM2', 'MIROC6', 'MRI-ESM2-0']

# # piControl
# allIn[5] = ['ACCESS-CM2', 'AWI-CM-1-1-MR', 'BCC-CSM2-MR', 'BCC-ESM1', 'CanESM5', 'CESM2', 'CESM2-FV2', 'CESM2-WACCM', 'CESM2-WACCM-FV2', 'FGOALS-g3', 'FIO-ESM-2-0', 'GISS-E2-1-G', 'GISS-E2-1-G-CC', 'GISS-E2-1-H', 'GISS-E2-2-G', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'MIROC6', 'MPI-ESM-1-2-HAM', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NESM3', 'NorCPM1', 'NorESM1-F', 'SAM0-UNICON']
# ##########################################################################################################
