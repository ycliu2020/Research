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
            if '{0}'.format(myIn)+'_' in name and '{0}'.format(mdl_nsm) in name and '{0}'.format(mdl_rlm) in name:
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
 

mdl_xpt = ['abrupt-4xCO2', 'ssp370', 'ssp245',
           'amip/CMIP', 'amip-hist', 'piControl']  # experiment

mdl_nsm = 'r1i1p1f1'  # ensemble member
mdl_rlm = 'Amon'  # model realm
data_path = '/data1/liuyincheng/CMIP6-mirror/'

# ################### modify this first!##################################################################
allIn = [[] for _ in range(len(mdl_xpt))]
# 4co2
allIn[0] = ['AWI-CM-1-1-MR', 'BCC-CSM2-MR', 'BCC-ESM1', 'CanESM5', 'CESM2', 'CESM2-WACCM', 'GISS-E2-1-G', 'GISS-E2-1-H', 'GISS-E2-2-G', 'MIROC6', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NESM3', 'NorCPM1', 'SAM0-UNICON']

# ssp370
allIn[1] = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'AWI-CM-1-1-MR', 'BCC-CSM2-MR', 'BCC-ESM1', 'CanESM5', 'CESM2', 'CESM2-WACCM', 'FGOALS-g3', 'MIROC6', 'MPI-ESM-1-2-HAM', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0']

# ssp245
allIn[2] = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'AWI-CM-1-1-MR', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CESM2-WACCM', 'FGOALS-g3', 'MIROC6', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NESM3']

# amip-cmip
allIn[3] = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'BCC-ESM1', 'CanESM5', 'CESM2', 'CESM2-WACCM', 'FGOALS-g3', 'GISS-E2-1-G', 'MIROC6', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NESM3', 'NorCPM1', 'SAM0-UNICON'] 

# amip-hist
allIn[4] = ['BCC-CSM2-MR', 'CESM2', 'MIROC6', 'MRI-ESM2-0']

# piControl
allIn[5] = ['ACCESS-CM2', 'AWI-CM-1-1-MR', 'BCC-CSM2-MR', 'BCC-ESM1', 'CanESM5', 'CESM2', 'CESM2-FV2', 'CESM2-WACCM', 'CESM2-WACCM-FV2', 'FGOALS-g3', 'FIO-ESM-2-0', 'GISS-E2-1-G', 'GISS-E2-1-G-CC', 'GISS-E2-1-H', 'GISS-E2-2-G', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'MIROC6', 'MPI-ESM-1-2-HAM', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NESM3', 'NorCPM1', 'NorESM1-F', 'SAM0-UNICON']
##########################################################################################################


if __name__ == '__main__':
    p_1 = input("开始分类还是反分类?(Y/N)：")
    if p_1 == 'y' or p_1 == 'Y':
        for idx, my_xpt in enumerate(mdl_xpt):
            print('====================================================')
            file_process(main_path=MainPath(my_xpt), allIn1=allIn[idx])
            print('{0}'.format(my_xpt)+' Done')
    elif p_1 == 'n' or p_1 == 'N':
        for idx, my_xpt in enumerate(mdl_xpt):
            antifile_process(main_path=MainPath(my_xpt))
            print('anti-classify Done')


