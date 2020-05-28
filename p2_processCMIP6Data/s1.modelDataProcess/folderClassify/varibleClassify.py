# 原作者说�?
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File  : 文件整理.py
# Author: roohom
# Date  : 2018/8/14
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

# 批量将每个模式中的数据按变量名分

import os
import shutil
import os.path as op

# 获取批量文件处理的总路�?
modelnames = ['BCC-CSM2-MR', 'CAMS-CSM1-0',
              'MIROC6', 'MRI-ESM2-0', 'CESM2']
for modelname in modelnames:

    def MainPath():
        main_path = '/data1/liuyincheng/CMIP6_mirror/amip-hist/{0}'.format(modelname)# 请手动查阅文件夹的路径：
        return main_path

    def file_process(main_path):
        os.chdir('{0}'.format(str(main_path)))# 将解释器的工作路径切换到要处理的文件夹的路径
        names = os.listdir('{0}'.format(main_path))# 获取当前目录下所有要批量处理的文件名names
        allIn = ["hus_", "ta_", "ts_",
                 "rlus_", "rsus_", "rsds_", "rlds_", "rsuscs_", "rsdscs_", "rldscs_",
                 "rlut_", "rsut_", "rsdt_", "rlutcs_", "rsutcs_"]
        for myIn in allIn:  # 你所要进行归类的关键字：
            # myIn = 'Amon_BCC'
            temp=myIn.split('_')
            myDst = '{0}/{1}'.format(main_path,temp[0]) # 请输入你所要放置的目标文件夹路径：
            # 如果不存在目录则创建该目�?
            if not op.exists(myDst):
                os.mkdir(myDst)
            for name in names:  # 遍历所有的文件�?
                if '{0}'.format(myIn) in name:
                    myScr = op.join(main_path, name)# 将上一级路径与文件名组合，得到文件的绝对路径，os.path.join(path,path)
                    shutil.move(myScr, myDst)# 进行文件移动 原来的路�?--> 目标路径
                    print('Done...')

    if __name__ == '__main__':
        file_process(main_path=MainPath())
