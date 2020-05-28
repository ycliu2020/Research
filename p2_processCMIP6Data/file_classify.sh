#!/bin/bash

ProjPath=`cd $(dirname "$0");pwd`

dirPath="$ProjPath/cultivar_process"
#�Զ�������ļ�����
FileName=$1

mkdir $FileName"_t"
mkdir $FileName"_t_b"
mkdir $FileName
mkdir $FileName"_b"

t_filePath="$ProjPath/$FileName"_t""
t_b_filePath="$ProjPath/$FileName"_t_b""
filePath="$ProjPath/$FileName"
b_filePath="$ProjPath/$FileName"_b""

#���д��Ҷ�ӵ��ļ���
dirList=`ls $dirPath`

#��������ҶƬ���ļ���
for DirName in $dirList
do
    cd $dirPath/$DirName
    file_1_1_t_list=`ls *_[1-3]_[1-3]_t.tif`
    for file_1_1 in $file_1_1_t_list
    do
        cd $t_filePath

        if [ ! -d "$DirName" ]
        then
            mkdir $DirName
        fi

        cp $dirPath/$DirName/$file_1_1 $t_filePath/$DirName 
    done

    cd $dirPath/$DirName
    file_1_1_t_b_list=`ls *_[1-3]_[1-3]_t_b.tif`
    for file_2_2 in $file_1_1_t_b_list
    do
        cd $t_b_filePath

        if [ ! -d "$DirName" ]
        then
            mkdir $DirName
        fi

        cp $dirPath/$DirName/$file_2_2 $t_b_filePath/$DirName
    done

    cd $dirPath/$DirName
    file_1_1_list=`ls *_[1-3]_[1-3].tif`
    for file_3_3 in $file_1_1_list
    do
        cd $filePath

        if [ ! -d "$DirName" ]
        then
            mkdir $DirName
        fi

        cp $dirPath/$DirName/$file_3_3 $filePath/$DirName
    done

    cd $dirPath/$DirName
    file_1_1_b_list=`ls *_[1-3]_[1-3]_b.tif`
    for file_4_4 in $file_1_1_b_list
    do
        cd $b_filePath

        if [ ! -d "$DirName" ]
        then
            mkdir $DirName
        fi

        cp $dirPath/$DirName/$file_4_4 $b_filePath/$DirName
    done
done
