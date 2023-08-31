# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import scanpy as sc
import os
import rpy2
#from rpy2.robjects import r
import rpy2.robjects as ro
r=ro.r
from io import BytesIO
import PIL.Image as Image
from rpy2 import robjects
import matplotlib.pyplot as pl
from matplotlib import rcParams
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt




r.source('ATAC/ATAC_all.R')
'''
r.read_1()
r.TSS_2()
r.QCVlnPlot_3()
r.Dimpro_4()
r.select_table_5()
r.show_select_6()
r.show_cluster_7()
'''

path=''

def up_dir(path):
    updir = os.path.abspath(os.path.join(path, "../"))
    updir = updir.replace('\\', '/')
    return updir

def show_figure(img_path):
    plt.figure(figsize=(7, 6))
    img = plt.imread(img_path)
    plt.imshow(img)
    plt.xticks([])  # 去掉x轴
    plt.yticks([])  # 去掉y轴
    plt.axis('off')  # 去掉坐标轴
    plt.show()

def tran_1(path,min_cells,min_features):
    upp = up_dir(path)
    print(path)
    print(upp)
    r.read_1(path,up_dir(path),min_cells,min_features)

def tran_2(path,TSS_enrichment):
    upp = up_dir(path)
    r.TSS_2(path,upp,TSS_enrichment)

    img_path = upp+'/figure/2_TSSPlot.png'
    show_figure(img_path)

def tran_3(path):
    upp = up_dir(path)
    r.QCVlnPlot_3(path,upp)

    img_path = upp + '/figure/3_QCVlnPlot.png'
    show_figure(img_path)


def tran_4(path,k_param,resolution,dims):
    upp = up_dir(path,)
    r.Dimpro_4(path, upp,k_param,resolution,dims)

    p_list=['4_DimPlot','4_DimLoad']
    for i in range(len(p_list)):
        img_path = upp + '/figure/'+p_list[i]+'.png'
        show_figure(img_path)

def tran_5(path):
    upp = up_dir(path)
    top10_list=list(r.select_table_5(path))
    top10 = '\n'.join(top10_list)
    return top10,top10_list



def tran_6(path,gene):
    global upp
    upp = up_dir(path)
    r.show_select_6(path, upp,gene)
    print(gene)

    p_list = ['6_VlnPlot', '6_FeaturePlot']
    for i in range(len(p_list)):
        img_path = upp + '/figure/'+p_list[i]+'.png'
        show_figure(img_path)

def tran_7(path):
    global upp
    upp = up_dir(path)
    r.show_cluster_7(path, upp)

    img_path = upp + '/figure/7_cluster.png'
    show_figure(img_path)


#rpy2.robjects.r('plot(t)')

#t2=rpy2.robjects.r('plot(t)')
#print(t2)


