#!/usr/bin/env python
# coding: utf-8



import pandas as pd
import math
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from skimage import data
from skimage.transform import resize
import argparse

parser = argparse.ArgumentParser(description= 'Preprocess .vcf files')
parser.add_argument('number', type=int, help= 'Number of observations')
args = parser.parse_args()
num_ = args.number


for i in range(num_):

    df = pd.read_csv("./Data/VCF/CSV/"+"output_"+str(i+1)+".csv",header=None)
    Array_df=np.array(df)
    mode_full = stats.mode(Array_df)
    Mode_sliced=list(mode_full[0][0])
    upd_df = np.matrix(df)

    for j in range(len(Mode_sliced)):
        if Mode_sliced[j] == 1 or Mode_sliced[j] == 2:
            upd_df[:, j][upd_df[:, j] == 1] = 3
            upd_df[:, j][upd_df[:, j] == 2] = 4
            upd_df[:, j][upd_df[:, j] == 0] = 1
            upd_df[:, j][upd_df[:, j] == 3] = 0
            upd_df[:, j][upd_df[:, j] == 4] = 0

    A=np.array(upd_df)
    indexlist=np.argsort(np.linalg.norm(A,axis=1, ord=1))
    sym_sortedA = A[indexlist]

    pd.DataFrame(sym_sortedA).to_csv("./Data/VCF/CSV/"+"output_"+str(i+1)+".csv")
    print(i)
  


window=100
stride=10



# generalized
for i in range(num_):
    df=pd.read_csv("./Data/VCF/CSV/"+"output_"+str(i+1)+".csv", index_col=0)
    df=np.asmatrix(df)
    s=np.shape(df)
    sortd=np.zeros(s)
    K=[]
    for k in range(0,s[1]-window+1,stride):
        A=np.array(df[:,k:k+window])
        df[:,k:k+window]=0
        indexlist=np.argsort(np.linalg.norm(A,axis=1, ord=1))
        sortedA = A[indexlist]
        D=np.pad(sortedA,((0,0),(k,s[1]-window-k)), 'constant')
        sortd=sortd+D
        K.append(k)

    sortd=sortd+df
    u= list(range(k+stride,s[1],stride))
    K =K+u
    g=1
    for j in K:
        if j+stride<window:
            sortd[:,j:j+stride]*= (1/g)
            g=g+1
        elif window <= j+stride <=s[1]-window:
            sortd[:,j:j+stride]*= (1/g)
        elif s[1]-window<= j+stride<=s[1]:
            sortd[:,j:j+stride]*= (1/g)
            g=g-1
    if s[1]>=190:
        sortd=np.delete(sortd, [list(range(92))+list(range(s[1]-92,s[1]))], 1)
    else:
        print('NA')
    
    d=resize(sortd, (64, 64))
    pd.DataFrame(d).to_csv("./Data/VCF/CSV/c/"+"output_"+str(i+1)+".csv")
#     print(i)



