#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description= 'Preprocess .vcf files')
parser.add_argument('filename', type=str, help= 'Name of the vcf file that is placed on the VCF folder')

args = parser.parse_args()

fname = args.filename

with open(f'./Data/VCF/{fname}.vcf') as f:
    lines = f.readlines()

with open('./Data/VCF/MissSWP_test.txt') as f:
    lines1 = f.readlines()

lines1[4] = 'segsites: 0\n'

lines1[5]='positions: 0\n'

w=[]
for i in range(1, len(lines)):
    k = int(lines[i].split("\t")[1])
    w.append(k)

w = np.asarray(w)

w_=range(len(w))
w_




M=[]
P=[]
k=i=p=0
while w[k]+1100000<w[-1]:
    U=[]
    while True:
        if w[p]> w[k]+1100000:
            #print(p)
            break
        elif w[p]>=w[k] and w[p]<=w[k]+1100000:
            U.append(p)
            #print(k,p,w[k]+1100000)
            p+=1
    if len(U)>=550:
        M.append(U)
        P.append(w[k])
        print(i,len(U))
        i+=1
    m = min(j for j in w if j > w[k]+1000)
    k=p=list(np.where(w==m))[0][0]





AA=[]
for i in range(99):
    A = ""
    B = ""
    for j in range(1,len(lines)):
        t = lines[j].rstrip('\n')
        s = t.split('\t')
        a=s[i+9].split('|')[0]
        b=s[i+9].split('|')[1]
        A+=a
        B+=b
    A+='\n'
    B+='\n'
    AA.append(A)
    AA.append(B)

np.shape(P)


for i in range(len(M)):
    B=[]
    for j in range(len(AA)):
        H=list(AA[j]);
        H=np.asarray(H)
        A = H[M[i]]
        C= ''.join(A)
        C+='\n'
        B.append(C)
    with open('./Data/VCF/MS_files/output_'+str(i)+'.ms','w') as f:
        f.write(lines1[0])
        f.write(lines1[1])
        f.write(lines1[2])
        f.write(lines1[3])
        f.write(lines1[4])
        f.write(lines1[5])
        f.writelines(B)




pd.DataFrame(P).to_csv("./Data/VCF/Positions 10000 window.csv")




