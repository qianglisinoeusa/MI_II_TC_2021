#!usr/bin/env python3
# -*- coding: utf-8 -*-

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# QiangLI
#
# University of Valencia
# Copyright (c) 2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import os
import sys
os.environ['THEANO_FLAGS'] = 'device=gpu'
import numpy as np
import matplotlib.pyplot as plt 
from scipy.io import savemat, loadmat
import time
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from scipy.stats import rv_histogram, norm
import pandas as pd
import seaborn as sns
plt.style.use('ggplot')

#sys.path.insert(0, 'rbig/')
#from rbig.rbig import RBIG, RBIGMI, RBIGKLD
#----------------------------------------------------------------------------------
# Option2:
sys.path.insert(0, '../rbig/')
from rbig.rbig import rbig

def plot_2d_joint(data, layer):
    fig  = plt.figure(figsize=(12, 5))
    g = sns.jointplot(x=data[0], y=data[1], color="#4CB391", alpha=1) #kind='kde'
    g.ax_joint.set_xticks([])
    g.ax_joint.set_yticks([])
    plt.tight_layout()
    if savename:
        g.savefig(f"rbig_"+str(layer)+"_data.png", transparent=False)
    #plt.show()
    return None


num_samples = 10000
x = np.abs(-5*np.random.randn(1, num_samples))
y = np.cos(x) + 0.45*np.random.randn(1, num_samples)
data = np.vstack((x, y))

#fig  = plt.figure(figsize=(12, 5))
#g = sns.jointplot(x=data[0], y=data[1], kind='scatter', color="#4CB391", alpha=0.1)
#g.ax_joint.set_xticks([])
#g.ax_joint.set_yticks([])
#plt.tight_layout()
#plt.savefig('rotation_a' + '.png')

#fig = plt.figure(figsize=(6, 6))
#plt.scatter(data[0], data[1], s=1)
#plt.xlabel('X')
#plt.xlabel('Y')
#plt.title('Original Data')
#plt.savefig('rotation' + '.png')
savename=True
plot_2d_joint(data,layer=00)

# Obtain the RBIG transform for the data using *PCA* rotation and 50 iterations
for layer in range(60):
    rbig_transformed_data, trans_params = rbig(data, layer, 'PCA')
    plot_2d_joint(rbig_transformed_data, layer)

    #fig  = plt.figure(figsize=(12, 5))
    #g = sns.jointplot(x=rbig_transformed_data[0], y=rbig_transformed_data[1], kind='scatter', color="#4CB391", alpha=0.1)
    #g.ax_joint.set_xticks([])
    #g.ax_joint.set_yticks([])
    #plt.tight_layout()
    #plt.savefig('rotation_a'+str(layer)+'.png')
    
    #fig = plt.figure(figsize=(6, 6))
    #plt.scatter(rbig_transformed_data[0], rbig_transformed_data[1], s=1)
    #plt.xlabel('X')
    #plt.xlabel('Y')
    #plt.axis('off')
    #plt.title('layer:%i' %layer)
    #plt.savefig('rotation'+str(layer)+'.png')
    #plt.show()
    
    
