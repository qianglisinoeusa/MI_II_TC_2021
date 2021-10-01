#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# QiangLI
#
# University of Valencia
# Copyright (c) 2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np 
from scipy.io import loadmat
import dit
from pyinform.mutualinfo import mutual_info
from pyitlib import discrete_random_variable
import matplotlib.pyplot as plt
from matplotlib.dates import date2num



SharedActivity = loadmat('increaseSharedActivity.mat')

before_A = SharedActivity['sAct1_nsAct1_nsN1_A_one']
before_B = SharedActivity['sAct1_nsAct1_nsN1_B_one']
before_C = SharedActivity['sAct1_nsAct1_nsN1_C_one']
After_A = SharedActivity['sAct2_nsAct1_nsN1_A_one']
After_B = SharedActivity['sAct2_nsAct1_nsN1_B_one'] 
After_C = SharedActivity['sAct2_nsAct1_nsN1_C_one']

#Mutual Information (Shannon mutual information 1948. PyInform)
bI_ab =  discrete_random_variable.information_mutual((before_A.T), (before_B.T))
bI_bc =  discrete_random_variable.information_mutual((before_A.T), (before_C.T))
bI_ac =  discrete_random_variable.information_mutual((before_A.T), (before_C.T))
aI_ab =  discrete_random_variable.information_mutual((After_A.T),  (After_B.T))
aI_bc =  discrete_random_variable.information_mutual((After_A.T),  (After_C.T))
aI_ac =  discrete_random_variable.information_mutual((After_A.T),  (After_C.T))

#Interaction Information (Bell03,https://pafoster.github.io/pyitlib/#discrete_random_variable.information_multi)
X_before = np.concatenate((before_A, before_B, before_C), axis=1)
X_after =  np.concatenate((After_A,  After_B,  After_C),  axis=1)
SI_before=discrete_random_variable.information_co(np.array((X_before.T))) 
SI_after= discrete_random_variable.information_co(np.array((X_after.T)))

#Total Correlation (Watanable 1960,https://pafoster.github.io/pyitlib/#discrete_random_variable.information_multi)
X_before = np.concatenate((before_A, before_B, before_C), axis=1)
X_after =  np.concatenate((After_A,  After_B,  After_C),  axis=1)
TC_before=discrete_random_variable.information_multi(X_before.T)
TC_after= discrete_random_variable.information_multi(X_after.T)

MI_before = np.mean((bI_ab,bI_bc,bI_ac))
MI_after = np.mean((aI_ab,aI_bc,aI_ac))

print([MI_before,MI_after])
print([SI_before,SI_after])
print([TC_before,TC_after])





'''
NonSharedActivity = loadmat('increaseNonSharedActivity.mat')

before_A = SharedActivity['sAct1_nsAct1_nsN1_A_one']
before_B = SharedActivity['sAct1_nsAct1_nsN1_B_one']
before_C = SharedActivity['sAct1_nsAct1_nsN1_C_one']
After_A = SharedActivity['sAct2_nsAct1_nsN1_A_one']
After_B = SharedActivity['sAct2_nsAct1_nsN1_B_one'] 
After_C = SharedActivity['sAct2_nsAct1_nsN1_C_one']

SharedNosharedActivity = loadmat('increaseSharedNonSharedActivity.mat')

before_A = SharedActivity['sAct1_nsAct1_nsN1_A_one']
before_B = SharedActivity['sAct1_nsAct1_nsN1_B_one']
before_C = SharedActivity['sAct1_nsAct1_nsN1_C_one']
After_A = SharedActivity['sAct2_nsAct1_nsN1_A_one']
After_B = SharedActivity['sAct2_nsAct1_nsN1_B_one'] 
After_C = SharedActivity['sAct2_nsAct1_nsN1_C_one']


NoiseActivity = loadmat('increaseNoiseActivity.mat')

before_A = SharedActivity['sAct1_nsAct1_nsN1_A_one']
before_B = SharedActivity['sAct1_nsAct1_nsN1_B_one']
before_C = SharedActivity['sAct1_nsAct1_nsN1_C_one']
After_A = SharedActivity['sAct2_nsAct1_nsN1_A_one']
After_B = SharedActivity['sAct2_nsAct1_nsN1_B_one'] 
After_C = SharedActivity['sAct2_nsAct1_nsN1_C_one']
'''
