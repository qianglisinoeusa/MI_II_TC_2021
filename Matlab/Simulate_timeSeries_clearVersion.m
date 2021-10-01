%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QiangLI
%
% University of Valencia
% Copyright (c) 2021
%
% Version 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
addpath(genpath('InteractionInformation'));
addpath(genpath('TC_code_RBIG_ITE'));
addpath(genpath('2017_RBIG'));
addpath(genpath('gcmi'));

numTrials=2000;
numSubjs=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                create sources of noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsNoiseA=.45*randn(10000,numTrials,numSubjs);
nsNoiseB=.45*randn(10000,numTrials,numSubjs);
nsNoiseC=.45*randn(10000,numTrials,numSubjs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                create sources of non-shared activity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsActivityA=randn(10000,numTrials,numSubjs);
nsActivityB=randn(10000,numTrials,numSubjs);
nsActivityC=randn(10000,numTrials,numSubjs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                create source of shared activity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sharedActivity=randn(10000,numTrials,numSubjs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                create time series in the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sAct1_nsAct1_nsN1_A=sharedActivity+nsActivityA+nsNoiseA;
sAct1_nsAct1_nsN1_B=sharedActivity+nsActivityB+nsNoiseB;
expo=2.5;
sAct1_nsAct1_nsN1_B = abs(sAct1_nsAct1_nsN1_B).^expo;
sAct1_nsAct1_nsN1_C=sharedActivity+nsActivityC+nsNoiseC;
expo=2.5;
sAct1_nsAct1_nsN1_C = abs(sAct1_nsAct1_nsN1_C).^expo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 1: Increased shared activity amplitude (0.8x) in both regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Increased shared activity amplitude (2x)')

sAct1_nsAct1_nsN1_A_one = sAct1_nsAct1_nsN1_A(:,1:60,1);
sAct1_nsAct1_nsN1_B_one = sAct1_nsAct1_nsN1_B(:,1:60,1);
expo=2.5;
sAct1_nsAct1_nsN1_B = abs(sAct1_nsAct1_nsN1_B).^expo;
sAct1_nsAct1_nsN1_C_one = sAct1_nsAct1_nsN1_C(:,1:60,1);
expo=2.5;
sAct1_nsAct1_nsN1_C = abs(sAct1_nsAct1_nsN1_C).^expo;


sAct2_nsAct1_nsN1_A=(0.8*sharedActivity)+nsActivityA+nsNoiseA;
sAct2_nsAct1_nsN1_B=(0.8*sharedActivity)+nsActivityB+nsNoiseB;
expo=2.5;
sAct2_nsAct1_nsN1_B = abs(sAct2_nsAct1_nsN1_B).^expo;
sAct2_nsAct1_nsN1_C=(0.8*sharedActivity)+nsActivityC+nsNoiseC;
expo=2.5;
sAct2_nsAct1_nsN1_C = abs(sAct2_nsAct1_nsN1_C).^expo;

sAct2_nsAct1_nsN1_A_one = sAct2_nsAct1_nsN1_A(:,1:60,1);
sAct2_nsAct1_nsN1_B_one = sAct2_nsAct1_nsN1_B(:,1:60,1);
sAct2_nsAct1_nsN1_C_one = sAct2_nsAct1_nsN1_C(:,1:60,1);
% save('increaseSharedActivity.mat', 'sAct1_nsAct1_nsN1_A_one', 'sAct1_nsAct1_nsN1_B_one', 'sAct1_nsAct1_nsN1_C_one', ...
%                                'sAct2_nsAct1_nsN1_A_one', 'sAct2_nsAct1_nsN1_B_one', 'sAct2_nsAct1_nsN1_C_one')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MI
PARAMS.lay=1000;
I_before_AB = RBIG_MSMI(sAct1_nsAct1_nsN1_A_one',sAct1_nsAct1_nsN1_B_one',PARAMS);
I_before_AC = RBIG_MSMI(sAct1_nsAct1_nsN1_A_one',sAct1_nsAct1_nsN1_C_one',PARAMS);
I_before_BC = RBIG_MSMI(sAct1_nsAct1_nsN1_B_one',sAct1_nsAct1_nsN1_C_one',PARAMS);
I_before = mean([I_before_AB,I_before_AC,I_before_BC]);

I_after_AB = RBIG_MSMI(sAct2_nsAct1_nsN1_A_one',sAct2_nsAct1_nsN1_B_one',PARAMS);
I_after_AC = RBIG_MSMI(sAct2_nsAct1_nsN1_A_one',sAct2_nsAct1_nsN1_C_one',PARAMS);
I_after_BC = RBIG_MSMI(sAct2_nsAct1_nsN1_B_one',sAct2_nsAct1_nsN1_C_one',PARAMS);
I_after = mean([I_after_AB,I_after_AC,I_after_BC]);

%TC
data_before = [sAct1_nsAct1_nsN1_A_one; sAct1_nsAct1_nsN1_B_one; sAct1_nsAct1_nsN1_C_one];
T_before = RBIG_TC(data_before',PARAMS);

data_after = [sAct2_nsAct1_nsN1_A_one; sAct2_nsAct1_nsN1_B_one; sAct2_nsAct1_nsN1_C_one];
T_after = RBIG_TC(data_after',PARAMS);

%II
data_before_AB = cmi_ggg(sAct1_nsAct1_nsN1_A_one,sAct1_nsAct1_nsN1_B_one,sAct1_nsAct1_nsN1_C_one);
data_before_AC = cmi_ggg(sAct1_nsAct1_nsN1_A_one,sAct1_nsAct1_nsN1_C_one,sAct1_nsAct1_nsN1_B_one);
data_before_BC = cmi_ggg(sAct1_nsAct1_nsN1_B_one,sAct1_nsAct1_nsN1_C_one,sAct1_nsAct1_nsN1_A_one);
data_before_SI = T_before - (data_before_AB+data_before_AC+data_before_BC);

data_after_AB = cmi_ggg(sAct2_nsAct1_nsN1_A_one,sAct2_nsAct1_nsN1_B_one, sAct2_nsAct1_nsN1_C_one);
data_after_AC = cmi_ggg(sAct2_nsAct1_nsN1_A_one,sAct2_nsAct1_nsN1_C_one, sAct2_nsAct1_nsN1_B_one);
data_after_BC = cmi_ggg(sAct2_nsAct1_nsN1_B_one,sAct2_nsAct1_nsN1_C_one, sAct2_nsAct1_nsN1_A_one);
data_after_SI = T_after - (data_after_AB+data_after_AC+data_after_BC);


[I_before_AB,I_before_AC,I_before_BC, data_before_SI, T_before]
[I_after_AB, I_after_AC,  I_after_BC, data_after_SI,  T_after]

[abs(I_after_AB-I_before_AB)/I_before_AB,  abs(I_after_AC-I_before_AC)/I_before_AC, abs(I_after_BC-I_before_BC)/I_before_BC, ...
    abs(data_after_SI-data_before_SI)/abs(data_before_SI), abs(T_after-T_before)/T_before]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 2: Increased non-shared activity amplitude (0.8x) in both regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Increased non-shared activity amplitude (2x) in both regions')
sAct1_nsAct1_nsN1_A_one = sAct1_nsAct1_nsN1_A(:,1:60,1);
sAct1_nsAct1_nsN1_B_one = sAct1_nsAct1_nsN1_B(:,1:60,1);
sAct1_nsAct1_nsN1_C_one = sAct1_nsAct1_nsN1_C(:,1:60,1);
 
sAct2_nsAct1_nsN1_A=sharedActivity+(1.5*nsActivityA)+nsNoiseA;
sAct2_nsAct1_nsN1_B=sharedActivity+(1.5*nsActivityB)+nsNoiseB;
expo=2.5;
sAct2_nsAct1_nsN1_B = abs(sAct2_nsAct1_nsN1_B).^expo;
sAct2_nsAct1_nsN1_C=sharedActivity+(1.5*nsActivityC)+nsNoiseC;
expo=2.5;
sAct2_nsAct1_nsN1_C = abs(sAct2_nsAct1_nsN1_C).^expo;

sAct2_nsAct1_nsN1_A_one = sAct2_nsAct1_nsN1_A(:,1:60,1);
sAct2_nsAct1_nsN1_B_one = sAct2_nsAct1_nsN1_B(:,1:60,1);
sAct2_nsAct1_nsN1_C_one = sAct2_nsAct1_nsN1_C(:,1:60,1);
% save('increaseNonSharedActivity.mat', 'sAct1_nsAct1_nsN1_A_one', 'sAct1_nsAct1_nsN1_B_one', 'sAct1_nsAct1_nsN1_C_one', ...
%                                'sAct2_nsAct1_nsN1_A_one', 'sAct2_nsAct1_nsN1_B_one', 'sAct2_nsAct1_nsN1_C_one')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MI
I_before_AB = RBIG_MSMI(sAct1_nsAct1_nsN1_A_one',sAct1_nsAct1_nsN1_B_one');
I_before_AC = RBIG_MSMI(sAct1_nsAct1_nsN1_A_one',sAct1_nsAct1_nsN1_C_one');
I_before_BC = RBIG_MSMI(sAct1_nsAct1_nsN1_B_one',sAct1_nsAct1_nsN1_C_one');
I_before = mean([I_before_AB,I_before_AC,I_before_BC]);

I_after_AB = RBIG_MSMI(sAct2_nsAct1_nsN1_A_one',sAct2_nsAct1_nsN1_B_one');
I_after_AC = RBIG_MSMI(sAct2_nsAct1_nsN1_A_one',sAct2_nsAct1_nsN1_C_one');
I_after_BC = RBIG_MSMI(sAct2_nsAct1_nsN1_B_one',sAct2_nsAct1_nsN1_C_one');
I_after = mean([I_after_AB,I_after_AC,I_after_BC]);

%TC
data_before = [sAct1_nsAct1_nsN1_A_one; sAct1_nsAct1_nsN1_B_one; sAct1_nsAct1_nsN1_C_one];
T_before = RBIG_TC(data_before');

data_after = [sAct2_nsAct1_nsN1_A_one; sAct2_nsAct1_nsN1_B_one; sAct2_nsAct1_nsN1_C_one];
T_after = RBIG_TC(data_after');

%SI
data_before = [sAct1_nsAct1_nsN1_A_one(:), sAct1_nsAct1_nsN1_B_one(:), sAct1_nsAct1_nsN1_C_one(:)];
data_before_AB = cmi_ggg(sAct1_nsAct1_nsN1_A_one,sAct1_nsAct1_nsN1_B_one,sAct1_nsAct1_nsN1_C_one);
data_before_AC = cmi_ggg(sAct1_nsAct1_nsN1_A_one,sAct1_nsAct1_nsN1_C_one,sAct1_nsAct1_nsN1_B_one);
data_before_BC = cmi_ggg(sAct1_nsAct1_nsN1_B_one,sAct1_nsAct1_nsN1_C_one,sAct1_nsAct1_nsN1_A_one);
data_before_SI = T_before - (data_before_AB+data_before_AC+data_before_BC);


data_after_AB = cmi_ggg(sAct2_nsAct1_nsN1_A_one,sAct2_nsAct1_nsN1_B_one, sAct2_nsAct1_nsN1_C_one);
data_after_AC = cmi_ggg(sAct2_nsAct1_nsN1_A_one,sAct2_nsAct1_nsN1_C_one, sAct2_nsAct1_nsN1_B_one);
data_after_BC = cmi_ggg(sAct2_nsAct1_nsN1_B_one,sAct2_nsAct1_nsN1_C_one, sAct2_nsAct1_nsN1_A_one);
data_after_SI = T_after - (data_after_AB+data_after_AC+data_after_BC);


[I_after_AB, I_after_AC,  I_after_BC, data_after_SI,  T_after]

[abs(I_after_AB-I_before_AB)/I_before_AB,  abs(I_after_AC-I_before_AC)/I_before_AC, abs(I_after_BC-I_before_BC)/I_before_BC, abs(data_after_SI-data_before_SI)/abs(data_before_SI), abs(T_after-T_before)/T_before]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 3: Increased shared (2x) & non-shared (2x) activity amplitude in both regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Increased shared & non-shared activity amplitude in both regions')
sAct1_nsAct1_nsN1_A_one = sAct1_nsAct1_nsN1_A(:,1:60,1);
sAct1_nsAct1_nsN1_B_one = sAct1_nsAct1_nsN1_B(:,1:60,1);
sAct1_nsAct1_nsN1_C_one = sAct1_nsAct1_nsN1_C(:,1:60,1);
 
sAct2_nsAct1_nsN1_A=(0.8*sharedActivity)+(0.8*nsActivityA)+nsNoiseA;
sAct2_nsAct1_nsN1_B=(0.8*sharedActivity)+(0.8*nsActivityB)+nsNoiseB;
expo=2.5;
sAct2_nsAct1_nsN1_B = abs(sAct2_nsAct1_nsN1_B).^expo;
sAct2_nsAct1_nsN1_C=(0.8*sharedActivity)+(0.8*nsActivityC)+nsNoiseC;
expo=2.5;
sAct2_nsAct1_nsN1_C = abs(sAct2_nsAct1_nsN1_C).^expo;

sAct2_nsAct1_nsN1_A_one = sAct2_nsAct1_nsN1_A(:,1:60,1);
sAct2_nsAct1_nsN1_B_one = sAct2_nsAct1_nsN1_B(:,1:60,1);
sAct2_nsAct1_nsN1_C_one = sAct2_nsAct1_nsN1_C(:,1:60,1);
% save('increaseSharedNonSharedActivity.mat', 'sAct1_nsAct1_nsN1_A_one', 'sAct1_nsAct1_nsN1_B_one', 'sAct1_nsAct1_nsN1_C_one', ...
%                                'sAct2_nsAct1_nsN1_A_one', 'sAct2_nsAct1_nsN1_B_one', 'sAct2_nsAct1_nsN1_C_one')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MI
I_before_AB = RBIG_MSMI(sAct1_nsAct1_nsN1_A_one',sAct1_nsAct1_nsN1_B_one');
I_before_AC = RBIG_MSMI(sAct1_nsAct1_nsN1_A_one',sAct1_nsAct1_nsN1_C_one');
I_before_BC = RBIG_MSMI(sAct1_nsAct1_nsN1_B_one',sAct1_nsAct1_nsN1_C_one');
I_before = mean([I_before_AB,I_before_AC,I_before_BC]);

I_after_AB = RBIG_MSMI(sAct2_nsAct1_nsN1_A_one',sAct2_nsAct1_nsN1_B_one');
I_after_AC = RBIG_MSMI(sAct2_nsAct1_nsN1_A_one',sAct2_nsAct1_nsN1_C_one');
I_after_BC = RBIG_MSMI(sAct2_nsAct1_nsN1_B_one',sAct2_nsAct1_nsN1_C_one');
I_after = mean([I_after_AB,I_after_AC,I_after_BC]);


%TC
data_before = [sAct1_nsAct1_nsN1_A_one; sAct1_nsAct1_nsN1_B_one; sAct1_nsAct1_nsN1_C_one;];
T_before = RBIG_TC(data_before');

data_after = [sAct2_nsAct1_nsN1_A_one; sAct2_nsAct1_nsN1_B_one; sAct2_nsAct1_nsN1_C_one;];
T_after = RBIG_TC(data_after');


%SI
data_before = [sAct1_nsAct1_nsN1_A_one(:), sAct1_nsAct1_nsN1_B_one(:), sAct1_nsAct1_nsN1_C_one(:)];
data_before_AB = cmi_ggg(sAct1_nsAct1_nsN1_A_one,sAct1_nsAct1_nsN1_B_one,sAct1_nsAct1_nsN1_C_one);
data_before_AC = cmi_ggg(sAct1_nsAct1_nsN1_A_one,sAct1_nsAct1_nsN1_C_one,sAct1_nsAct1_nsN1_B_one);
data_before_BC = cmi_ggg(sAct1_nsAct1_nsN1_B_one,sAct1_nsAct1_nsN1_C_one,sAct1_nsAct1_nsN1_A_one);
data_before_SI = T_before - (data_before_AB+data_before_AC+data_before_BC);

data_after_AB = cmi_ggg(sAct2_nsAct1_nsN1_A_one,sAct2_nsAct1_nsN1_B_one, sAct2_nsAct1_nsN1_C_one);
data_after_AC = cmi_ggg(sAct2_nsAct1_nsN1_A_one,sAct2_nsAct1_nsN1_C_one, sAct2_nsAct1_nsN1_B_one);
data_after_AB = cmi_ggg(sAct2_nsAct1_nsN1_B_one,sAct2_nsAct1_nsN1_C_one, sAct2_nsAct1_nsN1_A_one);
data_after_SI = T_after - (data_after_AB+data_after_AC+data_after_BC);


[I_after_AB, I_after_AC,  I_after_BC, data_after_SI,  T_after]

[abs(I_after_AB-I_before_AB)/I_before_AB,  abs(I_after_AC-I_before_AC)/I_before_AC, abs(I_after_BC-I_before_BC)/I_before_BC, abs(data_after_SI-data_before_SI)/abs(data_before_SI), abs(T_after-T_before)/T_before]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 4: Increased nosie activity amplitude in both regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Increased nosie activity amplitude')

sAct1_nsAct1_nsN1_A_one = sAct1_nsAct1_nsN1_A(:,1:60,1);
sAct1_nsAct1_nsN1_B_one = sAct1_nsAct1_nsN1_B(:,1:60,1);
sAct1_nsAct1_nsN1_C_one = sAct1_nsAct1_nsN1_C(:,1:60,1);
 
sAct2_nsAct1_nsN1_A=sharedActivity+nsActivityA+(1.8*nsNoiseA);
sAct2_nsAct1_nsN1_B=sharedActivity+nsActivityB+(1.8*nsNoiseB);
expo=2.5;
sAct2_nsAct1_nsN1_B = abs(sAct2_nsAct1_nsN1_B).^expo;
sAct2_nsAct1_nsN1_C=sharedActivity+nsActivityC+(1.8*nsNoiseC);
expo=2.5;
sAct2_nsAct1_nsN1_C = abs(sAct2_nsAct1_nsN1_C).^expo;

sAct2_nsAct1_nsN1_A_one = sAct2_nsAct1_nsN1_A(:,1:60,1);
sAct2_nsAct1_nsN1_B_one = sAct2_nsAct1_nsN1_B(:,1:60,1);
sAct2_nsAct1_nsN1_C_one = sAct2_nsAct1_nsN1_C(:,1:60,1);
% save('increaseNoiseActivity.mat', 'sAct1_nsAct1_nsN1_A_one', 'sAct1_nsAct1_nsN1_B_one', 'sAct1_nsAct1_nsN1_C_one', ...
%                                   'sAct2_nsAct1_nsN1_A_one', 'sAct2_nsAct1_nsN1_B_one', 'sAct2_nsAct1_nsN1_C_one')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MI
I_before_AB = RBIG_MSMI(sAct1_nsAct1_nsN1_A_one',sAct1_nsAct1_nsN1_B_one');
I_before_AC = RBIG_MSMI(sAct1_nsAct1_nsN1_A_one',sAct1_nsAct1_nsN1_C_one');
I_before_BC = RBIG_MSMI(sAct1_nsAct1_nsN1_B_one',sAct1_nsAct1_nsN1_C_one');
I_before = mean([I_before_AB,I_before_AC,I_before_BC]);

I_after_AB = RBIG_MSMI(sAct2_nsAct1_nsN1_A_one',sAct2_nsAct1_nsN1_B_one');
I_after_AC = RBIG_MSMI(sAct2_nsAct1_nsN1_A_one',sAct2_nsAct1_nsN1_C_one');
I_after_BC = RBIG_MSMI(sAct2_nsAct1_nsN1_B_one',sAct2_nsAct1_nsN1_C_one');
I_after = mean([I_after_AB,I_after_AC,I_after_BC]);

%TC
data_before = [sAct1_nsAct1_nsN1_A_one; sAct1_nsAct1_nsN1_B_one; sAct1_nsAct1_nsN1_C_one];
T_before = RBIG_TC(data_before');

data_after = [sAct2_nsAct1_nsN1_A_one; sAct2_nsAct1_nsN1_B_one; sAct2_nsAct1_nsN1_C_one];
T_after = RBIG_TC(data_after');

%SI
data_before = [sAct1_nsAct1_nsN1_A_one(:), sAct1_nsAct1_nsN1_B_one(:), sAct1_nsAct1_nsN1_C_one(:)];
data_before_AB = cmi_ggg(sAct1_nsAct1_nsN1_A_one,sAct1_nsAct1_nsN1_B_one,sAct1_nsAct1_nsN1_C_one);
data_before_AC = cmi_ggg(sAct1_nsAct1_nsN1_A_one,sAct1_nsAct1_nsN1_C_one,sAct1_nsAct1_nsN1_B_one);
data_before_BC = cmi_ggg(sAct1_nsAct1_nsN1_B_one,sAct1_nsAct1_nsN1_C_one,sAct1_nsAct1_nsN1_A_one);
data_before_SI = T_before - (data_before_AB+data_before_AC+data_before_BC);

data_after_AB = cmi_ggg(sAct2_nsAct1_nsN1_A_one,sAct2_nsAct1_nsN1_B_one, sAct2_nsAct1_nsN1_C_one);
data_after_AC = cmi_ggg(sAct2_nsAct1_nsN1_A_one,sAct2_nsAct1_nsN1_C_one, sAct2_nsAct1_nsN1_B_one);
data_after_BC = cmi_ggg(sAct2_nsAct1_nsN1_B_one,sAct2_nsAct1_nsN1_C_one, sAct2_nsAct1_nsN1_A_one);
data_after_SI = T_after - (data_after_AB+data_after_AC+data_after_BC);

[I_after_AB, I_after_AC,  I_after_BC, data_after_SI,  T_after]
[abs(I_after_AB-I_before_AB)/I_before_AB,  abs(I_after_AC-I_before_AC)/I_before_AC, abs(I_after_BC-I_before_BC)/I_before_BC, abs(data_after_SI-data_before_SI)/abs(data_before_SI), abs(T_after-T_before)/T_before]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




 