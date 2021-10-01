#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# QiangLI
#
# University of Valencia
# Copyright (c) 2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Step1: download the dataset
from nilearn import datasets
from scipy.io import savemat

development_dataset = datasets.fetch_development_fmri(n_subjects=122, age_group='child')

msdl_data = datasets.fetch_atlas_msdl()
msdl_coords = msdl_data.region_coords
n_regions = len(msdl_coords)
labels=msdl_data.labels
print('MSDL has {0} ROIs, part of the following networks :\n{1}.'.format(n_regions, msdl_data.networks))

# Visualization
from nilearn import plotting
msdl = datasets.fetch_atlas_msdl()
atlas_types = {'MSDL': msdl.maps}

for name, atlas in sorted(atlas_types.items()):
    plotting.plot_prob_atlas(atlas, title=name)

print('ready')
plotting.show()


# Step2:  Region signals extraction
from nilearn import input_data

masker = input_data.NiftiMapsMasker(
    msdl_data.maps, resampling_target="data", t_r=2, detrend=True,
    low_pass=.1, high_pass=.01, memory='nilearn_cache', memory_level=1).fit()


children = []
adult  = []
pooled_subjects = []
groups = []  # child or adult
for func_file, confound_file, phenotypic in zip(
        development_dataset.func,
        development_dataset.confounds,
        development_dataset.phenotypic):
    time_series = masker.transform(func_file, confounds=confound_file)
    pooled_subjects.append(time_series)
    if phenotypic['Child_Adult'] == 'child':
        children.append(time_series)
    if phenotypic['Child_Adult'] == 'adult':
        adult.append(time_series)
    groups.append(phenotypic['Child_Adult'])
savemat("child.mat", {"child":children})
savemat("adult.mat", {"adult":adult})

print('Data has {0} children.'.format(len(children)))
print('Data has {0} aadult.'.format(len(adult)))
