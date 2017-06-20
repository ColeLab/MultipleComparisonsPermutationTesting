# Taku Ito
# 6/16/2017

import numpy as np
import sys
sys.path.append('../pythonCode/')
import permutationTesting as pt


## Create data set 
numVoxels = 1000
numSubjs = 20
nCond = 2
sigEffects = 10
effectAmp = 3

print 'Running simulation...'
print 'Running contrasts on', str(numVoxels), 'voxels ('+ str(numVoxels), 'independent tests) with', str(numSubjs), 'subjects'
print nCond, 'conditions, with ', sigEffects, 'significant real effets'

dataSet = np.zeros((numVoxels,numSubjs,nCond))

for vox in range(numVoxels):
    if vox < sigEffects:
        dataSet[vox,:,0] = np.random.normal(0,1,(numSubjs,)) + effectAmp
        dataSet[vox,:,1] = np.random.normal(0,1,(numSubjs,)) 
    else:
        dataSet[vox,:,0] = np.random.normal(0,1,(numSubjs,))
        dataSet[vox,:,1] = np.random.normal(0,1,(numSubjs,)) 

contrastSet = dataSet[:,:,0] - dataSet[:,:,1]

# Run permutation test
alpha = .05
t, maxT_thresh = pt.maxT(contrastSet,nullmean=0, alpha=alpha,tail=1,permutations=10000, nproc=10)

print 'Number of true effects:', sigEffects
print 'Number of statistically significant effects (p < .05):',np.sum(t>maxT_thresh)
