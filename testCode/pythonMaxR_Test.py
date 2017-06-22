# Taku Ito
# 6/16/2017

import numpy as np
import sys
sys.path.append('../pythonCode/')
import permutationTesting as pt


## Create data set 
numVoxels = 1000
numSubjs = 100
sigEffects = 10
behavArr = np.random.normal(0,1,(numSubjs,))

print 'Running simulation for upper-tailed test...'
print 'Running correlations on', str(numVoxels), 'voxels ('+ str(numVoxels), 'independent tests) with', str(numSubjs), 'subjects'
print sigEffects, 'significant real effects'

dataSet = np.zeros((numVoxels,numSubjs))

for vox in range(numVoxels):
    if vox < sigEffects:
        dataSet[vox,:] = np.random.normal(0,1,(numSubjs,)) + behavArr
    else:
        dataSet[vox,:] = np.random.normal(0,1,(numSubjs,))

# Run permutation test
alpha = .05
r, maxR_thresh = pt.maxR(dataSet, behavArr,alpha=alpha,tail=1,permutations=10000, nproc=10)

print 'Number of true effects:', sigEffects
print 'Number of statistically significant effects (p < .05):',np.sum(r>maxR_thresh)


######################################################

print '\n**********************\n'
print 'Running simulation for lower-tailed test...'
print 'Running correlations on', str(numVoxels), 'voxels ('+ str(numVoxels), 'independent tests) with', str(numSubjs), 'subjects'
print  sigEffects, 'significant real effects'

dataSet = np.zeros((numVoxels,numSubjs))

for vox in range(numVoxels):
    if vox < sigEffects:
        dataSet[vox,:] = np.random.normal(0,1,(numSubjs,)) - behavArr
    else:
        dataSet[vox,:] = np.random.normal(0,1,(numSubjs,))

# Run permutation test
alpha = .05
r, maxR_thresh = pt.maxR(dataSet, behavArr,alpha=alpha,tail=-1,permutations=10000, nproc=10)

print 'Number of true effects:', sigEffects
print 'Number of statistically significant effects (p < .05):',np.sum(r<maxR_thresh)


######################################################

print '\n**********************\n'
print 'Running simulation for two-tailed test...'
print 'Running correlations on', str(numVoxels), 'voxels ('+ str(numVoxels), 'independent tests) with', str(numSubjs), 'subjects'
print sigEffects, 'significant real effects'

dataSet = np.zeros((numVoxels,numSubjs))

for vox in range(numVoxels):
    if vox < sigEffects:
        # Odd effects are lower-tailed
        if vox%2==1:
            dataSet[vox,:] = np.random.normal(0,1,(numSubjs,)) - behavArr 

        # Even ieffects a upper-tailed
        elif vox%2==0:
            dataSet[vox,:] = np.random.normal(0,1,(numSubjs,)) + behavArr
    else:
        dataSet[vox,:] = np.random.normal(0,1,(numSubjs,))


# Run permutation test
alpha = .05
r, maxR_thresh = pt.maxR(dataSet, behavArr,alpha=alpha,tail=0,permutations=10000, nproc=10)

print 'Number of true effects:', sigEffects
print 'Number of statistically significant effects (p < .05):',np.sum(r>maxR_thresh) + np.sum(r<-maxR_thresh)
