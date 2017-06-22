% Taku Ito
% 6/16/2017

% Test script to test code for permutation testing in matlab
clear all
clc

addpath('../matlabCode/')

%% Create data set
numVoxels = 100;
numSubjs = 20;
nCond = 2;
sigEffects = 20;
effectAmp = 3;

disp('Running simulation - upper-tailed test...')
disp(['Running contrasts on ' num2str(numVoxels) '  voxels (' num2str(numVoxels) ' independent tests) with ' num2str(numSubjs) ' subjects'])
disp([num2str(nCond) ' conditions, with ' num2str(sigEffects) ' significant real effects'])

tic
% Construct data set (already contrasted)
dataSet = zeros(numVoxels,numSubjs,nCond);
for vox=1:numVoxels
    if vox <= sigEffects
        dataSet(vox,:, 1) = randn(1,numSubjs) + effectAmp;
        dataSet(vox,:, 2) = randn(1,numSubjs);
    else
        dataSet(vox,:, 1) = randn(1,numSubjs);
        dataSet(vox,:, 2) = randn(1,numSubjs);
    end
end

% Create contrast map
contrastSet = dataSet(:,:,1) - dataSet(:,:,2);


%% Run permutation test
alpha = .05;
[t, maxT_thresh] = maxT(contrastSet,'nullmean', 0, ...
                              'alpha', alpha, ...
                              'tail', 1, ...
                              'permutations',10000, ...
                              'nproc', 10);
toc

disp(['Number of true effects: ' num2str(sigEffects)])
disp(['Number of statistically significant effects (p < .05): ' num2str(sum(t>maxT_thresh))])

%%%%%%%%%%%%%%%%%
disp('########################################')
disp('Running simulation - lower-tailed test...')
disp(['Running contrasts on ' num2str(numVoxels) '  voxels (' num2str(numVoxels) ' independent tests) with ' num2str(numSubjs) ' subjects'])
disp([num2str(nCond) ' conditions, with ' num2str(sigEffects) ' significant real effects'])

tic
% Construct data set (already contrasted)
dataSet = zeros(numVoxels,numSubjs,nCond);
for vox=1:numVoxels
    if vox <= sigEffects
        dataSet(vox,:, 1) = randn(1,numSubjs) - effectAmp;
        dataSet(vox,:, 2) = randn(1,numSubjs);
    else
        dataSet(vox,:, 1) = randn(1,numSubjs);
        dataSet(vox,:, 2) = randn(1,numSubjs);
    end
end

% Create contrast map
contrastSet = dataSet(:,:,1) - dataSet(:,:,2);


%% Run permutation test
alpha = .05;
[t, maxT_thresh] = maxT(contrastSet,'nullmean', 0, ...
                              'alpha', alpha, ...
                              'tail', -1, ...
                              'permutations',10000, ...
                              'nproc', 10);
toc

disp(['Number of true effects: ' num2str(sigEffects)])
disp(['Number of statistically significant effects (p < .05): ' num2str(sum(t<maxT_thresh))])

%%%%%%%%%%%%%%%%%
disp('########################################')
disp('Running simulation - two-tailed test...')
disp(['Running contrasts on ' num2str(numVoxels) '  voxels (' num2str(numVoxels) ' independent tests) with ' num2str(numSubjs) ' subjects'])
disp([num2str(nCond) ' conditions, with ' num2str(sigEffects) ' significant real effects'])

tic
% Construct data set (already contrasted)
dataSet = zeros(numVoxels,numSubjs,nCond);
for vox=1:numVoxels
    if vox <= sigEffects
        % if voxel is odd, mave it be a negative contrast effect
        if mod(vox,2) == 1
            dataSet(vox,:, 1) = randn(1,numSubjs) - effectAmp;
            dataSet(vox,:, 2) = randn(1,numSubjs);
        % if voxel is even, mave it be a positive contrast effect
        elseif mod(vox,2)==0
            dataSet(vox,:, 1) = randn(1,numSubjs) + effectAmp;
            dataSet(vox,:, 2) = randn(1,numSubjs);
        end
    else
        dataSet(vox,:, 1) = randn(1,numSubjs);
        dataSet(vox,:, 2) = randn(1,numSubjs);
    end
end

% Create contrast map
contrastSet = dataSet(:,:,1) - dataSet(:,:,2);


%% Run permutation test
alpha = .05;
[t, maxT_thresh] = maxT(contrastSet,'nullmean', 0, ...
                              'alpha', alpha, ...
                              'tail', 0, ...
                              'permutations',10000, ...
                              'nproc', 10);
toc

disp(['Number of true effects: ' num2str(sigEffects)])
disp(['Number of statistically significant effects (p < .05): ' num2str(sum(abs(t)>maxT_thresh))])
