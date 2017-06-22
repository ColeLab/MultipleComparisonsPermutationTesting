% Taku Ito
% 6/16/2017

% Test script to test code for permutation testing in matlab
clear all
clc

addpath('../matlabCode/')

%% Create data set
numVoxels = 100;
numSubjs = 100;
sigEffects = 20;
behavArr = randn(numSubjs,1);

disp('Running simulation - upper-tailed test...')
disp(['Running correlations on ' num2str(numVoxels) '  voxels (' num2str(numVoxels) ' independent tests) with ' num2str(numSubjs) ' subjects'])
disp([num2str(sigEffects) ' significant real effects'])

tic
% Construct data set
dataSet = zeros(numVoxels,numSubjs);
for vox=1:numVoxels
    if vox <= sigEffects
        dataSet(vox,:) = randn(1,numSubjs) + behavArr';
    else
        dataSet(vox,:) = randn(1,numSubjs);
    end
end


%% Run permutation test
alpha = .05;
[r, maxR_thresh] = maxR(dataSet,behavArr,...
                              'alpha', alpha, ...
                              'tail', 1, ...
                              'permutations',10000, ...
                              'nproc', 10);
toc

disp(['Number of true effects: ' num2str(sigEffects)])
disp(['Number of statistically significant effects (p < .05): ' num2str(sum(r>maxR_thresh))])

%%%%%%%%%%%%%%%%%
disp('########################################')
disp('Running simulation - lower-tailed test...')
disp(['Running correlations on ' num2str(numVoxels) '  voxels (' num2str(numVoxels) ' independent tests) with ' num2str(numSubjs) ' subjects'])
disp([num2str(sigEffects) ' significant real effects'])

tic
% Construct data set (already contrasted)
dataSet = zeros(numVoxels,numSubjs);
for vox=1:numVoxels
    if vox <= sigEffects
        dataSet(vox,:) = randn(1,numSubjs) - behavArr';
    else
        dataSet(vox,:) = randn(1,numSubjs);
    end
end


%% Run permutation test
alpha = .05;
[r, maxR_thresh] = maxR(dataSet, behavArr, ...
                              'alpha', alpha, ...
                              'tail', -1, ...
                              'permutations',10000, ...
                              'nproc', 10);
toc

disp(['Number of true effects: ' num2str(sigEffects)])
disp(['Number of statistically significant effects (p < .05): ' num2str(sum(r<maxR_thresh))])

%%%%%%%%%%%%%%%%%
disp('########################################')
disp('Running simulation - two-tailed test...')
disp(['Running correlations on ' num2str(numVoxels) '  voxels (' num2str(numVoxels) ' independent tests) with ' num2str(numSubjs) ' subjects'])
disp([num2str(sigEffects) ' significant real effects'])

tic
% Construct data set (already contrasted)
dataSet = zeros(numVoxels,numSubjs);
for vox=1:numVoxels
    if vox <= sigEffects
        % if voxel is odd, mave it be a negative contrast effect
        if mod(vox,2) == 1
            dataSet(vox,:) = randn(1,numSubjs) - behavArr';
        % if voxel is even, mave it be a positive contrast effect
        elseif mod(vox,2)==0
            dataSet(vox,:) = randn(1,numSubjs) + behavArr';
        end
    else
        dataSet(vox,:) = randn(1,numSubjs);
    end
end

%% Run permutation test
alpha = .05;
[r, maxR_thresh] = maxR(dataSet, behavArr, ...
                              'alpha', alpha, ...
                              'tail', 0, ...
                              'permutations',10000, ...
                              'nproc', 10);
toc

disp(['Number of true effects: ' num2str(sigEffects)])
disp(['Number of statistically significant effects (p < .05): ' num2str(sum(abs(r)>maxR_thresh))])
