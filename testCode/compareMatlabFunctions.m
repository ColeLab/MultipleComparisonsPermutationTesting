% Taku Ito
% 6/16/2017

% Test script to test code for permutation testing in matlab
clear all
clc

addpath('../matlabCode/')
addpath('../matlabCode/mult_comp_perm_t1/')

permutations=10000;

%% Create data set
numVoxels = 10000;
numSubjs = 20;
nCond = 2;
sigEffects = 20;
effectAmp = 3;

disp(['Running simulation analysis with ' num2str(permutations) ' permutations...'])
disp(['Running contrasts on ' num2str(numVoxels) '  voxels (' num2str(numVoxels) ' independent tests) with ' num2str(numSubjs) ' subjects'])
disp([num2str(nCond) ' conditions, with ' num2str(sigEffects) ' significant real effets'])

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
disp('Timing permutation testing using House code')
tic
alpha = .05;
[t, maxT_thresh] = maxT(contrastSet,'nullmean', 0, ...
                  'alpha',alpha,...
                  'tail', 1, ...
                  'permutations',permutations, ...
                  'nproc', 20);
toc

%disp(['Number of true effects: ' num2str(sigEffects)])
%disp(['Number of statistically significant effects (p < .05): ' num2str(sum(p_fwe>.95))])


%% Compare with Groppe's Matlab function
% Format data in his required format
tmpdat = contrastSet';
disp('Running permutations using Groppe''s Matlab function')
tic
[pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t1(tmpdat,permutations,1,alpha,0,0);
toc
