function [trueR, maxR_thresh] = maxR(data_arr, behav_arr, varargin)
    % Taku Ito
    % 07/15/2017
    %
    % Code to perform permutation testing to control for family-wise error (FWE)
    % Using max-R approach (variation of maxR approach described in Nichols & Holmes (2002))
    % Nichols TE, Holmes AP. (2002). Nonparametric permutation tests for functional neuroimaging: A primer with Examples. Hum. Brain Mapp., 15: 1-25. doi:10.1002/hbm.1058
    % 
    % MATLAB version
    %
    % Required Parameters:
    %     data_arr    =   MxN matrix of set of M independent measurements (e.g., FC-values) across N subjects
    %     behav_arr   =   Nx1 array of behavioral measures for N subjects
    % Optional Parameters:
    %     alpha       =   alpha value to return the maxR threshold {default = .05}
    %     tail        =   [0,1, or -1] 
    %                     If tail = 1, reject the null hypothesis if the correlation is greater than the null dist (upper tailed test).  
    %                     If tail = -1, reject the null hypothesis if the correlation is less than the null dist (lower tailed test). 
    %                     If tail = 0, reject the null hypothesis for a two-tailed test
    %                     {default : 0} 
    %     permutations =  Number of permutations to perform {default = 1000}
    %     nproc       =   number of processes to run in parallel {default = 1}
    %     pvals       =   if True, returns equivalent p-value distribution for all t-values {default = True}
    %
    % Returns:
    %     trueR           : Array of Pearson-r values of the true correlations map (Mx1 vector, for M tests)
    %     maxR_thresh     : The Pearson-r value corresponding to the corrected alpha value. If a two-tailed test is specified, the absolute value of the maxR threshold is provided.
    %     p (optional)    : Array of FWE-corrected p-values (Mx1 vector, for M tests); 
    %
    %
    % EXAMPLE USAGE:
    %     Data is in a 2D matrix, i.e., variable X observation (e.g., voxels X subjects or connections X subjects)
    %     We want to test the significance of resting-state FC values against an array of behavioral values, run 1000 permutations, and use 10 processors
    %     [realR, maxR_thresh] = permutationTesting(data, behav_array, 'alpha', .05, 'tail', 0, 'permutations', 1000, 'nproc', 10);

    % Instantiate input parser
    p = inputParser;
    % Specify default parameters (if keyword arguments not provided)
    default_alpha = .05;
    default_tail = 0;
    default_permutations = 1000;
    default_nproc = 1;
    addRequired(p, 'data_arr');
    addRequired(p, 'behav_arr');
    addOptional(p, 'alpha', default_alpha, @isnumeric);
    addOptional(p, 'tail', default_tail, @isnumeric);
    addOptional(p, 'permutations', default_permutations, @isnumeric);
    addOptional(p, 'nproc', default_nproc, @isnumeric);

    % Parse inputs
    parse(p, data_arr, behav_arr, varargin{:});
    alpha = p.Results.alpha;
    tail = p.Results.tail;
    permutations = p.Results.permutations;
    nproc = p.Results.nproc;

    % Normalize data
    data_arr = zscore(data_arr,0,2);
    behav_arr = zscore(behav_arr);
    % Obtain true correlation r-values
    tmp = repmat(behav_arr', size(data_arr,1),1);  
    trueR = mean(tmp.*data_arr,2);

    % Prepare inputs for parallel processing
    seeds = randi(10000000,permutations,1);
    maxR_dist = zeros(permutations,1);
    parfor (i=1:permutations, nproc)
        seed = seeds(i); 
        maxR_dist(i) = runPermutation(data_arr, behav_arr, tail, seed);
    end

    % Find threshold for alpha
    maxR_dist_sorted = sort(maxR_dist);
    % Specify which tail we want
    if tail == 1
        topPercVal_maxR_inx = length(maxR_dist_sorted)*(1-alpha);
    elseif tail == -1
        topPercVal_maxR_inx = length(maxR_dist_sorted)*(alpha);
    elseif tail == 0
        topPercVal_maxR_inx = length(maxR_dist_sorted)*(1-alpha);
    end

    maxR_thresh = maxR_dist_sorted(topPercVal_maxR_inx+1);

    % removed old functionality
%    % Construct ECDF from maxR_dist
%    %[Fe Xe]  = ecdf(maxR_dist);
%    sortedT = sort(maxR_dist);
%    p_fwe = zeros(size(realT));
%    % Compute p-values for each real T
%    for test=1:size(realT,1)
%        t = realT(test);
%        p_fwe(test) = (sum(sortedT < t))/permutations;
%        
%    end
end


function maxR = runPermutation(data_arr,behav_arr,tail,seed)
    % Helper function to perform a single permutation

    % Set random seed
    rand('seed',seed);
    

    % Randomly permute behavioral data
    indices = randperm(length(behav_arr));
    behav_arr = behav_arr(indices);
    
    % Obtain permutation correlation r-values
    tmp = repmat(behav_arr', size(data_arr,1),1);  
    permR = mean(tmp.*data_arr,2);

    if tail == 1
        maxR = max(permR);
    elseif tail == -1
        maxR = min(permR);
    elseif tail == 0
        maxR = max(abs(permR));
    end

end
