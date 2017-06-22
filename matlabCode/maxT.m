function [realT, maxT_thresh] = maxT(diff_arr, varargin)
    % Taku Ito
    % 07/15/2017
    %
    % Code to perform permutation testing to control for family-wise error (FWE)
    % Using max-T approach as described in Nichols & Holmes (2002)
    % Nichols TE, Holmes AP. (2002). Nonparametric permutation tests for functional neuroimaging: A primer with Examples. Hum. Brain Mapp., 15: 1-25. doi:10.1002/hbm.1058
    % 
    % MATLAB version
    %
    % Required Parameters:
    %    diff_arr    =   MxN matrix of set of M independent tests for condition 1 minus condition 2 across N subjects
    %                    diff_arr can also be an array of multiple values (or tests) compared against the nullmean (or null mean)
    % Optional Parameters:
    %    nullmean    =   Expected value of the null hypothesis {default = 0, for a t-test against 0}
    %    alpha       =   alpha value to return the maxT threshold {default = .05}
    %    tail        =   [0, 1, or -1] 
    %                    If tail = 1, reject the null hypothesis if the statistic is greater than the null dist (upper tailed test).  
    %                    If tail = -1, reject the null hypothesis if the statistic is less than the null dist (lower tailed test). 
    %                    If tail = 0, reject the null hypothesis for a two-tailed test
    %                    {default : 0} 
    %    permutations =  Number of permutations to perform {default = 1000}
    %    nproc       =   number of processes to run in parallel {default = 1}
    %    pvals       =   if True, returns equivalent p-value distribution for all t-values {default = True}
    %
    % Returns:
    %    t: Array of T-values of correct contrast map (Mx1 vector, for M tests)
    %    maxTThreshold   : The t-value threshold corresponding to the corrected alpha value. If a two-tailed test is specified, the maxR is provided as an absolute value
    %
    %
    % EXAMPLE USAGE:
    %     Data is in a 2D matrix, i.e., variable X observation (e.g., voxels X subjects or connections X subjects)
    %     We want to test the significance of resting-state FC values against 0, run 1000 permutations, and use 10 processors
    %     [realT, maxT_thresh] = permutationTesting(data, 'nullmean', 0, 'alpha', .05, 'tail', 1, 'permutations', 1000, 'nproc', 10);

    % Instantiate input parser
    p = inputParser;
    % Specify default parameters (if keyword arguments not provided)
    default_nullmean = 0;
    default_alpha = .05;
    default_tail = 0;
    default_permutations = 1000;
    default_nproc = 1;
    addRequired(p, 'diff_arr');
    addOptional(p, 'nullmean', default_nullmean, @isnumeric);
    addOptional(p, 'alpha', default_alpha, @isnumeric);
    addOptional(p, 'tail', default_tail, @isnumeric);
    addOptional(p, 'permutations', default_permutations, @isnumeric);
    addOptional(p, 'nproc', default_nproc, @isnumeric);

    % Parse inputs
    parse(p, diff_arr, varargin{:});
    nullmean = p.Results.nullmean;
    alpha = p.Results.alpha;
    tail = p.Results.tail;
    permutations = p.Results.permutations;
    nproc = p.Results.nproc;


    % Prepare inputs for parallel processing
    seeds = randi(10000000,permutations,1);
    maxT_dist = zeros(permutations,1);
    parfor (i=1:permutations, nproc)
        seed = seeds(i); 
        maxT_dist(i) = runPermutation(diff_arr, nullmean, tail, seed);
    end

    % Obtain real t-values 
    [H, P, CI, STATS] = ttest(diff_arr,nullmean,'dim',2);
    realT = STATS.tstat;

    % Find threshold for alpha
    maxT_dist_sorted = sort(maxT_dist);
    % Specify which tail we want
    if tail == 1
        topPercVal_maxT_inx = length(maxT_dist_sorted)*(1-alpha);
    elseif tail == -1
        topPercVal_maxT_inx = length(maxT_dist_sorted)*(alpha);
    elseif tail == 0
        topPercVal_maxT_inx = length(maxT_dist_sorted)*(1-alpha);
    end

    maxT_thresh = maxT_dist_sorted(topPercVal_maxT_inx+1);

    % removed old functionality
%    % Construct ECDF from maxT_dist
%    %[Fe Xe]  = ecdf(maxT_dist);
%    sortedT = sort(maxT_dist);
%    p_fwe = zeros(size(realT));
%    % Compute p-values for each real T
%    for test=1:size(realT,1)
%        t = realT(test);
%        p_fwe(test) = (sum(sortedT < t))/permutations;
%        
%    end
end


function maxT = runPermutation(diff_arr,nullmean,tail,seed)
    % Helper function to perform a single permutation

    % Set random seed
    rand('seed',seed)
    

    % Create a random matrix to shuffle conditions (randomly multiply contrasts by 1 or -1)
    shufflemat = randn(size(diff_arr)); % Sample values from a gaussian dist.
    pos = shufflemat > 0; % binarize into pos values 
    neg = shufflemat < 0; % binarize into neg values
    % matrix of 1 and -1
    shufflemat = pos + neg*(-1); % add pos and neg values into a 'shuffled matrix'

    % Shuffle original raw values
    diff_arr = diff_arr.*shufflemat;

    % Take t-test against 0 for each independent test 
    [H, P, CI, STATS] = ttest(diff_arr,nullmean,'dim',2);
    if tail == 1
        maxT = max(STATS.tstat);
    elseif tail == -1
        maxT = min(STATS.tstat);
    elseif tail == 0
        maxT = max(abs(STATS.tstat));
    end

end
