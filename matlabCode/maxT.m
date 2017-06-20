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
    % Parameters:
    %     diff_arr = MxN matrix of set of M independent tests for condition 1 minus condition 2 across N subjects
    %                diff_arr can also be an array of multiple values (or tests) compared against the nullmean (or null mean)
    % Optional parameters:
    %     nullmean      :       Expected value of the null hypothesis {default : 0, for a t-test against 0}
    %     alpha         :       Corrected alpha level {default : .05)
    %     tail          :       [1 or -1] If tail = 1, reject the null hypothesis if the mean of the data is greater than 0 (upper tailed test).  
    %                           If tail = -1, reject the null hypothesis if the mean of the data is less than nullmean {default = 1}
    %     permutations  :       Number of permutations to perform {default 1000}
    %     nproc         :       Number of processes to run in parallel {default 1}
    % 
    % Returns:
    %     t             :   Array of T-values of correct contrast map (Mx1 vector, for M tests)
    %     maxT_thresh   :   maximum t-value for associated alpha-level 
    % 
    % N.B.: Only works for paired one-sample t-tests
    %
    % EXAMPLE USAGE:
    %     Data is in a 2D matrix, i.e., variable X observation (e.g., voxels X subjects or connections X subjects)
    %     We want to test the significance of resting-state FC values against 0, run 1000 permutations, and use 10 processors
    %     [realT, p_fwe] = permutationTesting(data, 'nullmean', 0, 'permutations', 1000, 'nproc', 10);

    % Instantiate input parser
    p = inputParser;
    % Specify default parameters (if keyword arguments not provided)
    default_nullmean = 0;
    default_alpha = .05;
    default_tail = 1;
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
        maxT_dist(i) = runPermutation(diff_arr, nullmean, seed);
    end

    % Obtain real t-values 
    [H, P, CI, STATS] = ttest(diff_arr,nullmean,'dim',2);
    realT = STATS.tstat;

    % Find threshold for alpha
    maxT_dist_sorted = sort(maxT_dist)
    % Specify which tail we want
    if tail == 1
        topPercVal_maxT_inx = length(maxT_dist_sorted)*(1-alpha);
    elseif tail == -1
        topPercVal_maxT_inx = length(maxT_dist_sorted)*(alpha);
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


function maxT = runPermutation(diff_arr,nullmean,seed)
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
    maxT = max(STATS.tstat);

end
