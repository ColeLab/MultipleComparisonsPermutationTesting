function [realT, p_fwe, maxT_dist] = permutationTesting(diff_arr, varargin)
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
    %     nullmean = Expected value of the null hypothesis (default 0, for a t-test against 0)
    %     permutations = Number of permutations to perform (default 1000)
    %     nproc = number of processes to run in parallel (default 1)
    % 
    % Returns:
    %     t: Array of T-values of correct contrast map (Mx1 vector, for M tests)
    %     p: Array of FWE-corrected p-values (Mx1 vector, for M tests); 
    %        Note, p-values correspond to values on the CDF. One-sided or or two-sided p-values can be computed accordingly.
    % 
    % N.B.: Only works for paired one-sample t-tests

    % Instantiate input parser
    p = inputParser;
    % Specify default parameters (if keyword arguments not provided)
    default_nullmean = 0;
    default_permutations = 1000;
    default_nproc = 1;
    addRequired(p, 'diff_arr');
    addOptional(p, 'nullmean', default_nullmean, @isnumeric);
    addOptional(p, 'permutations', default_permutations, @isnumeric);
    addOptional(p, 'nproc', default_nproc, @isnumeric);
    % Parse inputs
    parse(p, diff_arr, varargin{:});
    nullmean = p.Results.nullmean;
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

    % Construct ECDF from maxT_dist
    %[Fe Xe]  = ecdf(maxT_dist);
    sortedT = sort(maxT_dist);
    p_fwe = zeros(size(realT));
    % Compute p-values for each real T
    for test=1:size(realT,1)
        t = realT(test);
        p_fwe(test) = (sum(sortedT < t))/permutations;
        
    end
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
