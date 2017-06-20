# Taku Ito
# 07/14/2017
##Modified by MWC June 20, 2017
##Modified by TI June 20, 2017

# Code to perform permutation testing to control for family-wise error (FWE)
# Using max-T approach as described in Nichols & Holmes (2002)
# Nichols TE, Holmes AP. (2002). Nonparametric permutation tests for functional neuroimaging: A primer with Examples. Hum. Brain Mapp., 15: 1-25. doi:10.1002/hbm.1058

import numpy as np
import scipy.stats as stats
import multiprocessing as mp
from statsmodels.distributions.empirical_distribution import ECDF


def maxT(diff_arr, nullmean=0, alpha=.05, tail=1, permutations=1000, nproc=1, pvals=False):
    """
    Performs family-wise error correction using permutation testing (Nichols & Holmes 2002)
    Note! Assumes a one-sided t-test (specify tail of test by tail parameter) 

    Citation: 
        Nichols TE, Holmes AP. (2002). Nonparametric permutation tests for functional neuroimaging: A primer with Examples. Hum. Brain Mapp., 15: 1-25. doi:10.1002/hbm.1058
    Required Parameters:
        diff_arr    =   MxN matrix of set of M independent tests for condition 1 minus condition 2 across N subjects
                        diff_arr can also be an array of multiple values (or tests) compared against the nullmean (or null mean)
    Optional Parameters:
        nullmean    =   Expected value of the null hypothesis {default = 0, for a t-test against 0}
        alpha       =   alpha value to return the maxT threshold {default = .05}
        tail        =   [1 or -1] If tail = 1, reject the null hypothesis if the mean of the data is greater than 0 (upper tailed test).  
                        If tail = -1, reject the null hypothesis if the mean of the data is less than nullmean {default = 1}
        permutations =  Number of permutations to perform {default = 1000}
        nproc       =   number of processes to run in parallel {default = 1}
        pvals       =   if True, returns equivalent p-value distribution for all t-values {default = True}

    Returns:
        t: Array of T-values of correct contrast map (Mx1 vector, for M tests)
        p: Array of FWE-corrected p-values (Mx1 vector, for M tests); 
           Note, p-values correspond to values on the CDF. One-sided or or two-sided p-values can be computed accordingly.

    N.B.: Only works for paired one-sample t-tests
    """
    # Focus on difference matrix -- more computationally feasible (and less data to feed into parallel processing)

    # Prepare inputs for multiprocessing
    inputs = []
    for i in range(permutations):
        seed = np.random.randint(0,100000,1)[0]
        inputs.append((diff_arr,nullmean,seed))

    pool = mp.Pool(processes=nproc)
    result = pool.map_async(_maxTpermutation,inputs).get()
    pool.close()
    pool.join()

    # Returns an array of T-values distributions (our null distribution of "max-T" values)
    maxT_dist = np.asarray(result)

    #Find threshold for alpha
    maxT_dist_sorted = np.sort(maxT_dist)
    # Specify which tail we want
    if tail == 1:
        topPercVal_maxT_inx = int(len(maxT_dist_sorted)*(1-alpha))
    elif tail == -1:
        topPercVal_maxT_inx = int(len(maxT_dist_sorted)*(alpha))

    maxT_thresh = maxT_dist_sorted[topPercVal_maxT_inx]

    # Obtain real t-values 
    t = stats.ttest_1samp(diff_arr, nullmean, axis=1)[0]

    if pvals:
        # Construct ECDF from maxT_dist
        ecdf = ECDF(maxT_dist)

        # Return p-values from maxT_dist using our empirical CDF (FWE-corrected p-values)
        p_fwe = ecdf(t)

        if tail == 1:
            p_fwe = 1.0 - p_fwe
        
        return t, maxT_thresh, p_fwe
    else:
        return t, maxT_thresh


def _maxTpermutation((diff_arr,nullmean,seed)):
    """
    Helper function to perform a single permutation
    """

    np.random.seed(seed)

    # Create a random matrix to shuffle conditions (randomly multiply contrasts by 1 or -1)
    shufflemat = np.random.normal(0,1,diff_arr.shape)
    pos = shufflemat > 0
    neg = shufflemat < 0
    # matrix of 1 and -1
    shufflemat = pos + neg*(-1)

    # Shuffle raw values
    diff_arr = np.multiply(diff_arr, shufflemat)

    # Take t-test against 0 for each independent test 
    t_matrix = stats.ttest_1samp(diff_arr,nullmean,axis=1)[0] 

    maxT = np.max(t_matrix)
    
    return maxT



