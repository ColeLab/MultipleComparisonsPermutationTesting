%function [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t2(data1,data2,n_perm,tail,alpha_level,mu,t_stat,reports,seed_state)
%
% mult_comp_perm_t2-Independent sample permutation test based on a
% t-statistic.  This function can perform the test on one variable or 
% simultaneously on multiple variables.  When applying the test to multiple
% variables, the "tmax" method is used for adjusting the p-values of each
% variable for multiple comparisons (Blair & Karnisky, 1993; Westfall & 
% Young, 1993).  Like Bonferroni correction, this method adjusts p-values 
% in a way that strongly controls the family-wise error rate.  However, the 
% permutation method will be more powerful than Bonferroni correction when 
% different variables in the test are correlated.
%
% Required Inputs:
%  data1   - 2D matrix of data (Observation x Variable)
%  data2   - 2D matrix of data (Observation x Variable)
%
% Optional Inputs:
%  n_perm      - Number of random permutations used to estimate the 
%                distribution of the null hypothesis.  Manly (1997) 
%                suggests using at least 1000 permutations for an alpha level 
%                of 0.05 and at least 5000 permutations for an alpha level of
%                0.01. {default=5000}
%  tail        - [1, 0, or -1] If tail=1, the alternative hypothesis is that the
%                mean of the data1 is greater than that of data2 (upper tailed test).
%                If tail=0, the alternative hypothesis is that the group
%                means differ (two tailed test).  If tail=-1, the alternative hypothesis
%                is that the mean of data1 is less than that of data2 (lower tailed test).
%                {default: 0}
%  alpha_level - Desired family-wise alpha level. Note, because of the finite
%                number of possible permutations, the exact desired family-wise
%                alpha may not be possible. Thus, the closest approximation 
%                is used and output as est_alpha. {default=.05}
%  mu          - The mean of the null hypothesis.  Must be a scalar or
%                1 x n vector (where n=the number of variables). {default: 0}
%  t_stat      - ['t','w', or 'td'] If 't', the standard independent samples
%                t-statistic is used that pools the variances between groups.
%                If 'w', Welch's t-statistic is used, which is somewhat
%                robust to differences in group variance and almost as 
%                powerful as the standard t statistic (Groppe et al., 2011a).
%                If 'td', each group is simply normalized by its estimated 
%                standard deviation and the difference between groups is 
%                taken. This is the t_dif statistic in Groppe et al. (2011a), 
%                which is very insensitive to differences in variances between 
%                groups but is much less powerful than the standard 
%                t-statistic. {default: 0}
%  reports     - [0 or 1] If 0, function proceeds with no command line
%                reports. Otherwise, function reports what it is doing to
%                the command line. {default: 1}
%  seed_state  - The initial state of the random number generating stream
%                (see MATLAB documentation for "randstream"). If you pass
%                a value from a previous run of this function, it should
%                reproduce the exact same values. Not supported for all
%                versions of MATLAB.
%
% Outputs:
%  pval       - p-value (adjusted for multiple comparisons) of each
%               variable
%  t_orig     - t-score for each variable
%  crit_t     - Lower and upper critical t-scores for given alpha level. 
%               t-scores that exceed these values significantly deviate from 
%               the null hypothesis.  For upper tailed tests, the lower
%               critical t-score is NaN. The opposite is true of lower
%               tailed tests.
%  est_alpha  - The estimated family-wise alpha level of the test.  With 
%               permutation tests, a finite number of p-values are possible.
%               This function tries to use an alpha level that is as close 
%               as possible to the desired alpha level.  However, if the 
%               sample size is small, a very limited number of p-values are 
%               possible and the desired family-wise alpha level may be 
%               impossible to approximately achieve.
%  seed_state - The initial state of the random number generating stream
%               (see MATLAB documentation for "randstream") used to 
%               generate the permutations. You can use this to reproduce
%               the output of this function.
%
% Note:
% -Unlike a parametric test (e.g., an ANOVA), a discrete set of p-values
% are possible (at most the number of possible permutations).  Since the
% number of possible permutations grows combinatorially with the number of
% participants, this is only an issue for small sample sizes (e.g., 6
% participants). When you have such a small sample size, the
% limited number of p-values may make the test overly conservative (e.g., 
% you might be forced to use an alpha level of .0286 since it is the biggest
% possible alpha level less than .05).
%
% -The null hypothesis of the permutation test is that the data from the
% two groups are independent samples from the same population. The
% permutation test will be sensitive to ANY difference in the distribution
% between the two groups (e.g., differences in means, variance, or
% skewness). This can be a problem there are uninteresting differences 
% between groups (e.g., the data from patients is more variable than the 
% data from control subjects). Choice of permutation statistic and 
% differences in the number of observations in each group will effect how 
% sensitive the permutation test is to various possible differences in 
% group population distributions (e.g., Groppe, et al., 2011a).
%
%
% Example:
% >> data1=randn(16,5);  %5 variables, 16 observations
% >> data1(:,1:2)=data1(:,1:2)+1.25; %mean of first two variables is 1
% >> data2=randn(15,5);  %5 variables, 15 observations
% >> [pval, t_orig, crit_t, est_alpha]=mult_comp_perm_t2(data1,data2);
% >> fprintf('Multiple comparison adjusted p-values:\n');
% >> disp(pval);
%
%
% References:
% Blair, R.C. & Karniski, W. (1993) An alternative method for significance
% testing of waveform difference potentials. Psychophysiology.
%
% Groppe, D.M., Urbach, T.P., & Kutas, M. (2011a) Mass univariate analysis 
% of event-related brain potentials/fields II: Simulation studies. 
% Psychophysiology, 48(12) pp. 1726-1737, DOI: 10.1111/j.1469-8986.2011.01272.x 
% http://www.cogsci.ucsd.edu/~dgroppe/PUBLICATIONS/mass_uni_preprint2.pdf
%
% Manly, B.F.J. (1997) Randomization, Bootstrap, and Monte Carlo Methods in
% Biology. 2nd ed. Chapman and Hall, London.
%
% Westfall, P.H. & Young, S.S. (1993) Resampling-Based Multiple Testing: 
% Examples and Methods for p-values Adjustment. Wiley, New York.
%
%
% For a review on permutation tests and other contemporary techniques for 
% correcting for multiple comparisons see:
%
%   Groppe, D.M., Urbach, T.P., & Kutas, M. (2011b) Mass univariate analysis 
% of event-related brain potentials/fields I: A critical tutorial review. 
% Psychophysiology, 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x 
% http://www.cogsci.ucsd.edu/~dgroppe/PUBLICATIONS/mass_uni_preprint1.pdf
%
%
% Author:
% David Groppe
% Dec, 2015
% Toronto
%

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% 3/17/2013-Randomization is now compatible with Matlab v13. Thanks to
% Aaron Newman for the fix.
%
%

function [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t2(data1,data2,n_perm,tail,alpha_level,mu,t_stat,reports,seed_state)

if nargin<2,
    error('You need to provide two sets of data.');
end

if nargin<3,
    n_perm=5000;
end

if nargin<4,
    tail=0;
elseif (tail~=0) && (tail~=1) && (tail~=-1),
    error('Argument ''tail'' needs to be 0,1, or -1.');
end

if nargin<5,
    alpha_level=0.05;
elseif (alpha_level>=1) || (alpha_level<=0),
   error('Argument ''alpha_level'' needs to be a number between 0 and 1.'); 
end

if nargin<6,
    mu=0;
end

if nargin<7,
    t_stat='t';
else
   if ~(strcmpi(t_stat,'t') || strcmpi(t_stat,'w') || strcmpi(t_stat,'td'))
      error('%s is an invalid value for t_stat argument',t_stat); 
   end
end

if nargin<8,
    reports=1;
end

if reports,
    if strcmpi(t_stat,'td')
        fprintf('Variance normalized (i.e., t_dif) t-statistic will be used.\n');
    elseif strcmpi(t_stat,'w')
        fprintf('Welch''s t-statistic will be used.\n');
    else
        fprintf('Pooled variance (i.e., standard independent samples) t-statistic will be used.\n');
    end
end

%Get random # generator state
try
    if verLessThan('matlab','8.1')
        defaultStream=RandStream.getDefaultStream;
    else
        defaultStream=RandStream.getGlobalStream;
    end
    if (nargin<9) || isempty(seed_state),
        %Store state of random number generator
        seed_state=defaultStream.State;
    else
        defaultStream.State=seed_state; %reset random number generator to saved state
    end
catch
    fprintf('Sorry, unable to set or return random number generator state with your version of Matlab.\n');
    seed_state=NaN;
end

if isempty(data1),
    error('You have no observations in data1.')
end
if isempty(data2),
    error('You have no observations in data2.')
end
[n_obs1, n_var]=size(data1);
[n_obs2, n_var2]=size(data2);
if n_var~=n_var2
   error('The number of variables in data1 and data2 are not equal.'); 
end

totalObs=n_obs1+n_obs2;
if reports,
    warning('off','all'); %for large # of subjects, nchoosek warns that its result is approximate
    n_psbl_prms=nchoosek(totalObs,n_obs1);
    if n_psbl_prms<100,
        watchit(sprintf(['Due to the very limited number of participants in each group,' ...
            ' the total number of possible permutations is small.\nThus only a limited number of p-values (at most %d) are possible and the test might be overly conservative.'], ...
            n_psbl_prms));
    end
    warning('on','all');
end

if alpha_level<=.01 && n_perm<5000,
    watchit(sprintf('You are probably using too few permutations for an alpha level of %f. Type ">>help mult_comp_perm_t2.m" for more info.', ...
        alpha_level));
elseif alpha_level<=.05 && n_perm<1000,
    watchit(sprintf('You are probably using too few permutations for an alpha level of %f. Type ">>help mult_comp_perm_t2.m" for more info.', ...
        alpha_level));
end


%% Remove null hypothesis mean from data
if isscalar(mu),
    data1=data1-mu;
elseif isvector(mu)
    s_mu=size(mu);
    if s_mu(1)>1,
        mu=mu';
        s_mu=size(mu);
    end
    if s_mu(2)~=n_var,
        error('mu needs to be a scalar or a 1 x %d vector (%d is the number of variables).',n_var,n_var);
    end
    data1=data1-repmat(mu,n_obs1,1);
else
    error('mu needs to be a scalar or a 1 x %d vector (%d is the number of variables).',n_var,n_var);
end

df=totalObs-2;
if reports,
    fprintf('mult_comp_perm_t2: Number of variables: %d\n',n_var);
    fprintf('mult_comp_perm_t2: Number of observations in Data Set 1: %d\n',n_obs1);
    fprintf('mult_comp_perm_t2: Number of observations in Data Set 2: %d\n',n_obs2);
    if ~strcmpi(t_stat,'td'),
        fprintf('t-score degrees of freedom: %d\n',df);
    end
end


%% Set up permutation test
if reports,
    fprintf('Executing permutation test with %d permutations...\n',n_perm);
    fprintf('Permutations completed: ');
end


%% Compute permutations
% Factor that is used to compute t-scores.  Saves time to compute it
% now rather than to compute it anew for each permutation.
mult_fact=(totalObs)/(n_obs1*n_obs2);

mx_t=zeros(1,n_perm);
for perm=1:n_perm
    if ~rem(perm,100)
        if (reports),
            if ~rem(perm-100,1000)
                fprintf('%d',perm);
            else
                fprintf(', %d',perm);
            end
            if ~rem(perm,1000)
                fprintf('\n');
            end
        end
    end
    %randomly assign participants to conditions
    r=randperm(totalObs);
    grp1=r(1:n_obs1);
    grp2=r((n_obs1+1):totalObs);
    %compute most extreme t-score
    if strcmpi(t_stat,'td'),
        mx_t(perm)=tdif([data1; data2],grp1,grp2,n_obs1,n_obs2);
    elseif strcmpi(t_stat,'w'),
        mx_t(perm)=tmax_welch([data1; data2],grp1,grp2,n_obs1,n_obs2); % ??
    else
        mx_t(perm)=tmax2([data1; data2],grp1,grp2,n_obs1,n_obs2,df,mult_fact);
    end
end
    
%End of permutations, print carriage return if it hasn't already been done
%(i.e., perm is NOT a multiple of 1000)
if (reports) && rem(perm,1000)
    fprintf('\n');
end

%% Compute critical t's
if tail==0,
    %two-tailed, test statistic is biggest absolute value of all t's
    mx_t=abs(mx_t);
    crit_t(2)=prctile(mx_t,100-100*alpha_level);
    crit_t(1)=-crit_t(2);
    est_alpha=mean(mx_t>=crit_t(2));
elseif tail==1,
    %upper tailed
    crit_t=prctile(mx_t,100-100*alpha_level);
    est_alpha=mean(mx_t>=crit_t);
else
    %tail=-1, lower tailed
    crit_t=prctile(mx_t,alpha_level*100);
    est_alpha=mean(mx_t<=crit_t);
end
if reports,
    fprintf('Desired family-wise alpha level: %f\n',alpha_level);
    fprintf('Estimated actual family-wise alpha level: %f\n',est_alpha);
end


%% Compute t-scores and p-values of actual observations
if strcmpi(t_stat,'td'),
    [dummy, t_orig]=tdif([data1; data2],grp1,grp2,n_obs1,n_obs2);
elseif strcmpi(t_stat,'w'),
    [dummy, t_orig]=tmax_welch([data1; data2],1:n_obs1,(n_obs1+1):totalObs,n_obs1,n_obs2); % ??
else
    [dummy, t_orig]=tmax2([data1; data2],1:n_obs1,(n_obs1+1):totalObs,n_obs1,n_obs2,df,mult_fact);
end
pval=zeros(1,n_var);
%compute p-values
for t=1:n_var,
    if tail==0,
        pval(t)=mean(mx_t>=abs(t_orig(t))); %note mx_t are now all positive due to abs command above
    elseif tail==1,
        pval(t)=mean(mx_t>=t_orig(t));
    elseif tail==-1,
        pval(t)=mean(mx_t<=t_orig(t));
    end
end




function [mx_t, all_t]=tmax2(dat,grp1,grp2,n_obs1,n_obs2,df,mult_fact)
% would it be faster if included this in the main function?

x1=dat(grp1,:);
x2=dat(grp2,:);

sm1=sum(x1);
mn1=sm1/n_obs1;
ss1=sum(x1.^2)-(sm1.^2)/n_obs1;

sm2=sum(x2);
mn2=sm2/n_obs2;
ss2=sum(x2.^2)-(sm2.^2)/n_obs2;

pooled_var=(ss1+ss2)/df;
stder=sqrt(pooled_var*mult_fact);

all_t=(mn1-mn2)./stder;
[~, mxId]=max(abs(all_t));
mx_t=all_t(mxId);


function [mx_t, all_t]=tdif(dat,grp1,grp2,n_obs1,n_obs2)
% Computes t_dif
% would it be faster if included this in the main function?

x1=dat(grp1,:);
x2=dat(grp2,:);

sm1=sum(x1);
mn1=sm1/n_obs1;
var1=(sum(x1.^2)-(sm1.^2)/n_obs1)/(n_obs1-1);

sm2=sum(x2);
mn2=sm2/n_obs2;
var2=(sum(x2.^2)-(sm2.^2)/n_obs2)/(n_obs2-1);

all_t=mn1./sqrt(var1)-mn2./sqrt(var2); % This effectively converts the data from each group into t-scores and then takes the difference
[~, mxId]=max(abs(all_t));
mx_t=all_t(mxId);


function [mx_t, all_t]=tmax_welch(dat,grp1,grp2,n_obs1,n_obs2)
% Computes Welch's t
% would it be faster if included this in the main function?

x1=dat(grp1,:);
x2=dat(grp2,:);

sm1=sum(x1);
mn1=sm1/n_obs1;
var1=(sum(x1.^2)-(sm1.^2)/n_obs1)/(n_obs1-1);

sm2=sum(x2);
mn2=sm2/n_obs2;
var2=(sum(x2.^2)-(sm2.^2)/n_obs2)/(n_obs2-1);

all_t=(mn1-mn2)./sqrt(var1/n_obs1+var2/n_obs2); 
[~, mxId]=max(abs(all_t));
mx_t=all_t(mxId);


function watchit(msg)
%function watchit(msg)
%
% Displays a warning message on the Matlab command line.  Used by 
% several Mass Univariate ERP Toolbox functions.
%

disp(' ');
disp('****************** Warning ******************');
disp(msg);
disp(' ');

%% Simple Code for Testing Family-Wise Error Rate or Power
% n_sim=1000;
% n_sig=zeros(n_sim,1);
% for a=1:n_sim,
%     if ~mod(a,10)
%        fprintf('Simulation %d/%d\n',a,n_sim); 
%     end
%     data1=randn(16,5);  %5 variables, 16 observations
%     %data1(:,1:2)=data1(:,1:2)+1.25; %mean of first two variables is 1
%     data2=randn(15,5);  %5 variables, 15 observations
%     [pval, t_orig, crit_t, est_alpha]=mult_comp_perm_t2(data1,data2,3000,0,0.05,0,'w',0);
%     n_sig(a)=sum(pval<0.05);
% end
% fprintf('Estimated FWER: %f\n',mean(n_sig>0));