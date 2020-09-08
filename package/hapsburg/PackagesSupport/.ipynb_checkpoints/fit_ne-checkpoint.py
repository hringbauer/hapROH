"""
Class the fits Ne with a MLE framework.
Standard error either from bootstrap over individuals,
or from likelihood curve (kind of assumping ROH are independent,
which is not fully correct)
@Author: Harald Ringbauer, 2020
"""

import numpy as np
import os as os
import sys as sys
import multiprocessing as mp
from statsmodels.base.model import GenericLikelihoodModel
from bisect import bisect_left, bisect_right
import pandas as pd
from hapsburg.PackagesSupport.roh_expectations import Expected_Roh
from hapsburg.PackagesSupport.pp_individual_roh_csvs import give_iid_paths


class MLE_ROH_Ne(GenericLikelihoodModel):
    '''
    Class for MLE estimation of block sharing
    between populations with modeled error
    Bins into length blocks and models block-sharing
    between populations as indep. Poisson
    Operates in cM unit space (not Morgan!)
    '''
    # Bins for mle_multi_run
    min_b, max_b = 0, 60.1  # Minimum/Maximum bin for calculations before error
    bin_width = 0.1  # Bin width for mle_multi_run
    min_len, max_len = 4.0, 20.0  # Minimum/maximum bin length actually analyzed    #4-15
    min_ind, max_ind = 0, 0  # Indices for start stop of bins of interest
    mid_bins = []  # Array for the bins
    chr_lgts=[] # The Lengths of the Chromosome [in centimorgan!]
    
    fp_rate = []  # Array for false positives
    theoretical_shr = []  # Array for theoretical expected sharing per bin
    trans_mat = np.zeros((2, 2))  # Transition matrix for theor. to expected block-sharing
    full_shr_pr = []  # Array for the full bin sharing 
      
    density_fun = 0  # function used to calculate the block sharing density; is required to be per cM!! 
    start_params = []  # List of parameters for the starting array
    error_model = True  # Parameter whether to use error model
    estimates = []  # The last parameter which has been fit
    summary = 0  # Full summary of the last fitted results is saved here
    output = True
    
    def __init__(self, endog=[], exog=[], start_params=[1000,],
                 min_len=4, max_len = 20.0,
                 error_model=False, output=False, 
                 chr_lgts=[], **kwds):
        '''endog: List of n individual ROH lists [roh_list1, ..., roh_listn]
        output: Whether to plot detailled output strings
        start_params: Where to start the fit
        min_len, max_len: Bin edges of ROH to fit.
        '''
        if len(exog)==0:
            exog = np.zeros(len(endog), dtype="int8") # Just empty placeholder
        self.initialize_ll_model(endog=endog, exog=exog, **kwds)
        self.set_params(start_params=start_params, error_model=error_model, output=output,
                        min_len=min_len, max_len=max_len, chr_lgts=chr_lgts)
        
        self.create_bins()  # Create the Mid Bin vector

        if self.error_model == True:  # In case required: 
            self.fp_rate = self.fp_rate(self.mid_bins) * self.bin_width  # Calculate the false positives per bin
            self.calculate_trans_mat()  # Calculate the Transformation matrix
        
        if self.output:
            print("Initial Minimum Length analyzed [cM]: %.4f" % self.min_len)
            print("Initial Maximum Length analyzed [cM]: %.4f" % self.max_len)
            print("Successfully initialized likelihood model")
        
    def initialize_ll_model(self, **kwds):
        '''Function to initialize LL model'''
        super(MLE_ROH_Ne, self).__init__(**kwds)  # Create the full object.
            
    def set_params(self, recalculate_bins=False, **kwargs):
        """Set the Parameters.
        Takes keyworded arguments"""
        for key, value in kwargs.items():
            setattr(self, key, value)
            
        if recalculate_bins:
            self.create_bins()  # Actually create the Bins
        
            # If needed calculate the false positive Rate:
            if self.error_model:
                self.fp_rate = fp_rate(self.mid_bins) * self.bin_width  # Calculate the false positives per bin
                self.calculate_trans_mat()  # Calculate the Transformation matrix
    
    ######################################################################################

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        """Fit Parameters via fit method of statsmodels"""
        if start_params == None:
            start_params = self.start_params  # Set the starting parameters for the fit
        fit = super(MLE_ROH_Ne, self).fit(start_params=start_params,
                                     maxiter=maxiter, maxfun=maxfun,
                                     **kwds)
        self.estimates = fit.params ### Save the estimates in field
        #self.summary = fit.summary()
        self.summary = self.get_summ_as_df(summary=fit.summary())
        return fit
    
    def fit_ll_profile(self, ns=[],level=1.92,num=1000):
        """Fit Parameters via likelihood profile.
        Only works for 1D search.
        Return lower, upper Ne supported by 
        ns: N values to test.
        level: loglikelihood difference to use."""

        if len(ns)==0:
            ns = np.logspace(2,5,num=num) # Default parameters
            
        ### Calculate Log Likelihoods
        lls = [self.loglike(params=[n]) for n in ns]
        i = np.argmax(lls)  # Get Maximum
        m = ns[i]   # Maximum N
        l = lls[i] - level  # The ll line for thresholding
        high = np.where(lls>=l)[0] # All values above threshold
        c0, c1 = ns[high[0]], ns[high[-1]]
        
        ### Save the result as dataframe
        df_res = get_default_res()
        df_res["coef"]=m
        df_res["0.025"]=c0
        df_res["0.975"]=c1
        df_res["n"]= len(self.endog)
        self.summary = df_res  # Set the result dataframe
        return df_res
    
    ######################################################################################
    
    def loglikeobs(self, params):
        '''Return vector of log likelihoods for every observation. (here pairs of pops)'''
        if self.output:
            for i in range(len(params)):
                print("Parameter %.0f : %.8f" % (i, params[i]))
        n_e = params[0]  # Absolute Parameter

        if n_e <= 0:  # If Parameters outside bound return infinitely negative likelihood
            return -np.ones(len(self.endog)) * (np.inf)
        
        ### Reccalculate the sharing probabilities for new Parameters
        self.calculate_thr_shr(params)  # Calculate theoretical sharing PER PAIR
        self.calculate_full_bin_prob()  # Calculate total sharing PER PAIR /TM Matrix and FP-rate dont need update
        
        ll = [self.ll_individual(self.endog[i], params) for i in range(len(self.endog))]
        if self.output:
            print("Total log likelihood: %.4f" % np.sum(ll))
        return np.array(ll).astype('float')  # Return negative log likelihood
    
    def ll_individual(self, roh, params):
        '''Log likelihood function for every raw of data (sharing between countries).
        Return log likelihood.'''
        roh = np.array(roh)  # Make l an Numpy vector for better handling
        
        bins = self.mid_bins[self.min_ind:self.max_ind + 1] - 0.5 * self.bin_width  # Rel. bin edges
        roh = roh[(roh >= bins[0]) * (roh <= bins[-1])]  # Cut out only blocks of interest
        
        ### Extract the sharing Probability
        shr_pr = self.full_shr_pr[self.min_ind:self.max_ind]
        
        log_pr_no_shr = -np.sum(shr_pr)   # The negative sum of all total sharing probabilities
        if len(roh) > 0:
            indices = np.array([(bisect_left(bins, l) - 1) for l in roh])  # Get indices of all shared blocks
            l1 = np.sum(np.log(shr_pr[indices]))
        else: l1 = 0
        ll = l1 + log_pr_no_shr
        return(ll)    
    
    ##################################################################################
    ##################################################################################
    ### Calculation bins
    
    def create_bins(self):
        '''Creates the bins according to parameters in Class.
        Also set the false positive rate'''
        bins = np.arange(self.min_b, self.max_b, self.bin_width)  # Create the actual bins
        
        self.min_ind = bisect_left(bins, self.min_len)  # Find the indices of the relevant points
        self.max_ind = bisect_left(bins, self.max_len)
        self.mid_bins = bins + 0.5 * self.bin_width
        k = len(self.mid_bins)
        
        if self.error_model == True:
            self.trans_mat = np.zeros((k, k)).astype(float)  # Create empty transition matrix
            self.fp_rate = fp_rate(self.mid_bins) * self.bin_width  # Calculate the false positives per bin
        
    def calculate_thr_shr(self, params):
        '''Calculates the expected Bessel-Decay per bin''' 
        bd = self.block_shr_density(self.mid_bins, params)
        self.theoretical_shr = bd * self.bin_width  # Integrate over bin width (Morgan tranform)
        
    def calculate_trans_mat(self):
        '''Calculate the transition matrix from true estimated to
        observed values for block sharing.'''
        k = len(self.trans_mat)
        for i in range(k):  # Iterate over all starting values
            x = self.mid_bins[i]
            pr_detect = (1 - censor_prob((x)))  # Probability of detecting block
            for j in range(0, i):
                y = self.mid_bins[j]
                # Down probability conditional on bigger than cut off:
                trans_pr = prob_down(x) * down_rate(x) * np.exp(-down_rate(x) * (x - y)) / (1 - np.exp(-down_rate(x) * (x - 1)))
                self.trans_mat[j, i] = pr_detect * trans_pr * self.bin_width 
                              
            for j in range(i + 1, k):
                y = self.mid_bins[j]
                trans_pr = (1 - prob_down(x)) * up_rate(x) * np.exp(-up_rate(x) * (y - max(x, 1)))     
                self.trans_mat[j, i] = pr_detect * trans_pr * self.bin_width
            
            y = x  # Now do the i,i case
            trans_pr_d = prob_down(x) * down_rate(x) * np.exp(-down_rate(x) * (x - y)) / (1 - np.exp(-down_rate(x) * (x - 1)))
            # trans_pr_u = (1 - prob_down(x)) * up_rate(x) * np.exp(-up_rate(x) * (y - max(x, 1)))
            trans_pr_u = (1 - prob_down(x)) * up_rate(x) * np.exp(-up_rate(x) * (y - x))
            self.trans_mat[i, i] = pr_detect * 1 / 2.0 * (trans_pr_d + trans_pr_u) * self.bin_width  # Prob of not going anywhere         
        
    def calculate_full_bin_prob(self):
        '''Calculate the full probablities per bin'''
        # Transition matrix times theoretically expected + false positives.
        if self.error_model:
            self.full_shr_pr = (np.dot(self.trans_mat, self.theoretical_shr) + self.fp_rate)  # Model with full error
        else:
            self.full_shr_pr = self.theoretical_shr  # Model without any error in detection
    
    ##################################################################################
    ### Output Functions
    
    def get_bl_shr_interval(self, interval, r, params=[]):
        '''Return the estimated block-sharing under the model in interval given distance r.
        For this use all the bins intersecting the interval and average.
        Assumes r is array and return array'''
        if len(params) == 0:  # If no parameters passed use last ones fit
            params = self.estimates
        ### Find the indices of the right ultimate bins:
        bins = self.mid_bins - 0.5 * self.bin_width
        ind = bisect_right(bins, interval[0])
        ind1 = bisect_left(bins, interval[1])
        
        estims = np.array([0.0 for i in r])
        for i in range(len(r)):
            self.calculate_thr_shr(r[i], params)  # Calculate the theoretical sharing
            self.calculate_full_bin_prob()  # Calculate the bin probability of sharing a block                
            ### Do numerical "Integral":
            mean_value = np.mean([self.full_shr_pr[j] for j in range(ind - 1, ind1 + 1)])
            estims[i] = mean_value * (interval[1] - interval[0]) / self.bin_width  # Normalize
        return(estims) 
    
    def print_block_nr(self):
        '''Calculate and print the Nr of analyzed Blocks of right Length.'''
        min_l, max_l = self.min_len, self.max_len
        fit_sharing = [np.sum((np.array(l) >= min_l) * (np.array(l) <= max_l)) for l in self.endog]
        bl_nr = np.sum(fit_sharing)
        if self.output:
            print("Min Length: %s" % min_l)
            print("Max. Length: %s" % max_l)
            print("Total Nr. of Blocks within Min and Max Length: %i" % bl_nr)
        return bl_nr
    
    def ci_li(self, ns=[],level=1.92, update=True):
        """Get confidence Intervall based on likelihood ration test.
        Return lower, upper Ne supported by 
        ns: N values to test.
        level: loglikelihood difference to use.
        update: Whether to update Confidence Intervals in the fit"""
        if len(ns)==0:
            ns = np.logspace(2,5,num=100)
        lls = [self.loglikeobs(params=[n]) for n in ns]
        l = np.max(lls) - level  # The line to find
        high = np.where(lls>=l)[0] # All values above threshold
        c0, c1 = ns[high[0]], ns[high[-1]]
        
        if update:
            self.summary["0.025"] = c0
            self.summary["0.975"] = c1
        return c0, c1
    
    ##################################################################################
    ### Block Density Function (Overwrite that for other Models)
    
    def block_shr_density(self, l, params):
        '''Calculate and print the Nr of analyzed Blocks of right Length.'''
        n_e = params[0]  # Population size (n_e = 2Ne)
        l = l / 100 # Switch to Morgan (for formula)
        
        es = Expected_Roh() # Load expected ROH class containing function
        chr_lgts=self.chr_lgts
        if len(chr_lgts)>0:
            chr_lgts=np.array(chr_lgts)/100 # Transform to centimorgan
        y = es.roh_pdf_allchr_N(x=l, N=n_e, chr_lgts=chr_lgts) / 100 # Calculate density per Centimorgan
        return y
    
    ##################################################################################
    ### Output function
    
    def get_summ_as_df(self, summary):
        """Transforms statsmodels summary to pandas dataframe.
        If no summary given use the last saved one
        Return that dataframe"""
        #if not summary:
        #    summary = self.summary
        results_as_html = summary.tables[1].as_html()
        df = pd.read_html(results_as_html, header=0, index_col=0)[0]
        df.rename(columns={'[0.025':'0.025', "0.975]":"0.975"}, inplace=True)
        df.drop(columns=["z","P>|z|"], inplace=True)
        return df
    
############################################################    
############################################################
### Method to simulate ROH vecs under the model

def ind_roh_from_pdf(ne=500, bin_range=[0.04, 0.5], 
                     nbins=100000, output=False):
    """Create ROH for single individual. Everything measured in Morgan.
    Return array of ROH lengths, in Morgan"""
    es = Expected_Roh()
    x = np.linspace(bin_range[0], bin_range[1], nbins)
    bw = x[1] - x[0]
    y = es.roh_pdf_allchr_N(x=x, N=ne)
    p = y * bw  # Probability Vector

    realized = np.random.random(len(p))<p   # Boolean Vector
    l = x[realized]
    
    if output:
        print(f"Created Blocks: {np.sum(realized)}")
        print(f"Maximum Block Length: {np.max(l):.4f}")
        print(f"Mean Block Length: {np.mean(l):.4f}")
        print(f"Total Block Length: {np.sum(l):.4f}")
    return l

def inds_roh_from_pdf(n_ind=5, ne=500, bin_range=[0.04, 0.5], 
                      nbins=100000, output=False, cm=True):
    """Create ROH for single individual. Everything measured in Morgan.
    Return array of ROH lengths, in Morgan. 
    cm: If true, output in centimorgan"""
    roh_vec=[]
    for _ in range(n_ind):
        l = ind_roh_from_pdf(ne=ne, bin_range=bin_range, 
                             nbins=nbins, output=output)
        
        if cm:
            l=l*100
        roh_vec.append(l)
    return roh_vec

############################################################
### Helper Functions related to Ne estimation.

def load_roh_vec(iids=[], base_path = "./output/roh/", suffix="_roh_full.csv"):
    """Load and return ROH length vector """
    paths = give_iid_paths(iids, base_folder=base_path, suffix=suffix)
    roh_dfs = [pd.read_csv(p) for p in paths]
    roh_vec = [df["lengthM"].values*100 for df in roh_dfs]
    return roh_vec

def get_default_res():
    """REturn default line of results dataframe.
    Used for unfittable scenarios"""
    dct= {"coef":[np.nan],
          "std err":[np.nan],
          "0.025":[np.nan],
          "0.975":[np.nan],
          "n":[np.nan]}
    df = pd.DataFrame(dct)
    return df