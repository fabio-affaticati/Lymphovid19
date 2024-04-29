import warnings
from typing import Union, Callable

import pandas as pd
import numpy as np
import parmap

from statsmodels.stats.multitest import multipletests
from scipy.stats import hypergeom, fisher_exact

def hypergeometric_test(x:int, n:int, M:int, N:int) -> float:
    """
    X: Number of hits in the sample
    M: population size
    n: number of hits in the population
    N: sample size

    Returns: probability of drawing X-1 or more hits from the population
    """
    return hypergeom.sf(x-1, M, n, N)


def fisher_exact_test(x:int, n:int, M:int, N:int) -> float:
    """
    x: Number of hits in the sample
    M: Population size
    n: Number of hits in the population
    N: Sample size

    Returns: Two-tailed p-value for Fisher's exact test
    """
    # Construct the 2x2 contingency table
    table = [[x, n - x], [N - x, M - N - n + x]]
    _, p_value = fisher_exact(table)
    return p_value


class Enrichment():
    """
    A class to perform enrichment analysis using a specified statistical test
    and multiple testing correction method.
    """

    def __init__(self, test: Union[str, Callable] = 'hypergeometric', mtc_method="fdr_bh"):
        """
        Initializes the Enrichment class.

        Parameters:
        -----------
        test : Callable, optional
            The statistical test to be used for enrichment analysis. Default is
            hypergeometric_test.
        mtc_method : ["fdr_bh", "bonferroni", ...], optional
            The method used for multiple testing correction. Default is
            "fdr_bh".
        """
        self.test = self._get_test_function(test)
        self.mtc_method = mtc_method
        return None
    
    def fit(self, y):
        """
        Fits the model with the given data.

        Parameters:
        -----------
        y : array-like
            The input counts, where the first column represents the positive
            class counts, and the second column the counts in the negative
            class.

        Returns:
        --------
        self : Enrichment
            The fitted model.
        """
        self.xs_ = y[:, 0]
        self.ns_ = y.sum(axis=1)
        self.M_ = y.sum()
        self.N_ = self.xs_.sum()
        self._validate()
        return self

    def fit_df(self, df, positive_class: Union[str, list]):
        """
        Fits the model with the given data.

        Parameters:
        -----------
        df : DataFrame
            The input counts in dataframe format. Provide only the subset of the
            columns with counts.
        positive_class : Union[str, list]
            The column(s) in the DataFrame that represent the positive class.
            All other columns will be considered part of the negative group,
            hence no additional columns should be present.

        Returns:
        --------
        self : Enrichment
            The fitted model.
        """
        if isinstance(positive_class, str):
            positive_class = [positive_class]

        self.xs_ = df[positive_class].sum(axis=1).to_numpy()
        self.ns_ = df.sum(axis=1).to_numpy()
        self.M_ = self.ns_.sum()
        self.N_ = self.xs_.sum()
        self.ids = df.index.to_numpy()
        self._validate()
        return self
    
    def _multiple_testing_correction(self, pvalues):
        """
        Applies multiple testing correction to the given p-values.

        Parameters:
        -----------
        pvalues : array-like
            The p-values to be corrected.

        Returns:
        --------
        corrected_pvalues : array-like
            The corrected p-values.
        """
        corrected_pvalues = multipletests(pvalues, method=self.mtc_method)[1]
        return corrected_pvalues


    def transform(self, y=None):
        """
        Transforms the data to calculate p-values, using multiprocessing
        implemented in parmap.

        Returns:
        --------
        pvals : array-like
            The calculated p-values.
        """
        pvals = parmap.starmap(self.test, zip(self.xs_, self.ns_), self.M_, self.N_)
        return np.array(pvals)

    def transform_df(self, multiple_testing_correction:bool=True):
        """
        Transforms the data to a DataFrame with calculated p-values.

        Parameters:
        -----------
        multiple_testing_correction : bool, optional
            Whether to apply multiple testing correction to the p-values.
            Default is True.

        Returns:
        --------
        df : DataFrame
            The DataFrame with calculated p-values.
        """
        pvals = self.transform()
        df = pd.DataFrame({
            "positive": self.xs_,
            "negative": self.ns_ - self.xs_,
            "p_value": pvals,
        })
        df.index = self.ids if self.ids is not None else df.index
        if multiple_testing_correction:
            df[f"corrected_p_value_{self.mtc_method}"] = self._multiple_testing_correction(df["p_value"])
        return df

    def _validate(self):
        if self.M_ >= 2**32:
            warnings.warn(f"""Population size exceeds {2**32-1} (32 bits).
                          p-values may be miscalculated.""")
            
    def _get_test_function(self, test: Union[str, Callable]) -> Callable:
        """
        Returns the corresponding statistical test function based on the input string or callable.

        Parameters:
        -----------
        test : Union[str, Callable]
            The statistical test to be used for enrichment analysis.

        Returns:
        --------
        function : Callable
            The corresponding statistical test function.
        """
        if callable(test):
            return test
        elif test == 'hypergeometric':
            return hypergeometric_test
        elif test == 'fisher':
            return fisher_exact_test
        # Add more elif conditions here if you have other tests
        else:
            raise ValueError(f"Invalid test: {test}. Allowed values are 'hypergeometric', or a callable function.")
        
        
def calculate_enrichment_features(enr):
    
    pseudocount = 0.01
    
    enr["log10_p_value"] = -enr["p_value"].apply(np.log10)
    enr["log10_corrected_p_value_fdr_bh"] = - enr["corrected_p_value_fdr_bh"].apply(np.log10)

    enr["positive_frac"] = (enr["positive"] + pseudocount) / enr["positive"].sum()
    enr["negative_frac"] = (enr["negative"] + pseudocount) / enr["negative"].sum()

    enr["fold_change"] = enr["positive_frac"]/enr["negative_frac"]
    enr["log2_fold_change"] = enr["fold_change"].apply(np.log2)
    return enr
        
        
def enrichment_test(res, cluster_count_df : pd.DataFrame, processed_data: pd.DataFrame, test_type: str,
                    minlog2fc: float, condition_col: str, pos_class: str, results_dir: str, file_name: str) -> None:
    
    
    positive_patients = list(set(processed_data[processed_data[condition_col]==pos_class]['sample_id']))
    print(f"Unique conditions {processed_data[condition_col].unique()}")
    print(f"Positive patients {positive_patients}")
    
    enr = Enrichment(test=test_type) # hypergeometric (one-sided) or fisher (two-sided)
    enr = enr.fit_df(cluster_count_df, positive_class=positive_patients) # all other columns are seen as the negative group
    enr = enr.transform_df()
    
    enr = calculate_enrichment_features(enr)
    
    
    if minlog2fc > 0:
        ### Establish a log fold change treshold
        enr = enr.query('log2_fold_change > @minlog2fc or log2_fold_change < -@minlog2fc')
        cluster_count_df_filtered = cluster_count_df.loc[cluster_count_df.index.astype(str).isin(list(enr.index))]

        enr = Enrichment(test=test_type)
        enr = enr.fit_df(cluster_count_df_filtered, positive_class=positive_patients) # all other columns are seen as the negative group
        enr = enr.transform_df()

        enr = calculate_enrichment_features(enr)
        
    
    enr.reset_index(names='cluster', inplace=True)
    enr['cluster'] = enr['cluster'].astype('int')

    # remove not clustered
    enr = enr.query('cluster != -1')
    
    # merge to obtain the motifs
    motifs_df = pd.merge(res.clusters_df, res.summary().reset_index().rename(columns={'index': 'cluster'}), on = 'cluster')

    pd.merge(motifs_df, enr, on='cluster').to_csv(results_dir+'total_'+file_name+'_results.tsv', sep ='\t', index = False)
    
    corrected_res = enr.query('corrected_p_value_fdr_bh < 0.05')
    pd.merge(motifs_df, corrected_res, on='cluster').to_csv(results_dir+file_name+'_results.tsv', sep ='\t', index = False)
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
