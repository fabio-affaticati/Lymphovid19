# Copyright © 2024 fabio-affaticati

import pandas as pd
import os
import pathlib

from src.utils import plotting
from src.utils import tools

DATADIR = os.path.join(pathlib.Path().absolute(), 'data/')
RESULTSDIR = os.path.join(pathlib.Path().absolute(), 'results/')
TASKSDIR = os.path.join(pathlib.Path().absolute(), 'tcrex_tasks/')
PLOTSDIR = os.path.join(pathlib.Path().absolute(), RESULTSDIR+'plots/')
TCREXDIR = '/Users/fabioaffaticati/Desktop/Work/TCRex/data/task'


if __name__ == "__main__":
    
    ### DANGER: This is a data cleaning step that removes files in the directory
    ####################################################################################################################################
    tools.clean_directory(TASKSDIR)
    tools.copy_most_recent_directory(TCREXDIR, TASKSDIR)
    ####################################################################################################################################
    ####################################################################################################################################
    for preds in os.listdir(TASKSDIR):
        preds_path = os.path.join(TASKSDIR, preds)
    tcrex_preds = pd.read_csv(os.path.join(TASKSDIR, preds_path+'/predictions.tsv'), sep = '\t', comment='#')
    
    # fixed the notation of the TRBV genes
    tcrex_preds['TRBV_gene'] = tcrex_preds['TRBV_gene'].str.replace('BV0', 'BV')
    tcrex_preds['TRBV_gene'] = tcrex_preds['TRBV_gene'].str.replace('-0', '-')
    tcrex_preds['clonotype'] = tcrex_preds['TRBV_gene'] + '_' + tcrex_preds['CDR3_beta']
    ####################################################################################################################################
    
    
    
    data = pd.read_csv(RESULTSDIR + 'preprocessed_data.csv', index_col=0)
    metadata = pd.read_excel(DATADIR+'metadata.xlsx')
    metadata.rename(columns={'SAMPLE_ID': 'sample_id'}, inplace=True)
    metadata['sample_id'] = metadata['sample_id'].astype(str)
    
    data['length'] = data['junction_aa'].str.len()
    data['clonotype'] = data['v_call'] + '_' + data['junction_aa']
    
    
    
    tcrex_merged = pd.merge(tcrex_preds, data.query('CONDITION == "Lymphomas"'), on='clonotype', how='inner')[data.columns]
    # plot in seaborn the distribution of the cloneFraction based on TIMEPOINTS
    import seaborn as sns
    import matplotlib.pyplot as plt
    # sort the x axis labels
    tcrex_merged['TIMEPOINTS'] = pd.Categorical(tcrex_merged['TIMEPOINTS'], categories=['baseline', 'V1', 'V3',], ordered=True)
    sns.boxplot(data=tcrex_merged, y = 'cloneFraction', x = 'TIMEPOINTS', hue = 'CONDITION', dodge= False, width=0.4, showfliers=False)
    # change y axis scale to log
    #plt.yscale('log')
    plt.tight_layout()
    plt.title('TCRex predicted COVID-19 specific ratios at different timepoints')
    plt.savefig(PLOTSDIR + 'tcrexpreds.png', dpi=600, bbox_inches='tight')
    exit()
    
    
    
    
    
    
    
    
    ####################################################################################################################################
    #plotting.barplot_unique_sequences(data, metadata, PLOTSDIR)
    ####################################################################################################################################
    
    ####################################################################################################################################
    clone_fractions_pivoted = data.pivot_table(index=['SAMPLE', 'TIMEPOINTS'], columns='junction_aa', values='cloneFraction', fill_value=0.0000001, aggfunc='sum')    
    #plotting.scatter_clone_fractions(clone_fractions_pivoted, PLOTSDIR)
    ####################################################################################################################################
    
    ####################################################################################################################################
    clustcr_data = pd.read_csv(RESULTSDIR + 'total_lymphoma_v_healthy_results.tsv', sep = '\t', index_col=0)
    # clustcr_data = clustcr_data.query('log10_corrected_p_value_fdr_bh > 0')
    clustcr_data.reset_index(drop=True, inplace=True)
    clustcr_data['chain'] = clustcr_data['v_call'].str[:3]
    clustcr_data['dot_size'] = clustcr_data.apply(plotting.calculate_dot_size, axis=1)
    #plotting.volcano_plot(clustcr_data, 'Lymphoma against Healthy', 'lvh', PLOTSDIR)
    ####################################################################################################################################
    
    ####################################################################################################################################
    pairs=[("Healthy", "Lymphomas")]
    clone_fractions_baseline_pivoted = data[data['TIMEPOINTS'] == 'baseline']
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted.pivot_table(index=['sample_id', 'CONDITION'], columns='junction_aa', values='cloneFraction', fill_value=0, aggfunc='sum')
    plotting.alpha_diversity_tcr(clone_fractions_baseline_pivoted, pairs, PLOTSDIR)
    ####################################################################################################################################
    
    
    
    def vgene_enrich_test(data, poscol, negcol):

        for chain in ['TRA', 'TRB', 'TRG', 'TRD']:
            
            # select columns that start with chain
            gene_fractions_pivoted_chain = data[[col for col in data.columns if col.startswith(chain)]]
            res = [mannwhitneyu(gene_fractions_pivoted_chain[gene_fractions_pivoted_chain.index.get_level_values(1) == poscol][vgene],
                                gene_fractions_pivoted_chain[gene_fractions_pivoted_chain.index.get_level_values(1) == negcol][vgene])[1]
                                for vgene in gene_fractions_pivoted_chain.columns]
            print(res)
            print(multipletests(res, method='fdr_bh')[1])
            
    
    
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    gene_fractions_pivoted =  data[data['TIMEPOINTS'] == 'baseline'].pivot_table(index=['sample_id', 'CONDITION'], columns='v_call', values='cloneFraction', fill_value=0, aggfunc='sum')
    #vgene_enrich_test(gene_fractions_pivoted, 'Lymphomas', 'Healthy')
