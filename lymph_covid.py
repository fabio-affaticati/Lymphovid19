# Copyright © 2024 fabio-affaticati

import pandas as pd
import os
import pathlib
import seaborn as sns

from src.utils import plotting
from src.utils import tools

DATADIR = os.path.join(pathlib.Path().absolute(), 'data/')
RESULTSDIR = os.path.join(pathlib.Path().absolute(), 'results/')
TASKSDIR = os.path.join(pathlib.Path().absolute(), 'tcrex_tasks/')
PLOTSDIR = os.path.join(pathlib.Path().absolute(), RESULTSDIR+'plots/')
TCREXDIR = '/Users/fabioaffaticati/Desktop/Work/TCRex/data/task'


if __name__ == "__main__":
    
    for directory in [DATADIR, RESULTSDIR, TASKSDIR, PLOTSDIR]:
        if not os.path.exists(directory):
            os.makedirs(directory)
            
    # Load preprocess data and metadata excel files
    data = pd.read_csv(RESULTSDIR + 'preprocessed_data.csv', index_col=0, low_memory=False)
    metadata = pd.read_excel(DATADIR+'metadata.xlsx')
    metadata.rename(columns={'SAMPLE_ID': 'sample_id'}, inplace=True)
    metadata['sample_id'] = metadata['sample_id'].astype(str)
    
    data['length'] = data['junction_aa'].str.len()
    data['clonotype'] = data['v_call'] + '_' + data['junction_aa'] + '_' + data['j_call'] 
    
    
    ####################################################################################################################################
    #plotting.barplot_unique_sequences(data, metadata, PLOTSDIR)
    ####################################################################################################################################

    
    ### DANGER: This is a data cleaning step that removes files in the directory
    ####################################################################################################################################
    tools.clean_directory(TASKSDIR)
    tools.copy_most_recent_directory(TCREXDIR, TASKSDIR)
    ####################################################################################################################################
    ####################################################################################################################################
    for preds in os.listdir(TASKSDIR):
        preds_path = os.path.join(TASKSDIR, preds)
    tcrex_preds = pd.read_csv(os.path.join(TASKSDIR, preds_path+'/predictions.tsv'), sep = '\t', comment='#')
    
    
    # fixed the notation of the TRBV and TRBJ genes
    tcrex_preds['TRBV_gene'] = tcrex_preds['TRBV_gene'].str.replace('BV0', 'BV')
    tcrex_preds['TRBV_gene'] = tcrex_preds['TRBV_gene'].str.replace('-0', '-')
    tcrex_preds['TRBJ_gene'] = tcrex_preds['TRBJ_gene'].str.replace('BJ0', 'BJ')
    tcrex_preds['TRBJ_gene'] = tcrex_preds['TRBJ_gene'].str.replace('-0', '-')
    tcrex_preds['clonotype'] = tcrex_preds['TRBV_gene'] + '_' + tcrex_preds['CDR3_beta'] + '_' + tcrex_preds['TRBJ_gene']

    
    
    
    
    
    
    
    
    
    
    
    ####################################################################################################################################
    tcrex_merged = pd.merge(tcrex_preds.drop_duplicates(subset=['clonotype']), data, on='clonotype', how='inner')[data.columns] 
    tcrex_merged['TIMEPOINTS'] = pd.Categorical(tcrex_merged['TIMEPOINTS'], categories=['baseline', 'V1', 'V3',], ordered=True)
    tcrex_merged.drop_duplicates(subset=['sample_id', 'clonotype', 'TIMEPOINTS'], inplace=True)
    tcrex_merged.reset_index(drop=True, inplace=True)
    tcrex_merged.dropna(inplace=True)
    #plotting.plot_tcrex_predictions(tcrex_merged, metadata, 'TIMEPOINTS', 'CONDITION', PLOTSDIR)
    ####################################################################################################################################
    
    ####################################################################################################################################
    # per epitope analysis
    tcrex_merged = pd.merge(tcrex_preds, data, on='clonotype', how='inner')#[data.columns] 
    tcrex_merged['TIMEPOINTS'] = pd.Categorical(tcrex_merged['TIMEPOINTS'], categories=['baseline', 'V1', 'V3',], ordered=True)
    tcrex_merged.dropna(inplace=True)
    #plotting.plot_tcrex_predictions_epitopes(tcrex_merged, metadata, 'sample_id', 'epitope', PLOTSDIR)
    ####################################################################################################################################
    
    
    
    ####################################################################################################################################
    # keep only beta chains
    beta_chain_data = data[data['TCR_Chain'].str.contains('TRB')]
    print(beta_chain_data)
    print(beta_chain_data.drop_duplicates(subset=['clonotype', 'sample_id']).shape)
    scatteratio_data = pd.merge(beta_chain_data, tcrex_preds.drop_duplicates(subset=['clonotype']), on='clonotype', how='left', indicator=True)
    
    #scatteratio_data = scatteratio_data.query('TIMEPOINTS == "V1"')
    scatteratio_data = scatteratio_data[['_merge', 'clonotype', 'sample_id', 'CONDITION', 'TIMEPOINTS', 'cloneFraction']]
    # groupby both sample id and CONDITION
    scatteratio_data = scatteratio_data.groupby(['sample_id', 'CONDITION', 'TIMEPOINTS'], as_index=False).apply(lambda x: x['_merge'].value_counts())
    scatteratio_data['totalBeta'] = scatteratio_data['left_only'] + scatteratio_data['both']
    #plotting.plot_scatteratio_tcr_specific(scatteratio_data, PLOTSDIR)
    
    scatteratio_data['fraction_betas'] = scatteratio_data['both']/scatteratio_data['totalBeta']
    plotting.plot_scatteratio_breadth(scatteratio_data, metadata, 'TIMEPOINTS', 'CONDITION', PLOTSDIR)
    ####################################################################################################################################

    exit(0)
    
    
    
    ####################################################################################################################################
    #covid_spec_data = beta_chain_data.query('_merge == "both"')
    covid_spec_data = beta_chain_data.copy()
    covid_spec_data.reset_index(drop=True, inplace=True)
    covid_spec_data['presence'] = 1
    # group by clonotype and condition and sum presence values
    covid_spec_data = covid_spec_data.groupby(['clonotype', 'CONDITION', 'TIMEPOINTS'], as_index=False).agg({'presence': 'sum'})
    covid_spec_data['C_T'] = covid_spec_data['CONDITION'] + '_' + covid_spec_data['TIMEPOINTS']
    
    covid_spec_data = covid_spec_data.pivot_table(index='clonotype', columns='C_T', fill_value=0, values='presence', aggfunc='sum')
    covid_spec_data.reset_index(inplace=True)
    
    covid_spec_data = pd.merge(covid_spec_data, tcrex_preds, on='clonotype', how='left', indicator=True)
    covid_spec_data.rename(columns={'_merge': 'tcrexpred'}, inplace=True)
    covid_spec_data['tcrexpred'] = covid_spec_data['tcrexpred'].map({'left_only': 'NotSpecific', 'both': 'CovidSpecific'})
    covid_spec_data = covid_spec_data[[col for col in covid_spec_data.columns 
                                    if col.startswith('Healthy') or col.startswith('Lymphomas') or col == 'clonotype' or col == 'tcrexpred']]    
    #plotting.plot_scatter_presence(covid_spec_data, PLOTSDIR)
    ####################################################################################################################################

    

    ####################################################################################################################################
    clone_fractions_pivoted = data.pivot_table(index=['SAMPLE', 'TIMEPOINTS'], columns='clonotype', values='cloneFraction', fill_value=0.0000001, aggfunc='sum')
    # create a pd.Series starting from clone_fractions_pivoted.columns that has value 'black' if the value is in test['clonotype'] and 'gray' otherwise
    tcrex_clonotypes = set(tcrex_preds['clonotype'])
    colors = ['red' if col in tcrex_clonotypes else sns.color_palette(palette='Greys')[1] for col in clone_fractions_pivoted.columns]
    plotting.scatter_clone_fractions(clone_fractions_pivoted, colors, PLOTSDIR)
    ####################################################################################################################################
    

    ####################################################################################################################################
    clustcr_data = pd.read_csv(RESULTSDIR + 'total_lymphoma_v_healthy_results.tsv', sep = '\t', index_col=0)
    # clustcr_data = clustcr_data.query('log10_corrected_p_value_fdr_bh > 0')
    clustcr_data.reset_index(drop=True, inplace=True)
    clustcr_data['chain'] = clustcr_data['v_call'].str[:3]
    clustcr_data['dot_size'] = clustcr_data.apply(plotting.calculate_dot_size, axis=1)
    #plotting.volcano_plot(clustcr_data, 'Lymphoma - Healthy clusters contrast at baseline', 'lvh', PLOTSDIR)
    ####################################################################################################################################
    
    ####################################################################################################################################
    pairs=[("Healthy", "Lymphomas")]
    clone_fractions_baseline_pivoted = data[data['TIMEPOINTS'] == 'baseline']
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted[clone_fractions_baseline_pivoted['TCR_Chain'].str.contains('TRB')]
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted.pivot_table(index=['sample_id', 'CONDITION'], columns='junction_aa', values='cloneFraction', fill_value=0, aggfunc='sum')
    plotting.alpha_diversity_tcr(clone_fractions_baseline_pivoted, pairs, PLOTSDIR)
    ####################################################################################################################################
    
    exit(0)
    
    
    
    
    
    
    
    
    
    
    def vgene_enrich_test(data, poscol, negcol):

        for chain in ['TRA', 'TRB', 'TRG', 'TRD']:
            
            # select columns that start with chain
            gene_fractions_pivoted_chain = data[[col for col in data.columns if col.startswith(chain)]]
            res = [mannwhitneyu(gene_fractions_pivoted_chain[gene_fractions_pivoted_chain.index.get_level_values(1) == poscol][vgene],
                                gene_fractions_pivoted_chain[gene_fractions_pivoted_chain.index.get_level_values(1) == negcol][vgene])[1]
                                for vgene in gene_fractions_pivoted_chain.columns]
            
            corrected = multipletests(res, method='fdr_bh')[1]
            corrected = list(zip(corrected, gene_fractions_pivoted_chain.columns))
            print(corrected)
            corrected = [x for x in corrected if x[0] < 0.75]
            
    
    
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    gene_fractions_pivoted =  data[data['TIMEPOINTS'] == 'baseline'].pivot_table(index=['sample_id', 'CONDITION'], columns='v_call', values='cloneFraction', fill_value=0, aggfunc='sum')
    print(gene_fractions_pivoted)
    vgene_enrich_test(gene_fractions_pivoted, 'Lymphomas', 'Healthy')
