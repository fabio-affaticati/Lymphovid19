# Copyright © 2024 fabio-affaticati

import pandas as pd
import os
import pathlib
import seaborn as sns

from src.utils import plotting
from src.utils import tools

DATADIR = os.path.join(pathlib.Path().absolute(), 'data/')
RESULTSDIR = os.path.join(pathlib.Path().absolute(), 'results/')
PLOTSDIR = os.path.join(pathlib.Path().absolute(), RESULTSDIR+'plots/')


if __name__ == "__main__":
    
    for directory in [DATADIR, RESULTSDIR, PLOTSDIR]:
        if not os.path.exists(directory):
            os.makedirs(directory)
            
    # Load preprocess data and metadata excel files
    data = pd.read_csv(RESULTSDIR + 'preprocessed_data.csv', index_col=0, low_memory=False)
    metadata = pd.read_excel(DATADIR+'metadata.xlsx')
    metadata.rename(columns={'SAMPLE_ID': 'sample_id'}, inplace=True)
    metadata['sample_id'] = metadata['sample_id'].astype(str)
    

    metadata['AB_TITER'] = metadata['AB_TITER'].astype(str)
    metadata['AB_TITER'] = metadata['AB_TITER'].str.replace('>400', '400')
    metadata['AB_TITER'] = metadata['AB_TITER'].str.replace('<3,80', '3.80')
    metadata['AB_TITER'] = metadata['AB_TITER'].str.replace('<3.80', '3.80')
    metadata['AB_TITER'] = metadata['AB_TITER'].str.replace('V2 ', '')
    metadata['AB_TITER'] = metadata['AB_TITER'].str.replace(',', '.')
    metadata['AB_TITER'] = metadata['AB_TITER'].astype(float)
    

    data['length'] = data['junction_aa'].str.len()
    data['clonotype'] = data['v_call'] + '_' + data['junction_aa'] + '_' + data['j_call'] 
                                                                                                                  
    data = data.drop_duplicates(subset=['sample_id', 'clonotype', 'TIMEPOINTS', 'CONDITION', 'junction_aa'])
    data.reset_index(drop=True, inplace=True)
    
    DETECT_predictions = pd.read_csv(RESULTSDIR + 'DETECT_predictions.tsv', sep = '\t', low_memory=False)
    DETECT_predictions.drop_duplicates(inplace=True)
    DETECT_predictions = DETECT_predictions.query('Score > 0.2 and Antigen == "Spike/surface glycoprotein (S)" and Species == "SARS-CoV-2"')
    DETECT_predictions.reset_index(drop=True, inplace=True)
    DETECT_predictions.rename(columns={'Epitope': 'epitope'}, inplace=True)
    DETECT_predictions['clonotype'] = DETECT_predictions['v_call'] + '_' + DETECT_predictions['junction_aa'] + '_' + DETECT_predictions['j_call'] 
    ####################################################################################################################################
    #plotting.barplot_unique_sequences(data, metadata, PLOTSDIR)
    ####################################################################################################################################

    
    
    
    ####################################################################################################################################
    scatteratio_data = pd.merge(data, DETECT_predictions, on=['clonotype', 'junction_aa', 'j_call', 'v_call'], how='left', indicator=True)
    
    scatteratio_data = scatteratio_data[scatteratio_data['TCR_Chain'].str.contains('TRB')]
    
    scatteratio_data.drop_duplicates(subset=['sample_id', 'clonotype', 'TIMEPOINTS', 'cloneFraction'], inplace=True)   
 
    scatteratio_data = scatteratio_data[['_merge', 'clonotype', 'sample_id', 'CONDITION', 'TIMEPOINTS']]
    scatteratio_data = scatteratio_data.groupby(['sample_id', 'CONDITION', 'TIMEPOINTS'], as_index=False).apply(lambda x: x['_merge'].value_counts())

    scatteratio_data['total'] = scatteratio_data['left_only'] + scatteratio_data['both']
    scatteratio_data['fraction_sequences'] = scatteratio_data['both']/scatteratio_data['total']
    
    # save data for mixed effect analysis
    pd.merge(scatteratio_data, metadata[['sample_id', 'SAMPLE']], on='sample_id').to_csv('scatteratio_data.csv')
    plotting.plot_scatteratio_breadth(scatteratio_data, metadata, 'TIMEPOINTS', 'CONDITION', PLOTSDIR)
    ####################################################################################################################################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ####################################################################################################################################
    scatteratio_data = pd.merge(data, DETECT_predictions, on=['clonotype', 'junction_aa', 'j_call', 'v_call'], how='outer', indicator=True)
    
    
    scatteratio_data = scatteratio_data[scatteratio_data['TCR_Chain'].str.contains('TRB')]
    
    scatteratio_data.drop_duplicates(inplace=True)    
    scatteratio_data = scatteratio_data[['_merge', 'clonotype', 'sample_id', 'CONDITION', 'TIMEPOINTS', 'cloneFraction']]
    scatteratio_data = scatteratio_data.query('_merge == "both"')
    scatteratio_data.reset_index(drop=True, inplace=True)
    
    # save data for mixed effect analysis
    plotting.plot_predictions(scatteratio_data, metadata, 'TIMEPOINTS', 'CONDITION', PLOTSDIR)
    scatteratio_data = scatteratio_data.groupby(['sample_id', 'CONDITION', 'TIMEPOINTS'], as_index=False).agg({'cloneFraction': 'sum'})
    pd.merge(scatteratio_data, metadata[['sample_id', 'SAMPLE']], on='sample_id').to_csv('clonal_fractions_scatteratio_data.csv')
    ####################################################################################################################################

    
    
    



    
    ####################################################################################################################################
    pairs=[("Healthy", "Lymphomas")]
    # Select only the baseline timepoint
    clone_fractions_baseline_pivoted = pd.merge(data, DETECT_predictions, on=['clonotype', 'junction_aa', 'j_call', 'v_call'], how='left', indicator=True)
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted[clone_fractions_baseline_pivoted['_merge'] == 'both']
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted[clone_fractions_baseline_pivoted['TIMEPOINTS'] == 'baseline']
    
    
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted[clone_fractions_baseline_pivoted['TCR_Chain'].str.contains('TRB')]
    # for each value of 'v_call' remove the part of the string after the '-' or '/' character
    #clone_fractions_baseline_pivoted['v_call'] = clone_fractions_baseline_pivoted['v_call'].str.split('-').str[0]
    #clone_fractions_baseline_pivoted['v_call'] = clone_fractions_baseline_pivoted['v_call'].str.split('/').str[0]
    

    clone_fractions_baseline_pivoted.drop_duplicates(subset=['sample_id', 'clonotype', 'TIMEPOINTS', 'CONDITION', 'cloneFraction'], inplace=True)
    clone_fractions_baseline_pivoted.reset_index(drop=True, inplace=True)

    #clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted.pivot_table(index=['sample_id', 'CONDITION'], columns='clonotype', values='cloneFraction', fill_value=0, aggfunc='sum')
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted.pivot_table(index=['sample_id', 'CONDITION', 'TIMEPOINTS'], columns='v_call', values='cloneFraction', fill_value=0, aggfunc='sum')
    clone_fractions_baseline_pivoted.to_csv('vgeneusage.csv')
    plotting.alpha_diversity_tcr(clone_fractions_baseline_pivoted, pairs, PLOTSDIR, 'baseline')
    ####################################################################################################################################
    ####################################################################################################################################
    pairs=[("Healthy", "Lymphomas")]
    # Select only the baseline timepoint
    clone_fractions_baseline_pivoted = pd.merge(data, DETECT_predictions, on=['clonotype', 'junction_aa', 'j_call', 'v_call'], how='left', indicator=True)
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted[clone_fractions_baseline_pivoted['_merge'] == 'both']
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted[clone_fractions_baseline_pivoted['TIMEPOINTS'] == 'V1']
    
    
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted[clone_fractions_baseline_pivoted['TCR_Chain'].str.contains('TRB')]
    # for each value of 'v_call' remove the part of the string after the '-' or '/' character
    #clone_fractions_baseline_pivoted['v_call'] = clone_fractions_baseline_pivoted['v_call'].str.split('-').str[0]
    #clone_fractions_baseline_pivoted['v_call'] = clone_fractions_baseline_pivoted['v_call'].str.split('/').str[0]
    

    clone_fractions_baseline_pivoted.drop_duplicates(subset=['sample_id', 'clonotype', 'TIMEPOINTS', 'CONDITION', 'cloneFraction'], inplace=True)
    clone_fractions_baseline_pivoted.reset_index(drop=True, inplace=True)

    #clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted.pivot_table(index=['sample_id', 'CONDITION'], columns='clonotype', values='cloneFraction', fill_value=0, aggfunc='sum')
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted.pivot_table(index=['sample_id', 'CONDITION', 'TIMEPOINTS'], columns='v_call', values='cloneFraction', fill_value=0, aggfunc='sum')
    clone_fractions_baseline_pivoted.to_csv('vgeneusage.csv')
    plotting.alpha_diversity_tcr(clone_fractions_baseline_pivoted, pairs, PLOTSDIR, 'V1')
    ####################################################################################################################################
    ####################################################################################################################################
    pairs=[("Healthy", "Lymphomas")]
    # Select only the baseline timepoint
    clone_fractions_baseline_pivoted = pd.merge(data, DETECT_predictions, on=['clonotype', 'junction_aa', 'j_call', 'v_call'], how='left', indicator=True)
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted[clone_fractions_baseline_pivoted['_merge'] == 'both']
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted[clone_fractions_baseline_pivoted['TIMEPOINTS'] == 'V3']
    
    
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted[clone_fractions_baseline_pivoted['TCR_Chain'].str.contains('TRB')]
    # for each value of 'v_call' remove the part of the string after the '-' or '/' character
    #clone_fractions_baseline_pivoted['v_call'] = clone_fractions_baseline_pivoted['v_call'].str.split('-').str[0]
    #clone_fractions_baseline_pivoted['v_call'] = clone_fractions_baseline_pivoted['v_call'].str.split('/').str[0]
    

    clone_fractions_baseline_pivoted.drop_duplicates(subset=['sample_id', 'clonotype', 'TIMEPOINTS', 'CONDITION', 'cloneFraction'], inplace=True)
    clone_fractions_baseline_pivoted.reset_index(drop=True, inplace=True)

    #clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted.pivot_table(index=['sample_id', 'CONDITION'], columns='clonotype', values='cloneFraction', fill_value=0, aggfunc='sum')
    clone_fractions_baseline_pivoted = clone_fractions_baseline_pivoted.pivot_table(index=['sample_id', 'CONDITION', 'TIMEPOINTS'], columns='v_call', values='cloneFraction', fill_value=0, aggfunc='sum')
    clone_fractions_baseline_pivoted.to_csv('vgeneusage.csv')
    plotting.alpha_diversity_tcr(clone_fractions_baseline_pivoted, pairs, PLOTSDIR, 'V3')
    ####################################################################################################################################
     
    
    exit(0)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ####################################################################################################################################
    vgenedata = data[data['TCR_Chain'].str.contains('TRB')]
    unique_clonotypes = vgenedata.groupby(['sample_id', 'v_call', 'TIMEPOINTS', 'CONDITION'])['clonotype'].nunique().reset_index()
    total_clonotypes = vgenedata.groupby('sample_id', 'TIMEPOINTS', 'CONDITION')['clonotype'].nunique().reset_index()
    total_clonotypes.rename(columns={'clonotype': 'total_clonotypes'}, inplace=True)

    vgenedata = pd.merge(unique_clonotypes, total_clonotypes, on='sample_id')
    vgenedata['percentage'] = vgenedata['clonotype'] / vgenedata['total_clonotypes']
    
    print(vgenedata)    
    
    
    vgenedata = vgenedata[['_merge', 'clonotype', 'sample_id', 'CONDITION', 'TIMEPOINTS', 'v_call']]
    vgenedata = vgenedata.groupby(['sample_id', 'CONDITION', 'TIMEPOINTS', 'v_call'], as_index=False).apply(lambda x: x['_merge'].value_counts())
    scatteratio_data['total'] = scatteratio_data['left_only'] + scatteratio_data['both']
    scatteratio_data['fraction_sequences'] = scatteratio_data['both']/scatteratio_data['total']
    clone_fractions_baseline_pivoted.drop_duplicates(subset=['sample_id', 'clonotype', 'TIMEPOINTS'], inplace=True)  
    clone_fractions_baseline_pivoted['v_call'] = clone_fractions_baseline_pivoted['v_call'].str.split('-').str[0]
    clone_fractions_baseline_pivoted['v_call'] = clone_fractions_baseline_pivoted['v_call'].str.split('/').str[0] 
    clone_fractions_baseline_pivoted.to_csv('vgeneusage.csv')
    ####################################################################################################################################
    exit(0)
    ####################################################################################################################################
    clone_fractions_pivoted = data.pivot_table(index=['SAMPLE', 'TIMEPOINTS'], columns='clonotype', values='cloneFraction', fill_value=0.0000001, aggfunc='sum')
    # create a pd.Series starting from clone_fractions_pivoted.columns that has value 'black' if the value is in test['clonotype'] and 'gray' otherwise
    clonotypes = set(DETECT_predictions['clonotype'])
    colors = ['red' if col in clonotypes else sns.color_palette(palette='Greys')[1] for col in clone_fractions_pivoted.columns]
    opacities = [1 if col in clonotypes else 0.1 for col in clone_fractions_pivoted.columns]
    plotting.scatter_clone_fractions(clone_fractions_pivoted, colors, opacities, PLOTSDIR)
    ####################################################################################################################################