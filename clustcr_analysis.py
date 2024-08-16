# Copyright © 2024 fabio-affaticati

import os
import pandas as pd
from clustcr.clustering.clustering import Clustering
from src.utils.Enrichment import *


RESULTSDIR = os.path.join(os.path.dirname(__file__), 'results/')


if __name__ == "__main__":
    
    for directory in [RESULTSDIR]:
        if not os.path.exists(directory):
            os.makedirs(directory)
    
    
    processed_data = pd.read_csv(RESULTSDIR + 'preprocessed_data.csv', index_col=0)
    
    # keep baseline samples only
    #processed_data = processed_data[processed_data['TIMEPOINTS'] == 'baseline']
    

    processed_data.rename(columns={'IMGT_VGene_Name': 'v_call', 'IMGT_JGene_Name': 'j_call', 'SAMPLE_ID':'sample_id', 'aaSeqCDR3' : 'junction_aa'}, inplace=True)
    processed_data['clonotype'] = processed_data['v_call'] + '_' + processed_data['junction_aa']
    processed_data.reset_index(drop=True, inplace=True)

    ### ClusTCR
    res = Clustering(method='mcl').fit(data=processed_data.drop_duplicates(subset=["clonotype"]),
                                    include_vgene=True, 
                                    cdr3_col="junction_aa",
                                    v_gene_col="v_call"
                                    )
    print(f'Percentage of sequences that were successfully clustered: {round(res.metrics(processed_data).retention(),3) * 100}%')
    res.clusters_df["clonotype"] = res.clusters_df.apply(lambda c: c["v_call"]+"_"+c["junction_aa"], axis=1)


    # or cluster everything and then remove the baseline samples at the end
    processed_data = processed_data[processed_data['TIMEPOINTS'] == 'baseline']
    
    
    data = pd.merge(left=processed_data, right=res.clusters_df[["clonotype", "cluster"]], on="clonotype", how="left")
    data["cluster"] = data["cluster"].fillna(-1).astype(int).astype(str)
    cluster_count_df = pd.crosstab(index=data["cluster"], columns=data["sample_id"])


    ### Remove private clusters
    cluster_count_df = cluster_count_df[cluster_count_df.ne(0).sum(axis=1) > 1]

    ### TO-DO set a threshold for the minimum number of sequences in a cluster

    ### Enrichment analysis
    enrichment_test(res, cluster_count_df, processed_data, 'fisher', 1.0, 'CONDITION', 'Lymphomas', RESULTSDIR, 'lymphoma_v_healthy')