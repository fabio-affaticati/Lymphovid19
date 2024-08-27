# 2024 fabio-affaticati

import os
import pandas as pd
import numpy as np
import re
from src.utils import tools


DATADIR = os.path.join(os.path.dirname(__file__), 'data/')
RESULTSDIR = os.path.join(os.path.dirname(__file__), 'results/')
PROCESSEDDIR = os.path.join(RESULTSDIR, 'processed_data/')
MIXCRDIR = os.path.join(DATADIR, 'mixcr_results/')


aminoacids = "ACDEFGHIKLMNPQRSTVWY"
_aminoacids_set = set(aminoacids)


def _is_aaseq(seq: str):
    """
    Check if string contains non-amino acid characters.
    Returns True if string only contains standard amino acid characters.
    """
    try:
        return all(c in _aminoacids_set for c in seq)
    except TypeError:
        return False
    
    
def _is_cdr3(seq: str):
    """
    Checks if string is a valid CDR3 amino acid sequence,
    according to the following defenitions:
        - First amino acid character is C.
        - Last amino acid is F or W.
        - Sequence exclusively contains valid amino acid characters.
    """
    try:
        return (
            _is_aaseq(seq)
            and (seq[0] == "C")
            and (seq[-1] in ["F", "W", "C"])
            and (len(seq) <= 30)
            and (len(seq) >= 4)
        )
    # Exclude non-string type input
    except TypeError:
        return False


if __name__ == "__main__":
    
    for directory in [DATADIR, RESULTSDIR, PROCESSEDDIR]:
        if not os.path.exists(directory):
            os.makedirs(directory)
    
    metadata = pd.read_excel(DATADIR+'metadata.xlsx')
    metadata['SAMPLE_ID'] = metadata['SAMPLE_ID'].astype(str)
    
    
    raw_data = tools.read_raw_data(MIXCRDIR, metadata)
    raw_data = pd.concat(raw_data, ignore_index=True)
    raw_data.reset_index(drop=True,inplace=True)
    
    # remove CD3s that are shorter than 3 amino acids, longer that 25 and do not start with C and end with F
    #raw_data = raw_data[raw_data['aaSeqCDR3'].str.contains('^C[A-Z]{3,25}F$', regex=True, na=False)]
    raw_data = raw_data[raw_data.aaSeqCDR3.apply(lambda x: _is_cdr3(x))]
    #raw_data = raw_data[~raw_data['aaSeqCDR3'].str.contains('_') & ~raw_data['aaSeqCDR3'].str.contains('\*')]
    
    raw_data.reset_index(drop=True, inplace=True)

    ### Remove IGs
    raw_data = raw_data[~raw_data['TCR_Chain'].str.contains('IGL') & ~raw_data['TCR_Chain'].str.contains('IGH') & ~raw_data['TCR_Chain'].str.contains('IGK')]
    raw_data.reset_index(drop=True, inplace=True)
    
    ### Normalize notation of V genes
    r = re.compile(r'[\S]{4,}[D]+[V]+[0-9]{1,2}')
    regmatch = np.vectorize(lambda x: bool(r.match(x)))
    raw_data.loc[regmatch(raw_data.IMGT_VGene_Name.values),'IMGT_VGene_Name'] = raw_data.loc[regmatch(raw_data.IMGT_VGene_Name.values),:]['IMGT_VGene_Name'].str.replace('DV', '/DV')


    # Keep only functional genes
    raw_data = tools.keep_functional_genes('IMGT_VGene_Name', raw_data)
    raw_data = tools.keep_functional_genes('IMGT_JGene_Name', raw_data)
    
    
    raw_data.rename(columns={'IMGT_VGene_Name': 'v_call', 'IMGT_JGene_Name': 'j_call', 'SAMPLE_ID':'sample_id', 'aaSeqCDR3' : 'junction_aa'}, inplace=True)
    raw_data = raw_data[['cloneCount', 'cloneFraction','junction_aa', 'sample_id', 'SAMPLE', 'TIMEPOINTS', 'CONDITION', 'j_call', 'v_call', 'TCR_Chain']]
    raw_data.to_csv(PROCESSEDDIR+'preprocessed_data.csv', sep=',')
    
    # prep data for DETECT
    raw_data[['junction_aa', 'v_call', 'j_call']].drop_duplicates().to_csv(PROCESSEDDIR+'detect_data.tsv', sep='\t', index=False)
    
    
    # prep data for TCRex
    raw_data = raw_data[raw_data['TCR_Chain'].str.contains('TRB')]
    raw_data.rename(columns={'v_call': 'TRBV_gene', 'j_call': 'TRBJ_gene', 'junction_aa' : 'CDR3_beta'}, inplace=True)
    raw_data = raw_data[['CDR3_beta', 'TRBJ_gene', 'TRBV_gene']]
    raw_data.drop_duplicates(inplace=True)
    raw_data.reset_index(drop=True, inplace=True)
    raw_data.to_csv(PROCESSEDDIR+'tcrex_data.tsv', sep='\t', index=False)
    
