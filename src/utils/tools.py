# Copyright © 2024 fabio-affaticati

import os
import shutil
import pandas as pd
import requests
from bs4 import BeautifulSoup


def keep_functional_genes(gene_col, data):
    
    if gene_col == 'IMGT_VGene_Name':
        url = (
        "https://www.imgt.org/genedb/resultPage.action?"
        "gene.id.species=Homo+sapiens&"
        "molComponent=TR&"
        "geneTypeLike=variable&"
        "allele.fcode=functional&"
        "cloneName=&"
        "locusLike=any&"
        "mainLocusLike=any&"
        "cosLocusLike=any&"
        "groupLike=any&"
        "subgroup=-1&"
        "geneLike=&"
        "selection="
    )
    elif gene_col == 'IMGT_JGene_Name':
        url = (
        "https://www.imgt.org/genedb/resultPage.action?gene.id.species=Homo+sapiens&molComponent=TR&geneTypeLike=joining&"
        "allele.fcode=functional&"
        "cloneName=&"
        "locusLike=any&"
        "mainLocusLike=any&"
        "cosLocusLike=any&"
        "groupLike=any&"
        "subgroup=-1&"
        "geneLike=&"
        "selection="
    )
    else: raise ValueError("Wrong column name given.")
        
    page = requests.get(url)

    soup = BeautifulSoup(page.content, "html.parser")
    results = soup.find("table",{"class":"mainresult"})


    list_genes = []
    for a in results.find_all('a', href=True):
        if a['href'].startswith('individualEntry?name='):
            list_genes.append(a.get_text())
    data = data[data[gene_col].isin(list_genes)]
    data.reset_index(drop=True, inplace=True)
    
    print(f"Unique {gene_col} genes: {data[gene_col].unique()}")
    
    return data


def read_raw_data(mixcrdir, metadata):
    all_files = [f for f in os.listdir(mixcrdir) if "_clones.txt" in f]

    print(f"Files: \n{all_files}")
    print(f'\nNumber of files: {len(all_files)}\n')
    raw_data = []

    for sample_name in metadata['SAMPLE_ID'].unique():
        
        matching = [s for s in all_files if str(sample_name.split('_')[0]) == str(s.split('_')[0]) and not pd.read_csv(mixcrdir+s, sep='\t').empty]
        
        if(matching):
            raw_tcr = pd.concat([pd.read_csv(mixcrdir+match, sep='\t') for match in matching], ignore_index=True)
            raw_tcr['SAMPLE_ID'] = sample_name
            raw_tcr['CONDITION'] = metadata.loc[metadata['SAMPLE_ID'] == sample_name]['CONDITION'].iloc[0]
            raw_tcr['SAMPLE'] = metadata.loc[metadata['SAMPLE_ID'] == sample_name]['SAMPLE'].iloc[0]
            raw_tcr['TIMEPOINTS'] = metadata.loc[metadata['SAMPLE_ID'] == sample_name]['TIMEPOINTS'].iloc[0]
            raw_tcr['IMGT_VGene_Name'] = raw_tcr['allVHitsWithScore'].str.split('\*00').str[0]
            raw_tcr['IMGT_JGene_Name'] = raw_tcr['allJHitsWithScore'].str.split('\*00').str[0]
            raw_tcr['TCR_Chain'] = raw_tcr['allVHitsWithScore'].str[:3]
            
            raw_tcr = raw_tcr.groupby(['aaSeqCDR3', 'SAMPLE_ID', 'SAMPLE', 'TIMEPOINTS', 'CONDITION', 'IMGT_JGene_Name', 'IMGT_VGene_Name', 'TCR_Chain'], as_index=False).agg({'cloneFraction': 'sum', 'cloneCount': 'sum'})
            raw_tcr.reset_index(drop=True,inplace=True)
            
            raw_tcr['cloneFraction'] = raw_tcr['cloneFraction']/len(matching)
            
            print(f"Sample name: {sample_name} \t\t N_unique CDR3: {len(raw_tcr['aaSeqCDR3'].unique())} \t\t Number of files: {len(matching)}")
            
            if not raw_tcr.empty:
                raw_data.append(pd.DataFrame(raw_tcr))
        else: pass

    return raw_data



### Function that removes all files in a directory but leaves the directory itself
def clean_directory(directory:str) -> None:
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        try:
            if os.path.isfile(item_path):
                os.unlink(item_path)
        except Exception as e:
            print(e)
        try:
            if os.path.isdir(item_path):
                shutil.rmtree(item_path)
        except Exception as e:
            print(e)



def copy_directories(source_path, destination_path):
    for item in os.listdir(source_path):
        item_path = os.path.join(source_path, item)
        if os.path.isdir(item_path):
            # copy the directory to the destination path, do not move
            shutil.copytree(item_path, os.path.join(destination_path, item))



def copy_most_recent_directory(source_path, destination_path):
    directories = [item for item in os.listdir(source_path) if os.path.isdir(os.path.join(source_path, item))]
    directories.sort(key=lambda x: os.path.getctime(os.path.join(source_path, x)), reverse=True)
    most_recent_directory = directories[0]
    item_path = os.path.join(source_path, most_recent_directory)
    # if the directory exists, remove it
    if os.path.exists(os.path.join(destination_path, most_recent_directory)):
        shutil.rmtree(os.path.join(destination_path, most_recent_directory))
    shutil.copytree(item_path, os.path.join(destination_path, most_recent_directory))
    
    