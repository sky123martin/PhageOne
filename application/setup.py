import sys
import pickle
import pandas as pd
import numpy as np 
import random
import subprocess
import math
import re
import requests 
import os.path
import multiprocessing 
from multiprocessing.pool import Pool
from os import path

#################################
#### DOWNLOAD PHAGESDB DATA #####
#################################
def download_all_fastas():
    """Download all fasta from phagesDB"""
    if not path.isfile("data/fasta_files/ZygoTaiga.fasta"):
        proc = subprocess.check_call("mkdir -p data/fasta_files", shell=True)

        url = "https://phagesdb.org/media/Actinobacteriophages-All.fasta"
        r = requests.get(url, allow_redirects=True)

        if r.status_code == 404:
            return "unable to find multifasta"

        # split fasta
        phage_fastas = str(r.content).split(">")

        # iterate through and download
        for fasta in phage_fastas:
            fasta = ">" + fasta
            if len(fasta)>1000:
                    phage_name = str(fasta).split('\\n',1)[0].split(" ")[2].replace(",", "")
                    f = open('data/fasta_files/{}.fasta'.format(phage_name), 'w')
                    f.writelines(fasta.replace("\\n"," \n"))
                    f.close()

def download_phage_metadata():
    """Download sequenced phage metadata from phagesDB"""
    url = "https://phagesdb.org/data/?set=seq&type=full"

    r = requests.get(url, allow_redirects=True)

    if r.status_code == 404:
        return "unable to find multifasta"
    
    lines = str(r.content).split("\\n")
    lines = [l.split("\\t") for l in lines]
    df = pd.DataFrame(columns = lines[0], data = lines[1:])

    # rename to fit with research code
    df.rename(columns={'b\'Phage Name':'phage',
                        'Temperate?':'temperate',
                        'Cluster':'cluster',
                        'Subcluster':'subcluster',
                        'Morphotype':'morphotype',
                        'Host':'host',
                        'Genome Length(bp)':'genome length',
                        'GC%':'gcpercent',
                        'Phamerated?':'is phamerated',
                        'Annotation Status':'is annotated'
                       }, inplace = True)
    
    # only retrieve specific columns
    df = df[['phage','temperate','cluster','subcluster','morphotype','host','genome length','is annotated','is phamerated','gcpercent']]
   
    # set types
    df = df.astype({'phage':'str',
                'temperate':'bool',
                'cluster':'str',
                'subcluster':'str',
                'morphotype':'str',
                'host':'str',
                'genome length':'float',
                'is annotated':'bool',
                'is phamerated':'bool',
                'gcpercent':'float'})

    # export to csv
    df.to_csv("data/phage_metadata.csv")

def download_phage_genes(phage):
    query_url = "https://phagesdb.org/api/genesbyphage/{}/".format(phage)

    try:
        response = requests.get(url = query_url).json()['results']
        df = pd.DataFrame(response)

        a_file = open("static_data/conversion_table.pkl", "rb")
        conversion_table = pickle.load(a_file)
        df["phage"] = df["PhageID"].apply(lambda x: x['PhageID'])
        df["gene number"] = df['GeneID'].apply(lambda x: x.split("_")[-1])
        df["pham"] = df['phams'].apply(lambda x: x[0])
        df["Notes"] = df["Notes"].apply(lambda x: str(x).lower()[2:-1])
        df["function"] = df["Notes"].apply(lambda x: conversion_table[x] if x in conversion_table.keys() else "NKF")

        df.rename(columns={'GeneID':'gene ID',
                            'translation':'translation',
                            'Orientation':'orientation',
                            'Start':'start',
                            'Stop':'stop',
                            'Notes':'uncleaned function'
                            }, inplace = True)
                            
        df = df[['gene ID','pham','function','translation','orientation','phage','gene number','start','stop','uncleaned function']]
        df.to_csv("data/genes_by_phage/{}_genes.csv".format(phage))
        print(phage)
        return df
    except: # if phage has no genes df then edit the metadata csv and return empty DF
        phage_df = pd.read_csv("data/phage_metadata.csv")
        phage_df[phage_df["phage"]!=phage].to_csv("data/phage_metadata.csv")
        return pd.DataFrame(columns=['gene ID','pham','function','translation','orientation','phage','gene number','start','stop','uncleaned function'])

def download_all_phage_genes():
    """Download list of phage genes for each phage that is sequenced in phagesDB"""

    if not path.isfile("data/cleaned_gene_list.csv"):

        meta_file_path = "data/phage_metadata.csv"
        if not path.isfile(meta_file_path):
            download_phage_metadata()
        phage_df = pd.read_csv(file_path)
        phages = phage_df["phage"].unique()

        proc = subprocess.check_call("mkdir -p data/genes_by_phage", shell=True)

        with Pool(multiprocessing.cpu_count()-2) as p:  # to check multiprocessing.cpu_count()
            output_dfs = p.map(download_phage_genes, phages)

        all_phages_genes_df = pd.concat(output_dfs)
        all_phages_genes_df.to_csv("data/cleaned_gene_list.csv")