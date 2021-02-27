from flask import session
from application import app
import sys
from bs4 import BeautifulSoup as bs
from contextlib import contextmanager
import pandas as pd
import numpy as np 
import random
import subprocess
import time
import math
import re
import requests 
import os.path
from os import path

# define Python user-defined exceptions
class Error(Exception):
    """Base class for other exceptions"""
    pass


class PrimerNotFound(Error):
    """Raised when primer is not found in DNA sequence"""
    pass

class GeneNotFound(Error):
    """Raised when gene is unable to be found in phagesdb"""
    pass

class PhageNotFound(Error):
    """Raised when phage is unable to be found in phagesdb"""
    pass


def download_phage_fasta(phage):
    """ Download phage fasta from phagesDB """
    file_path = 'fasta_files/{}'.format(phage)
    # download file if file is not in directory
    if not path.isfile(file_path):
        url = "https://phagesdb.org/media/fastas/{}.fasta".format(phage)
        r = requests.get(url, allow_redirects=True)
        if r.status_code == 404:
            return "unable to find phage"
        open('fasta_files/{}'.format(phage), 'wb').write(r.content)
        return ""

def fasta_to_DNA(phage):
    """ Extract DNA into string form from a fasta file """
    file_path = 'fasta_files/{}'.format(phage)
    f = open(file_path, "r")
    DNA = ''.join(f.read().split("\n")[1:])
    return DNA

def collect_gene_info(phage, gp_num):
    """ Collect gene info using phagesdb API """
    file_path = "genes_by_phage/{}.csv".format(phage)

    if path.isfile(file_path):
        # if file is in directory do not query API
        genes_df = pd.read_csv(file_path)
    else:
        # if file is not in directory query API and save
        query_url = "https://phagesdb.org/api/genesbyphage/{}/".format(phage)
        out = requests.get(url = query_url).json()["results"]
        if out == {}: # phage not found
            return genes_df.shape[0]
        genes_df = pd.DataFrame(out)

        genes_df.to_csv()

    # extract gene info
    target_gene_df = genes_df[genes_df["Name"] == str(gp_num)]

    # is gene number is not in list 
    if target_gene_df.shape[0] == 0:
        return genes_df.shape[0]

    start_bp = int(target_gene_df["Start"].iloc[0])
    stop_bp = int(target_gene_df["Stop"].iloc[0])
    function = target_gene_df["Notes"].iloc[0]
    pham = target_gene_df["phams"].iloc[0][0]
    
    return start_bp, stop_bp, pham, function

def find_primers(DNA, bp_position_lower, bp_position_upper, template_DNA = "", melting_temp = 60):
    """ Finds forward, middle, reverse primers given DNA sequence 
    
    Insertion Primers:
        - FORWARD upstream of lower bp
        - REVERSE downstream of upper bp
        - EXPERIMENTAL of region unique to insertion DNA

    Swap Primers:
        - FORWARD upstream of lower bp
        - REVERSE downstream of upper bp
        - EXPERIMENTAL of region unique to insertion DNA

    Deletion Primers:
        - FORWARD upstream of lower bp
        - REVERSE downstream of upper bp
        - EXPERIMENTAL of region spanning across deleted region
                ie sequene is 123456789 we want to delete 56 then a EXPERIMENTAL primer could be 3479
    """
    if bp_position_lower == bp_position_upper: # insertion
        pass

    elif len(template_DNA) > 0 and bp_position_lower < bp_position_upper: # swap
        pass
    elif len(template_DNA) == 0 and bp_position_lower < bp_position_upper: # deletion
        pass
    
    return primers

def find_primer(DNA, region, is_experiment = False):
    """ Greedy search to find unique primer to minimize hairpin score and homodimer"""
    # https://libnano.github.io/primer3-py/quickstart.html#calcHairpin


def find_amplicon(DNA, forward_primer, reverse_primer):
    """ Using primers find amplicon DNA 
    
        ex. 
            params
                DNA = ATGGAGGGGCGATG
                forward_primer = TGG
                reverse_primer = ATC
            return 
                AGGGGC
                
    """

    # find primer F location
    F_start_bp = DNA.find(primer_forward)+1
    F_stop_bp = loc_F_start + len(primer_forward)

    # find primer R location
    R_stop_bp = DNA.find(reverse(complement(primer_reverse)))+1
    R_start_bp = R_stop_bp + len(primer_reverse)

    if F_start_bp == 0 or R_stop_bp == 0: # primers are not found in DNA
        raise PrimerNotFound

    amplicon = DNA[F_stop_bp-1:R_stop_bp-1]

    return amplicon

def find_editing_substrate(DNA, bp_position_lower, bp_position_upper, template_DNA = None, homo_length = app.config["HOMOLOGOUS_LENGTH"]):
    """ Find DNA sequence for editing substrate this includes homologous arms
    
    Insertion:
        ex.
            params
                DNA = abcdefghi
                bp_position_lower = 4
                bp_position_upper = 4
                template_DNA = xyz

            return
                cdxyzef

    Swap:
        ex. 
            params
                DNA = abcdefghi
                bp_position_lower = 4
                bp_position_upper = 6
                template_DNA = xyz

            return
                bcxyzgh

    Deletion:
        ex. 
            params
                DNA = abcdefghi
                bp_position_lower = 4
                bp_position_upper = 6
                template_DNA = ""

            return
                bcgh

    """
    template_DNA = ("" if template_DNA == None else template_DNA)
    if bp_position_upper == bp_position_lower: # for insertion
        homo_pre = DNA[bp_position_lower-homo_length:bp_position_lower]
        homo_post = DNA[bp_position_upper:bp_position_upper+homo_length]
        edited_DNA = DNA[:bp_position_lower] + template_DNA + DNA[bp_position_upper:]
    else: # for swap or deletion
        homo_pre = DNA[bp_position_lower-homo_length-1:bp_position_lower-1]
        homo_post = DNA[bp_position_upper:bp_position_upper+homo_length]
        edited_DNA = DNA[:bp_position_lower-1] + template_DNA + DNA[bp_position_upper:]

    substrate = homo_pre + template_DNA + homo_post

    return substrate, edited_DNA

def complement(seq):
    comp_table = {"A":"T", "T":"A", "C":"G", "G":"C"}

    comp_seq = ""
    for bp in seq:
        comp_seq += comp_table[bp]
    return comp_seq
    
def reverse(seq):
    return seq[::-1]