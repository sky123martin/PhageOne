from flask import session
from application import app
import sys
from bs4 import BeautifulSoup as bs
from contextlib import contextmanager
import primer3
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


def primer3_calculate_hairpin(seq):
    return primer3.calcHairpin(seq).dg

def primer3_calculate_tm(seq):
    tm = primer3.calcTm(seq)
    return tm if isinstance(tm, float) else None

def primer3_calculate_heterodimer(seq1, seq2):
    return primer3.calcHeterodimer(seq1, seq2).dg

def find_start_position(primer, DNA, orientation):
    if orientation == "F":
        start = DNA.find(primer)
    elif orientation == "R":
        start = DNA.find(reverse(complement(primer))) + len(primer)
    return start + 1

def find_primers(DNA, edited_DNA, bp_position_start, bp_position_stop, edit_type, template_DNA = "", melting_temp = 60, primer_length = 20):
    """ Finds forward, middle, reverse primers given DNA sequence 
    
    Insertion Primers:
        - FORWARD upstream of start bp
        - REVERSE downstream of stop bp
        - EXPERIMENTAL of region unique to insertion DNA

    Swap Primers:
        - FORWARD upstream of start bp
        - REVERSE downstream of stop bp
        - EXPERIMENTAL of region unique to insertion DNA

    Deletion Primers:
        - FORWARD upstream of start bp
        - REVERSE downstream of stop bp
        - EXPERIMENTAL of region spanning across deleted region
                ie sequene is 123456789 we want to delete 56 then a EXPERIMENTAL primer could be 3479
    """
    possible_primers = find_possible_primers(DNA, edited_DNA, bp_position_start, bp_position_stop, edit_type, template_DNA)
    possible_primers["hairpin"] = possible_primers["primer_seq"].apply(primer3_calculate_hairpin)
    possible_primers["tm"] = possible_primers["primer_seq"].apply(primer3_calculate_tm)
    possible_primers.dropna(inplace=True)
    # set hairpin to zero, below zero indicates binding
    possible_primers = possible_primers[possible_primers["hairpin"]>=-1000]

    # sort by target melting temp
    possible_primers = possible_primers.sort_values(by=["tm"]).reset_index(drop=True)

    max_tm_range = 0.1

    primer_sets = pd.DataFrame(columns=["primer f","primer r","primer e","primer c", "tm f", "tm r", "tm c", "tm e", "min heterodimer", "mean tm", "tm range", "min hairpin"])

    while primer_sets.empty:
        for index, row in possible_primers.iterrows(): 
            current_tm = row["tm"]
            primers_in_range = possible_primers.iloc[index:][abs(possible_primers.iloc[index:]["tm"]-current_tm) < max_tm_range]
            if len(primers_in_range["primer"].unique()) == 4: # all needed primers
                for _, f_row in primers_in_range[primers_in_range["primer"]=="f"].iterrows():
                    for _, r_row in primers_in_range[primers_in_range["primer"]=="r"].iterrows():
                        for _, c_row in primers_in_range[primers_in_range["primer"]=="c"].iterrows():
                            for _, e_row in primers_in_range[primers_in_range["primer"]=="e"].iterrows():
                                # temperature calculations
                                mean_tm = (f_row["tm"]+r_row["tm"]+c_row["tm"]+e_row["tm"])/4
                                min_tm = min([f_row["tm"],r_row["tm"],c_row["tm"],e_row["tm"]])
                                max_tm = max([f_row["tm"],r_row["tm"],c_row["tm"],e_row["tm"]])

                                # primer compatibilty (heterodimer score)
                                dg_f = primer3_calculate_heterodimer(f_row["primer_seq"], r_row["primer_seq"])
                                dg_c = primer3_calculate_heterodimer(c_row["primer_seq"], r_row["primer_seq"])
                                dg_e = primer3_calculate_heterodimer(e_row["primer_seq"], r_row["primer_seq"])
                                min_heterodimer = min([dg_f, dg_c, dg_e])

                                # min hairpin
                                min_hairpin = min([f_row["hairpin"], r_row["hairpin"], c_row["hairpin"], e_row["hairpin"]])
                                
                                primer_sets = primer_sets.append({  "primer f": f_row["primer_seq"],
                                                                    "primer r": r_row["primer_seq"],
                                                                    "primer e": e_row["primer_seq"],
                                                                    "primer c": c_row["primer_seq"],
                                                                    "tm f": round(f_row["tm"],2),
                                                                    "tm r": round(r_row["tm"],2),
                                                                    "tm c": round(c_row["tm"],2),
                                                                    "tm r": round(e_row["tm"],2),
                                                                    "hetero f dg": round(dg_f,2),
                                                                    "hetero c dg": round(dg_c,2),
                                                                    "hetero e dg": round(dg_e,2),
                                                                    "hairpin f": round(f_row["hairpin"],2),
                                                                    "hairpin r": round(r_row["hairpin"],2),
                                                                    "hairpin c": round(c_row["hairpin"],2),
                                                                    "hairpin e": round(e_row["hairpin"],2),
                                                                    "min heterodimer": -round(min_heterodimer), # sign change and rounding are for final sort
                                                                    "mean tm": mean_tm,
                                                                    "diff tm target": round(abs(melting_temp-mean_tm),1),
                                                                    "tm range": round(max_tm-min_tm, 1),
                                                                    "min hairpin":min_hairpin}, ignore_index=True)
        max_tm_range += 0.1
        if max_tm_range > 1:
            return primer_sets

    primer_sets.drop_duplicates(inplace=True)

    # confirm starts in each primer
    primer_sets["f start"] = primer_sets["primer f"].apply(find_start_position, args=(DNA,"F"))
    primer_sets["f edited start"] = primer_sets["primer f"].apply(find_start_position, args=(edited_DNA,"F"))

    primer_sets["r start"] = primer_sets["primer r"].apply(find_start_position, args=(DNA,"R"))
    primer_sets["r edited start"] = primer_sets["primer r"].apply(find_start_position, args=(edited_DNA,"R"))

    primer_sets["c start"] = primer_sets["primer c"].apply(find_start_position, args=(DNA,"F"))
    primer_sets["c edited start"] = primer_sets["primer c"].apply(find_start_position, args=(edited_DNA,"F")) # this should be zero becuase it doesn't bind!

    primer_sets["e start"] = primer_sets["primer e"].apply(find_start_position, args=(DNA,"F"))
    primer_sets["e edited start"] = primer_sets["primer e"].apply(find_start_position, args=(edited_DNA,"F"))  # this should be zero becuase it doesn't bind!


    return primer_sets.sort_values(by = ["tm range","min heterodimer", "diff tm target"]).reset_index(drop=True)


def all_possible_subsets(seq, primer_length, primer_type):
    """ Takes sequence and find all possible string subsets of primer_length
    """
    return pd.DataFrame(np.array([[seq[i:i+primer_length], primer_type] for i in range(len(seq)-primer_length+1)]), columns=['primer_seq','primer'])

def find_possible_primers(DNA, edited_DNA, bp_position_start, bp_position_stop, edit_type, template_DNA, primer_length = 20, search_size = 200, testing = False):
    """ Finds region to extract primers from then get ever possible subset
    """
    start_index = bp_position_start-1
    stop_index = bp_position_stop-1 
    if edit_type == "replacement": # insertion
        # possible forward primers
        forward_region = DNA[start_index-primer_length-search_size:start_index-primer_length]
        f_df = all_possible_subsets(forward_region, primer_length, "f")
        # possible reverse primers
        if testing:
            reverse_region = DNA[stop_index+primer_length:stop_index+primer_length+search_size]
        else:    
            reverse_region = reverse(complement(DNA[stop_index+primer_length:stop_index+primer_length+search_size]))
        r_df = all_possible_subsets(reverse_region, primer_length, "r")

        # possible control primers
        control_region = DNA[start_index:stop_index]
        c_df = all_possible_subsets(control_region, primer_length, "c")

        # possible experimental primers
        experimental_region = edited_DNA[start_index:start_index+len(template_DNA)]
        e_df = all_possible_subsets(experimental_region, primer_length, "e")

    elif edit_type == "insertion":
        # possible forward primers
        forward_region = DNA[start_index-primer_length-search_size:start_index-primer_length]
        f_df = all_possible_subsets(forward_region, primer_length, "f")

        # possible reverse primers
        if testing:
            reverse_region = DNA[stop_index+primer_length:stop_index+primer_length+search_size]
        else:    
            reverse_region = reverse(complement(DNA[stop_index+primer_length:stop_index+primer_length+search_size]))

        r_df = all_possible_subsets(reverse_region, primer_length, "r")

        # possible control primers
        control_region = DNA[start_index-primer_length+1:stop_index+primer_length-1]
        c_df = all_possible_subsets(control_region, primer_length, "c")

        # possible experimental primers
        experimental_region = edited_DNA[start_index:start_index+len(template_DNA)]

        e_df = all_possible_subsets(experimental_region, primer_length, "e")

    elif  edit_type == "deletion": # deletion
        # possible forward primers
        forward_region = DNA[start_index-primer_length-search_size:start_index-primer_length]
        f_df = all_possible_subsets(forward_region, primer_length, "f")
        
        # possible reverse primers
        if testing:
            reverse_region = DNA[stop_index+primer_length:stop_index+primer_length+search_size]
        else:    
            reverse_region = reverse(complement(DNA[stop_index+primer_length:stop_index+primer_length+search_size]))
        r_df = all_possible_subsets(reverse_region, primer_length, "r")

        # possible control primers
        control_region = DNA[start_index:stop_index]
        c_df = all_possible_subsets(control_region, primer_length, "c")
        # possible experimental primers
        experimental_region = edited_DNA[start_index-primer_length+1:start_index+primer_length-1]
        e_df = all_possible_subsets(experimental_region, primer_length, "e")
    
    possible_primers = pd.concat([f_df, r_df, c_df, e_df])

    # check for duplicate strands
    if not possible_primers["primer_seq"].is_unique:
        raise TypeError

    return possible_primers

    

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
    F_start_bp = DNA.find(forward_primer)+1
    F_stop_bp = F_start_bp + len(forward_primer)

    # find primer R location
    R_stop_bp = DNA.find(reverse(complement(reverse_primer)))+1
    R_start_bp = R_stop_bp + len(reverse_primer)

    if F_start_bp == 0 or R_stop_bp == 0: # primers are not found in DNA
        raise PrimerNotFound

    amplicon = DNA[F_stop_bp-1:R_stop_bp-1]

    return amplicon

def find_editing_substrate(DNA, bp_position_start, bp_position_stop, template_DNA = None, homo_length = app.config["HOMOLOGOUS_LENGTH"]):
    """ Find DNA sequence for editing substrate this includes homologous arms
    
    Insertion:
        ex.
            params
                DNA = abcdefghi
                bp_position_start = 4
                bp_position_stop = 4
                template_DNA = xyz

            return
                abcdxyzefghi

    Swap:
        ex. 
            params
                DNA = abcdefghi
                bp_position_start = 4
                bp_position_stop = 6
                template_DNA = xyz

            return
                bcxyzgh

    Deletion:
        ex. 
            params
                DNA = abcdefghi
                bp_position_start = 4
                bp_position_stop = 6
                template_DNA = ""

            return
                bcgh

    """
    template_DNA = ("" if template_DNA == None else template_DNA)
    if bp_position_stop == bp_position_start: # for insertion
        homo_pre = DNA[bp_position_start-homo_length:bp_position_start]
        homo_post = DNA[bp_position_stop:bp_position_stop+homo_length]
        edited_DNA = DNA[:bp_position_start] + template_DNA + DNA[bp_position_stop:]
    else: # for swap or deletion
        homo_pre = DNA[bp_position_start-homo_length-1:bp_position_start-1]
        homo_post = DNA[bp_position_stop:bp_position_stop+homo_length]
        edited_DNA = DNA[:bp_position_start-1] + template_DNA + DNA[bp_position_stop:]

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