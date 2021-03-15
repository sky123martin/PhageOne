from __future__ import print_function
from application import app

from flask import Flask, render_template, flash, request, redirect, session
from application.forms import BRED_edit_form, editing_guide_form
from application.utility import *
import pandas as pd 
import numpy as np 
import random
import sys
import math
import os
import json

from application import utility

@app.before_first_request
def setup_data():
    # make data directory
    proc = subprocess.check_call("mkdir -p data", shell=True)

    # download phagesDB info
    download_phage_metadata()
    download_all_phage_genes()
    download_all_fastas()

@app.after_request
def after_request_cleaning(response):
    print("after_request executing!", print(response))
    return response

@app.route('/')
def home():
    return render_template("home.html")

@app.route('/BRED', methods=['GET', 'POST'])
def BRED(error = ""):

    # <h3>Subcluster: <a href="https://phagesdb.org/subclusters/{{phage_info['subcluster']}}/">{{phage_info['subcluster']}}</a></h3>
    # <h3>Cluster: <a href="https://phagesdb.org/clusters/{{phage_info['cluster']}}/">{{phage_info['subcluster']}}</a></h3>
    # <h3>Cluster Lifecycle: <a href="https://phagesdb.org/glossary/#{{phage_info['Cluster Life']}}/">{{phage_info['subcluster']}}</a></h3>
    # initialize forms
    edit_form = BRED_edit_form()

    colors = {"homo":"#ffddd2",
              "old_region":"#006d77",
              "new_region":"#e29578",
              "primer":"#83c5be",
              "background": "#edf6f9"} # four colors

    results = {}

    phage = ""
    bp_position_start = ""
    bp_position_stop = ""
    gp_number = ""
    template_DNA = ""

    # receive inputs
    if error != "":
        print(error)
        return render_template("BRED.html", results = results, phage = phage, bp_position_start = bp_position_start, bp_position_stop = bp_position_stop, gp_number = gp_number, template_DNA = template_DNA, error = error, edit_form = edit_form, colors = colors)

    elif request.method == 'POST' and edit_form.validate_on_submit(): # receive deletion inputs
        # extract inputs out of form submission
        phage = edit_form.phage.data
        edit_type = edit_form.edit_type.data
        location_type = edit_form.location_type.data
        bp_position_start = edit_form.bp_position_start.data
        bp_position_stop = edit_form.bp_position_stop.data
        gp_number = edit_form.gp_number.data
        template_DNA = edit_form.template_DNA.data


        if location_type == "base pair":
            if bp_position_start == None and bp_position_stop == None and edit_type == "insertion":
                error = "insertion requires bp location"
            elif bp_position_start <= 0 and bp_position_stop == None and edit_type == "insertion":
                error = "base pair must be non-negative"
            elif bp_position_start > 0 and bp_position_stop == None and edit_type == "insertion":
                pass
            elif bp_position_start == None or bp_position_stop == None:
                error = "Missing start or stop bound on location"
            elif bp_position_start<=0 or bp_position_stop<=0:
                error = "base pair must be non-negative"
            elif bp_position_start>bp_position_stop:
                error = "start bp larger than stop bp"
        elif location_type == "gene product number":
            if gp_number == None:
                error = "missing gene product number"
            
        if edit_type == "insertion":
            bp_position_stop = bp_position_start

        if error != "":
            return render_template("BRED.html", results = results, phage = phage, primer = {}, bp_position_start = bp_position_start, bp_position_stop = bp_position_stop, gp_number = gp_number, template_DNA = template_DNA, error = error, edit_form = edit_form, colors = colors)
            

        # if gene product is used to location convert to bps
        if location_type == "gene product number":
            out = collect_gene_info(phage, gp_number)

            if isinstance(out, str): # if phage is not in DB 
                error = "Inputted phage '{}' is not in the phagesDB, check name spelling and be sure it is sequenced".format(phage)
            elif isinstance(out, int): # if geneproducts is out of bounds
                if out < gp_number:
                    error = "Gene product {} is out of bounds, {} has a total of {} genes".format(gp_number, phage, out)
                else:
                    error = "Gene product {} is not a valid gene product (could be labeled as mRNA)".format(gp_number) 
            else:
                # if gene is found unpack object
                bp_position_start, bp_position_stop, pham, function = out
                results["pham"] = pham
                results["function"] = function[2:-1]
                results["gene number"] = gp_number
        
        # check if template DNA is valid DNA sequence
        if template_DNA.replace("A","").replace("T","").replace("G","").replace("C","") != "":
            error = "Unkown character found in input template DNA"

        if error != "":
            return render_template("BRED.html", results = results, primer = {}, phage = phage, bp_position_start = bp_position_start, bp_position_stop = bp_position_stop, gp_number = gp_number, template_DNA = template_DNA, error = error, edit_form = edit_form, colors = colors)

        # retrieve phage DNA
        DNA = fasta_to_DNA(phage)

        # find substrate and edit DNA sequence
        substrate, edited_DNA = find_editing_substrate(DNA, bp_position_start, bp_position_stop, template_DNA)

        # find primers
        primer_sets = find_primers(DNA, edited_DNA, bp_position_start, bp_position_stop, edit_type, template_DNA)
        if primer_sets.empty:
            error = "unable to find primers"
            return render_template("BRED.html", results = results, primer = {}, phage = phage, bp_position_start = bp_position_start, bp_position_stop = bp_position_stop, gp_number = gp_number, template_DNA = template_DNA, error = error, edit_form = edit_form, colors = colors)

        primer = primer_sets.iloc[1].to_dict()

        # find primers
        results["region"] = DNA[bp_position_start-1: bp_position_stop-1]
        results["substrate"] = substrate
        results["edited DNA"] = edited_DNA
        results["DNA"] = DNA
        results["template DNA"] = template_DNA if template_DNA != None else ""
        results["edit type"] = edit_type
        
        return render_template("BRED.html", results = results, primer = primer, phage = phage, bp_position_start = bp_position_start, bp_position_stop = bp_position_stop, gp_number = "", template_DNA = template_DNA, error = error, edit_form = edit_form, colors = colors)
        
    return render_template("BRED.html", results = results, primer = {}, phage = phage, bp_position_start = bp_position_start, bp_position_stop = bp_position_stop, gp_number = gp_number, template_DNA = template_DNA, error = error, edit_form = edit_form, colors = colors)


@app.route('/EditingGuide', methods=['GET', 'POST'])
def EditingGuide(error = ""):
    form = editing_guide_form()
    return render_template("EditingGuide.html", error=error)

@app.route('/GenomeSynteny', methods=['GET', 'POST'])
def GenomeSynteny(error = ""):
        print(error)
        return render_template("GenomeSynteny.html", error=error)

