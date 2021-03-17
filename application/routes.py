from __future__ import print_function
from application import app

from flask import Flask, render_template, flash, request, redirect, session, send_from_directory
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
    proc = subprocess.check_call("mkdir -p data/fasta_files", shell=True)
    proc = subprocess.check_call("mkdir -p data/genes_by_phage", shell=True)
    proc = subprocess.check_call("mkdir -p output", shell=True)

    # download phagesDB info
    download_phage_metadata()
    # download_all_phage_genes()
    # download_all_fastas()

@app.before_request
def clean_output_data():
    proc = subprocess.check_call("find output -type f -mmin +30 -delete", shell=True) # delete files in output that are >30 minutes old

@app.route('/')
def home():
    return render_template("home.html")

@app.route('/file/<path:file_name>', methods=['GET', 'POST'])
def download_file(file_name):
    print(file_name)
    try:
        return send_from_directory("../output/", filename=file_name, as_attachment=True)
    except:
        return BRED(">30 minutes has passed since original request, please rerun search to have access to files")

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

    if request.method == 'POST' and edit_form.validate_on_submit(): # receive deletion inputs
        # extract inputs out of form submission
        phage = edit_form.phage.data
        edit_type = edit_form.edit_type.data
        location_type = edit_form.location_type.data
        bp_position_start = edit_form.bp_position_start.data
        bp_position_stop = edit_form.bp_position_stop.data
        gp_number = edit_form.gp_number.data
        template_DNA = edit_form.template_DNA.data
        orientation = edit_form.orientation.data

        print(phage, edit_type, location_type, bp_position_start, bp_position_stop, gp_number, orientation)
        if edit_form.EGFP.data != None and edit_form.EGFP.data:
            template_DNA = eGFP_DNA()

        if edit_type == "insertion":
            bp_position_stop = bp_position_start

        if location_type == "base pair":
            if bp_position_start == None or bp_position_stop == None:
                Bred("Missing start or stop bound on location")
            elif bp_position_start <= 0 or  bp_position_stop <= 0:
                Bred("base pair must be non-negative")
            elif bp_position_start>bp_position_stop:
                Bred("start bp larger than stop bp")
        
        if edit_type != "deletion" and template_DNA=="":
            Bred("Empty template DNA input")

        phage_info = collect_phage_info(phage)

        if isinstance(phage_info, str):
            Bred("Input phage '{}' is not in the phagesDB, check name spelling and sequenced status".format(phage))

        # if gene product is used to location convert to bps
        if location_type == "gene product number":
            if gp_number == None:
               Bred("missing gene product number")
            out = collect_gene_info(phage, gp_number)
            if isinstance(out, str): # if phage is not in DB 
                BRED("Input phage '{}' is not in the phagesDB, check name spelling and sequenced status".format(phage))
            elif isinstance(out, int): # if gene products is out of bounds
                if out < gp_number:
                    BRED("Gene product {} is out of bounds, {} has a total of {} genes".format(gp_number, phage, out))
                else:
                    BRED("Gene product {} is not a valid gene product (could be labeled as mRNA)".format(gp_number))
            else:
                # if gene is found unpack object
                bp_position_start = out["start"]
                bp_position_stop = out["stop"]
                results["pham"] = out["pham"]
                results["function"] = out["function"]
                results["gene number"] = out["gene number"]
                results["region orientation"] = out["orientation"]

        bp_position_start = int(bp_position_start)
        bp_position_stop = int(bp_position_stop)

        # check if position has buffer upstream and downstream
        if 200 > bp_position_start:
             BRED("Start position of {} does not leave 200bp upstream buffer needed in order to find suitable primers".format(bp_position_start))
        elif phage_info["genome length"]-200 < bp_position_stop:
             BRED("Stop position of {} does not leave downstream 200bp buffer needed in order to find suitable primers, genome length is {}".format(bp_position_stop, phage["genome length"]))

        # check if template DNA is valid DNA sequence
        if template_DNA.replace("A","").replace("T","").replace("G","").replace("C","") != "":
            BRED("Unknown character found in input template DNA")

        # Check if deletion
        if edit_type != "deletion" and orientation=="":
            error = "orientation not specified"
        elif orientation=="R" or (orientation=="Same" and results["region orientation"]=="R"):
            template_DNA = reverse(complement(template_DNA))

        # retrieve phage DNA
        DNA = fasta_to_DNA(phage)

        # find substrate and edit DNA sequence
        substrate, edited_DNA = find_editing_substrate(DNA, bp_position_start, bp_position_stop, template_DNA)

        # find primers
        results["ID"] = process_id = str(random.randint(0,sys.maxsize))
        results["region"] = DNA[bp_position_start-1: bp_position_stop-1]
        results["substrate"] = substrate
        results["edited DNA"] = edited_DNA
        results["DNA"] = DNA
        results["template DNA"] = template_DNA if template_DNA != None else ""
        results["edit type"] = edit_type
        results["homo length"] = int((len(substrate)-len(results["template DNA"]))/2)

        # find primers
        primer_sets = find_primers(DNA, edited_DNA, bp_position_start, bp_position_stop, edit_type, template_DNA)
        if primer_sets.empty:
            error = "unable to find primers"
            return render_template("BRED.html", results = results, primer = {}, phage = phage, bp_position_start = bp_position_start, bp_position_stop = bp_position_stop, gp_number = gp_number, template_DNA = template_DNA, error = error, edit_form = edit_form, colors = colors)

        primer = primer_sets.iloc[0].to_dict()

        # make results downloadable
        primer_sets.drop(["tm range","min heterodimer", "diff tm target"], axis=1, inplace=True)
        primer_sets.to_csv("output/{}_{}_PRIMERS_{}.csv".format(phage, edit_type, results["ID"]))

        BRED_to_fasta(results["ID"], phage, results["edit type"], results["substrate"], primer)
        
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

