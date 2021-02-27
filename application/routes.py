from __future__ import print_function
from application import app

from flask import Flask, render_template, flash, request, redirect, session
from application.forms import BRED_edit_form
from application.utility import *
import pandas as pd 
import numpy as np 
import random
import sys
import math
import json

from application import utility

@app.route('/', methods=['GET', 'POST'])
def home(error = ""):
    return render_template("home.html")

@app.route('/BRED', methods=['GET', 'POST'])
def BRED(error = ""):
    # initialize forms
    edit_form = BRED_edit_form()

    results = {}

    phage = ""
    bp_position_lower = ""
    bp_position_upper = ""
    gp_number = ""
    template_DNA = ""

    # receive inputs
    if error != "":
        print(error)
        return render_template("BRED.html", results = results, phage = phage, bp_position_lower = bp_position_lower, bp_position_upper = bp_position_upper, gp_number = gp_number, template_DNA = template_DNA, error = error, edit_form = edit_form)

    elif request.method == 'POST' and edit_form.validate_on_submit(): # receive deletion inputs
        # extract inputs out of form submission
        phage = edit_form.phage.data
        bp_position_lower = edit_form.bp_position_lower.data
        bp_position_upper = edit_form.bp_position_upper.data
        gp_number = edit_form.gp_number.data
        template_DNA = edit_form.template_DNA.data

        # error handling on location
        if bp_position_lower == None and bp_position_upper == None and gp_number == None:
            error = "Missing location! Must input either start/stop bp position or gene production number"
        elif (bp_position_lower != None or bp_position_upper != None) == (gp_number != None):
            error = "Must input either start/stop bp position or gene production number, both not both"
        elif ((bp_position_lower != None) != (bp_position_upper != None)):
            error = "Missing start or stop bound on location"
        elif bp_position_lower != None and bp_position_upper != None and bp_position_lower > bp_position_upper:
            error = "start bp larger than stop bp"

        if error != "":
            return render_template("BRED.html", results = results, phage = phage, bp_position_lower = bp_position_lower, bp_position_upper = bp_position_upper, gp_number = gp_number, template_DNA = template_DNA, error = error, edit_form = edit_form)
            

        # if gene product is used to location convert to bps
        if bp_position_lower == None and bp_position_upper == None and gp_number != None:
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
                bp_position_lower, bp_position_upper, pham, function = out
                results["pham"] = pham
                results["function"] = function[2:-1]
                results["gene number"] = gp_number
        
        # check if template DNA is valid DNA sequence
        if template_DNA.replace("A","").replace("T","").replace("G","").replace("C","") != "":
            error = "Unkown character found in input template DNA"

        if error != "":
            return render_template("BRED.html", results = results, phage = phage, bp_position_lower = bp_position_lower, bp_position_upper = bp_position_upper, gp_number = gp_number, template_DNA = template_DNA, error = error, edit_form = edit_form)
          

        # Determine edit type
        if bp_position_upper == bp_position_lower:
            edit_type = "Insertion"
        elif template_DNA == "":
            edit_type = "Deletion"
        elif template_DNA != "":
            edit_type = "Replacement"
        else:
            edit_type = ""

        # retrieve phage DNA
        download_phage_fasta(phage)

        DNA = fasta_to_DNA(phage)

        # find substrate and edit DNA sequence
        substrate, edited_DNA = find_editing_substrate(DNA, bp_position_lower, bp_position_upper, template_DNA)

        # find primers

        results["substrate"] = substrate
        results["edited DNA"] = edited_DNA
        results["DNA"] = DNA
        results["edit type"] = edit_type
        

        return render_template("BRED.html", results = results, phage = phage, bp_position_lower = bp_position_lower, bp_position_upper = bp_position_upper, gp_number = "", template_DNA = template_DNA, error = error, edit_form = edit_form)
        
    return render_template("BRED.html", results = results, phage = phage, bp_position_lower = bp_position_lower, bp_position_upper = bp_position_upper, gp_number = gp_number, template_DNA = template_DNA, error = error, edit_form = edit_form)
