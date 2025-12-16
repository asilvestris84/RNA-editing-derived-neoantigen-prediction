#!/lustrehome/pietrolucamazzacuva/AnaConda/envs/nndev/bin/python


"""\
This script is used to spot potential A-to-I RNA editing sites in GTEx samples REDItools output tables, using fixed empirical filters on base depth, number of A-to-G substitutions and A-to-G rate, as well as absence of genomic source of variation. This script can only be applied to unstranded RNA-seq REDItools output tables. The name of each output table must contain the SRA accession number (SRR) of the corresponding RNA-seq experiment, followed by an optional suffix (euqals for all experiments). Each output table must be contained in a folder named after the SRA accession number, containing any additional data related to that experiment. All folders containing data from individual experiments on the same tissue must be contained in a single folder named after the tissue. There must be a separate folder for each tissue.

Usage: python3 AtoI_Candidate_Sites_Search.py inputpath outputpath tablesuffix tissue_1 srr_1 tissue_2 srr_2 ... tissue_n srr_n

inputpath: path to the folder containing the tissue folders
outputpath: path to the folder used store the script output
tablesuffix: REDItools output tables optional suffix
tissue_x: name of the tissue folder 
srr_x: accession number of the experiment

note: srr_x must be an accession number present in the tissue_x folder  

Author: Pietro Luca Mazzacuva
Ccopyright: Copyright 2025, RNA-editing-derived-neoantigen-prediction
Credits: Alessandro Silvestris, Miriam Rahina Semente
License: GPL
Version: 1.0.1
Maintainer: Pietro Luca Mazzacuva
Email: pietroluca.mazzacuva@unicampus.it
Status: Production
"""


#importing the required libraries
import pandas as pd
import gzip, os, sys
from tqdm import tqdm
from subprocess import Popen, PIPE
from multiprocessing import Pool


#function to search candidate A-to-I sites 
def samplecheck(tissue, srr):

    #define global variables and local variables
    global inputpath
    global columns
    global outputpath
    global teblesuffix
    savepath = os.path.join(outputpath, tissue)
    if not os.path.isdir(savepath):
        os.mkdir(savepath)
    #srr is the accession number of the 
    srrpath = os.path.join(inputpath, f"{srr}/{srr}{teblesuffix}")
    
    #search of cancidate A-to-I sites
    if os.path.isfile(srrpath):
        editing = []
        p = Popen(["zgrep", "-c", "$", srrpath], stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()

        #progressbar setting to monitor the execution progress
        with tqdm(total=int(str(out.strip().split()[0]).replace("'", " ")[2:]), desc=srr) as progressbar:
            with gzip.open(srrpath, "r") as File:
                for Row in File:
                    row = Row.decode("utf-8").strip().split("\t")

                    #candidate A-to-I sites conditions to satisfy 
                    if int(row[3])==2 and row[12] == "-" and int(row[4]) >= 10:
                        if row[9] != "-":
                            if int(row[9]) >= 10:
                                counts = eval(row[6])
                                if row[7]=="AG":
                                    if counts[2] >=3 and counts[2]/sum(counts)>=0.01:
                                        row[3] = "1"
                                        editing.append(row + [srr+","])
                                elif row[7]=="TC":
                                    if counts[1] >=3 and counts[1]/sum(counts)>=0.01:
                                        row[3] = "0"
                                        editing.append(row + [srr+","])
                        else:
                            counts = eval(row[6])
                            if row[7]=="AG":
                                if counts[2] >=3 and counts[2]/sum(counts)>=0.01:
                                    row[3] = "1"
                                    editing.append(row + [srr+","])
                            elif row[7]=="TC":
                                if counts[1] >=3 and counts[1]/sum(counts)>=0.01:
                                    row[3] = "0"
                                    editing.append(row + [srr+","])
                    progressbar.update(1)   

        #search result saving            
        editing = pd.DataFrame(editing, columns=columns+["SRR"])
        editing.sort_values(by=["Region", "Position"], inplace=True)
        editing.to_csv(os.path.join(savepath, f"{srr}.txt"), sep="\t", index=None)
    else:
        editing = pd.DataFrame()
        editing.to_csv(os.path.join(savepath, f"{srr}.txt"), sep="\t", index=None)


#execution parallelization on multiple cores
inputpath, outputpath, teblesuffix = sys.argv[1], sys.argv[2], sys.argv[3]
columns = ["Region", "Position", "Reference", "Strand" , "Coverage-q30", "MeanQ", "BaseCount[A,C,G,T]", 
           "AllSubs", "Frequency", "gCoverage-q30", "gMeanQ", "gBaseCount[A,C,G,T]", "gAllSubs", "gFrequency"]

counter = 0
inputs = []
temp = []
for i in sys.argv[4:]:
    temp.append(i)
    counter +=1
    if counter == 2:
        inputs.append(temp)
        temp = []
        counter = 0

with Pool(32) as pool:
    pool.starmap(samplecheck, inputs)
