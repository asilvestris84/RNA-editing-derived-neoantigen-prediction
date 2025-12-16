#!/lustrehome/pietrolucamazzacuva/AnaConda/envs/nndev/bin/python


"""\
This script is used to aggregate in a single table the the AtoI_Candidate_Sites_Search.py outputs and compute additional data on the candidates A-to-I RNA editing sites previously identified.

Usage: python3 AtoI_Candidate_Sites_Data_Aggregation.py path

path: path to the folder containing the AtoI_Candidate_Sites_Search.py outputs, used also to store this script output

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
import os, glob, sys
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool

#functions to calculate A-to-I rates
def editingAG(x):
    return round(100*eval(x)[2]/sum(eval(x)), 2)

def editingTC(x):
    return round(100*eval(x)[1]/sum(eval(x)), 2)

#function to calculate the number of samples affected by specific A-to-I sites
def nsamples(x):
    return float(len(x.split(",")))

#function to change the dtype from int to string
def tostr(x):
    return str(x)+","

#function to aggregate the candidate A-to-I sites data
def returntissue(tissue):

    #global and local variables assignment
    global outputpath
    global cols
    savepath = os.path.join(outputpath, tissue)
    df = pd.DataFrame()
    paths = glob.glob(savepath+"/SR*")

    #progressbar setting to monitor the execution 
    with tqdm(total = len(paths), desc=tissue, leave=True) as progressbar:
        for path in paths:
            try:
                df=pd.concat([df, pd.read_table(path, usecols=cols)], axis=0)
            except:
                continue
            progressbar.update(1)

    #candidate A-to-I sites data aggregation and mining
    df.reset_index(drop=True, inplace=True)
    df.loc[:, "Editing"] = 0.0
    A = df[df.loc[:, "Reference"]=="A"].index.tolist()
    T = df[df.loc[:, "Reference"]=="T"].index.tolist()
    df.loc[A, "Editing"] = df.loc[A, "BaseCount[A,C,G,T]"].apply(editingAG)
    df.loc[T, "Editing"] = df.loc[T, "BaseCount[A,C,G,T]"].apply(editingTC)
    del A, T
    df.drop(["Reference", "BaseCount[A,C,G,T]"], axis=1, inplace=True)
    srr = df.loc[:, ["Region","Position", "SRR"]]
    ed = df.loc[:, ["Region","Position", "Editing"]]
    del df
    srr = srr.groupby(by=["Region","Position"]).sum()
    srr.reset_index(drop=False, inplace=True)
    ed = ed.groupby(by=["Region","Position"]).median()
    ed.reset_index(drop=False, inplace=True)
    df = ed.merge(srr, how="inner", on=["Region","Position"])
    del srr, ed
    df.loc[:, "N_Tissue_GTEx_Samples"] = len(paths)
    df.loc[:, "N_Editing_GTEx_Samples"] = df.loc[:, "SRR"].apply(nsamples)
    df.loc[:, "Samples"] = df.loc[:, "N_Editing_GTEx_Samples"]/df.loc[:, "N_Tissue_GTEx_Samples"]
    df.loc[:, "Samples"] = 100*df.loc[:, "Samples"]
    df = df.round({"Samples": 2})
    df.loc[:, "SRR_Samples_%"] = df.loc[:, "Samples"].apply(tostr)
    df.loc[:, "Median_Editing_Frequency_%"] = df.loc[:, "Editing"].apply(tostr)
    df.loc[:, "Tissue"] = tissue+","
    df.drop(["Editing", "Samples"], axis=1, inplace=True)
    df.loc[:, "SRR"] = "[" + df.loc[:, "SRR"] + "],"
    print(f"{tissue} completed", flush=True)
    return df


#global variables assignment
outputpath = sys.argv[1]

tissues = ['ovary','breast','nerve','testis','brain','fallopian_tube','lung','liver','adrenal_gland','spleen',
           'colon','adipose_tissue','blood_vessel','pituitary','salivary_gland','kidney','stomach','pancreas',
           'esophagus','skin','small_intestine','heart','uterus','thyroid','muscle','cervix_uteri','vagina',
           'prostate','blood','bladder']

cols = ["Region","Position", "Reference", "BaseCount[A,C,G,T]", "SRR"]

table = []

#execution parallelization on multiple cores
with Pool(30) as pool:
    for result in pool.map(returntissue, tissues):
        table.append(result)

#result saving
table = pd.concat(table, axis=0)
table = table.groupby(by=["Region","Position"]).sum()
table.reset_index(drop=False, inplace=True)
table.to_csv(os.path.join(outputpath, "Candidate_Sites_FinalTable.txt"), sep="\t", index=None)
