#!/lustrehome/pietrolucamazzacuva/AnaConda/envs/nndev/bin/python


"""\
This script is used to interpolate the candidates A-to-I RNA editing sites data with the Strong Binders data by means of genomic coordinates and store the result in a single table for each original strong binders table.

Usage: python3 Strong_Binders_AtoI_Candidate_sites_Aggregation.py inputpath outputpath tableprefix 

inputpath: path to the folder containing the AtoI_Candidate_Sites_Data_Aggregation.py output table
outputpath: path to the folder used store this script output and the Strong_Binders_Genomic_Coordinates_Liftover_Formatting.py output
tableprefix: same table prefix used in the Strong_Binders_Genomic_Coordinates_Liftover_Formatting.py script
 

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
import os, sys
import pandas as pd

#functions to retrieve genomic coordinates from in house format
def returnregion(x):
    return x.split("_")[0]

def returnposition(x):
    return int(x.split("_")[1])

#functions to retrive A-to-I editing sites presence over the GTEx tissues
def returnhealthy(x):
    y = x.split(",")
    y = [i for i in y if i != "pancreas"]
    return "yes" if len(y) > 0 else "no"

def returnpancreas(x):
    return "yes" if "pancreas" in x else "no"

#global variables assigment
inputpath, outputpath, tableprefix = sys.argv[1], sys.argv[2], sys.argv[3]

table = pd.read_table(os.path.join(inputpath, "Candidate_Sites_FinalTable.txt"), sep="\t")
table.drop(["N_Tissue_GTEx_Samples", "N_Editing_GTEx_Samples"], axis=1, inplace=True)

#strong binders data and candidate A-to-I editing sites data interpolation
for i in ["I", "II"]:
    df = pd.read_table(os.path.join(outputpath, f"{tableprefix}{i}_hg19_cohordinates.txt"))
    df.loc[:, "Region"] = df.loc[:, "Mutation_Position_Hg19"].apply(returnregion)
    df.loc[:, "Position"] = df.loc[:, "Mutation_Position_Hg19"].apply(returnposition)
    df2 = pd.read_table(os.path.join(outputpath, f"{tableprefix}{i}.tsv"))
    cols = df2.columns.tolist()
    cols[4] = "Mutation_Position_Hg38"
    df2.columns = cols
    df2 = pd.concat([df2, df], axis=1)
    df2 = df2.merge(table, how="left", on=["Region", "Position"])
    cols = ['Evaluated_Mutant_Peptide', 'Gene', 'NM_ID', 'Change', 'Mutation_Position_Hg38', 'Mutation_Position_Hg19', 
            'Peptide_Sequence', 'Samples', 'N_Samples', 'SRR', 'SRR_Samples_%', 'Median_Editing_Frequency_%', 'Tissue']
    df2 = df2.loc[:, cols]
    df2.loc[:, "Healthy_Tissues"] = df2.loc[:, 'Tissue'].astype("str")
    df2.drop(["Tissue"], axis=1, inplace=True)
    df2.loc[:, "Also_In_Healthy_Pancreas"] = df2.loc[:, "Healthy_Tissues"].apply(returnpancreas)
    df2.loc[:, "Also_In_Others_Healthy_Tissues"] = df2.loc[:, "Healthy_Tissues"].apply(returnhealthy) 
    cols = ['Evaluated_Mutant_Peptide', 'Gene', 'NM_ID', 'Change', 'Mutation_Position_Hg38', 'Mutation_Position_Hg19', 
            'Peptide_Sequence', 'Samples', 'N_Samples', 'SRR', 'SRR_Samples_%', 'Median_Editing_Frequency_%', 
            "Also_In_Healthy_Pancreas", "Also_In_Others_Healthy_Tissues", "Healthy_Tissues"]
    df2 = df2.loc[:, cols]
    df2.to_csv(os.path.join(outputpath, f"{tableprefix}{i}_with_GTEx_data.tsv"), sep="\t", index=None)
