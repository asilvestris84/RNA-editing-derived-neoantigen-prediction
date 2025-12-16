#!/lustrehome/pietrolucamazzacuva/AnaConda/envs/nndev/bin/python


"""\
This script is used to format the inputs and the outputs of the UCSC liftover online tool used to interpolate the candidates A-to-I RNA editing sites data with the Strong Binders data by means of genomic coordinates. 

Usage: Strong_Binders_Genomic_Coordinates_Liftover_Formatting.py path mode tableprefix tablesuffix

path: path to the folder containing the AtoI_Candidate_Sites_Search.py outputs, used also to store this script output
mode: output (if the script is used to format UCSC liftover output) or input (if the script is used to format UCSC liftover input)
tableprefix: common suffix used in the strong binders table names  
tablesuffix: file extention (es: .tsv)

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
import os, sys

#functions to retrive genomic coordinates from an in house format
def returnregion(x):
    return x.split("_")[0]

def returnposition(x):
    return int(x.split("_")[1])

#arguments retrieving
outputpath, mode, tableprefix, tablesuffix  = sys.arg[1], sys.argv[2], sys.argv[3], sys.argv[4]

if mode == "output":
    #UCSC liftover results formating
    for i in ["I", "II"]:
        df = pd.read_table(os.path.join(outputpath, f"liftover_{tableprefix}{i}.bed"), header=None)
        df.iloc[:, 1] = df.iloc[:, 1].astype("str")
        df.loc[:, "Mutation_Position_Hg19"] = df.iloc[:, 0] + "_" + df.iloc[:, 1]
        df.loc[:, "Mutation_Position_Hg19"].to_csv(os.path.join(outputpath, f"{tableprefix}{i}_hg19_cohordinates.txt"), index=None)
else:
    #UCSC liftover inputs formatting
    for i in ["I", "II"]:
        df = pd.read_table(os.path.join(outputpath, f"{tableprefix}{i}{tablesuffix}"), usecols=["Mutation_Position"])
        df.loc[:, "1"] = df.loc[:, "Mutation_Position"].apply(returnregion)
        df.loc[:, "2"] = df.loc[:, "Mutation_Position"].apply(returnposition)
        df.loc[:, "3"] = df.loc[:, "Mutation_Position"].apply(returnposition)
        df.loc[:, ["1", "2", "3"]].to_csv(os.path.join(outputpath, f"liftover_{tableprefix}{i}.txt"), sep="\t", index=None, header=False)
