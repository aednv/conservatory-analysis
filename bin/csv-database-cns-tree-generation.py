#!/usr/bin/env python3

#csv-database-cns-tree-generation.py, Amber, 2/6/2023
#This script maps cns conservation onto different gene trees using the final cns csv files generated from Conservatory.

import os
import sys
import pandas as pd
import numpy as np
import subprocess

#read arguments
if (len(sys.argv) != 3):
    print("Expected 2 input, got ", len(sys.argv)-1)
myGeneID = sys.argv[1] # ex. "Zm003452253"
cnsFilesPath = sys.argv[2] # All.merged.map.csv path

#check that required files are present
mapPath = cnsFilesPath + '/conservatoryV9.5.3d.final.map.csv'
if (os.path.exists(mapPath) == False):
    print("Missing 'conservatoryV9.5.3d.final.map.csv'")
    exit()
if (os.path.exists('./' + myGeneID + '_TreeGenes.txt') == False): #the list of genes from your gene tree of interest
    print("Missing " + myGeneID + "_TreeGenes.txt")
    exit()
allOrthologsPath = cnsFilesPath + '/angiospermsOrthologs.csv'
if (os.path.exists(allOrthologsPath) == False): #all the ortholog.csv files from conservatory/genomes/ concatenated together
    print("Missing 'angiospermsOrthologs.csv'")
    exit()

#functions
print("checkpoint1")
#generateDataframe(gene) makes a new pandas dataframe with the genes in <gene>_TreeGenes.txt as the first column. Input is the geneID. It returns the new dataframe.
def generateDataframe(gene):
    df = pd.read_csv(gene + '_TreeGenes.txt', header=None).drop_duplicates() #read in and remove duplicates
    df.rename(columns={0:'tree_genes'}, inplace=True)
    df = df.reset_index(drop=True)
    return df

#lines_that_contain() is a helper function for extractCnsIDs(). It only returns lines in the filepath (fp) that contain the search string (string).
def lines_that_contain(string, fp):
    return [line for line in fp if string in line]

#extractCnsIDs(gene) takes a gene ID string as input. It looks up the gene in the All.merged.map.csv, and finds all of the CNS regions for that gene.
def extractCnsIDs(gene):
    cnsIDs = []
    total_count = 0
    with open(mapPath, "r") as fp:
        for line in lines_that_contain(gene, fp):
            lineArray = line.split(",")
            cnsIDs.append(str(lineArray[0]))
    cnsIDs = np.unique(cnsIDs)
    print(str(len(cnsIDs)) + " CNS regions found for " + str(gene))
    return cnsIDs

#cnsMap(cnsID, df) makes a new column in the specified dataframe (df) for a cnsID, and marks if that region is conserved for each of the ortholog tree genes.
#The new df column will be named the same as cnsID.
def cnsMap(cnsID, df):
    orthoCNSGeneIDs = []
    with open(mapPath, "r") as fp:
        for line in lines_that_contain(cnsID, fp):
            lineArray = line.split(",")
            orthoCnsGeneIDRaw = str(lineArray[2]) #Cgrandiflora-Cgrandiflora-Cagra.1968s0046
            orthoCNSGeneIDs.append(orthoCnsGeneIDRaw.split("-")[2]) #Cagra.1968s0046
    #orthoCNSGeneIDs = filter(None, orthoCNSGeneIDs) #make sure there are no empty strings
    print(str(len(orthoCNSGeneIDs)) + " orthologous CNS genes found for " + str(cnsID))
    print(orthoCNSGeneIDs)
    df[str(cnsID)]=df['tree_genes'].str.contains('|'.join(orthoCNSGeneIDs))
    updatedDf = df
    print(str(df[cnsID].sum()) + " / " + str(len(orthoCNSGeneIDs)) + " CNS regions mapped to tree")
    return updatedDf

#falsePosCheck(gene) extracts all of the orthologs for (gene) included in the conservatory analysis and saves them in a new file <gene>.orthologs.txt in the tmp directory. It gets all of the orthologs from the all.orthologs.csv file.
#it returns the orthoFileName, to be used with dataframeMerge().
def falsePosCheck(gene, df):
    subprocess.call("grep " + gene + " " + allOrthologsPath + " | awk -F'-' '{print $3}' > " + gene + ".orthologs.txt", shell=True)
    print("check1")
    with open(str(gene + ".orthologs.txt")) as data:
        newGenes = data.read().split('\n')
    print("check2")
    newGenes = filter(None, newGenes) #make sure there are no empty strings
    df[str(falsePosCheck)]=df['tree_genes'].str.contains('|'.join(newGenes))
    updatedDf = df
    print("check3")
    return updatedDf
print("checkpint2")
#~~~~~~~~~~~~~~~~~~
#step 1: make a dataframe with all the tree genes
myTreeDataframe = generateDataframe(myGeneID)
print("checkppoint3")
#step 2: check that all tree genes were included in the conservatory analysis
myTreeDataframe = falsePosCheck(myGeneID, myTreeDataframe)
print("check4")
#step 3: get all cns regions for reference gene
myCnsIDs = extractCnsIDs(myGeneID)
print("checkpoint5")
print(len(myCnsIDs))
#step 4: cycle through each cns region, checking which CNS regions are conserved
for id in myCnsIDs:
    print("Processing " + id)
    myTreeDataframe = cnsMap(id, myTreeDataframe) #merge all checked

#step 5: save final dataframe as csv file for graphing in R
myTreeDataframe.to_csv(myGeneID + "_conservedCNSTable.csv")

