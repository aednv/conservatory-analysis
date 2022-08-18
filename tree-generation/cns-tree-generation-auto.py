#cns-tree-generation-auto.py, Amber, 8/18/2022
#this script combines all of the python and bash scripts in Part2 of ramosa3-tree-generation.ipynb into a single python file.
#this script should be reproducable for mapping different cns conservation onto different gene trees

#set up
import os
import sys
import pandas as pd
from numpy import loadtxt
import re

#read arguments
if (len(sys.argv) != 2):
    print("Expected 1 input, got ", len(sys.argv)-1)
gene = sys.argv[1]

#load samtools from command line
os.system("module load samtools/1.9")

#check that required files are present
os.path.exists('./combinedCNS.csv') #all of the maize CNSs from conservatory/CNS/ concatenated together with cat
os.path.exists('./Poaceae.bam') #the final conserved ortholog alignment file from conservatory/alignments/
os.path.exists('./Poaceae.bam.bai') #the index file for the above from conservatory/alignments/
os.path.exists('./all.orthologs.txt') #all the ortholog.csv files from conservatory/genomes/ concatenated together
os.path.exists('./'+ gene + 'TreeGenes.txt') #TODO the list of genes from your gene tree of interest, change to be based on geneID


def generateDataframe(myGeneID):
    df = pd.read_csv('ra3TreeGeneList2.txt', header=None).drop_duplicates() #read in and remove duplicates
    df.rename(columns={0:'ra3_tree_genes'}, inplace=True)
    df = df.reset_index(drop=True)
    return df

def lines_that_contain(string, fp):
    return [line for line in fp if string in line]

#look up myGeneID in the combinedCNS.csv, and find all of the CNSs for that gene.
def extractAlignments(myGeneID):
    cnsIDs = []
    total_count = 0
    with open("combinedCNS.csv", "r") as fp:
        for line in lines_that_contain(myGeneID, fp):
            lineArray = line.split(",")
            cnsID = str(lineArray[0]) + "_" + str(lineArray[1]) + "_" + str(lineArray[2]) + "_" + str(lineArray[5]) + "_" + str(lineArray[6]) + "_" + str(lineArray[7]) #convert CNS line into cnsID
            print(cnsID)
            cnsIDs.append(cnsID)
            #use samtools to extract alignments to a temp file with cnsID as name
            os.system("samtools view Poaceae.bam " + str(lineArray[5]) + ":" + str(lineArray[6]) + "-" + str(lineArray[7]) + " | awk '{print $1}' | awk -F: '{print $3}' > " + cnsID + ".txt")
            total_count += 1
    print(str(len(cnsIDs)) + " CNS ortholog alignment files generated for " + str(myGeneID))
    return cnsIDs


def dataframeMerge(myGeneID, newColName):
    geneList = loadtxt(str(newColName + ".txt"), comments="#", delimiter="/n", unpack=False, dtype="str")
    df[str(myGeneID + '_includedInAnalysis')]=df['ra3_tree_genes'].str.contains('|'.join(map(re.escape, geneList)))

def falsePosCheck(myGeneID):
    os.system("grep " + myGeneID + " all.orthologs.csv | awk -F, '{print $1}' > " + myGeneID + ".orthologs.txt")
    orthoFileName = myGeneID + ".orthologs"
    dataframeMerge(myGeneID, orthoFileName)

def importCNSs(myGeneID, cnsIDs):
    for cns in cnsIDs:
        dataframeMerge(myGeneID, cns)

def extractSameRefGenes(myGeneID):
    sameRefGenes = df['ra3_tree_genes'].str.contains(myGeneID[0:5])
    return sameRefGenes


#step 5: identify other orthologs from the same refrence genome on the tree. Repeat and add that data to dataframe as well.


#step 6: save as a csv file for import into R for graphing
def saveCSV():
    df.to_csv(myGeneID + "_conservedCNS.csv")


#~~~
#step 1: make a dataframe with all the tree genes
df = generateDataframe(gene)

#step 5: identify other orthologs from the same refrence genome on the tree. Repeat and add that data to dataframe as well.
sameRefGeneList = extractSameRefGenes(gene)

for gene in sameRefGeneList:
    cnsIDList = extractAlignments(gene)
    falsePosCheck(gene)
    importCNSs(gene, cnsIDList)
saveCSV(gene)















