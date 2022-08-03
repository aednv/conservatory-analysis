#import libraries
import networkx as nx
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import sys

#read in the cns node graph with blast weights from the adjacency list
G = nx.read_weighted_edgelist("adjacencyList.txt")
sys.stdout.write("CNS nodes: " + str(G.order()) + "\n") #nodes
sys.stdout.write("connections: " + str(G.size()) + "\n") #edges

#add species attribute to all cns nodes (first 2 letters of GeneID)
for n in G.nodes():
    G.nodes[n]['species']=n[0:2]

#add gene attributes to all cns nodes, gene length based on species
for n in G.nodes():
    if G.nodes[n]['species'] == 'Zm': #maize
        G.nodes[n]['gene'] = n[0:20]
    elif G.nodes[n]['species'] == 'Os': #rice
        G.nodes[n]['gene'] = n[0:16]
    elif G.nodes[n]['species'] == 'Br': #brachy
        G.nodes[n]['gene'] = n[0:12]
    elif G.nodes[n]['species'] == 'Se': #setaria
        G.nodes[n]['gene'] = n[0:14]
        
directory="orthologs" #set otholog.csv directory. Files are generated from conservatory and in conservatory/genomes. Modified the filesnames to add the ref genome to the beginning.
import os

for filename in os.listdir(directory):
    sys.stdout.write("Current file" + str(filename) + "\n")
    with open(str(directory+"/"+filename), "r") as fp: #open ortholog.csv file
        numLines = sum(1 for line in fp)
        #print(numLines)
        quarterDone = int(numLines/4)
        #print(quarterDone)
        lineCount = 1
    with open(str(directory+"/"+filename), "r") as fp:
        for line in fp:
            lineArray = line.split(",")
            targetGene = str(lineArray[0])
            refGene = str(lineArray[1])
            #get CNS nodes from graph if they exist
            refCNSs = [n for n,atr in G.nodes(data=True) if atr['gene'] == refGene]
            targetCNSs = [n for n,atr in G.nodes(data=True) if atr['gene'] == targetGene]
            if targetCNSs != [] and refCNSs != []:
                for cns_R in refCNSs:
                    for cns_T in targetCNSs:
                        if G.has_edge(cns_R, cns_T): #check for cns edge
                            #print("Connecting", cns_R, cns_T)
                            G.edges[cns_R, cns_T]["orthologous"] = "True"
            
            if lineCount == quarterDone:
                sys.stdout.write("25%\n")
            if lineCount == quarterDone*2:
                sys.stdout.write("50%\n")
            if lineCount == quarterDone*3:
                sys.stdout.write("75%\n")
            if lineCount == numLines:
                sys.stdout.write("100%\n")
            lineCount+=1
    nx.write_gml(G, "cns_graph_orthologs.gml") #save after processing each file

