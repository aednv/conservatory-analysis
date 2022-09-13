#!/usr/bin/env nextflow

//syntax version
nextflow.enable.dsl=2

//Inputs

//params.mainGene(string) - Required user input - Gene ID of your main gene of interest, as annotated in the protein database.
params.mainGene = ""

//params.outgroup(string) - Required user input - Name of the outgroup gene, as annotated in the protein database.
params.outgroup = ""

//params.noSearch(string = "true"/"false") - If false (default), search for new genes. If true, the gene list provided in params.startGenes will be used directly to build the tree.
params.noSearch = "false"

//params.startGenes(string) - Path to fasta file of starting genes with their protein sequences. One of these sequences should be the outgroup gene. This parameter is only used if params.noSearch is set to true.
params.startGenes = ""

//params.dbPath(string) - Path to the protein database directory.
params.dbPath = file("./gene_data/protein_db/")

//params.cnsPath(string) - Path to the cns mapping resource folder. Contains the CNS map files (All.merged.cns.csv and All.merged.cns.map).
params.cnsPath = file("./gene_data/cns_mapping_resources/")

//params.phytoolsEnv(string) - Path to phytools conda enviornment.
params.phytoolsEnv  = file("./phytoolsConda")

//params.colorful(string = "true"/"false") - If false (default), makes a greyscale graph. If true, makes a colorful graph.
params.colorful = "false"

/*
 * Processes
 */

//scaleGenomeInput searches the CNS map files and makes a list of each homologous gene with conserved CNS regions. It records which family the species is from and adds that to
//the list of family genomes to be included for the downstream hmm search. The inputs are params.mainGene, and the outputs are (1) a text file with the list of homologous genes, newline seperated and 
//(2) a list of genome directories to be used in the hmm analysis.  
process scaleGenomeInput {
	tag {"scaleGenomeInput $mainGene"}
	executor 'lsf'
	queue 'short'
	cpus 1
	time '1h'
	
	publishDir "gene_data", mode: 'copy'
	
	input:
		val mainGene
	
	output:
		path(".txt."), emit : startGenes
		val famList, emit : familyList
	
	"""
	#verify that input gene has conserved CNS regions (that it is in the csv file)
	#check all merged cns. for each matching line in ColG
		#get all matching ortholog genes for every CNS region in borth the merged.map and family.map
		#for any CNS that are not family level, use the genome db csv to link the species to its Family
		#add that family to the family list
	#remove any duplicates in the gene list and the family list
	
	"""
}

//generateFasta builds a fasta file of protein sequences. Inputs: (1) path to a text file list of genes, newline seperated. (2) list of directories making up the protein database.
process generateFasta {
	tag {"generate Fasta $startGenes"}
	executor 'lsf'
	queue 'short'
	cpus 1
	time '1h'
	
	publishDir "gene_data", mode: 'copy'
	
	input:
		path startGenes
		val familyList
	
	output:
		path("*.fa"), emit : startProteinSeqs
	
	"""
	#build a fasta from the gene ids
	#modify from code below
	"""
}

process findGenesRoundOne {
	tag {"findGenesRoundOne $startGenes"}
	executor 'lsf'
	queue 'short'
	cpus 1
	time '2h'
	
	publishDir "gene_data", mode: 'copy'
	
	input:
		path startGenes
	
	output:
		path("*.aln")
		path("*.hmm")
		path("*.roundOne.fa"), emit: genes
		path("*.roundOne.txt")
	
	"""
	#load cluster modules
	module load MAFFT/7.313
	module load hmmer/3.1b2
	module load samtools/1.9
	#align startGenes
	dos2unix $startGenes
	mafft-linsi $startGenes > ${startGenes}.aln
	#make profile hidden markov model
	hmmbuild -o ${params.mainGene}.summary.roundOne.txt ${params.mainGene}.roundOne.hmm ${startGenes}.aln
	#use model to pull top hits out of protein database
	hmmsearch -o ${params.mainGene}.hmm.roundOne.txt --noali --tblout ${params.mainGene}.search.output.roundOne.txt ${params.mainGene}.roundOne.hmm $params.dbPath
	#extract all search matches with a higher score than the outgroup
	sed '1,3d' ${params.mainGene}.search.output.roundOne.txt > ${params.mainGene}.search.output.trimmed.roundOne.txt
	cat ${params.mainGene}.search.output.trimmed.roundOne.txt | while read line
	do
		currentGene=\$(echo \$line | awk '{print \$1}')
		echo \$currentGene
		samtools faidx $params.dbPath \$currentGene >> ${params.mainGene}.roundOne.fa
		if [ \$currentGene == $params.outgroup ]; then break; fi
	done
	"""
}

process findGenesRoundTwo {
	tag {"findGenesRoundTwo $roundOneGenes"}
	executor 'lsf'
	queue 'short'
	cpus 1
	time '2h'
	
	publishDir "gene_data", mode: 'copy'
	
	input:
		path roundOneGenes
	
	output:
		path("*.aln")
		path("*.hmm")
		path("*.roundTwo.fa"), emit: genesFinal
		path("*.roundTwo.txt")
	
	"""
	#load cluster modules
	module load MAFFT/7.313
	module load hmmer/3.1b2
	module load samtools/1.9
	#align round one genes
	mafft-linsi $roundOneGenes > ${roundOneGenes}.aln
	#make profile hidden markov model
	hmmbuild -o ${params.mainGene}.summary.roundTwo.txt ${params.mainGene}.roundTwo.hmm ${roundOneGenes}.aln
	#use model to pull top hits out of protein database
	hmmsearch -o ${params.mainGene}.hmm.roundTwo.txt --noali --tblout ${params.mainGene}.search.output.roundTwo.txt ${params.mainGene}.roundTwo.hmm $params.dbPath
	#extract all search matches with a higher score than the outgroup
	sed '1,3d' ${params.mainGene}.search.output.roundTwo.txt > ${params.mainGene}.search.output.trimmed.roundTwo.txt
	cat ${params.mainGene}.search.output.trimmed.roundTwo.txt | while read line
	do
		currentGene=\$(echo \$line | awk '{print \$1}')
		echo \$currentGene
		samtools faidx $params.dbPath \$currentGene >> ${params.mainGene}.roundTwo.fa
		if [ \$currentGene == $params.outgroup ]; then break; fi
	done
	"""
}

process buildTree {
	tag {"buildTree $treeGenes"}
	executor 'lsf'
	queue 'long'
	clusterOptions '-R "rusage[mem=2000]" "span[hosts=1]"'
	cpus 4
	time '10h'
	
	publishDir "results", mode: 'copy'
	
	input:
		path treeGenes
	
	output:
		path("*.trimmed"), emit: treeGenesTrimmed
		path("*.raxml.supportFBP"), emit: finalTree
		path("*.aln")
		path("*raxml*")
	
	"""
	#remove any duplicate genes from the fasta input (awk command from https://bioinformatics.stackexchange.com/questions/15647/removing-duplicate-fasta-sequences-based-on-headers-with-bash)
	awk '/^>/ { f = !(\$0 in a); a[\$0]++ } f' $treeGenes > ${treeGenes}.trimmed
	#load cluster modules
	module load raxml-ng/0.9.0
	module load MAFFT/7.313
	mafft-linsi ${treeGenes}.trimmed > ${treeGenes}.trimmed.aln
	singularity exec \$RAXMLNGIMG raxml-ng-mpi --all --msa ${treeGenes}.trimmed.aln --model JTT+G --threads 4 --bs-metric fbp,tbe
	"""
}

process mapCnsData {
	tag {"mapCnsData $treeGenesTrimmed"}
	executor 'lsf'
	queue 'short'
	cpus 1
	time '1h'
	
	publishDir "results", mode: 'copy'
	
	input:
		path treeGenesTrimmed
	
	output:
		path("*.csv"), emit: cnsTable
		path("*_TreeGenes.txt")
	
	"""
	module load python3/3.5.0
	module load samtools/1.9
	module load python3/3.5.0_packages/pandas/0.18.0
	#extract just gene names from list
	cat $treeGenesTrimmed | grep '>' | cut -b 2- > ${params.mainGene}_TreeGenes.txt
	cns-tree-generation-auto.py  $params.mainGene $params.cnsPath
	"""
}

process graphCnsTree {
	tag {"graphCnsTree $params.mainGene"}
	executor 'lsf'
	queue 'short'
	cpus 1
	clusterOptions '-R "rusage[mem=10000]" "span[hosts=1]"'
	time '1h'
	
	publishDir "results", mode: 'copy'
	
	input:
		path tree
		path cnsTable
	
	output:
		path("*tree.pdf")
	
	script:
	"""
	module load anaconda3/2019.03
	source activate $params.phytoolsEnv
	graph-cns-tree.R $tree $cnsTable $params.outgroup $params.mainGene $params.colorful
	"""
}
	
/*
 * Workflow
 */
workflow {
	//!!modifications for deep cns searching
	
	//create channel from main reference gene input
	if( params.mainGene == "" ) {
		error "No reference gene input. <params.mainGene> is empty."
	}
	mainGene_ch = Channel.from( params.mainGene )
	 
	//automatically search for genes, build tree, and map CNSs. (default)
	if( params.noSearch == false ) {
		//scale genomes to be included in the hmm search, and generate homologous gene list
		scaleGenomeInput ( mainGene_ch )
		//build a protein fasta file for homologous gene list
		generateFasta ( scaleGenomeInput.out.startGenes, scaleGenomeInput.out.familyList )
		//discover more homologous genes to expand tree
		findGenesRoundOne( generateFasta.out.startGenes, scaleGenomeInput.out.familyList )
		findGenesRoundTwo( findGenesRoundOne.out.genes, scaleGenomeInput.out.familyList )
		//build gene tree
		buildTree( findGenesRoundTwo.out.genesFinal )
		//!!map CNS data chould use the map file now, not the bam file
		mapCnsData( buildTree.out.treeGenesTrimmed )
	}
	//build tree from user input protein fasta and map CNSs.
	else if( params.noSearch == true ) {
		//build gene tree
		buildTree ( startGenes_ch )
		//!!use map file, not bam file
		mapCnsData( buildTree.out.treeGenesTrimmed )
	}
	else {
		error "Invalid noSearch parameter. Use false or true (lowercase)."
	}
	
	//graph in R
	graphCnsTree( buildTree.out.finalTree, mapCnsData.out.cnsTable )
}
