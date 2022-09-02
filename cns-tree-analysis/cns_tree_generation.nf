#!/usr/bin/env nextflow

//syntax version
nextflow.enable.dsl=2

//inputs

//discover new homologous genes and expand the tree. If false, the gene list provided in params.startGenes will be used directly to build the tree
params.discoverGenes = "TRUE"

//path to fasta file of starting genes with their protein sequences. One of these should be an outgroup gene
params.startGenes = "./gene_data/ra3_start_genes.fa"

//name of the outgroup gene, as annotated in the protein database
params.outgroup = "Aco009575"

//path to the protein database file
params.dbPath = "/project/uma_madelaine_bartlett/AmberDeneve/conservatory/conservatory/genomes/Poaceae/Poaceae.proteins.combined.fa"

//path to the cns mapping resouce folder
params.cnsPath = file("./gene_data/cns_mapping_resources/")

//gene id of your main gene of interest, as annotated in the protein database
//maize ramosa3 (ra3)
params.mainGene = "Zm00001eb327910_P002"

//path to phytools conda enviornment
params.phytoolsEnv  = file("./phytoolsConda")

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
	graph-cns-tree.R $tree $cnsTable $params.outgroup $params.mainGene
	"""
}
	
/*
 * Workflow
 */
workflow {
	//create channel
	startGenes_ch = Channel.fromPath( params.startGenes, checkIfExists: true )
	
	//build tree and map CNS data
	if( params.discoverGenes == "TRUE" ) {
		//discover homologous genes
		findGenesRoundOne( startGenes_ch )
		findGenesRoundTwo( findGenesRoundOne.out.genes )
		//build gene tree
		buildTree( findGenesRoundTwo.out.genesFinal )
		mapCnsData( buildTree.out.treeGenesTrimmed )
	}
	else if( params.discoverGenes == "FALSE" ) {
		//build gene tree
		buildTree ( startGenes_ch )
		mapCnsData( buildTree.out.treeGenesTrimmed )
	}
	else {
		error "Invalid discoverGenes parameter."
	}
	
	//graph in R
	graphCnsTree( buildTree.out.finalTree, mapCnsData.out.cnsTable )
}
