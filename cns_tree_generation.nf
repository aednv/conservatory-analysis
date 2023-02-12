#!/usr/bin/env nextflow

//syntax version
nextflow.enable.dsl=2

//inputs
//don't discover new homologous genes and expand the tree. If true, the gene list provided in params.startGenes will be used directly to build the tree
params.noSearch = "false"
//path to fasta file of starting genes with their protein sequences.
params.startGenes = ""
//name of the outgroup gene, as annotated in the protein database
params.outgroup = ""
//path to the protein database file
params.dbPath = file("./gene_data/protein_database/all.proteins.fasta")
//path to the cns mapping resouce folder
params.cnsPath = file("./gene_data")
//gene id of your main gene of interest, as annotated in the protein database
params.mainGene = ""
//path to phytools/ cns mapping conda environment
params.phytoolsEnv = ""
//path to bioinformatics tools conda
params.bioinfEnv = ""
//use colorful graph option (default is white and grey)
params.colorful = "false"

/*
 * Processes
 */
process findGenesRoundOne {
	tag {"findGenesRoundOne $startGenes"}
	executor 'slurm'
	queue 'cpu'
	cpus 1
	time '2h'
	memory '1 GB'
	
	publishDir "gene_data", mode: 'copy'
	
	input:
		path startGenes
	
	output:
		path("*.aln")
		path("*.hmm")
		path("*.roundOne.fa"), emit: genes
		path("*.roundOne.txt")
	
	"""
	#load conda env
	module load anaconda/2022.10
	conda activate $params.bioinfEnv
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
	executor 'slurm'
	queue 'cpu'
	cpus 1
	time '2h'
	memory '1 GB'

	publishDir "gene_data", mode: 'copy'
	
	input:
		path roundOneGenes
	
	output:
		path("*.aln")
		path("*.hmm")
		path("*.roundTwo.fa"), emit: genesFinal
		path("*.roundTwo.txt")
	
	"""
	#load conda env
	module load anaconda/2022.10
	conda activate $params.bioinfEnv
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
	executor 'slurm'
	queue 'cpu'
	memory '4 GB'
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
	#load conda env
	module load anaconda/2022.10
	conda activate $params.bioinfEnv
	mafft-linsi ${treeGenes}.trimmed > ${treeGenes}.trimmed.aln
	raxml-ng-mpi --all --msa ${treeGenes}.trimmed.aln --model JTT+G --threads 4 --bs-metric fbp,tbe
	"""
}

process mapCnsData {
	tag {"mapCnsData $treeGenesTrimmed"}
	executor 'slurm'
	queue 'cpu'
	cpus 1
	time '1h'
	memory '1 GB'

	publishDir "results", mode: 'copy'
	
	input:
		path treeGenesTrimmed
	
	output:
		path("*.csv"), emit: cnsTable
		path("*_TreeGenes.txt")
	
	"""
	#load conda env
	module load anaconda/2022.10
	conda activate $params.bioinfEnv
	#extract just gene names from list
	cat $treeGenesTrimmed | grep '>' | cut -b 2- > ${params.mainGene}_TreeGenes.txt
	csv-database-cns-tree-generation.py $params.mainGene $params.cnsPath
	"""
}

process graphCnsTree {
	tag {"graphCnsTree $params.mainGene"}
	executor 'slurm'
	queue 'cpu'
	cpus 1
	memory '1 GB'
	time '1h'
	
	publishDir "results", mode: 'copy'
	
	input:
		path tree
		path cnsTable
	
	output:
		path("*tree.pdf")
	
	script:
	"""
	#load conda env
	module load anaconda/2022.10
	conda activate $params.phytoolsEnv
	graph-cns-tree.R $tree $cnsTable $params.outgroup $params.mainGene $params.colorful
	"""
}
	
/*
 * Workflow
 */
workflow {
	//create channel
	startGenes_ch = Channel.fromPath( params.startGenes, checkIfExists: true )
	 
	//build tree and map CNS data
	if( params.noSearch == false ) {
		//discover homologous genes
		findGenesRoundOne( startGenes_ch )
		findGenesRoundTwo( findGenesRoundOne.out.genes )
		//build gene tree
		buildTree( findGenesRoundTwo.out.genesFinal )
		mapCnsData( buildTree.out.treeGenesTrimmed )
	}
	else if( params.noSearch == true ) {
		//build gene tree
		buildTree ( startGenes_ch )
		mapCnsData( buildTree.out.treeGenesTrimmed )
	}
	else {
		error "Invalid noSearch parameter. Use false or true (lowercase)."
	}
	
	//graph in R
	graphCnsTree( buildTree.out.finalTree, mapCnsData.out.cnsTable )
}
