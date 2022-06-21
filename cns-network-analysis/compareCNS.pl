#!/usr/bin/perl

# compareCNS 5/12/2022
# by Amber de Neve
# This script (compareCNS.pl) builds a blast database of every CNS (Conserved Non-coding Sequence) of every species in compareSpecies.csv. It then blasts each single CNS against every other CNS in the database. 
# It generates an adjacency list, AdjacencyList.txt, which is used to build the graph in the analyzeCNSNetwork.py script.
#
# Needs the following files:
# 1) compareSpecies.csv - a text file with a comma seperating each input species CNS file. The species name must match the name of the reference genome in the genome_database.csv. For example, the maize Poaceae CNSs would be specified with ZmB73NAM5.
# 2) input CNS csv files (multiple) - Conservatory generates a final csv file in the Conservatory/CNS folder called "<Family>.csv". Rename this to the refrence genome name (for example, ZmB73NAM5.csv). Put every input CNS file into the ConservatoryCompare/CNS_input folder.
# 3) Genome fasta files for all refrence species, unzipped. Should be placed in the ConservatoryCompare/genomes directory.

# perl setup
use POSIX;
use strict;
use Getopt::Long;
use Cwd 'abs_path';

# Step 1: Directory and file setup, check that files exist
print localtime() . ": Step 1 ~ Checking file structure and building indicies.\n";
my $conservatoryDir = abs_path(".");
my $compareSpeciesFile = $conservatoryDir . "/compareSpecies.csv";
die "ERROR: Cannot find species input file ($compareSpeciesFile). See compareSpeciesExample.csv file for input specifications.\n" unless -e $compareSpeciesFile;

# get species and check that each CNS input exists, that genomes exist and that they are indexed.
open(my $compareSpecies, "<", $compareSpeciesFile);
chomp(my $compareSpeciesLine = <$compareSpecies>);
my @species = split /,/, $compareSpeciesLine;

foreach (@species){
    my $cnsDatabaseFile = $conservatoryDir . "/CNS_input/" . $_ . ".csv";
    die "ERROR: Cannot find CNS Database file ($cnsDatabaseFile) in the CNS_input folder.\n" unless -e $cnsDatabaseFile;
    my $speciesGenome = $conservatoryDir . "/genomes/" . $_ . ".fasta";
    die "ERROR: Cannot find genome file ($speciesGenome) in the genomes folder. Is it still zipped?\n" unless -e $speciesGenome;
    my $speciesGenomeIndex = $conservatoryDir . "/genomes/" . $_ . ".fasta.fai";
    # build genome index if it doesn't exist
    if(! -e $speciesGenomeIndex){
        system("samtools faidx $speciesGenome");
    }
}
close($compareSpecies);

my $cnsAllFasta = $conservatoryDir . "/CNS_input/combinedCNS.fasta";
my $allCNSDatabaseFile = $conservatoryDir . "/CNS_input/combinedCNS.csv";

# check if Step 2 has already been completed. If so, skip.
if(! -e $cnsAllFasta && ! -e $allCNSDatabaseFile){

    # Step 2: Set up CNS blast database. Goes through each species and reads in each line of its CNS csv file.
    print localtime() . ": Step 2 ~ Building CNS sequence databases.\n";
    foreach (@species){
        print localtime() . ": Building CNS sequence database for $_.\n";
        my $cnsDatabaseFile = $conservatoryDir . "/CNS_input/" . $_ . ".csv";

        open(my $cnsDatabase, "<", $cnsDatabaseFile) or die "Can't open file: $cnsDatabaseFile\n";
        while(my $curLine = <$cnsDatabase>){
            chomp $curLine;
            my ($cnsGene, $cnsRelativeStart, $cnsRelativeEnd, $score, $relativePosition, $chromosome, $cnsStart, $cnsEnd, $cnsType) = split /,/, $curLine;
            
            # generate the cnsID by combining the gene name plus the relative start and stop positions
            my $cnsID = "$cnsGene" . "_" . "$cnsRelativeStart" . "_" . "$cnsRelativeEnd";
            # specify the file name for the CNS fasta file
            my $cnsFasta = $conservatoryDir . "/CNS_input/" . $_ . ".cns.fasta";
            my $region = $chromosome . ":" . $cnsStart . "-" . $cnsEnd;
            # go to the refrence genome and extract the CNS region to the new fasta file
            # note: There might be an error here if the chromosome annotation in the cnsDatabaseFile doesn't match the chromosome annotation in the index.
            my $speciesGenome = $conservatoryDir . "/genomes/" . $_ . ".fasta";
            system("samtools faidx $speciesGenome $region >> $cnsFasta 2>&1 ");
            # replace the chromosome position with the cnsID in the new fasta file
            system("sed -i '/>/s/$region/$cnsID/' $cnsFasta");
        }
        close($cnsDatabase);
    }

    # Step 3: Combine CNS fastas & CNS databases for all species, and build the blast database
    print localtime() . ": Step 3 ~ Combining CNS sequences for all species and building the blast database.\n";
    system("cat CNS_input/*.cns.fasta > CNS_input/combinedCNS.fasta");
    system("cat CNS_input/*.csv > CNS_input/combinedCNS.csv");

    # if the CNS_input/cns_blastdb directory doesn't exist, make it
    if(! -e "$conservatoryDir/CNS_input/cns_blastdb"){
        system("mkdir $conservatoryDir/CNS_input/cns_blastdb");
    }

    # make blast database
    system("makeblastdb -in $cnsAllFasta -input_type fasta -dbtype nucl -out \"$conservatoryDir/CNS_input/cns_blastdb/allCNS\"" );

}else{
    print localtime() . ": SKIPPED Step 2 ~ Building CNS sequence databases.\nSKIPPED Step 3 ~ Combining CNS sequences for all species and building the blast database.\nDelete combinedCNS.fasta and combinedCNS.csv to repeat steps 2 and 3.\n";
}


# Step 4: Blast every CNS against every other CNS. Keep those above a certain score. Write the cnsID1, cnsID2, and similarity score to the final adjacencyList.txt file.
print localtime() . ": Step 4 ~ Starting CNS reciprocal blast.\n";
my $adjacencyListFile = $conservatoryDir . "/adjacencyList.txt";

open(my $adjacenyList, ">", $adjacencyListFile) or die "Can't open file: $adjacencyListFile\n"; 
open(my $cnsDatabase, "<", $allCNSDatabaseFile) or die "Can't open file: $allCNSDatabaseFile\n";
while(my $curCNSline = <$cnsDatabase>){
    chomp $curCNSline;
    my ($cnsGene, $cnsRelativeStart, $cnsRelativeEnd, $score, $relativePosition, $chromosome, $cnsStart, $cnsEnd, $cnsType) = split /,/, $curCNSline;
    # generate the cnsID by combining the gene name plus the relative start and stop positions
    my $cnsID = "$cnsGene" . "_" . "$cnsRelativeStart" . "_" . "$cnsRelativeEnd";
    # make temporary files to hold the current cnsID fasta sequence and the current cnsID blast results.
    my $currentCNSseq = $conservatoryDir . "/curCNSseq.tmp";
    my $tmpBlastResults = "$conservatoryDir/CNS_input/" . "$cnsGene" . "blast.tmp";
    # find the cnsID sequence in the combined fasta
    system("samtools faidx $cnsAllFasta $cnsID >> $currentCNSseq");
    # blast the current CNS sequence against the CNS database. 35 max hits is set. E-value cut off is at 15.
    system("blastn -query $currentCNSseq -db $conservatoryDir/CNS_input/cns_blastdb/allCNS -out $tmpBlastResults -evalue 15 -outfmt 6 -max_target_seqs 35");

    # open the temp blast results file, remove any self hits and add the high similarity matches to the adjacency list
    open(my $tmpResults, "<", $tmpBlastResults);
    while(my $curLine = <$tmpResults>){
        # note: can remove downstream to upstream identical hits here in the future
        chomp $curLine;
        # note: _q stands for query, _s stands for subject
        my ($query, $subject, $identical_percentage, $overlap, $mismatches, $gaps, $start_q, $end_q, $start_s, $end_s, $evalue, $bitscore) = split /\t/, $curLine;
        #print localtime() . ": Filtering blast results for " . "$query" . " and " . "$subject" . "\n";
        if ($query eq $subject){
            #print "Skip self hit.\n";
            next;
        }
        # calculate the opposite bitscore. Larger bitscores mean higher similarity, but we need to find the shortest paths.
        my $oppBitscore = floor( (1 / $bitscore) * 100000);
        # format input for networkX graph read in
        my $str = "$query" . " " . "$subject" . " " . "$oppBitscore" . "\n";
        print "Writing to adjacency list: " . "$str";
        print $adjacenyList $str;
    }
    close($tmpResults) or die "Can't close file: $tmpResults";
    # delete temporary files
    unlink($currentCNSseq);
    unlink($tmpBlastResults);
}
close($cnsDatabase);
close($adjacenyList);

my $nodeNum = `wc -l < $allCNSDatabaseFile`;
my $edgeNum = `wc -l < $adjacencyListFile`;
chomp $nodeNum;
chomp $edgeNum;
print "Done. Graph has $nodeNum nodes and $edgeNum edges.\n";

