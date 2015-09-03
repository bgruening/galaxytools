#! /usr/bin/perl -w
#
#
# Copyright (C) 2009 by Pathogene Group, Sanger Center
#
# Author: JIT
# Description: 
#		a wrapper for all the walk stages
#		parameters designed for Bacteria - should be able to run in every stage
#		
#		tested on Salmonella, C. difficile
#		warning: need to be careful how many files are produced
#
#
#


use strict;
use warnings;
use Getopt::Long;

my $PI = `echo $$` ;	chomp($PI) ;

my $path = '' ; 
if ($0 =~ /(\S+)\/image.pl/ ) {
	#print "The path is $1\n" ; 	
	$path = $1 ; 
}

#my $IMAGE_ver = '1' ;
#my $aligner = 'smalt' ;




#required parameters
my $solexa_prefix = '';
#my $ram = 4;
my $gzipped = 0;
my $iteration = '';
my $iteration_all = '';
my $tmp_file = 'yes';
my $reference = '' ;
my $dir_prefix = '' ;

#predefined parameters
my $start = 'yes' ;

my $ref = 'image.contigs.fa' ;
#my $solexa_dir = 'Solexa_reads';
#my $reads_subset = 500000 ; # number of reads to partition

#my $bwa_index = 'is' ;

my $insertsize = 600 ;
my $range1 = 1   ;
my $range2 = $insertsize ;
my $kmer = 31 ;
#my $expcov = 60 ;
#my $covcut = 10 ;

my $toignore = 10 ;
my $overhang = 50 ;

my $velvet_insert_size = 500 ;

my $velvet_option = 0 ; 

#my $walk_mem = 4 ;
#my $walkqueue = 'long' ;

my $SMALT_minScore = 0 ;


GetOptions (	'prefix=s' => \$solexa_prefix, 
		'scaffolds=s' => \$velvet_option,
#		'ssaharam=i' => \$ram, 		
		'iteration=i'    => \$iteration,
		'all_iteration=i'     => \$iteration_all,
#    	'rm_tmp=s' => \$tmp_file,
		'kmer=i' => \$kmer,
#		'expcov=i' => \$expcov,
#		'covcut=i' => \$covcut,
		'toignore=i' => \$toignore,
		'vel_ins_len=i' => \$velvet_insert_size,
#		'walkram=i' => \$walk_mem,
#		'aligner_aligner=s' => \$aligner, 
#		'aligner_ram=i' => \$ram,
#		'aligner_reads_partition=i' => \$reads_subset,
#		'ssaha_queue=s' => \$ssahaqueue,
#		'walk_queue=s' => \$walkqueue,
		'overhang=i' => \$overhang,
		'reference=s' => \$reference,
		'dir_prefix=s' => \$dir_prefix,
#		'bwa_index=s' => \$bwa_index, 
#		'version=i' => \$IMAGE_ver , 

                'smalt_minscore=i' => \$SMALT_minScore, 


) or die "Incorrect usage!\n";


my $inappropriate = 0 ;

$inappropriate = 1 if !$solexa_prefix ;
$inappropriate = 1 if !$iteration ;
$inappropriate = 1 if !$iteration_all ;
$inappropriate = 1 if !$dir_prefix ;



#my $reads = $reads_subset * 4;


if ($inappropriate == 1) { 	

	print "\n\nIMAGE: The ultimate gap close pipeline\n" ;
	print "Author: jit\@sanger.ac.uk\n\n" ;
	
	print "Usage (velvet assembly, with the output contigs.fa ):\n \t/absolute_path_of_image/image.pl -scaffolds contigs.fa -prefix 76bp -iteration 1 -all_iteration 10 -dir_prefix ite\n\n" ; 
	print "Usage (scaffolded assembly):\n \t/absolute_path_of_image/image.pl -scaffolds scaffolds.fa -prefix 76bp -iteration 1 -all_iteration 10 -dir_prefix ite\n\n" ;
#	print "Usage (newbler assembly):\n \t/absolute_path_of_image/image.pl -454 \n\n" ; 
	print "Usage (preformatted): \n" ; 
	
	print "\tYou only need: scaffolds.fa\n" ; 
	print "\tFor customised scaffolds info need two files: a contigs.fa.original and read.placed.original (Please see manual)\n" ; 

	print "\tTo run provided example:\n\t\t/absolute_path_of_image/image.pl -prefix 76bp -iteration 1 -all_iteration 10 -dir_prefix ite\n\n" ; 


	print "\nCompulsory parameters:		-454 (for newbler assembly)\n" ;
	print "                                -velvet (for velvet assembly)\n" ; 
	print "                                -prefix 2540_5   	(Illumina prefix lane)\n" ;
	print "				-iteration 1     	(current iteration)\n" ;
	print "				-all_iteration 5 	(number of iterations to be run after 1 iteration)\n" ;
	print "                		-dir_prefix ite  	(prefix for the iterations)\n" ;
	print "\n\noptional parameters: 		-kmer 31	 	(kmer option in velvet)\n" ;
	print "				-toignore 10     	(exclude extension of less than n bases)\n" ;
	print "				-overhang 50 	 	(to not extend if overhang has > 50bp) \n" ;
	print "				-vel_ins_len 500 	(insert length specified in velvet)\n" ;
        print "                                -smalt_minScore 0       (used for -m option in smalt)\n" ;

#	print "				-reference blah.fasta 	 (optional; if a fasta is provided then nucmer will be used)\n" ;
#	print "                         	-bwa_index is 	 	(used for different parameters of bwa)\n\n" ;

exit 0; 

}









# check if the iteration is the first or not
if ( $iteration eq '1'  ) {
    $start = 'yes' ;
}
else {
    $start = 'no' ;
}

if ( $velvet_option && $start eq 'yes' ) {

    system("cp $velvet_option $velvet_option.bak") ;
    system("$path/get_fasta_containNs_toIMAGE_format.pl $velvet_option.bak 100 100") ;
    system("ln -s $velvet_option.bak.above100.read.placed image.read.placed") ;
    system("ln -s $velvet_option.bak.above100.fasta image.contigs.fa") ;
    system("rm tmp.fa") ; 

}
elsif ($start eq 'yes') {

    system("$path/01_prepare_new_read_placed.pl read.placed.original") ;
    system("$path/02_prepare_new_contigs_fa.pl image.read.placed contigs.fa.original") ;
    #system("ln -s read.placed.new image.read.placed") ;
    #system("ln -s contigs.fa.new image.contigs.fa") ;

}





print '------------------------------------------------' ;
print "\nIMAGE parameters\n\n" ;
print "\nProcess ID      = $PI\n" ;
print "directory prefix  = $dir_prefix\n\n" ;
print "reference         = $ref \n";
print "solexa_prefix     = $solexa_prefix \n" ;
print "insert size       = $insertsize\n" ;
print "toignore          = $toignore\n" ;
print "overhang          = $overhang\n" ;
#print "aligner           = $aligner\n" ;
print "velvet ins_len	 = $velvet_insert_size\n" ;
print "velvet kmer	 = $kmer\n" ;
print "\n" ;
print "current iteration = $iteration\n" ;
print "iteration left    = $iteration_all\n" ;
#print "version    	  = $IMAGE_ver\n" ;
#print "bwa_index          = $bwa_index\n" ;
print '------------------------------------------------' . "\n\n";


open my $out , '>>' , "submit.$dir_prefix.log" or die "cannot open.." ;
print $out '------------------------------------------------' ;
print $out "\nWalk parallelise scripts\n\n" ;
print $out "\nProcess ID      = $PI\n" ;
print $out "directory prefix  = $dir_prefix\n\n" ;
print $out "reference         = $ref \n";
print $out "solexa_prefix     = $solexa_prefix \n" ;
#print $out "read chunk size   = $reads lines \n" ;
print $out "insert size       = $insertsize\n" ;
print $out "toignore          = $toignore\n" ;
print $out "overhang          = $overhang\n" ;
#print $out "aligner           = $aligner\n" ;
print $out "velvet ins_len    = $velvet_insert_size\n" ;
print $out "velvet kmer	= $kmer\n" ;
print $out "\n" ;
print $out "current iteration = $iteration\n" ;
print $out "iteration left    = $iteration_all\n" ;
#print $out "bwa_index          = $bwa_index\n" ;
print $out '------------------------------------------------' . "\n\n";


print $out "WHAT DID I JUST DO (compulsory parameters) :\n" ;
print "WHAT DID I JUST DO (compulsory parameters) :\n" ;



######################################
# checking all the paths
print "checking all the paths..\n" ;
#checkpath("ssaha2Build") ;


checkpath("$path/ssaha2") ;
checkpath("$path/smalt_x86_64") ; 
#checkpath("bwa") ;


checkpath("$path/velveth") ;
checkpath("$path/velvetg") ;
checkpath("$path/image.pl") ;

checkpath("$path/10_20_smalt_and_gather.pl") ; 
#checkpath("$path/10_20_ssaha_and_gather.pl ") ;
#checkpath("$path/10_20_bwa_and_gather.pl ") ;

checkpath("$path/01_prepare_new_read_placed.pl ") ;
checkpath("$path/check_velvet_contig.pl") ;
checkpath("$path/60_new_iteration_raw.pl") ;


print "path checking successfully.. next!\n\n" ;


####                                                                                                                                                                                                                                                                          

#######################################
# some initial processing

if ($start eq 'yes') {


#    system("$path/01_prepare_new_read_placed.pl read.placed.original") ;
#    system("$path/02_prepare_new_contigs_fa.pl image.read.placed contigs.fa.original") ;
#    system("ln -s read.placed.new image.read.placed") ;
#    system("ln -s contigs.fa.new image.contigs.fa") ;
    

       
	my $initial_dir = "$dir_prefix$iteration" ;
	mkdir $initial_dir || die "can not create inital directory $initial_dir" ;

	print "FIRST RUN: $initial_dir created\n" ;

	my $command_tmp = 'sed \'s/>//p;D\' image.contigs.fa | awk \'{print $1 "\t0\t" $1 "\t0\t"}\' > image.insert_size' ;
	system ("$command_tmp") ;
	chdir $initial_dir ;

	system("ln -s ../image.read.placed") ;
	system("ln -s ../image.contigs.fa") ;
	system("ln -s ../image.insert_size") ;
}



system("$path/contig_length_fasta.pl image.contigs.fa > contigs.fa.len.txt") ;
system("mkdir bridge") ;



##########################################################################################################
#stage 1 and 2


#aligner component


mkdir "sam" ;
chdir "sam" ;


print "smalt built on reference.. submit\n" ;
    
#index
system("$path/smalt_x86_64 index -k 13 -s 2 ref.fa ../$ref") ; 


# start gene models
system("$path/10_20_smalt_and_gather.pl 0 ../../$solexa_prefix\_1.fastq ../../$solexa_prefix\_2.fastq $insertsize final.sam $SMALT_minScore")  ;
chdir("../") ;

# all obsolete, using only smalt now...
#if ($aligner eq 'ssaha') {
#	print "ssahaBuild on reference!\n" ;
#	system("ssaha2Build -solexa -save hs36k13s2 ../$ref") ;
#	if ( -e "../../$solexa_prefix\_1.fastq.gz") {
#	    system("gunzip ../../$solexa_prefix\_1.fastq.gz ../../$solexa_prefix\_1.fastq.gz") ;
#	}
#	system("10_20_ssaha_and_gather.pl ../../$solexa_prefix\_1.fastq ../../$solexa_prefix\_2.fastq $insertsize final.sam") ;
#	chdir("../") ; 
#}
#elsif ($aligner eq 'bwa') {
#	print "bwa index on reference!\n" ;
#	print "index parameter: $bwa_index\n" ;
#	system("bwa index -a $bwa_index -p reference ../$ref") ;
#	system("10_20_bwa_and_gather.pl ../../$solexa_prefix\_1.fastq ../../$solexa_prefix\_2.fastq $insertsize final.sam") ;
#	chdir("../") ; 
#}


##########################################################################################################
# stage 3
# the actual walking...



	# merge the list file!
	system("cat bridge/*.list > final.list") ;
	



	if (-e "Illumina_concensus.fa") {
		system("rm Illumina_concensus.fa") ;
	}
	
	mkdir "tmp" or die "can not create tmp\n" ;
	mkdir "velvet_$kmer.auto" or die "can not create velvet_$kmer.auto\n" ;
	mkdir "walks" or die "can not create walks\n" ;

	
	# walk 1
	system("$path/30_walk_stage1_v2.pl kmer$kmer $insertsize $toignore $overhang $velvet_insert_size") ;


	# walk2 stage																	         
	system("$path/40_walk_stage2_v2.pl $insertsize $toignore $overhang > walk2.auto") ;

	# walk3 	stage here?
	system("$path/50_walk_stage3_v2.pl > walk3.auto") ;





##########################################################################################################
# start next iteration
# if there are still some iterations left
$iteration_all-- ;

if ($iteration_all != 0 ) {

	
	chdir("../") ;

	unless ( -e "partitioned_1.fastq" ) {
	    my $path_partition = "$dir_prefix" . "1" ; 
	    system("ln -s $path_partition\/partitioned_1.fastq") ; 
            system("ln -s $path_partition\/partitioned_2.fastq") ;

	}
	    
	system("$path/60_new_iteration_raw.pl $dir_prefix $dir_prefix" . ($iteration+1) . " $dir_prefix" . "$iteration partitioned ". ($iteration+1) .  " $iteration_all $kmer $toignore $overhang $velvet_insert_size");

	


	
}

close($out) ;




print "all done! all done!\n" ;
















sub checkpath {
my $program = shift ;

	if ( `which $program` ) {
		print "found! $program path: " , `which $program` ;
	}
	else {
		print "$program path not found!\n" ;
		exit ;
	}
}

