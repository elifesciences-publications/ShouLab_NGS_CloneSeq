#!/usr/bin/perl
# Requirements:
# Have X11 installed
# log into lynx with: ssh -X -u username@lynx.fhcrc.org
# Required programs: samtools, vcftools, GATK, picard, R, pindel, java1.7, mrsfast, splitread, X11 locally	
# Currently installed versions
# samtools: 1.0
# VCFtools (v0.1.12b)
# GATK 3.2-2-gec30cee (also needs Java v1.7)
# Picard 1.102
# mrsfast 2.6.0.4 (Newer versions dont work for splitread!!)
# splitread v??? not the version of the website but obtained from Eichler Lab
# pindel 0.2.5a7
# R version 3.1.1: Libraries: Rsamtools, cn.mops, rtracklayer, BSgenome, zoo , plyr, GenomicRanges + all dependencies 

# OPTIONAL Improvements: 
# add timestamp to files
# Understand & summarize pindel	output
# explain that bam file need to exist already for all cnmops comparisons or change pipeline so it cycles through tools then strains instead of strains then tools.
# get useful splitread output, auto-detect read length and adjust trimming parameter

use strict; use warnings;  use File::Copy; use Getopt::Std; use List::Util qw( min max );

# check that the -s and -p flags were set
our($opt_s, $opt_p, $opt_n); getopts('ns:p:');
unless($opt_p){die "\n\nUSAGE: CloneSeqPipeline.pl -s Strains.txt -p params.txt\n\n"}
unless($opt_s){die "\n\nUSAGE: CloneSeqPipeline.pl -s Strains.txt -p params.txt\n\n"}

# Setting up timestamp and runtime start.
my $totalstart=time; my $strainstart; my $functionstart;
my@time=localtime();
my$stamp= ($time[5]+1900).sprintf("%02d",$time[4]+1).sprintf("%02d",$time[3]).sprintf("%02d",$time[2]).sprintf("%02d",$time[1]);

# saving original STDOUT & STDERR
open(SAVEDSTDOUT, ">&STDOUT") || warn "Can't save STDOUT\n"; #saves STDOUT so it can later be restored
open(SAVEDSTDERR, ">&STDERR") || warn "Can't save STDERR\n"; #saves STDERR so it can later be restored

# exporting paths to perl library and bin
system 'export PERL5LIB=$PERL5LIB:/fh/fast/shou_w/src/vcftools_0.1.12b/perl';
system 'export PERL5LIB=$PERL5LIB:/fh/fast/shou_w/src/perl5';
system 'export PATH=$PATH:/fh/fast/shou_w/bin';

##########################################			
#  Declaring Variables 				#
##########################################
open PARAMS, "<$opt_p" or die $!;
my$straindir; # directory path containing the folders for each strain
my@straindirs;  #subdirectories found in #straindir
my$refdir; # directory path where reference sequences are found
my$jardir; # directory path where java jars are found
my$bindir; # directory path where scripts arefound
my$logdir; #directory path where logs are saved
my%params; # program specific parameter settings
my%overwrite; #overwrite existing files settings 
my %rmint; # remove intermediate files settings
my $paramprint=""; #parameter print messsage

##########################################			
#  Default Parameters 					#
##########################################
$params{gzip}="";
$params{fastqc}="";
$params{align}= 'mem -M -R \'@RG\tID:Sample_WS\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit\'';
$params{chrcov}="-A";
$params{gatk}="-ploidy 1 -glm BOTH -stand\_emit\_conf 10 -stand\_call\_conf 30 %%% --filter Qual=200/MinDP=10";
$params{pindel}="250";
$params{splitread}="";
$params{cnmops}="";

##########################################			
#  Read/Check Parameter file from $opt_p #
##########################################

$paramprint.="\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\nRunning the following programs:\n\n";

while(<PARAMS>){
	chomp;
	s/""/ /g; s/"//g; 	#remove quotes around fields in param file, replace empty set quotes with single space
	if(/^#/){next}elsif(/^\s*$/){next}		#skip commented and empty lines
	elsif(/^strainDir\t(\S+)$/){$straindir=$1; $paramprint.="\tstrains are in folder $straindir\n"}	#saving strain directory
	elsif(/^refDir\t(\S+)$/){$refdir=$1; $paramprint.="\treference genomes are in folder $refdir\n"}	#saving strain directory
	elsif(/^jarDir\t(\S+)$/){$jardir=$1; $paramprint.="\tjava jar files are in folder $jardir\n"}		#saving reference genome directory
	elsif(/^binDir\t(\S+)$/){$bindir=$1; $paramprint.="\texecutable files are in folder $bindir\n\n"}		#saving reference genome directory
	elsif(/^(\S+)\t([01])\t([^\t]+)\t([01]+)\t([01])$/){		# finding programs, parameters, and settings
		if($2==1){	# if set to run (1)
			if($3 ne " "){$params{$1}=$3} # overwrite default parameters if $3 wasn't empty
			$overwrite{$1}=$4; #saving overwrite setting
			$rmint{$1}=$5; # saving remove intermediate files setting
			$paramprint.= "\t$1 has been set to RUN w/ options: $params{$1}\toverwrite:$overwrite{$1}\trmint:$rmint{$1}\n";
		}
		else{
			$params{$1}="dontrun";
			$paramprint.= "\t$1 has been set to NOT RUN \n";
		}
	}
	else{$paramprint.= "\nMalformated line in $opt_p:\n$_\n\n"}; # line in parameter file didn't match any of the above formats
}
print "$paramprint\n";
close PARAMS or die $!;
die "\nERROR:\t\t Reference directory $refdir doesn't exist\n" unless -e $refdir; #stop if reference directory can't be found
die "\nERROR:\t\t BIN directory $refdir doesn't exist\n" unless -e $bindir; #stop if binary directory can't be found
die "\nERROR:\t\t JAR directory $jardir doesn't exist\n" unless -e $jardir; #stop if reference directory can't be found
die "\nERROR:\t\t Strain directory $straindir doesn't exist\n" unless -e $straindir;  #stop if directory of strains can't be found
##########################################			
#  Prepare the Reference Directory	 #
##########################################
print"\nIndexing Reference Files\n";
opendir REFS, $refdir or die "Cannot open $!";
foreach(readdir REFS){	
	if(/\.fasta$/||/\.fa$/){ #for every .fasta or .fa file do the following 4 commands (indexing & Dictionary building)
		system "$bindir"."samtools faidx $refdir$_" unless -s "$refdir$_.fai";
		system "$bindir"."bwa index $refdir$_" unless -s "$refdir$_.ann";
		system "$bindir"."java -jar $jardir"."CreateSequenceDictionary.jar R=$refdir$_ O=$refdir$_.dict" unless -s "$refdir$_.dict";
		system "mrsfast --index $refdir$_" unless -s "$refdir$_.index";
	}
}	
##########################################			
### Read/Check Strain file from $opt_s ###
##########################################

opendir STRAINDIR, $straindir or die "Cannot open $!"; #open $straindir	
foreach (readdir STRAINDIR){push @straindirs, $_} 
closedir STRAINDIR;

open STRAINS, "<$opt_s" or die $!; #read in the file specified after -s
my@strains; my%refs;
while(<STRAINS>){
	chomp;
	if(/^#/){next}elsif(/^\s*$/){next}		#skip commented and empty lines
	elsif(/^(\S+)\t(\S+)\t?\S+?$/){		#get the strain and reference
		unless(-e $straindir.$1){warn "\nERROR:\t\tCouldn't find any directory named $1 in $straindir. Not running line $_ \n";next;} #check strain exists
		unless(-s $refdir.$2.".fasta"){print "\nERROR:\t\tCouldn't find reference genome $2.fasta in $refdir. Not running line $_\n"; next;} # check reference exists
	
		unless(grep(/^$1$/,@strains)){push @strains,$1} # if not yet in strains list, add it;
		$refs{$1}=$2;  #save its reference genome
	}
	else{print "\nSkipping malformated line in $opt_s:\n$_\n\n"} # report lines that aren't in the right format
}
close STRAINS or die $!;

print "\nRunning Analysis on the following strains:\n\n";
foreach my$strain (@strains){
	print "\tStrain $strain with reference $refs{$strain}\n";
}
print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";

##########################################################			
## Analysis of strains start with this loop ###
##########################################################
my$straincount=0;
foreach my$strain (@strains){ #for each strain in the list

	$strainstart=time;  $straincount++; #for progress reporting 
	print "\n#######	Starting analysis of strain $strain with reference $refs{$strain} #######\n\n";
	my$reffile = $refdir.$refs{$strain}.".fasta";	#make path to reference genome 
	my$subdir=$straindir.$strain."\/";	#path to strain folder

	unless( -e $subdir."VS_Ref"){ system "mkdir $subdir"."VS_Ref"} #summary folder for comparison to reference
	
	unless($opt_n){ #if -n is set don't generate log and print to screen
		unless( -e $subdir."logs"){ system "mkdir $subdir"."logs"} #log folder for this strain
		open STDERR, ">>", $subdir."logs/$strain"."_$stamp.log" or die "Can't open logfile.txt: $!\n";  #all STDERR is redirected to the log file, only things printed to SAVEDSTDERR are printed to the screen
		open STDOUT, ">>", $subdir."logs/$strain"."_$stamp.log" or die "Can't open logfile.txt: $!\n";  #all STDOUT is redirected to the log file, only things printed to SAVEDSTDOUT are printed to the screen
		print "$paramprint\n";
	}	

##########################################			
## GZIP & CONCATENATION OF FASTQ FILES ###
##########################################

	opendir STRAINDIR, $subdir or die "Cannot open $!"; #open $straindir	
	
	my @readfiles;
	unless($params{gzip} eq "dontrun"){
		$functionstart=time;							#used to determine how long it takes for the function to run
		foreach(sort (readdir STRAINDIR)){				#for each file in $straindir
			if(/^\S+(_R[12])\S+?(\.fastq)(\.gz)?$/){	# check if it starts with $strain and ends in .fastq or .fastq.gz 
				if(-s "$subdir$strain$1$2.gz"){ #if it has already been concatenated and zipped
						print  SAVEDSTDOUT  "GZIP		Looks like $subdir$strain$1$2.gz already exists. Skipping gzip & concatenation of read files\n";
						next;
				}
				else{	#if it hasn't
					push @readfiles, $subdir.$_; 
				}
			}
		}
		foreach my$readfile (@readfiles){ #for each file added to @readfiles combines the fastq files for the forward (R1) and reverse (R2) reads
			if($readfile =~ /^\S+?(_R[12])\S+?(\.fastq(\.gz)?)$/){
				if(defined$3){										# if already a .gz file
					print  "GZIP		Already zipped, Executing: cat $readfile >> $subdir$strain$1$2\n";
					system("cat $readfile >>$subdir$strain$1$2");			#append / initiate .gz file
				} 	
				else{					# otherwise gzip then concatenate
					print  "GZIP		Not yet zipped, Executing: gzip $subdir$readfile then cat $readfile.gz >> $subdir$strain$1$2.gz\n";
					system("gzip $readfile"); #zip it 
					system("cat $readfile.gz >> $subdir$strain$1$2.gz"); #append / initiate .gz file
				}	
			}			
		}
		print SAVEDSTDOUT "GZIP on $strain COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n GZIP runtime: ".((time-$functionstart)/60)." min\n\n";
	}
	closedir STRAINDIR;			

### finding names of the R*.fastq.gz files
	my@strainfiles;
	opendir STRAINDIR, $subdir or die "Cannot open $!"; #open $straindir	
	foreach (readdir STRAINDIR){
		if(/^(\Q$strain\E_R[12])\.fastq\.gz$/){push @strainfiles, $1}
	}	
	closedir STRAINDIR;
	
	foreach my $file (@strainfiles){			#start of foreach @strainfiles
##########################################			
####			FASTQC				  ####
##########################################	
		unless($params{fastqc} eq "dontrun"){
			$functionstart=time;
			my $infile="$subdir$file.fastq.gz"; #specify the file to be read
		
			if((-s "$subdir$file"."_fastqc.html")&&($overwrite{fastqc} == 0)){print SAVEDSTDOUT "FASTQC		Looks like fastQC was already run for $file.fastq.gz Skipping...\n"} #check if it's already been run and isn't set to overwrite
			else{
				system "rm $subdir$file"."_fastqc.html" if -e "$subdir$file"."_fastqc.html"; #remove any preexisting files used for later success check
				print SAVEDSTDOUT "\nRunning FASTQC on $strain\n";
				my @exec =("$bindir"."fastqc", $infile); #build command: fastq
				if($overwrite{fastqc}==1){print "FASTQC		Overwriting any existing files \n"}
				print"FASTQC		Executing: @exec\n";
				system(@exec); #run command
				@exec =("cp","$subdir$file"."_fastqc.html",$subdir."VS_Ref"); #build command: cp to summary folder
				system(@exec); #run command
			}
			print SAVEDSTDERR "\nERROR		FASTQC on $file FAILED!!!!\n\n" unless -s "$subdir$file"."_fastqc.html";
			print SAVEDSTDOUT "FASTQC on $file COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n FASTQC runtime: ".((time-$functionstart)/60)." min\n\n";
		}
		
##########################################			
## bwa | samtools view | samtools sort  ##
##########################################		
		if($file eq $strainfiles[-1]){ #if $file is the last file in $strainfiles
			$file =~ s/_R[12]$//; # drop R1/R2 from end of file name
		
			unless($params{align} eq "dontrun"){
				$functionstart=time;
				my $infile="$subdir$file"."_R1.fastq.gz"; #file to be read
				if(-s "$subdir$file"."_R2.fastq.gz"){$infile.=" $subdir$file"."_R2.fastq.gz"}			#if paired end run include read2 file
				my $outfile="$subdir$file.sorted"; #file to be generated
				my $params = $params{align}; # get parameters for ALIGN
		
				if((-s "$outfile.mdup.bam")&&($overwrite{align} == 0)){print SAVEDSTDOUT "ALIGN		Looks $outfile.mdup.bam already exists. Skipping....\n"}	
				else{
					system "rm $outfile.mdup.bam" if -e "$outfile.mdup.bam"; #remove any preexisting files used for later success check
					print SAVEDSTDOUT "\nRunning ALIGN on $strain\n";
					
					# Build command bwa | samtools view | samtools sort
					my @exec= ("$bindir"."bwa", $params{align}, $reffile, $infile, "|");
					push @exec, ("$bindir"."samtools view", "-Su","-","|");
					push @exec, ("$bindir"."samtools sort", "-", $outfile);
		
					if($overwrite{align}==1){print "ALIGN		Overwriting any existing files\n"}
					print "\nALIGN		Executing: @exec \n\n";	
					system ("@exec"); #run bwa | samtools view | samtools sort
				
					# Build picard markduplicates
					@exec = ("$bindir"."java -jar $jardir"."MarkDuplicates.jar","I=$outfile.bam","O=$outfile.mdup.bam","M=$outfile.mdup.metric");
				
					print "\nALIGN		Executing: @exec \n\n";
					system ("@exec"); #run picard markduplicates
				
					#samtools index
					@exec = ("$bindir"."samtools index", "$outfile.mdup.bam");
				
					print "\nALIGN		Executing: @exec \n\n";
					system ("@exec"); # run samtools index
				
					if($rmint{align}==1){
						print "\nALIGN		Removing intermediate file: $outfile.bam\n\n";
						system "rm $outfile.bam";
					}
				}
				print SAVEDSTDOUT "\nERROR		ALIGN on $strain FAILED!!!!\n\n" unless -s $outfile.".mdup.bam"; #success check
				print SAVEDSTDOUT "ALIGN on $strain COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n ALIGN runtime: ".((time-$functionstart)/60)." min\n\n";		
			} #end of ALIGN

####################################################		
######	ChrCov:	mpileup | pile2cov |chrcov.R #######
####################################################	
			unless($params{chrcov} eq "dontrun"){
				$functionstart=time;
			
				my $infile="$subdir$file.sorted.mdup.bam";
				my $outfile="$subdir$file.cov";
				my $params = $params{chrcov}."f$reffile";	

				if((-s "$subdir$strain"."_ChrCov.png")&&($overwrite{chrcov} == 0)){print SAVEDSTDOUT "CHRCOV		Looks $subdir$strain"."_ChrCov.png already exists. Skipping...\n"}
				else{
					system "rm $subdir$strain"."_ChrCov.png" if -e "$subdir$strain"."_ChrCov.png"; #remove any preexisting files used for later success check
					print SAVEDSTDOUT "\nRunning CHRCOV on $strain\n";
					if($overwrite{chrcov}==1){print "CHRCOV  	Overwriting any existing files \n"}
					open COV, ">$outfile"; #open .cov file to save to
					print "CHRCOV		$bindir"."samtools mpileup $params $infile|\n\n";
					open(PILEUP, "$bindir"."samtools mpileup $params $infile|"); #read the output from mpileup into PILEUP file handle
					while(<PILEUP>){
							chomp;
							/(^\S+)\t(\d+)\t\w\t(\d+)/; # find the relevant data
							print COV "$1~$2\t$1\t$2\t$3\n";	#print them to the .cov file
					}
					close COV;
				
					$infile=$outfile; #now the .cov file is the input
					my @exec = ("Rscript --vanilla $bindir"."ChrCov.R",$strain,$infile,$refs{$strain}); #build command to send .cov file to ChrCov.R
					print"\nCHRCOV 		Executing: @exec \n\n";
					system("@exec");		
				
					my$exec= "cp $subdir$strain"."_ChrCov.png $subdir"."VS_Ref"; #copy result to summary directory
					print"\nCHRCOV 		Executing: $exec \n\n";
					system($exec);
					$exec= "cp $subdir$strain"."_ChrCov.tab $subdir"."VS_Ref"; #copy result to summary directory
					print"\nCHRCOV 		Executing: $exec \n\n";
					system($exec);
					
					if($rmint{chrcov}==1){system "rm $subdir$strain.cov"}# remove .cov to save space (it's the largest file generated)
				}
				print SAVEDSTDOUT "\nERROR		CHRCOV on $strain FAILED!!!!\n\n" unless -s "$subdir$strain"."_ChrCov.png"; #success check
				print SAVEDSTDOUT "CHRCOV on $strain COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n CHRCOV runtime: ".((time-$functionstart)/60)." min\n\n";
				
			}#samtools mpileup | cov | ChrCov	

#######################################			
######	GATK SNP Calling	      #####
#######################################			
			unless($params{gatk} eq "dontrun"){
				$functionstart=time;
				if((-s "$subdir$file.vcf")&&($overwrite{gatk} == 0)){print SAVEDSTDOUT "GATK 		Looks $subdir$file.vcf already exists. Skipping...\n"}
				else{
					system "rm $subdir$file.vcf" if -e "$subdir$file.vcf"; #remove any preexisting files used for later success check
					print SAVEDSTDOUT "\nRunning GATK on $strain\n";				
					if($overwrite{gatk}==1){print "GATK		Overwriting any existing files \n"}
				
					#GATK RealignerTargetCreator				
					my $infile="$subdir$file.sorted.mdup.bam";
					my $outfile="$subdir$file.list";
					my @params = split /%%%/, $params{gatk}; # multiple GATK programs have parameter settings, they are separated by %%% in the parameter file
				
					my @exec = ("$bindir"."java -jar $jardir"."GenomeAnalysisTK.jar", "-T RealignerTargetCreator","-R $reffile","-I $infile", "-o $outfile");
				
					print"\nGATK	 	Executing: @exec \n\n";
					system ("@exec");
				
					#GATK IndelRealigner
					$infile="$subdir$file.sorted.mdup.bam";
					$outfile="$subdir$file.realigned.bam";
				
					@exec = ("$bindir"."java -jar $jardir"."GenomeAnalysisTK.jar", "-T IndelRealigner","-R $reffile","-targetIntervals $subdir$file.list","-I $infile", "-o $outfile");
				
					print"\nGATK 		Executing: @exec \n\n";
					system ("@exec");
				
					#GATK Unified Genotyper Round 1
					$infile=$outfile;
					$outfile="$subdir$file.raw_vcf";
					my $params = $params[0];
				
					@exec = ("$bindir"."java -jar $jardir"."GenomeAnalysisTK.jar", "-T UnifiedGenotyper","-R $reffile","$params","-I $infile", "-o $outfile");
				
					print"\nGATK 		Executing: @exec \n\n";
					system ("@exec");
				
					#vcftools vcf-annotate
					$infile=$outfile;
					$outfile="$subdir$file.filtered.raw_vcf";
					$params = $params[1];
				
					@exec = ("$bindir"."vcf-annotate",$params,$infile, ">".$outfile);
				
					print"\nGATK 		Executing: @exec \n\n";
					system("@exec");
				
					#BaseRecalibrator
					$infile="$subdir$file.realigned.bam";
					$outfile="$subdir$file.recal.table";
				
					@exec = ("$bindir"."java -jar $jardir"."GenomeAnalysisTK.jar", "-T BaseRecalibrator","-R $reffile","-knownSites $subdir$file.filtered.raw_vcf","-I $infile", "-o $outfile");
				
					print"\nGATK 		Executing: @exec \n\n";
					system("@exec");
				
					#PrintReads
					$infile="$subdir$file.realigned.bam";
					$outfile="$subdir$file.recal.bam";
				
					@exec = ("$bindir"."java -jar $jardir"."GenomeAnalysisTK.jar", "-T PrintReads","-R $reffile","-BQSR $subdir$file.recal.table","-I $infile", "-o $outfile");
				
					print"\nGATK 		Executing: @exec \n\n";
					system("@exec");
				
					#GATK Unified Genotyper Round 1#GATK Unified Genotyper Round 1
					$infile=$outfile;
					$outfile="$subdir$file.vcf";
					$params = $params[0];
				
					@exec = ("$bindir"."java -jar $jardir"."GenomeAnalysisTK.jar", "-T UnifiedGenotyper","-R $reffile","$params","-I $infile", "-o $outfile");
				
					print"\nGATK 		Executing: @exec \n\n";
					system("@exec");
				
					my $exec =("cp $subdir$strain.vcf $subdir"."VS_Ref"); #copy to summary folder
					system $exec;
					$exec =("cp $subdir$strain.vcf.idx $subdir"."VS_Ref"); #copy to summary folder
					system $exec;
				
					if($rmint{gatk}==1){ #removing intermediate files
						print "GATK			Removing intermediate file $subdir$file.list\n";
						system("rm $subdir$file.list");
						print "GATK			Removing intermediate file $subdir$file.realigned.bam\n";
						system("rm $subdir$file.realigned.bam");
						print "GATK			Removing intermediate file $subdir$file.raw_vcf\n";
						system("rm $subdir$file.raw_vcf");
						print "GATK			Removing intermediate file $subdir$file.recal.table\n";
						system("rm $subdir$file.recal.table");
						print "GATK			Removing intermediate file $subdir$file.recal.bam\n";
						system("rm $subdir$file.recal.bam");
					}			
				}
				print SAVEDSTDOUT "\nERROR		GATK on $strain FAILED!!!!\n\n" unless -s "$subdir$file.vcf"; #success check
				print SAVEDSTDOUT "GATK on $strain COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n GATK runtime: ".((time-$functionstart)/60)." min\n\n";
			}#end of GATK

######################################			
######	Pindel 			 		 #####
######################################			
			unless($params{pindel} eq "dontrun"){
				$functionstart=time;
				my $infile="$subdir$file.sorted.mdup.bam";
				my $outfile="$subdir$file.pindel";
				my $params = $params{pindel};
				
				my @exec = ("$bindir"."pindel", "-f $reffile","-i $subdir"."temp.pindel", "-o $outfile");
			
				opendir STRAINDIR, $subdir;
				if((-e $subdir.$file.".pindel_BP")&&($overwrite{pindel} == 0)){print SAVEDSTDOUT "PINDEL 		Looks $outfile already exists. Skipping...\n"}
				else{
					system "rm $subdir.$file".".pindel_D" if -e $subdir.$file.".pindel_D";  #remove any preexisting files used for later success check
					print SAVEDSTDOUT "\nRunning PINDEL on $strain\n";					
					if($overwrite{pindel}==1){print "PINDEL  	Overwriting any existing files \n"}
					
					open TEMP, ">$subdir"."temp.pindel"; #make a temporary pindel file
					print TEMP "$infile $params sample1"; #print the required line to it 
					close TEMP;
					
					print"\nPINDEL 		Executing: @exec \n\n";
					system("@exec");
					if($rmint{pindel}==1){
						print "PINDEL		Removing intermediate file $subdir"."temp.pindel\n";
						system" rm $subdir"."temp.pindel";
					}
				}
				close STRAINDIR;
			
				my $exec =("cp $subdir$strain*pindel* $subdir"."VS_Ref"); #copy results to summary folder
				system $exec;
				print SAVEDSTDOUT "\nERROR		PINDEL on $strain FAILED!!!! (tiny possibility that it worked but detected no deletions)\n\n" unless -s $subdir.$file.".pindel_D"; #success check
				print SAVEDSTDOUT "PINDEL on $strain COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n PINDEL runtime: ".((time-$functionstart)/60)." min\n\n";
			}#end of pindel
	
#######################################			
######	Splitread			 	  #####
#######################################			
			#file test to make sure its paired end
		
			unless($params{splitread} eq "dontrun"){
				$functionstart=time;
				if((-e $subdir.$file.".splitread.pair.discordand.txt")&&($overwrite{splitread} == 0)){print SAVEDSTDOUT "SPLITREAD 	Looks $subdir$file.splits already exists. Skipping...\n"}
				else{
					system "rm $subdir.$file".".splitread.pair.discordand.txt" if -e $subdir.$file.".splitread.pair.discordand.txt"; #remove any preexisting files used for later success check
					print SAVEDSTDOUT "\nRunning SPLITREAD on $strain\n";	
					if($overwrite{splitread}==1){print "SPLITREAD  	Overwriting any existing files \n"}
					
					my$params; my $trim;
					
					if($params{splitread} eq ""){ #calculate most common read length and use as trim length
						print "\nSPLITREAD 	Calculating mode read length\n";	
						my $infile="$subdir$file.sorted.mdup.bam";
						my $exec = ("$bindir"."java -jar $jardir"."GenomeAnalysisTK.jar -T ReadLengthDistribution -R $reffile -I $infile");
						print "SPLITREAD Executing: $exec\n\n";
						my @test=`$exec`;
						my %readlength;
						foreach my$a (@test){
							if($a =~ /\s+(\d+)\s+(\d+)/){$readlength{$2}=$1} 
						}		
						
						$trim=$readlength{max(keys(%readlength))}; #set params (trim length) to most frequent readlength
						if(($trim/2-int($trim/2))!= 0){$trim = $trim-1} # if value is odd subtract 1
						
						my $mindist=$trim+5;
						my $maxdist=$trim*5;
						
						my $hamming=int($trim/20); #set hamming distance to 5% of read length
						if(($hamming/2-int($hamming/2))!= 0){$hamming = $hamming+1} # if hamming distance is odd add 1						
						
						$params = "-len $trim -ed $hamming -ed2 ".($hamming/2)." -u $maxdist -l $mindist -m 10000000 ";	

						print "\nSPLITREAD 	Setting splitread parameters to $params\n";
					}
					else{  #or use pre-specified parameter from params file
						$params = $params{splitread};
						$params{splitread} =~ /-len(\d+)/; #extracting the trim length from params
						$trim = $1;
					} 
					
					my $infile = $subdir.$strain;
					
				
					# Unzipping the fastq files
					print "SPLITREAD 	Unzipping fastq files\n\n";			
					system "gunzip -k -f $infile"."_R1.fastq.gz";
					system "gunzip -k -f $infile"."_R2.fastq.gz";
				
					# mixing the read1 & read 2 fastq files and splitting it into smaller files of 200,000 lines
					my @exec = ("$bindir"."mixcatnewHeader", "-f $infile"."_R1.fastq");
					push @exec, "-r $infile"."_R2.fastq";
					push @exec , "| split -l 200000 - $strain/$strain.mix.split.";
					print"\nSPLITREAD 	Executing: @exec \n\n";
					system("@exec");
				
					# Rezipping the fastq files
					print "SPLITREAD 	Rezipping fastq files\n\n";
					system "rm $infile"."_R1.fastq";
					system "rm $infile"."_R2.fastq";

					print"\nSPLITREAD	Trimming split files\n\n";
					opendir STRAINDIR, "$subdir";
					foreach(readdir STRAINDIR){
						if( /\.mix\.split\.\w\w$/ ){
						
							#Running Readtrimmer
							@exec =("$bindir"."Readtrimmer -i $subdir$_", "-o $subdir$_.trimmed -l $trim");
							print"\nSPLITREAD 	Executing: @exec \n\n";
							system "@exec";
							system "rm $subdir$_";
						
							#Running Split read
							@exec = ("$bindir"."SplitReadAll_lite");
							push @exec, ("-ref $reffile","-refai $reffile.fai");
							push @exec, ("$params", "-seq $subdir$_.trimmed", "-out $subdir$_.trimmed.map","-outunmap $subdir$_.trimmed.unmapped", "-outdir $subdir$_.trimmed.match");
							push @exec, ("-path $bindir", "-mrsfast $bindir", "-samtools $bindir");
		
							print"\nSPLITREAD 		Executing: @exec \n\n";
							system("@exec");
						}
					}
					close STRAINDIR;
				
					opendir STRAINDIR, "$subdir";

					#merging BAM files and removing individual files
					print "\nSPLITREAD	Merging BAM Files\n";
					print "\nSPLITREAD  Executing: $bindir"."samtools merge -f $subdir$strain.splitread.single.sorted.bam  $subdir$strain.mix.split.*.trimmed.match.single.sorted.bam\n\n";
					
					system "$bindir"."samtools merge -f $subdir$strain.splitread.single.sorted.bam  $subdir$strain.mix.split.*.trimmed.match.single.sorted.bam";
					system "rm $subdir$strain.mix.split.*.trimmed.match.single.sorted.bam";
				
					system "$bindir"."samtools merge -f $subdir$strain.splitread.sorted.bam  $subdir$strain.mix.split.*.trimmed.match.sorted.bam";
					system "rm $subdir$strain.mix.split.*.trimmed.match.sorted.bam";
				
					system "$bindir"."samtools merge -f $subdir$strain.splitread.split.match.single.sorted.bam  $subdir$strain.mix.split.*.trimmed.match.split.match.single.sorted.bam";
					system "rm $subdir$strain.mix.split.*.trimmed.match.split.match.single.sorted.bam";
				
					system "$bindir"."samtools merge -f $subdir$strain.splitread.split.match.sorted.bam  $subdir$strain.mix.split.*.trimmed.match.split.match.sorted.bam";
					system "rm $subdir$strain.mix.split.*.trimmed.match.split.match.sorted.bam";
				
					system "$bindir"."samtools merge -f $subdir$strain.splitread.transchr.sorted.bam  $subdir$strain.mix.split.*.trimmed.match.transchr.sorted.bam";
					system "rm $subdir$strain.mix.split.*.trimmed.match.transchr.sorted.bam";
				
					system "$bindir"."samtools merge -f $subdir$strain.splitread.inv.evert.sorted.bam  $subdir$strain.mix.split.*.trimmed.match.inv.evert.sorted.bam";
					system "rm $subdir$strain.mix.split.*.trimmed.match.inv.evert.sorted.bam";
					
					system "rm $subdir$strain.mix.split.*.trimmed.unmapped";
				
					# merging text files
					foreach(readdir STRAINDIR){
						if(/\.bam$/){next}
						elsif(/^\Q$strain\E\.mix\.split\.\w\w\.trimmed\.match(\.\S+)$/){
							system "cat $subdir$_ >> $subdir$strain.splitread$1";
							system "rm $subdir$_";
						}
						elsif(/^\Q$strain\E\.mix\.split\.\w\w\.trimmed(\.unmapped)?/){system "rm $subdir$_"}
					}
					close STRAINDIR;
					print "cp $subdir$strain*.txt $subdir"."VS_Ref"; #copy results to summary folder
					
					#Robin's Python Code				
					#python /fh/fast/shou_w/bin/Robin/user_scripts/splitread_pipeline/bin_total.py sample.match.pair.discordant.txt > sample.bins
					#python /fh/fast/shou_w/bin/Robin/user_scripts/splitread_pipeline/get_raw_input_from_targets_Supercontig.py sample.bins sample.match.pair.discordant.txt >sample.binned.events
					#python /fh/fast/shou_w/bin/Robin/user_scripts/splitread_pipeline/get_unique_call_reads_Supercontig.py sample.binned.events > sample.unique.bins
				}
				print SAVEDSTDOUT "\nERROR		SPLITREAD on $strain FAILED!!!!\n\n" unless -s $subdir.$file.".splitread.pair.discordand.txt"; #success check
				print SAVEDSTDOUT "SPLITREAD on $strain COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n SPLITREAD runtime: ".((time-$functionstart)/60)." min\n\n";	
			}#end of splitread 	

#######################################			
######	cn.mops			 		  #####
#######################################		
			unless($params{cnmops} eq "dontrun"){
				$functionstart=time;
				if((-e $subdir.$strain."_CNVs.gff")&&($overwrite{cnmops} == 0)){print SAVEDSTDOUT "CNMOPS 	Looks $subdir$file.gff already exists. Skipping...\n"}
				else{
					system "rm $subdir$strain"."_CNVs.gff" if -e $subdir.$strain."_CNVs.gff"; #remove any preexisting files used for later success check
					print SAVEDSTDOUT "\nRunning CNMOPS on $strain\n";				
					my $infile="$subdir$file.sorted.mdup.bam";
					my $params=$params{cnmops}; #should be a list paths to *.sorted.mdup.bam files for each strains we want to compare
					open BAMS, "<$params" or print SAVEDSTDOUT "\nERROR	Couldn't find list of mdup.bam files $params used for cnmops comparison!!!!\n\n";
					
					my @links;
					while (<BAMS>){ # read list of bam files and make a symbolic link in the current folder
						chomp;
						if(/^#/){next}elsif(/^\s*$/){next} #skip empty & commented lines
						unless(/\/$strain\./){ #unless its the current strain being analyzed make a symbolic link to the other .bam files and .bai files
							unless(-s $_){print SAVEDSTDOUT "\nERROR	Couldn't find $_ for cnmops comparison!!!!\n\n";next}
							
							my $name=substr($_,0,length($_)-1)."*";
							print "cp -s $name $subdir\n";
							system "cp -s $name $subdir";
						}
						
						if(/[^\/]+bam$/){push @links, $&;push @links, $&.".bai"; } #save file name (not the path) to @links
					}
					close BAMS;
				
					my @exec = ("Rscript --vanilla $bindir"."cnmops.R",$strain,$refs{$strain}); #build cnmops command
								
					if($overwrite{cnmops}==1){print "CNMOPS  	Overwriting any existing files \n"}
					print"\nCNMOPS 		Executing: @exec \n\n";
					system("@exec");
					system "cp $subdir*CNV* $subdir"."VS_Ref"; #copy results to summary folder
					if($rmint{cnmops}==1){foreach my$link (@links){system "rm $subdir$link" unless $link=~ /$strain/}} #remove the links to the other bam files make earlier
				}
				print SAVEDSTDOUT "\nERROR		CNMOPS on $strain FAILED!!!!\n\n" unless -s $subdir.$strain."_CNVs.gff"; #success check
				print SAVEDSTDOUT "CNMOPS on $strain COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n CNMOPS runtime: ".((time-$functionstart)/60)." min\n\n";
			}#end of cn.mops
			
#######################################			
######	End of tools		  #####
######################################	
			if($rmint{gzip}==1){system "rm $subdir$strain"."_R?.fastq.gz"} #remove combined readfile archive to save space							
		}#end of if last strainfile
	}#end of foreach @strainfiles	
	open(STDOUT, ">&SAVEDSTDOUT") || warn "Can't restore STDOUT\n";
	open(STDERR, ">&SAVEDSTDERR") || warn "Can't restore STDERR\n";
	print "\n$strain COMPLETE!! ($straincount of ".(scalar@strains).")\n  Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n";
}#foreach @strains
#system 'chmod -Rf 2770 /fh/fast/shou_w/NextGenSeq/CloneSeq/';
print "\nRUN COMPLETE!!\n  Total runtime: ".((time-$totalstart)/60)." min\n\n";
__END__
