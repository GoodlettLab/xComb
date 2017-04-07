#!/usr/bin/perl-w

#############################################################################
# Copyright (C) 2017, William Fondrie.                                      #
#                                                                           #
# This script is part of xComb.		                                          #
#                                                                           #
# xComb is free software. You can redistribute and/or                       #
# modify it under the terms of the GNU General Public License               #
# as published by the Free Software Foundation; either version 3            #
# of the License or, (at your option), any later version.                   #
#                                                                           #
# These modules are distributed in the hope that they will be useful,       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              #
# GNU General Public License for more details.                              #
#############################################################################
#############################################################################
# Copyright (C) 2010, Alexandre Panchaud.                                   #
#                                                                           #
# This script is free software. You can redistribute and/or                 #
# modify it under the terms of the GNU General Public License               #
# as published by the Free Software Foundation; either version 2            #
# of the License or, (at your option), any later version.                   #
#                                                                           #
# These modules are distributed in the hope that they will be useful,       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              #
# GNU General Public License for more details.                              #
#############################################################################

use strict;
#use warnings;
#use File::Copy;

#################################################################################
#################################################################################
######-->	Main section						<--######
#################################################################################
#################################################################################

### declare arguments for this program, check arguments and store argument values
unless(@ARGV==7){
	print "\nWrong or missing arguments!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
# Inter/intra type
my $type = $ARGV[0];
unless($type =~ m/-i=[0-2]{1}/){
	print "\nWrong or missing arguments!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
# Cross-linker type
my $cx = $ARGV[1];
unless($cx =~ m/-x=[0-6]{1}/){
	print "\nWrong or missing arguments!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
# Amine cross-linker type
my $mc = $ARGV[2];
unless($mc =~ m/-mc=[0-2]{1}/){
	print "\nWrong or missing arguments!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
# Min peptide length
my $pep_length = $ARGV[3];
unless($pep_length =~ m/-pl=[1-6]{1}/){
	print "\nWrong or missing arguments!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
# Phenyx format
my $format = $ARGV[4];
unless($format =~ m/-f=[0-1]{1}/){
	print "\nWrong or missing arguments!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
# Test mode
my $test = $ARGV[5];
unless($test =~ m/-t=[0-1]{1}/){
	print "\nWrong or missing arguments!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
# directory containing file(s).digest
my $directory = $ARGV[6];
#$directory =~ s|/||;
### specify some variables, chdir, count *.digest files
chdir("$directory");
my (@files,$filenb,$sep);
@files = <*.digest>;
$filenb = @files;

### check if inter specified and only 1 protein
if($filenb == 1 && $type eq "-i=1"){
	print "\n WARNING: You cannot perform inter cross-linking (-i=1) with 1 protein sequence only!!! See instruction:\n\n";
	&instruction;
	exit 1;
}
### warning for inter&intra with only 1 protein
if($filenb == 1 && $type eq "-i=2"){
	print "\n WARNING: Only intra cross-links will be performed as you have only 1 protein sequence and specified -i=2!!!\n\n";
}
### check if no more than 50 proteins
if($filenb > 50){
	print "\n WARNING: The maximum protein sequences allowed is 50!!! See instruction:\n\n";
	&instruction;
	exit 1;
}
### check input test mode and set separator accordingly
if($test eq "-t=1"){
	$sep = "|";
}else{
	$sep = "";
}

### check Missed cleavage mode and set value accordingly
if($mc eq "-mc=2"){
	$mc = 2;
}elsif($mc eq "-mc=1"){
	$mc = 1;
}else{
	$mc = 0;
}
### retrieve Min peptide length value
$pep_length =~ s/-pl=//;
### create output file
my $outfile = $directory."_cxdb.fasta";
print "\n--> Output filename: $outfile\n\n";
### IMPORTANT !!! Use solution 1 for command line or change to solution 2 for CGI launch
open(OUT,">../$outfile"); # solution 1
#open(OUT,">$outfile"); # solution 2

### increment to open one file less at each loop for inter-cx and start with file 2
my $count = 1;
### print some STDOUT for user
print "\n--> Creating CXMS database based on following parameters:\n\n";
print "\t--> DIGEST directory = $directory\n";
print "\t--> Inter/intra mode = $type\n";
print "\t--> Type of crosslinker = $cx\n";
print "\t--> Missed cleavage = $mc\n";
print "\t--> Phenyx format = $format\n";
print "\t--> Test mode = $test\n\n";
print "--> Converting...\n\n";

#################################################################################
#################################################################################
######-->		DoNotCleave mode				<--######
#################################################################################
#################################################################################

for(my $h=0; $h<@files; ++$h){
	print "\tNow crosslinking file: $files[$h]\n";
	### declare all variables
	my (@header,@line,@pos,@seq,@shift,@shift_pos);
	my ($intra_last,$intra_ab,$intra_ba);
	### Open reference file
	open(IN, $files[$h]);
	### read line content and store into array
	while(<IN>){
		### match AC, DE and peptide lines / avoid table header line
		next if /^mass.*$/;
		if($_ =~ m/AC=.*/){
			chomp $_;
			$_ =~ s/AC=//;
			push @header, $_;
		}
		elsif($_ =~ m/DE=.*/){
			chomp $_;
			$_ =~ s/DE=//;
			push @header, $_;
		}
		elsif($_ =~ m/ID=.*/){
			chomp $_;
			$_ =~ s/ID=//;
			push @header, $_;
		}else{
		chomp $_;
		@line = split(/\t/, $_);
		push @pos, $line[1];
		push @seq, $line[3];
		}
	}
	close IN;
	### test print
	#print "#########@header\n@pos\n@seq\n";

	if($type eq "-i=0" || $type eq "-i=2"){
		#################################################################################
		######-->	first perform intra cross-linked peptides		<--######
		#################################################################################

		### create new array for shifting with loops
		@shift = @seq;
		@shift_pos = @pos;
		for(my $j=0; $j<@seq; ++$j){
			### test print
			#print "-->$seq[$j]\n";
			### reinitialize both variables for each loop
			$intra_ab = "";
			$intra_ba = "";
			### calculates the size of the array for each loop after shift of previous loop
			$intra_last = @shift-1;
			### print the sequences like: aa -> ab & ba -> etc
			for(my $i=0; $i<=$intra_last; ++$i){
				### test both peptide length according to input
				my $length1 = length($seq[$j]);
				my $length2 = length($shift[$i]);
				if($length1 < $pep_length || $length2 < $pep_length){
					next;
				}
				$intra_ab = $seq[$j].$sep.$shift[$i];
				$intra_ba = $shift[$i].$sep.$seq[$j];
				### test for whether to print or not based on crosslinker input
				my $print_intra = &crosslinker(\$seq[$j],\$shift[$i],\$cx,\$mc);

				### print FASTA output based on print value from SUB
				if($print_intra == 1){
					### print sequences into FASTA format
					if($format eq "-f=1"){
						printfasta_phenyx_intra(\$header[0],\$header[1],\$header[2],\$j,\$i,\$seq[$j],\$shift[$i],\$pos[$j],\$shift_pos[$i],\$intra_ab,\$intra_ba);
					}elsif($format eq "-f=0"){
						printfasta_normal_intra(\$header[0],\$header[1],\$header[2],\$j,\$i,\$seq[$j],\$shift[$i],\$pos[$j],\$shift_pos[$i],\$intra_ab,\$intra_ba);
					}
				}
			}
			### remove first item for array for next loop to perform only: bb -> bc & cb -> etc
			shift(@shift);
			shift(@shift_pos);
		}
	}

	#################################################################################
	######-->	perform inter cross-linked peptides			<--######
	#################################################################################

	### no inter-crosslink peptides created if last digest file or only a single digest file
	next if($h==@files-1);
	next if($filenb==1);

	if($type eq "-i=0" || $type eq "-i=1"){
		######-->	Build inter cross-linking sequences			<--######
		for(my $i=$count; $i<@files; ++$i){
			print "\t\twith: $files[$i]\n";
			my (@header1,@line1,@pos1,@seq1);
			my ($inter_aA,$inter_Aa);
			open(INTER, $files[$i]);
			### read line content and store into array
			while(<INTER>){
				### match AC, DE and peptide lines / avoid table header line
				next if /^mass.*$/;
				if($_ =~ m/AC=.*/){
					chomp $_;
					$_ =~ s/AC=//;
					push @header1, $_;
				}
				elsif($_ =~ m/DE=.*/){
					chomp $_;
					$_ =~ s/DE=//;
					push @header1, $_;
				}
				elsif($_ =~ m/ID=.*/){
					chomp $_;
					$_ =~ s/ID=//;
					push @header1, $_;
				}else{
				chomp $_;
				@line1 = split(/\t/, $_);
				push @pos1, $line1[1];
				push @seq1, $line1[3];
				}
			}
			close INTER;
			### test print
			#print "!!!@header\n@pos\n@seq\n";
			#print "!!!!!!@header1\n@pos1\n@seq1\n";
			for(my $j=0; $j<@seq; ++$j){
				### reinitialize both variables for each loop
				$inter_aA = "";
				$inter_Aa = "";
				### Build swap 1
				for(my $i=0; $i<@seq1; ++$i){
					### test both peptide length according to input
					my $length1 = length($seq[$j]);
					my $length2 = length($seq1[$i]);
					if($length1 < $pep_length || $length2 < $pep_length){
						next;
					}
					### assemble both sequences
					$inter_aA = $seq[$j].$sep.$seq1[$i];
					$inter_Aa = $seq1[$i].$sep.$seq[$j];
					### test for whether to print or not based on crosslinker input
					my $print_inter = &crosslinker(\$seq[$j],\$seq1[$i],\$cx,\$mc);
					### print FASTA output based on print value from SUB
					if($print_inter == 1){
						### print sequences into FASTA format
						if($format eq "-f=1"){
							printfasta_phenyx_inter(\$header[0],\$header[1],\$header1[0],\$header1[1],\$j,\$i,\$seq[$j],\$seq1[$i],\$pos[$j],\$pos1[$i],\$inter_aA,\$inter_Aa);
						}elsif($format eq "-f=0"){
							printfasta_normal_inter(\$header[0],\$header[1],\$header1[0],\$header1[1],\$j,\$i,\$seq[$j],\$seq1[$i],\$pos[$j],\$pos1[$i],\$inter_aA,\$inter_Aa);
						}
					}
				}
			}
		}
		### increment count for next loop to open one file less for inter-CX
		++$count;
	}
}

### print STDOUT for user
print "\n--> Finished :o\n\n";
### exit
exit;

#################################################################################
#################################################################################
######-->	Subroutine section					<--######
#################################################################################
#################################################################################

sub instruction {
print '
digest2cxdb.pl (xComb v1.4)

FUNCTION:
Creates inter & intra cross-linking databases for search with standard algorithm
such as Sequest, Mascot or Phenyx. Uses pre-digested protein sequences in DIGEST
format. If a single digest is present, performs only intra database. The maximum
protein sequences allowed per directory in digest format is 50!!!

USAGE:
digest2cxdb.pl -i=[value] -x=[value] -mc=[value] -pl=[value] -f=[value] -t=[value] [digest_dir]

	-i=[value]: Inter/Intra type mode

	  0= Both intra and inter crosslinks are generated
	  1= Only inter crosslinks are generated
	  2= Only intra crosslinks are generated

	-x=[value]: Specify type of cross-linker

	  0= Amine/Amine (e.g. BS2, DST, DFDNB)
	  1= Sulfhydryl/Sulfhydryl (e.g. BMOE, BMDB or S-S bridge)
	  2= Amine/carboxyl (e.g. EDC)
	  3= Amine/Non-selective (e.g. SDA)
	  4= Amine/Sulfhydryl (e.g. AMAS, SPDP)
	  5= Tyrosine/Tyrosine (e.g. dityrosine)
	  6= Lysine/Lysine

	-mc=[value]: Specify whether amine cross-linking needs 1 trypsin missed cleavage at least. on=1 / off=0

	-pl=[value]: Specify the minimum peptide length for building a permutation [1-6]

	-f=[value]: Phenyx FASTA formatted output. on=1 / off=0

	-t=[value]: Test mode which adds "|" between the two peptides for easier reading. on=1 / off=0

	[digest_dir]: Specifies directory containing files in DIGEST format

Ex: digest2cxdb.pl -i=1 -x=1 -mc=1 -pl=4 -f=0 -t=0 digest_dir


';
}

sub printfasta_phenyx_inter {
	my ($header1a,$header1b,$header2a,$header2b,$j,$i,$seq1,$seq2,$pos1,$pos2,$inter_aA,$inter_Aa) = @_;
	print OUT ">".$$header1a."_".$$header2a."_aA_".$$j."_".$$i." \\DE=a=".$$seq1." (".$$pos1.") ".$$header1b." cx A=".$$seq2." (".$$pos2.") ".$$header2b."\n";
	print OUT $$inter_aA."\n";
	print OUT ">".$$header2a."_".$$header1a."_aA_".$$i."_".$$j." \\DE=a=".$$seq2." (".$$pos2.") ".$$header2b." cx A=".$$seq1." (".$$pos1.") ".$$header1b."\n";
	print OUT $$inter_Aa."\n";

}

sub printfasta_normal_inter {
	my ($header1a,$header1b,$header2a,$header2b,$j,$i,$seq1,$seq2,$pos1,$pos2,$inter_aA,$inter_Aa) = @_;
	print OUT ">".$$header1a."_".$$header2a."_aA_".$$j."_".$$i." a=".$$seq1." (".$$pos1.") ".$$header1b." cx A=".$$seq2." (".$$pos2.") ".$$header2b."\n";
	print OUT $$inter_aA."\n";
	print OUT ">".$$header2a."_".$$header1a."_aA_".$$i."_".$$j." a=".$$seq2." (".$$pos2.") ".$$header2b." cx A=".$$seq1." (".$$pos1.") ".$$header1b."\n";
	print OUT $$inter_Aa."\n";

}

sub printfasta_phenyx_intra {
	my ($header1a,$header1b,$header1c,$j,$i,$seq1,$seq2,$pos1,$pos2,$intra_ab,$intra_ba) = @_;
	### test if aa or ab peptide. Write only once for aa and ab/ba for others
	if($$intra_ab eq $$intra_ba){
		print OUT ">".$$header1a."_aa_".$$j."_".$$i." \\AC=".$$header1a." \\ID=".$$header1c." \\DE=a=".$$seq1." (".$$pos1.") cx a=".$$seq2." (".$$pos2.") ".$$header1b."\n";
		print OUT $$intra_ab."\n";
	}else{
		print OUT ">".$$header1a."_ab_".$$j."_".$$i." \\AC=".$$header1a." \\ID=".$$header1c." \\DE=a=".$$seq1." (".$$pos1.") cx b=".$$seq2." (".$$pos2.") ".$$header1b."\n";
		print OUT $$intra_ab."\n";
		print OUT ">".$$header1a."_ba_".$$i."_".$$j." \\AC=".$$header1a." \\ID=".$$header1c." \\DE=b=".$$seq2." (".$$pos2.") cx a=".$$seq1." (".$$pos1.") ".$$header1b."\n";
		print OUT $$intra_ba."\n";
	}

}

sub printfasta_normal_intra {
	my ($header1a,$header1b,$header1c,$j,$i,$seq1,$seq2,$pos1,$pos2,$intra_ab,$intra_ba) = @_;
	### test if aa or ab peptide. Write only once for aa and ab/ba for others
	if($$intra_ab eq $$intra_ba){
		print OUT ">".$$header1a."_aa_".$$j."_".$$i." a=".$$seq1." (".$$pos1.") cx a=".$$seq2." (".$$pos2.") ".$$header1b."\n";
		print OUT $$intra_ab."\n";
	}else{
		print OUT ">".$$header1a."_ab_".$$j."_".$$i." a=".$$seq1." (".$$pos1.") cx b=".$$seq2." (".$$pos2.") ".$$header1b."\n";
		print OUT $$intra_ab."\n";
		print OUT ">".$$header1a."_ba_".$$i."_".$$j." b=".$$seq2." (".$$pos2.") cx a=".$$seq1." (".$$pos1.") ".$$header1b."\n";
		print OUT $$intra_ba."\n";
	}

}

sub crosslinker {
	my ($pep1,$pep2,$cx_type,$mc_type) = @_;
	### set print decision variables
	my $decision = 0;
	my $amine1_decision = 0;
	my $amine2_decision = 0;
	my ($pos1,$pos2);
	### calculates number of amine/carboxyl/sulfhydryl/tyrosyl/lysine and length of both peptides
	my $carboxyl1 = ($$pep1 =~ tr/DE//);
	my $amine1 = ($$pep1 =~ tr/K//);
	my $sulfhydryl1 = ($$pep1 =~ tr/C//);
	my $tyrosine1 = ($$pep1 =~ tr/Y//);
        my $lysine1 = ($$pep1 =~ tr/K//);
	my $carboxyl2 = ($$pep2 =~ tr/DE//);
	my $amine2 = ($$pep2 =~ tr/K//);
	my $sulfhydryl2 = ($$pep2 =~ tr/C//);
	my $tyrosine2 = ($$pep2 =~ tr/Y//);
        my $lysine2 = ($$pep2 =~ tr/K//);
	my $length1 = length($$pep1);
	my $length2 = length($$pep2);
	### checking position of K in peptide with K=1
	if($amine1 == 1){
		$$pep1 =~ m/K/;
		$pos1 = length($`)+1;
	}
	if($amine2 == 1){
		$$pep2 =~ m/K/;
		$pos2 = length($`)+1;
	}
	### test print
	#print "$$pep1 $$pep2 $$cx_type $$mc_type D1=$carboxyl1 K1=$amine1 D2=$carboxyl2 K2=$amine2\n";
	#print "$length1\t$pos1\n";
	#print "$length2\t$pos2\n";
	### set the amine decision for each peptide based on MC params allowed for subsequent test
	if($$mc_type == 0 && $amine1 > 0){
		$amine1_decision = 1;			# any peptide with K >0 as MC=0
	}
	if($$mc_type == 0 && $amine2 > 0){
		$amine2_decision = 1;			# any peptide with K >0 as MC=0
	}
	if($$mc_type == 1 && $amine1 > 1){
		$amine1_decision = 1;			# any peptide with K >1 like NVAEVKEYVSEK as MC=1
	}
	if($$mc_type == 1 && $amine2 > 1){
		$amine2_decision = 1;			# any peptide with K >1 like NVAEVKEYVSEK as MC=1
	}
	## Added by Michael/Bao
	if($$mc_type == 2 && $amine1 > 2){
		$amine1_decision = 2;			# any peptide with K >2 like NKVAEVKEYVSEK as MC=2
	}
	if($$mc_type == 2 && $amine2 > 2){
		$amine2_decision = 2;			# any peptide with K >2 like NKVAEVKEYVSEK as MC=2
	}
	##------------------------
	if($$mc_type == 1 && $amine1 == 1){
		if($pos1 < $length1){
			$amine1_decision = 1;		# any peptide with K=1 but not terminal like NVAEVKEYVSER as MC=1
		}
	}
	if($$mc_type == 1 && $amine2 == 1){
		if($pos2 < $length2){
			$amine2_decision = 1;		# any peptide with K=1 but not terminal like NVAEVKEYVSER as MC=1
		}
	}
	## Added by Michael/Bao
	if($$mc_type == 2 && $amine1 == 2){
		if($pos1 < $length1){
			$amine1_decision = 2;		# any peptide with K=2 but not terminal like NVAEVKEYVKSER as MC=2
		}
	}
	if($$mc_type == 2 && $amine2 == 2){
		if($pos2 < $length2){
			$amine2_decision = 2;		# any peptide with K=2 but not terminal like NVAEVKEYVKSER as MC=2
		}
	}
	##--------------------------
	### test print
	#print "$amine1_decision\t$amine2_decision\n";
	### test for cross-linkers
	if($$cx_type eq "-x=0"){
		if($amine1_decision == 1 && $amine2_decision == 1){
			$decision = 1;
		}
	}
	### test for Sulfhydryl/Sulfhydryl cross-linkers
	if($$cx_type eq "-x=1"){
		if($sulfhydryl1 > 0 && $sulfhydryl2 > 0){
			$decision = 1;
		}
	}
	### test for Amine/Carboxyl crosslinkers
	if($$cx_type eq "-x=2"){
		if($carboxyl1 > 0 && $amine2_decision == 1){
			$decision = 1;
		}
		elsif($amine1_decision == 1 && $carboxyl2 > 0){
			$decision = 1;
		}
	}
	### test for Amine/Non-selective crosslinkers
	if($$cx_type eq "-x=3"){
		if($amine1_decision == 1 || $amine2_decision == 1){
			$decision = 1;
		}
	}
	### test for Amine/Sulfhydryl crosslinkers
	if($$cx_type eq "-x=4"){
		if($sulfhydryl1 > 0 && $amine2_decision == 1){
			$decision = 1;
		}
		elsif($amine1_decision == 1 && $sulfhydryl2 > 0){
			$decision = 1;
		}
	}
	### test for Tyrosine/Tyrosine cross-linkers
	if($$cx_type eq "-x=5"){
		if($tyrosine1 > 0 && $tyrosine2 > 0){
			$decision = 1;
		}
	}
        ### test for Lysine/Lysine cross-linkers
        if($$cx_type eq "-x=6"){
		if($lysine1 > 0 && $lysine2 > 0){
			$decision = 1;
		}
	}
	### test print
	#print "Decision: $decision\n";
	return $decision;
}
