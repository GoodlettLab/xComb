#!/usr/bin/perl

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

### declare arguments for this program, check arguments and store argument values
unless(@ARGV==4){
	print "\nWrong or missing arguments!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
### Input format type
my $type = $ARGV[0];
unless($type =~ m/-in=[0-2]{1}/){
	print "\nWrong or missing arguments!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
### Enzyme
my $enzyme = $ARGV[1];
unless($enzyme =~ m/-ez=[0-9]{1}/){
	print "\nWrong or missing arguments!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
### Missed cleavage
my $MC = $ARGV[2];
unless($MC =~ m/-mc=[0-3]{1}/){
	print "\nWrong or missing arguments!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
### input file
my $infile = $ARGV[3];
### check for correct file extension
unless($infile =~ m/.*\.fasta$/ || $infile =~ m/.*\.dat$/){
	print "\nFile extension must be of type .fasta or .dat!!! See instructions:\n\n";
	&instruction;
	exit 1;
}
### retrieve MC value
$MC =~ s/-mc=//;
### pass input file to directory variable
my $directory = $infile;
### open any *.fasta file in the current directory
open(IN, "$infile");

### set std line break according to format
if($type eq "-in=1"){
	### change the standard line break input reader to > for FASTA format
	$/=">";
	### set the cursor to start of file + 1 character, (i.e. skip first ">")
	seek(IN, 1, 0);
	### create directory name based on input filename
	$directory =~ s/.fasta//i;
}
if($type eq "-in=2"){
	### change the standard line break input reader to > for FASTA format
	$/=">";
	### set the cursor to start of file + 1 character, (i.e. skip first ">")
	seek(IN, 1, 0);
	### create directory name based on input filename
	$directory =~ s/.fasta//i;
}
elsif($type eq "-in=0"){
	### change the standard line break input reader to // for DAT format
	$/="//\n";
	### create directory name based on input filename
	$directory =~ s/.dat//i;
}

### create directory for digest file based on filename
mkdir("$directory",0777) || die "Cannot mkdir $directory!!!\n\n";

### loops through the lines
my $prot_count = 1;
while(<IN>) {
  chomp $_;
  ### remove CR LF for Dos and replace with LF for Unix
  $_ =~ s/\015\012/\012/g;
  ### remove CR for Mac and replace with LF for Unix
  $_ =~ s/\015/\012/g;
	my ($header,$AC,$ID,$DE,$sequence,@header);
	my $FToffset = 0;
	if($type eq "-in=1"){
		### get header information
		$header = &getheader_fasta($_);
		$header =~ m/([A-Z]{1}[0-9]{1}[A-Z0-9]{1}[A-Z0-9]{1}[A-Z0-9]{1}[0-9]{1})\|([A-Z0-9]{1,}_[A-Z0-9]{1,}) (.*)/;
		$AC = $1;
		$ID = $2;
		$DE = $3;
		### get sequence information
		$sequence = &getsequence_fasta($_);
	}
	if($type eq "-in=2"){
		### get header information
		$header = &getheader_fasta($_);
		$AC = "PROT".$prot_count;
		$ID = "N/A";
		$DE = $header;
		### get sequence information
		$sequence = &getsequence_fasta($_);
	}
	elsif($type eq "-in=0"){
		### get AC, ID, DE and sequence
		@header = &getdat_info($_);
		$AC = $header[0];
		$ID = $header[1];
		$DE = $header[2];
		$sequence = $header[3];
		$FToffset = $header[4];
		### test print
		#print "$AC\t$ID\t$DE\n$sequence\n";

	}
	### create file nased on AC number with extension .digest
	open(OUT,">$directory/$AC.digest");
	### test print
	#print "$header\n$sequence\n";
	### print AC/DE/ID and header line
	print OUT "AC=$AC\nDE=$DE\nID=$ID\nmass\tposition\t#MC\tpeptide sequence\n";
	### calculate sequence length
	my $countAZ = ( $sequence =~ tr/[a-zA-Z]//);
	### declare an array which will store the position of all cleavage sites
	my @position;
	### loops through sequence, find cleavage sites according to enzyme specificity, return position and add to array
	### TRYPSIN
	if($enzyme eq "-ez=0"){
		while ($sequence =~ m/[KR]/ig){ # match any K or R followed by any aa
			my $start = pos($sequence); # returns the position of the match to $start
			my $match = substr($sequence,$start-1,2); # find the di-aa sequence of the cleavage site
			#### add cleavage site to array if not equals to "KP" or "RP"
			if ($start < $countAZ && $match ne "KP" && $match ne "RP"){ # if not equal to last aa + eq to KP or RP, then add to array
			push(@position,$start); # add position to the array
			}
		}
	}
	### LYS-C
	if($enzyme eq "-ez=1"){
		while ($sequence =~ m/[K]/ig){ # match any K followed by any aa
			my $start = pos($sequence); # returns the position of the match to $start
			if($start ne $countAZ){ # if not equal to last position, then add to array
				push(@position,$start); # add position to the array
			}
		}
	}
	### ARG-C
	if($enzyme eq "-ez=2"){
		while ($sequence =~ m/[R]/ig){ # match R followed by any aa
			my $start = pos($sequence); # returns the position of the match to $start
			my $match = substr($sequence,$start-1,2); # find the di-aa sequence of the cleavage site
			#### add cleavage site to array if not equals to "KP" or "RP"
			if ($start < $countAZ && $match ne "RP"){ # if not equal to last aa + eq to RP, then add to array
			push(@position,$start); # add position to the array
			}
		}
	}
	### GLU-C (E)
	if($enzyme eq "-ez=3"){
		while ($sequence =~ m/[E]/ig){ # match any E followed by any aa
			my $start = pos($sequence); # returns the position of the match to $start
			my $match = substr($sequence,$start-1,2); # find the di-aa sequence of the cleavage site
			#### add cleavage site to array if not equals to "KP" or "RP"
			if ($start < $countAZ && $match ne "EE" && $match ne "EP"){ # if not equal to last aa + eq to EE or EP, then add to array
			push(@position,$start); # add position to the array
			}
		}
	}
	### GLU-C (DE)
	if($enzyme eq "-ez=4"){
		while ($sequence =~ m/[DE]/ig){ # match any D or E followed by any aa
			my $start = pos($sequence); # returns the position of the match to $start
			my $match = substr($sequence,$start-1,2); # find the di-aa sequence of the cleavage site
			#### add cleavage site to array if not equals to "KP" or "RP"
			if ($start < $countAZ && $match ne "DE" && $match ne "DP" && $match ne "EE" && $match ne "EP"){ # if not equal to last aa + eq to DE, DP, EE or EP, then add to array
			push(@position,$start); # add position to the array
			}
		}
	}
	### ASP-N
	if($enzyme eq "-ez=5"){
		while ($sequence =~ m/[D]/ig){ # match any D followed by any aa
			my $start = pos($sequence)-1; # returns the position of the match to $start minus 1 for N-terminal cleavage position
			if($start ne $countAZ){ # if not equal to last position, then add to array
				push(@position,$start); # add position to the array
			}
		}
	}
	### TRYPSIN, NO P RULE
	if($enzyme eq "-ez=6"){
		while ($sequence =~ m/[KR]/ig){ # match any K or R followed by any aa
			my $start = pos($sequence); # returns the position of the match to $start
			if($start ne $countAZ){ # if not equal to last position, then add to array
				push(@position,$start); # add position to the array
			}
		}
	}
	### ARG-C, NO P RULE
	if($enzyme eq "-ez=7"){
		while ($sequence =~ m/[R]/ig){ # match any R followed by any aa
			my $start = pos($sequence); # returns the position of the match to $start
			if($start ne $countAZ){ # if not equal to last position, then add to array
				push(@position,$start); # add position to the array
			}
		}
	}
	### Chymotrypsin
	if($enzyme eq "-ez=8"){
		while ($sequence =~ m/[FLMWY]/ig){ # match any F,L,M,W or Y followed by any aa
			my $start = pos($sequence); # returns the position of the match to $start
			my $match = substr($sequence,$start-1,2); # find the di-aa sequence of the cleavage site
			#### add cleavage site to array if not equals to "KP" or "RP"
			if ($start < $countAZ && $match ne "FP" && $match ne "LP" && $match ne "MP" && $match ne "WP" && $match ne "YP"){ # if not equal to last aa + eq to FP, LP, MP, WP or yP, then add to array
			push(@position,$start); # add position to the array
			}
		}
	}
	### Thermolysin (WEF)
	#### Note that Thermolysin cuts at the N-term of L, V, I, A, M, F so the "-1" are needed
	if($enzyme eq "-ez=9"){
		while($sequence =~ m/[LVIAMF]/ig){ # match any L, V, I, A, M, F
			my $start = pos($sequence)-1; # returns the position of the match to $start
	  	if($start ne ($countAZ-1)){ # if not equal to last position, then add to array
	    	push(@position,$start); # add position to the array
			}
	  }
	}


	### declare a new array
	my @position1;
	@position1 = sort {$a <=> $b} @position;# sort first array in numeric order and store in new array
	push (@position1,$countAZ);# add position of last aa to array
	unshift (@position1,0);# add 1st position of aa which is 0 and not 1 for a computer array
	### test print
	#print "@position1\n$countAZ\n@position\n";
	### loops through until all position have been treated, store sequence and calculates mass of each peptide
	for (my $i=1; $i<@position1; ++$i){
		for (my $j=0; $j<=$MC; ++$j){
			### calculates length of peptide
			my $length = $position1[$i+$j] - $position1[$i-1];
			### test if length is negative meaning at the end of peptide and cannot write nb peptides proportional to MC
			if($length < 1){
				next;
			}
			my $seq = substr($sequence,$position1[$i-1],$length);
			### substract each aa and assign corresponding value in hash to a scalar. Then sum through out the loop each Mr
			my $sumMr1 = &getmass($seq);
			### Adds H2O to the peptide for exact mass (H+ for N-terminus and OH for C-terminus)
			my $sumMr1final = $sumMr1 + 18.01056;
			### calculates start and end positions
			my $start = $position1[$i-1]+1+$FToffset;
			my $end = $position1[$i+$j]+$FToffset;
			printf OUT "%.4f", $sumMr1final;
			print OUT "\t$start-$end\t$j\t$seq\n";
			### Reinitialize variables for next while loop
			$sumMr1 = 0;
		}
	}
	close OUT;
	++$prot_count;
}

close IN;
exit 0;

##########		Subroutine(s) section		##########

sub instruction {
print '
protein2digest.pl (xComb v1.4)

FUNCTION:
Transforms DAT or FASTA protein files into DIGEST format for further use with compatible
softwares.

USAGE:
protein2digest.pl -in=[value] -ez=[value] -mc=[value] [file_name]

	-in=[value]: Specify input format

	  0= UniProt DAT
	  1= UniProt FASTA
	  2= Any FASTA

	-ez=[value]: Specify enzyme

	  0= Trypsin
	  1= Lys-C
	  2= Arg-C
	  3= Glu-C (E)
	  4= Glu-C (DE)
	  5= Asp-N
	  6= Trypsin, no P rule
	  7= ARG-C, no P rule
	  8= Chymotrypsin, no P rule
		9= Thermolysin

	-mc=[value]: Specify nb of missed cleavage allowed (max. = 3)

	[file_name]: Specify file to be converted

Ex: protein2digest.pl -in=1 -ez=0 -mc=2 filename


';
}

sub getsequence_fasta {
  my @entry = split("\n",$_);
	my $a = @entry;
	my $seq = "";
	for (my $i=1; $i <= $a; $i +=1){
		$seq .= $entry[$i];
	}
	return $seq;
}

sub getheader_fasta {
  my @entry = split("\n",$_);
	my $header = $entry[0];
	return $header;
}

sub getdat_info {
	my ($AC,$ID,$DE,$SQ,$start,$end,$seq,$FToffset);
	### retrieves ID
	$_ =~ m/ID[\s]{1,}([A-Z0-9]{1,}_[A-Z0-9]{1,}).*/;
	$ID = $1;
	### retrieves AC
	$_ =~ m/AC[\s]{1,}([A-Z0-9]{6}).*/;
	$AC = $1;
	### retrieves DE
	$_ =~ m/DE[\s]{1,}RecName: Full=(.*);.*/;
	$_ =~ m/DE[\s]{1,}SubName: Full=(.*);.*/;
	$DE = $1;
	### retrieves SQ
	$_ =~ m/SQ[\s]{1,}.*;\n([\w\W]*)/;
	$SQ = $1;
	$SQ =~ s/\W//g;
	### retrieves FT CHAIN information
	if ($_ =~ m/FT[\s]{1,}CHAIN[\s]{1,}([0-9]{1,})[\s]{1,}([0-9]{1,}).*/){
		$start = $1;
		$end = $2;
	}else{
		$start = 1;
		$end = length($SQ);
	}
	$FToffset = $start-1;
	### test print
	print "$AC\t$ID\t$DE\t$start\t$end\n$SQ\n\n";
	$seq = substr($SQ,$start-1,($end-$start)+1);
	return ($AC,$ID,$DE,$seq,$FToffset);
}

sub getmass {
	my ($seq) = @_;
	my $aamass = 0;
	my $protmass = 0;
	my(%aa_Mr) = (
		'A' => '71.03711', 'B' => '114.53494', 'C' => '103.00919', 'D' => '115.02694', 'E' => '129.04259', 'F' => '147.06841', 'G' => '57.02146', 'H' => '137.05891', 'I' => '113.08406', 'K' => '128.09496', 'L' => '113.08406',
		'M' => '131.04049', 'N' => '114.04293', 'P' => '97.05276', 'Q' => '128.05858', 'R' => '156.10111', 'S' => '87.03203', 'T' => '101.04768', 'U' => '150.95364', 'V' => '99.06841', 'W' => '186.07931', 'X' => '110.13560', 'Y' => '163.06333',
		'Z' => '128.55059');
	for (my $i=0; $i < (length($seq)); $i += 1) {
		$aamass = $aa_Mr{substr($seq,$i,1)};
		$protmass += $aamass;
		}
	return $protmass;
}
