#!/usr/bin/perl

use Getopt::Std;
use List::Util qw( max);
getopt('hlqasi');


##########################################
# PooLD last updated 16 November 2016
# Information on how to run the program available at
# http://sourceforge.net/p/ldx/wiki/Home/
#
# Alison Feder
# afeder@stanford.edu
##########################################

##########################################
#Copyright (c) 2012, Alison Feder
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
#    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#############################################



##########################################
#Inputs
##########################################
#This is simple parameter setting
##########################################

#print "opt h is $opt_h. opt l is $opt_l\n";

if($opt_h == "" or $opt_h < 0){
    $rejectHigh = 100;
}
else{
    $rejectHigh = $opt_h;
}

#What is too low a read depth?
if($opt_l == "" | $opt_l < 0){
    $rejectLow = 10;
}
else{
    $rejectLow = $opt_l;
}

#What quality scores should we reject below?
if($opt_q == "" | $opt_q < 0){
    $qualityThres = 20;
}
else{
    $qualityThres = $opt_q;
}

#How fine should the maximum likelihood mesh be?
$mleIntervals = .01;

#Noisy?
$print= 0;

#Insert size and
if($opt_s == "" | $opt_s < 0){
    $searchlen = 500;
}
else{
    $searchlen = $opt_s;
}


#Allele frequency
if($opt_a == "" | $opt_a < 0 | $opt_a > 1){
    $allelefreq = .1;
}
else{
    $allelefreq = $opt_a;
}



#Allele frequency
if($opt_i == "" | $opt_i < 0 ){
    $mininter = 11;
}
else{
    $mininter = $opt_i;
}




##########################################
#Initialize some variables to keep track of things we'll want at the end
##########################################

#Number of locations checked as SNPs
$totalSNPs = 0;

#Number of polymorphic SNPs
$polySNPs = 0;

#Number of high read depth or lower minor allele frequency  filtering SNPs
$filtSNPs = 0;

$snpstot = 0;
$snps_rd_filt = 0;
$snps_triallelic = 0;
$snps_low_maf = 0;
$snps_raw = 0;

$pairs_tot = 0;
$pairs_comp = 0;
$pairs_11 = 0;
$pairs_noedge = 0;

$newentry = 0;
$nonpoly = 0;

@triallelic_snps = ();

#Number of total SNP pairs identified

$print = 0;
##########################################
#Convert/Parse the data
##########################################
if($print == 1){ print "Convert and parse stage (0/2)\n";}
#For now, I'll just take .sam files. It wouldn't be too hard to samtools .bam files


#Ensure that we have files provided
if(@ARGV[0] eq "" | @ARGV[1] eq ""){

    if(@ARGV[0] eq ""){
	print "Argv[0] matches null";
    }
    if(@ARGV[1] eq ""){
	print "Argv[1] matches null";
    }
    print "Must provide two input files: @ARGV[0] and @ARGV[1] \n";
    die;

}

#Grab a vcf file and cut out the SNPs
$sam_in = @ARGV[0];
$vcf_in = @ARGV[1];


if($print == 1){ print "Convert and parse stage (1/2)- file locations processed\n";}

##########################################
#SNP processing
##########################################
#Determine which read each SNP maps to
#Determine which SNPs are polymorphic
##########################################

#NOTE: Relies on sorted data


#Read in the vcf file
if($vcf_in =~ /\.gz$/) {

    open(SNPS, "gunzip -c $vcf_in |") or die "problem opening $vcf_in";

}else{
    open(SNPS, $vcf_in) or die;
}

#We want a list of all the polymorphic locations, so we'll parse out the second column
@snps = ();
while(<SNPS>){

	@toParse = split("\t", $_);

	if(@toParse[0] !~ /^#/){

		push(@snps, @toParse[1]);
	}

}

$snpnum= @snps;


for($i = 0; $i< @snps; $i++){
    chomp(@snps[$i]);

}
#@snps should now contain all the polymorphic sites


#initialize some hashes that we'll need
#
#%HoH is going to be indexed by SNPs. Each SNP will point to another hash
#This hash will be indexed by read ID numbers. Each ID will point to an allele.
%HoH =  ();

#freq will also be indexed by SNPs and will keep a running total of frequencies
%freq = ();


#ind represents the index before which each read shouldn't check SNPs.
$ind = 0;

open (FILE, $sam_in);

$start = 0;
$length = 0;
$len_unparsed = 0;


#foreach line in the file,
while(<FILE>){

    @in = split("\t", $_);

    #first parse out the start location and length

    $start = @in[3];


    #This section deals with handling indels
    #Probably the best thing to do here is to insert/delete in a copy of the sequence so that the indexing works out in a straightforward way

    $sequence = @in[9];
    $quality = @in[10];


    #For now, get the length of the sequence
    $length =  length($sequence);
    $stop = 0;

    #Don't go into the indel modification regime if the string matches either 93M or *

    #first, figue out how long the match should be:

    $lenmatch = $length . "M";

    if(@in[5] =~ $lenmatch){
	#print "No indel...\n";
    }
    elsif(@in[5] =~ "\\*"){
	#print "no match. Discard.\n";
	$stop = 1;
    }
    else{
	#print "Indel detected. Correcting.\n";

	#split the CIGAR string into alpha and numeric components
	@alp = split('[0-9]+', @in[5]);
	@num = split('[MIDSHP]', @in[5]);

	$sequence = "";
	$quality = "";

	#$i is the index in the alpha matrix
	$i = 1;
	#$printind is the corresponding index in the string being parsed
	$printind = 0;

	while($i < @alp){

	    #Examine the type of match: match, insertion, deletion
	    $parse = @alp[$i];

	    #If match, simply concatenate the matching portion of the string to the eventual string to be returned
	    if($parse =~ "M"){
		$sub_string =  substr(@in[9], $printind, $num[$i-1]);
		$qual_sub = substr(@in[10], $printind, $num[$i-1]);

		#print "In M loop. Trying to add the $sub_string\n";
		$sequence = $sequence  .  $sub_string;
		$quality = $quality . $qual_sub;

		#print "So far, newstring looks like: $newstring\n";
		$printind = $printind + @num[$i - 1];
	    }
	    #If there's an insertion, don't add anything to the reference string. Skip the indices.
	    elsif($parse =~ "I"){
		#print "So far, newstring looks like: $newstring\n";
		$printind = $printind + @num[$i - 1];
	    }
	    #If there's a deletion, add question marks where the deleted nucleotides would have been
	    elsif($parse =~ "D"){
		#print "In D loop. Trying to add @num[$i-1] ?s\n";
		$sequence = $sequence . ('?'x@num[$i-1]);
		$quality = $quality . ('?'x@num[$i-1]);
		#print "So far, newstring looks like: $newstring\n";
	    }
	    #If there's a soft clip, don't add anything to the reference string. skip the indices
	    elsif($parse =~ "S"){
		$printind = $printind + @num[$i - 1];
	    }
	    #Still need to deal with HP types
	    else{
		#print "I haven't dealt with this CIAGR type yet";
		$stop = 1;
	    }
	    $i++;
	}

    }

    #Now $sequence should contain the string we would like to calculate with.

    #examine each SNP
    #Start the search for the SNPs at the $ind variable


	$tmpind = $ind;

#The idea behind this part is that we can only handle the full length of processing
#if we do it on the fly, throwing out things after they are no longer useful
#After a SNP is out of intersection range with the current index, we'll forget it
    while($tmpind < @snps && $stop != 1){



	#for length, we would like $length to be 93 in most cases, unless there's an insertion or a soft clipping. Have I dealt with insertions yet? Yes

	$query = @snps[$tmpind];

	#print("query is $query\n");
	#print("tmpind is $tmpind\n");
	#print("snps = @snps\n");

	chomp($query);

#	print "We are currently looking at position $query and our read runs from $start to $start + $length\nThe current tmpind is $tmpind (@snps[$tmpind])";

	#if the SNP falls before the start of the read, we no longer need to examine this SNP
	#Therefore, we have all the information we're going to get about this SNP and we can see if we can throw it out because of read depth or polymorphism
	#reasons

	if($query < $start){


#	    print "\n\nProcessing: $query\n\n";

	    $tmpind++;
	    $ind = $tmpind;

#	    print "query was $query, but next query will be @snps[$tmpind] \n";
#	    print "tmpind is currently $tmpind ?\n";
	    #high read depth filter

	    $si = keys %{$HoH{$query}};
	    @si_tmp = keys %{$HoH{$query}};
	    $snpstot++;

		#print "\tChecking to see if anything is around. Si = $si (@si_tmp)\n Comparing to $rejectHigh and $rejectLow\n";

	    if(keys %{$HoH{$query}} >= $rejectHigh){

	    	@tmp22 = keys %{$HoH{$query}};
	#    	print "@tmp22\n";

			#print "Cannot proceed. Read depth too high (higher than $rejectHigh) \n";
	    }
	    if(keys %{$HoH{$query}} < $rejectHigh && keys %{$HoH{$query}} > $rejectLow){
			#print "low enough score\n\n";

		$snps_rd_filt++;

		$first = 1;
		$second = 1;
		$third = 1;
		$fourth = 1;
		$triallelic = 0;

		$Acount = 0;
		$Bcount = 0;
		$Ccount = 0;
		$Dcount = 0;
		$A = "";
		$B = "";
		$C = "";
		$D = "";

	#print "examining new location\n";
		for $read ( keys %{ $HoH{$query} } ) {

#	    print "\t$query: looking at read $HoH{$query}{$read}\n";
		    if($first == 1){
			$A = $HoH{$query}{$read};
			$Acount++;
			$first = 0;
		    }
		    else{

			$read = $HoH{$query}{$read};
			if($read eq $A){
			    $Acount++;
			}
			else{
			    if($second == 1){

				$second = 0;
				$B = $read;
				$Bcount++;

			    }
			    elsif($read eq $B){
					$Bcount++;
			    }
			    else{

			    	if($third == 1){

			    		#If the third, specify third parameters
			    		$third = 0;
			    		$C = $read;
			    		$Ccount++;
			    	}
			    	elsif($read eq $C){
			    		$Ccount++;

			    	}
			    	else{
			    		if($fourth == 1){
			    			$fourth = 0;
			    			$D = $read;
			    			$Dcount++;
			    		}
			    		elsif($read eq $D){
			    			$Dcount++;
			    		}
			    		else{
			    			$triallelic = 0;
			    		}

			    	}
			    }
			}

		    }

#	    print "$read=$HoH{$snploc}{$read} ";
		}

#		print "pre: Acount = $Acount, Bcount = $Bcount, Ccount = $Ccount, DCount = $Dcount. A = $A, B = $B, C = $C, D = $D, triallelic = $triallelic\n";
		#determine two most frequent SNPs
		$mostfreq = $Acount;
		$pos = $A;
		if($Bcount > $mostfreq){
			$mostfreq = $Bcount;
			$pos = $B;
		}
		if($Ccount > $mostfreq){
			$mostfreq = $Ccount;
			$pos = $C;
		}
		if($Dcount > $mostfreq){
			$mostfreq = $Dcount;
			$pos = $D;
		}


		$secfreq = $Acount;
		$secpos = $A;
		if($pos eq $A){
			$secfreq = $Bcount;
			$secpos = $B;
		}
		#determine second most frequent SNP
		if($Bcount > $secfreq && $B ne $pos){
			$secfreq = $Bcount;
			$secpos = $B;
		}
		if($Ccount > $secfreq && $C ne $pos){
			$secfreq = $Ccount;
			$secpos = $C;
		}
		if($Dcount > $secfreq && $D ne $pos){
			$secmostfreq = $Dcount;
			$secpos = $D;
		}

		$Acount = $mostfreq;
		$A = $pos;
		$Bcount = $secfreq;
		$B = $secpos;

	#	print "$query: post: Acount = $Acount, Bcount = $Bcount, Ccount = $Ccount, DCount = $Dcount. A = $A, B = $B, C = $C, D = $D, triallelic = $triallelic\n";

		#Toss triallelic SNPs
		#Edited out. No longer in use
		#Triallelic SNPs are now handled by looking at the two most common SNPs
		if($triallelic == 1){

		    $snps_triallelic++;

		#    print "$query Not included. Triallelic\n";
	#	    for $locus ( keys %{ $HoH{$query} } ) { print "$locus=$HoH{$query}{$locus}\n"; }

		    push(@triallelic_snps, $query);
		    delete $HoH{$query};


		}
		#Check allele frequencies
		elsif($Acount + $Bcount == 0){



###  for $locus ( keys %{ $HoH{$query} } ) { print "$locus=$HoH{$query}{$locus} "; }

#		    print "Deleting %HoH{$query} in acount+bcount = 0\n";

		    delete $HoH{$query};
		    #What does it mean if Acount + Bcount == 0? There's nothing that maps.

		    $nonpoly++;

		}
		else{



		    $PA = $Acount/ ($Acount + $Bcount);
		#	print "Examining $query. PA = $PA, A = $Acount, B = $Bcount\n";
		    #If the minor allele frequency isn't greater than var, toss it.
		    if($PA > 1 - $allelefreq || $PA < $allelefreq){

			$snps_low_maf++;
		#	print "Deleting %HoH{$query} in MAF < .1 \n";
			delete $HoH{$query};

		    }

		    else{


			$maf = $PA;
			if($PA < .5){
			    $maf = 1 - $PA;
			}

			$freq{$query} = $maf;

			#Here, we could potentially check and compute any locations before 300 before it.


#############################
			@snp_spots = sort {$a <=> $b } keys( %HoH);

			#This is the position we're checking everything against

			if(@snp_spots > 1){
			#We should iterate over all possible values of seq2.

			    for($l = 0; $l < @snp_spots - 1; $l++){

				$seq1 = $l;

				$sequence1 = $snp_spots[$seq1];
				$sequence2 = $query;


				if($query > $snp_spots[$seq1] + $searchlen){

				    @keys = sort{$a <=> $b} keys( %HoH);
#####				    print "l = $l, Here are the keys: @keys\n";
#				    delete($HoH{$snps_spots[$seq1]});

				    $tmp1 = $snp_spots[$seq1];
#####				    print "Deleting %HoH{$tmp1} in no longer in range \n";
				    delete $HoH{$tmp1} or die;

				    @keys = sort{$a <=> $b} keys( %HoH);
#####				    print "l = $l, Here are the keys: @keys\n";

				}
				elsif($query > @snp_spots[$seq1]){
				    @snp1 = ();
				    @snp2 = ();

#####				    print "Examinable.\n";

				    #Alright, grab the keys for the first SNP
				    @keys1 = keys %{$HoH{$snp_spots[$seq1]}};

				    @keys2 = keys %{$HoH{$query}};
###				    @keys2 = keys %{$HoH{$snp_spots[$seq2]}};


			#Iterate through these and try to match with the second pair
				    $i = 0;
				    foreach $matchkey (@keys1){

					if(exists $HoH{$query}{$matchkey}){

					    @snp1[$i] = $HoH{$snp_spots[$seq1]}{$matchkey};
					    @snp2[$i] = $HoH{$query}{$matchkey};
###					    @snp2[$i] = $HoH{$snp_spots[$seq2]}{$matchkey};

					    $i++;
					}

				    }

				    $x12 = 0;
				    $x21 = 0;
				    $x22 = 0;
				    $x11 = 1;

				    $a1 = @snp1[0];
				    $a2 = @snp2[0];
				    $b1 = "U";
				    $b2 = "U";

				    $len= @snp1;

#specify the first read in common as 11
#Put in some checking for only two alleles
#Put in some checking for

				    $first = 1;
				    $first_2 = 1;

				    $pairs_tot++;

				    for($i = 1; $i< @snp1; $i++){

					if(@snp1[$i] eq $a1){
					    if(@snp2[$i] eq $a2){

						$x11++;
					    }
					    else{
						if($first == 1){ $b2 = @snp2[$i]; $first = 0 ; }
						if($b2 eq @snp2[$i]){ $x12++;}# print "x12++ ";}

					    }
					}
					else{
					    if($first_2 == 1){ $b1 = @snp1[$i]; $first_2 = 0; }

					    if($b1 eq @snp1[$i]){
						if(@snp2[$i] eq $a2){

						    $x21++;
						}
						else{
						    if($first == 1){ $b2 = @snp2[$i]; $first = 0 ; }
						    if($b2 eq @snp2[$i]){ $x22++; }#print "x22+ "; }
						    else{

						    }
						}
					    }
					    else{

					    }

					}
				    }

				    $sum = $x11 + $x12 + $x21 + $x22;

				    #Let's sort for major and minor alleles!


				    if($x11 + $x12 < $x21 + $x22){

					#switch the nucleotides

					$tmp = $a1;
					$a1 = $b1;
					$b1 = $tmp;

					#switch the rows

					$tmp = $x11;
					$x11 = $x21;
					$x21 = $tmp;

					$tmp = $x12;
					$x12 = $x22;
					$x22 = $tmp;
				    }

				    if($x11 + $x21 < $x12 + $x22){

					#switch the nucleotides

					$tmp = $a2;
					$a2 = $b2;
					$b2 = $tmp;

					#switch the cols

					$tmp = $x11;
					$x11 = $x12;
					$x12 = $tmp;

					$tmp = $x21;
					$x21 = $x22;
					$x22 = $tmp;
				    }


				    $PA = $freq{@snp_spots[$seq1]};
				    $PB = $freq{$query};
			#    	print "What are we trying to find? @snp_spots[$seq1] \n";
			#    	print "looking at $PA and $PB\n";
###				    $PB = $freq{@snp_spots[$seq2]};

				    $x11f = $x11/$sum;
				    $x12f = $x12/$sum;
				    $x21f = $x21/$sum;
				    $x22f = $x22/$sum;

				    if($x11 + $x12 != 0 && $x11 + $x21 != 0 && $x22 + $x21 != 0 && $x22 + $x12 != 0){

					$pairs_comp++;

					$R2 = (($x11f - ($x11f + $x12f)*($x11f + $x21f)) ** 2)/(($x11f + $x12f) * ($x12f + $x22f) * ($x11f + $x21f) * ($x21f + $x22f));

					$readdepth1 = keys %{$HoH{$snp_spots[$seq1]}};
					$readdepth2 = keys %{$HoH{$query}};
###					$readdepth2 = keys %{$HoH{$snp_spots[$seq2]}};

					$shareddepth = @snp1;
					$shareddepth_true = $x11 + $x12 + $x21 + $x22;

					if($shareddepth_true > $mininter){

					    $pairs_11 ++;
					    $x = 0;
					    $prind = 0;
					    @r2 = ();
					    @vals = ();


					    while($x <1.0000000000001 & $exclude == 0){


						$D = sqrt($x * $PB * (1-$PB) * $PA * (1-$PA));

		    #Calculate the expected values in the 2x2 table
						@exp = ();

						@exp[0] = ($PA * $PB + $D);
						@exp[1] = ($PA * (1-$PB) - $D);
						@exp[2] = ((1-$PA) * $PB - $D);
						@exp[3] = ((1-$PA) * (1-$PB) + $D);

						if(@exp[1] <= 0 | @exp[2] <= 0){

						    @vals[$prind] = -10000;

						}
						else{

						    $match = $x11* log(@exp[0]) + $x12*log(@exp[1]) +$x21*log(@exp[2]) + $x22*log(@exp[3]) ;

						    @vals[$prind] = $match;
						    @r2[$prind] = $match;

						}

						$prind++;
						$x = $x+$mleIntervals;

					    }

					    $highestInd = 0;

					    for($i = 1; $i<@vals; $i++){

						if(@vals[$highestInd]< @vals[$i]){ $highestInd = $i;}# print "\tSwitching\n";}

					    }
#So, now we have our max
#Find a point that's higher than it.

					    $j = $highestInd;

					    while($j < @vals & @vals[$j] > @vals[$highestInd] - 1.92){
							$j++;
					    }

					    $size = @vals;

					    if(@vals[$j] == @vals[$size]){ $high = "2"; }
					    else{ $high = ($j-1) * $mleIntervals; }


					    $j = $highestInd;
					    while($j >= 0 & @vals[$j] > @vals[$highestInd] - 1.92){
						$j--;
					    }
					    if($j == 0){$low = 0; }
					    else{ $low = ($j+1) * $mleIntervals; }

					    $r2 = $mleIntervals*$highestInd;

		#We really only want certain lines:
		#those in which low != r2 and high != r2

					    #if($low != $r2 && $high != $r2 && $low > 0){

						$pairs_noedge++;

						$fPA = ($x11 + $x12) / ($x11 + $x12 + $x21 + $x22);
						$fPB = ($x11 + $x21) / ($x11 + $x12 + $x21 + $x22);
						$r22x2 = ($x11f - $fPA * $fPB)*($x11f - $fPA * $fPB) / ($fPA * $fPB * (1- $fPA) * (1-$fPB));



						print "@snp_spots[$seq1]\t$query\t$x11\t$x12\t$x21\t$x22\t$PA\t$PB\t$readdepth1\t$readdepth2\t$shareddepth_true\t$low\t$r2\t$high\t$r22x2\t$a1\t$b1\t$a2\t$b2\n";
###						print "@snp_spots[$seq1]\t@snp_spots[$seq2]\t$x11\t$x12\t$x21\t$x22\t$PA\t$PB\t\t$readdepth1\t$readdepth2\t$shareddepth_true\t$low\t$r2\t$high\n";
					   #}
					}

				    }
				}
				else{
#####				    print "Not yet examinable\n";
				}

			    }
			}

		    }

		}

	    }
	    else{

		#read deptht oo high. Throw it out
		delete $HoH{$query};

	    }

	}

	#if it falls within the range, do several things
	elsif($query >= $start && $query <= $start + $length){

#   	 print "\tStatus: in range. Incrementing tmpind.\n";
	    #find the allelic state


	    $localind = $query - $start;
	    $allele = substr($sequence, $localind, 1);
	    $score = substr($quality, $localind, 1);
	    $scorepr = ord($score)-33;

	    if($allele eq ""){
		#don't output, I guess
	    }
	    elsif($qualityThres < ord($score)-33 && $allele ne "?"){
#		print "Found allele: $allele\n";

		if(!defined $HoH{$query}){ $newentry++; }
		$HoH{$query}{@in[0]} = $allele;
#		print "HoH{$query}{@in[0]} = $allele - $HoH{$query}{@in[0]} \n";

	    }



	    $snps_raw++;

	    #write to the SNP file
	    $tmpind++;

	}
	#if it falls beyond the range, stop looking at SNPs on this read
	elsif($query > $start + $length){
	    $stop = 1;
#    	    print "\tStatus: after range. Starting new read\n";
	}

    }
}
