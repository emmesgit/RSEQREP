#############################################################################################################
# RSEQREP: RNA-Seq Reports, an open-source cloud-enabled framework for reproducible
# RNA-Seq data processing, analysis, and result reporting
# 
# https://github.com/emmesgit/RSEQREP
#
# Copyright (C) 2017 The Emmes Corporation 
# 
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
# 
# https://github.com/emmesgit/RSEQREP/SOFTWARE.xlsx  
# 
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# To cite this software, please reference doi:10.12688/f1000research.13049.1
#
# Program:  parse-rseqc-bam-jc-results.pl 
# Version:  RSEQREP 1.1.2
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:	Parse RSeqC junction_annotation.py results
# Input:	list of absolute file paths
# Output: 	N/A
#############################################################################################################

use strict;
use warnings;
use File::Basename;

my $isFirst=1;
while(<>) {
	chomp;
	my $file = $_;
	my $sampleId = basename($file);
	$sampleId =~ s/_bam_jc\.txt//;	
	open(DATA, "<$file") or die "Couldn't open file $file, $!";
	my $resLine=undef;
	my $hedLine=undef;
	
	while(<DATA>){
		chomp;		
		my $line = $_;		
		unless($line=~/Reading reference/ || $line=~ m/Load/ || $line =~ m/^$/ || $line=~ m/^\s/ || $line=~ m/^#/ || $line=~ m/^=/) {
			my @linePart = split('\:',$line,2);			
			my $key  = lc(trim($linePart[0]));
			my $value= trim($linePart[1]);
			$key =~ s/\s\s/ /g;
			$key =~ s/\s/_/g;
			
			
			if(! defined $hedLine) {
				$hedLine="sample_id\t$key";
			} else{
				$hedLine="$hedLine\t$key";
			}		
			
			if(! defined $resLine) {
				$resLine="$sampleId\t$value";
			} else{
				$resLine="$resLine\t$value";
			}			
		}
	}
	if($isFirst) {
		print $hedLine."\n";
	} 
	if(defined $resLine) {
		print $resLine."\n";
	}
	close DATA;
	$isFirst=0;	
}

sub trim {
	my $s = shift;
	$s =~ s/^\s+|\s+$//g ; 
	return $s;
}