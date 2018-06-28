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
# https://github.com/emmesgit/RSEQREP/blob/master/SOFTWARE.xlsx  
# 
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# To cite this software, please reference doi:10.12688/f1000research.13049.1
#
# Program:  parse-rseqc-bam-gc-results.pl 
# Version:  RSEQREP 1.1.3
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:	Parse RSeqC read_GC.py results
# Input:	list of absolute file paths
# Output: 	N/A
#############################################################################################################

use strict;
use warnings;
use File::Basename;

my $isFirst=1;

my $header = "sample_id\tmin\tq1\tmedian\tmean\tq3\tmax";

while(<>) {
	chomp;
	my $file = $_;
	
	open(DATA, "<$file") or die "Couldn't open file $file, $!";
	
	if($isFirst) {
		print $header."\n";
		$isFirst=0;
	} 
	
	while(<DATA>){
		print;
	}
	close DATA;	
}