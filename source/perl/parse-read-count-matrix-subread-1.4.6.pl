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
# Program:  parse-read-count-matrix-subread-1.4.6.pl
# Version:  RSEQREP 1.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:	Parse read counts from vertical format and outputs gene by sample matrix
# Input:	File list with tab delimited read count results
# Output: 	N/A
#############################################################################################################

use strict;
use warnings;
use File::Basename;

my $hashRef;
my @sampleIds=();

while(<>) {
	chomp;
	my $file = $_;
	
	## id read from file name
	$file=~/([^\/]*)_count.tab$/;
	my $id = $1;
	
	push(@sampleIds,$id);

	open(DATA, "<$file") or die "Couldn't open file $file, $!";
	
	## id read from within file
	my $sampleId = undef;
	
	while(<DATA>){
		unless(m/^# Program/) {
			chomp;
			
			if(m/^Geneid/) {
				(undef,undef,undef,undef,undef,undef,$sampleId) = split('\t',$_);
				$sampleId=~m/([^\/]*)\.bam$/;
				$sampleId=$1;
				#print("file name sample id $id, within sample id: $sampleId\n")
			} else {
				my ($geneId,undef,undef,undef,undef,undef,$geneCount) = split('\t',$_);
				
				## check to see if ids match between file name and in file
				if($sampleId ne $id) {
					die("IDs do not match ($sampleId:$id");
				}
				$hashRef->{$geneId}->{$id}=$geneCount;
			}
		}
	}
	close DATA;		
}

my $row = 'gene_id';
foreach my $id (@sampleIds) {
	$row.="\t$id";
}
$row.="\n";

foreach my $geneId (keys %$hashRef) {
	$row .= $geneId;
	
	foreach my $id (@sampleIds) {
		$row.= "\t$hashRef->{$geneId}->{$id}";
	}
	$row.="\n";
}
print $row;