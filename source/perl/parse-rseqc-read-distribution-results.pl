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
# Program:  parse-rseqc-read-distribution-results.pl 
# Version:  RSEQREP 2.1.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:	Parse RSeqC read_distribution.py results
# Input:	list of absolute file paths
# Output: 	N/A
#############################################################################################################

use strict;
use warnings;
use File::Basename;

my $isFirst=1;
my $header = "sample_id\ttotal_reads\ttotal_tags\ttotal_assigned_tags\tcds_exons_bases\tcds_exons_tags\tcds_exons_tags_kb\tperc_cds_exons_tags\tfp_utr_exons_bases\tfp_utr_exons_tags\tfp_utr_exons_tags_kb\tperc_fp_utr_exons_tags\ttp_utr_exons_bases\ttp_utr_exons_tags\ttp_utr_exons_tags_kb\tperc_tp_utr_exons_tags\tintron_bases\tintron_tags\tintron_tags_kb\tperc_intron_tags\tintergenic_up_bases\tintergenic_up_tags\tintergenic_up_tags_kb\tperc_intergenic_up_tags\tintergenic_down_bases\tintergenic_down_tags\tintergenic_down_tags_kb\tperc_intergenic_down_tags\tperc_exon_tags\tperc_intergenic_tags";

while(<>) {
	chomp;
	my $file = $_;
	my $sampleId = basename($file);
	$sampleId =~ s/_bam_rc\.txt//;
	
	open(DATA, "<$file") or die "Couldn't open file $file, $!";
	
	my ($totalReads,$totalTags,$totalAssignedTags,$cdsExonBases,$cdsExonTags,$cdsExonTagsKb,$percCdsExonTags,
	$fpUtrExonBases,$fpUtrExonTags,$fpUtrExonTagsKb,$percFpUtrExonTags,$tpUtrExonBases,
	$tpUtrExonTags,$tpUtrExonTagsKb,$percTpUtrExonTags,$intronBases,$intronTags,$intronTagsKb,$percIntronTags,$intergenicUpBases,
	$intergenicUpTags,$percIntergenicUpTags,$intergenicUpTagsKb,$intergenicDownBases,
	$intergenicDownTags,$intergenicDownTagsKb,$percIntergenicDownTags,$percIntergenicTags,$percExonTags) = undef; 
	
	while(<DATA>){
		chomp;
		my $line = $_;
		
		unless($line=~m/^=/ || $line=~ m/^Group/) {	
			if($line=~ m/^Total Reads/)  {
				my @linePart = split('\s+',$line,5);
				$totalReads =trim($linePart[2]);
			} elsif($line=~ m/^Total Tags/) {
				my @linePart = split('\s+',$line,5);
				$totalTags =trim($linePart[2]);
			} elsif($line=~ m/^Total Assigned Tags/) {
				my @linePart = split('\s+',$line,5);
				$totalAssignedTags =trim($linePart[3]);
			} elsif($line=~ m/^CDS_Exons/) {
				my @linePart = split('\s+',$line,5);
				$cdsExonBases=trim($linePart[1]);
				$cdsExonTags=trim($linePart[2]);
				$cdsExonTagsKb = trim($linePart[3]);
				$percCdsExonTags = sprintf("%.4f",$cdsExonTags/$totalAssignedTags*100);
			} elsif($line=~ m/^5'UTR_Exons/) {
				my @linePart = split('\s+',$line,5);
				$fpUtrExonBases=trim($linePart[1]);
				$fpUtrExonTags=trim($linePart[2]);
				$fpUtrExonTagsKb=trim($linePart[3]);
				$percFpUtrExonTags = sprintf("%.4f",$fpUtrExonTags/$totalAssignedTags*100);
			} elsif($line=~ m/3'UTR_Exons/) {
				my @linePart = split('\s+',$line,5);
				$tpUtrExonBases=trim($linePart[1]);
				$tpUtrExonTags=trim($linePart[2]);
				$tpUtrExonTagsKb=trim($linePart[3]);
				$percTpUtrExonTags = sprintf("%.4f",$tpUtrExonTags/$totalAssignedTags*100);
			} elsif($line=~ m/^Introns/) {
				my @linePart = split('\s+',$line,5);
				$intronBases=trim($linePart[1]);
				$intronTags=trim($linePart[2]);
				$intronTagsKb=trim($linePart[3]);
				$percIntronTags = sprintf("%.4f",$intronTags/$totalAssignedTags*100);
			} elsif($line=~ m/^TSS_up_10kb/) {
				my @linePart = split('\s+',$line,5);
				$intergenicUpBases=trim($linePart[1]);
				$intergenicUpTags=trim($linePart[2]);
				$intergenicUpTagsKb=trim($linePart[3]);
				$percIntergenicUpTags = sprintf("%.4f",$intergenicUpTags/$totalAssignedTags*100);
			} elsif($line=~ m/^TES_down_10kb/) {
				my @linePart = split('\s+',$line,5);
				$intergenicDownBases=trim($linePart[1]);
				$intergenicDownTags=trim($linePart[2]);
				$intergenicDownTagsKb=trim($linePart[3]);
				$percIntergenicDownTags = sprintf("%.4f",$intergenicDownTags/$totalAssignedTags*100);
			}
		}		
	}
	if($isFirst) {
		print "$header\n";
		$isFirst = 0;
	}
		
	$percExonTags 		= $percCdsExonTags + $percFpUtrExonTags + $percTpUtrExonTags;
	$percIntergenicTags = $percIntergenicUpTags + $percIntergenicDownTags; 
	
	print "$sampleId\t$totalReads\t$totalTags\t$totalAssignedTags\t$cdsExonBases\t$cdsExonTags\t$cdsExonTagsKb\t$percCdsExonTags\t$fpUtrExonBases\t$fpUtrExonTags\t$fpUtrExonTagsKb\t$percFpUtrExonTags\t$tpUtrExonBases\t$tpUtrExonTags\t$tpUtrExonTagsKb\t$percTpUtrExonTags\t$intronBases\t$intronTags\t$intronTagsKb\t$percIntronTags\t$intergenicUpBases\t$intergenicUpTags\t$intergenicUpTagsKb\t$percIntergenicUpTags\t$intergenicDownBases\t$intergenicDownTags\t$intergenicDownTagsKb\t$percIntergenicDownTags\t$percExonTags\t$percIntergenicTags\n";
}

sub trim {
	my $s = shift;
	$s =~ s/^\s+|\s+$//g ; 
	return $s;
}