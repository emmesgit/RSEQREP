#!/bin/bash
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
# Program:  run-pre-processing.sh
# Version:  RSEQREP 1.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  script to run pre-processing only steps (FASTQ download, decrypt, align, 
#				FASTQC, RSEQC, and featureCounts)
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Point to config file
CONFIG=config/config.xlsx

## Locate Source Directory from root
CUR=`pwd`;
SRCDIR="${CUR}/source";

#################
##
## Parse configuration file
##
#################

## Create CSV files for pre-processing, sample metadata, and report 
## configuration based on user input in rnaseq-configuration-template.xlsx 
## convert GTF to bed -- Unless already done
## Index genome -- Unless already done
Rscript $SRCDIR/r/parse-rnaseq-configuration.r $CONFIG $SRCDIR;

## get workflow directory, analysis directory, workflow configuration and metadata csv file locations
WCD=`head -1 dir.csv`;
ACD=`head -2 dir.csv | tail -1`;
WFC=`head -3 dir.csv | tail -1`;
MTA=`tail -1 dir.csv`;

#################
##
## Run pre-processing
##
#################

## Run sample pre-processing Perl script
echo "cd $WCD";
echo "perl $SRCDIR/perl/preprocess-rnaseq.pl $WFC $MTA";
cd $WCD;
perl $SRCDIR/perl/preprocess-rnaseq.pl $WFC $MTA;

#################
##
## Parse pre-processing data
##
#################

## Generate feature counts matrix
echo 'Creating Feature Counts Matrix...';
find $WCD | grep -P '_count.tab$' > feature_counts_outfiles.txt;
perl $SRCDIR/perl/parse-read-count-matrix-subread-1.4.6.pl feature_counts_outfiles.txt > $ACD/data/fragment_count_matrix.tab;
gzip $ACD/data/fragment_count_matrix.tab;
rm feature_counts_outfiles.txt;

## Generate read distribution matrix
echo 'Creating Read Distribution Matrix...';
find $WCD | grep -P 'bam_rc.txt$' > read-distribution_outfiles.txt;
perl $SRCDIR/perl/parse-rseqc-read-distribution-results.pl read-distribution_outfiles.txt > $WCD/rseqc/bam_rc_parsed.tab;
rm read-distribution_outfiles.txt;

## Generate BAM GC matrix
echo 'Creating BAM GC Matrix...';
find $WCD | grep -P 'bam_gc.txt$' > bam_gc_outfiles.txt;
perl $SRCDIR/perl/parse-rseqc-bam-gc-results.pl bam_gc_outfiles.txt > $WCD/rseqc/bam_gc_parsed.tab;
rm bam_gc_outfiles.txt;

## Generate BAM JC matrix
echo 'Creating BAM JC Matrix...';
find $WCD | grep -P 'bam_jc.txt$' > bam_jc_outfiles.txt;
perl $SRCDIR/perl/parse-rseqc-bam-jc-results.pl bam_jc_outfiles.txt > $WCD/rseqc/bam_jc_parsed.tab;
rm bam_jc_outfiles.txt;

## Generate BAM QC matrix
echo 'Creating BAM QC Matrix...';
find $WCD | grep -P 'bam_qc.txt$' > bam_qc_outfiles.txt;
perl $SRCDIR/perl/parse-rseqc-bam-qc-results.pl bam_qc_outfiles.txt > $WCD/rseqc/bam_qc_parsed.tab;
rm bam_qc_outfiles.txt;
