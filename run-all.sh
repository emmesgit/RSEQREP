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
#   example: sh run-all.sh --config config/config.xlsx --log /home/repuser/run_all_log.txt
#
# Program:  run-all.sh
# Version:  RSEQREP 2.1.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Master Shell script -- runs analysis from start to end.
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Locate Source Directory from root
SRCDIR="$(cd `dirname $0` && pwd)/source";


#################
##
## Run pre-processing
##
#################

## add config script for fastq-dump to throw the cache in /tmp.
## otherwise .sra cache files are stored on the home directory
mkdir $HOME/.ncbi -p;
echo '/repository/user/main/public/root = "/tmp"' > $HOME/.ncbi/user-settings.mkfg

## Run sample pre-processing Perl script
python3 $SRCDIR/python/rseqrep --run-preprocessing $@

## get workflow directory, analysis directory, workflow configuration and metadata csv file locations
WCD=`head -1 $SRCDIR/../dir.csv`;
ACD=`head -2 $SRCDIR/../dir.csv | tail -1`;
WFC=`head -3 $SRCDIR/../dir.csv | tail -1`;
MTA=`tail -1 $SRCDIR/../dir.csv`;

#################
##
## Initialize analysis
##
#################

## migrate to R directory
cd $SRCDIR/r

## update sample metadata -- add bam statistics
Rscript 00-qc-normalization/init-sample-meta-data.r $ACD/data $WCD/rseqc $WCD;

## import annotations from ensembl
Rscript 00-qc-normalization/init-export-annotations.r;

## perform TMM normalization
Rscript 00-qc-normalization/init-tmm-normalization-fragments.r;

## calculate log fold changes
Rscript 00-qc-normalization/init-log-cpm-fold-change-from-baseline.r

## calculate PCA, Eucledian, and Spearman distances (Run in parallel)
Rscript 01-bias-confounding-effects/init-euclid-dist.r &
Rscript 01-bias-confounding-effects/init-spearman-dist.r &
Rscript 01-bias-confounding-effects/init-pca.r &
wait;

## discover significantly differentially expressed genes (SDEG)
Rscript 02-sdeg-identification/init-edgeR-glm-model.r;

## initialize pvclusters
Rscript 03-sdeg-clusters/init-pvclusters.r;

## run gene set enrichment analysis (GSEA)
Rscript 04-sdeg-organization-known-modules/init-gsea-sampling.r;

## construct report
rm $ACD/report/cache/* $ACD/report/figs/* $ACD/report/tables/*
Rscript -e "library(knitr); setwd('$ACD/report'); knit('$SRCDIR/knitr/rseq-report.Rnw')";
cd $ACD/report;
pdflatex rseq-report.tex;
pdflatex rseq-report.tex;