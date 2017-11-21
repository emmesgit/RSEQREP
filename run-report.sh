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
# Program:  run-report.sh
# Version:  RSEQREP 1.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  script to run reporting only steps (knitR/latex based results reporting)
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

## construct report
rm $ACD/report/cache/* $ACD/report/figs/* $ACD/report/tables/*
Rscript -e "library(knitr); setwd('$ACD/report'); knit('$SRCDIR/knitr/rseq-report.Rnw')";
cd $ACD/report;
pdflatex rseq-report.tex;
pdflatex rseq-report.tex;

## convert PDF to PNG
cd $ACD/report/figs
for f in *.pdf; do
	gs -q -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r300x300 -sOutputFile=${f%.pdf}.png $f
done

## return to RSEQREP dir
cd $SRCDIR/../;