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
# Program:  setup.sh
# Version:  RSEQREP 2.3.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  script to run setup only steps (R script to run sanity checks on configuration, 
#				download genome/transcriptome, format GMT files, and parse configuration)
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Locate Source Directory from root
SRCDIR="$(cd `dirname $0` && pwd)/source";

## Create CSV files for pre-processing, sample metadata, and report 
## configuration based on user input in rnaseq-configuration-template.xlsx 
## convert GTF to bed -- Unless already done
## Index genome -- Unless already done
python3 $SRCDIR/python/rseqrep $@