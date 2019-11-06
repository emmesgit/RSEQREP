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
# Program:  run-pre-processing.sh
# Version:  RSEQREP 2.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  script to run pre-processing only steps (FASTQ download, decrypt, align, 
#				FASTQC, RSEQC, and featureCounts)
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Locate Source Directory from root
SRCDIR="$(cd `dirname $0` && pwd)/source";

## add config script for fastq-dump to throw the cache in /tmp.
## otherwise .sra cache files are stored on the home directory
mkdir $HOME/.ncbi -p;
echo '/repository/user/main/public/root = "/tmp"' > $HOME/.ncbi/user-settings.mkfg

## Pass command line arguments to python program
python3 $SRCDIR/python/rseqrep --run-preprocessing $@