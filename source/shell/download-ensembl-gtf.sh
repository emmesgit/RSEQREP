#!/bin/bash
##############################################################################################################################
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
# Program:  download-ensembl-gtf.sh
# Version:  RSEQREP 2.2.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Download GTF formatted Ensembl genome annotations.
# Input:    N/A
# Output:   N/A
##############################################################################################################################

## Command line arguments
## $1) ensembl version # ie: 87
## $2) download results to directory ie: /home/user01/annot
ensemblVer=$1
resDir=$2

## migrate to results directory
cd $resDir;

## get Ensembl GTF genome annotations
echo "downloading Ensembl Version $ensemblVer GTF annotations..."
echo "wget \"ftp://ftp.ensembl.org/pub/release-$ensemblVer/gtf/homo_sapiens/Homo_sapiens.*.chr.gtf.gz\""
wget -nv "ftp://ftp.ensembl.org/pub/release-$ensemblVer/gtf/homo_sapiens/Homo_sapiens.*.chr.gtf.gz";

## Unzip
gunzip Homo_sapiens.*.chr.gtf.gz;

## rename
mv Homo_sapiens.*.chr.gtf Homo_sapiens.ensembl.version$ensemblVer.chr.gtf

echo 'Done with GTF annotations...'