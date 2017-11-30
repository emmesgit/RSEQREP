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
# https://github.com/emmesgit/RSEQREP/SOFTWARE.xlsx  
# 
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# To cite this software, please reference doi:10.12688/f1000research.10464.1
#
# Program:  download-gene-sets.sh
# Version:  RSEQREP 1.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Download Blood transcription module gene sets, Reactome gene sets, and kegg pathway AMI based gene sets
# Input:    N/A
# Output:   N/A
##############################################################################################################################

## download and unzip blood transcript modules
wget 'http://www.nature.com/ni/journal/v15/n2/extref/ni.2789-S5.zip';
unzip ni.2789-S5.zip;
mv SupplementaryData_TutorialPackage/BTM_for_GSEA_20131008.gmt ./;
rm -r SupplementaryData_TutorialPackage ni.2789-S5.zip;

## Download, filter and reformat reactome data set
Rscript ../r/init-reactome-gmt.r

## Download Kegg genes, map to pathways and format for GMT
Rscript ../r/init-kegg-gmt.r