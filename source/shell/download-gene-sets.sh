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
# Program:  download-gene-sets.sh
# Version:  RSEQREP 2.1.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Download Blood transcription module gene sets, Reactome gene sets, and kegg pathway AMI based gene sets
# Input:    N/A
# Output:   N/A
##############################################################################################################################

## Locate Source Directory from root
SRCDIR="$(cd `dirname $0` && pwd)/../../source";

## identify command line parameters: can be btm, kegg, reactome, or all
if [ -z "$1" ] ; then
	echo 'Please use command line option btm for blood transcription modules gene sets, reactome for Reactome gene sets, kegg for KEGG gene sets or all for all 3.';
	exit 0;
fi
SETS=$1;

## download and unzip blood transcript modules
if [ $SETS = 'btm' ] || [ $SETS = 'all' ] ; then
	wget 'http://www.nature.com/ni/journal/v15/n2/extref/ni.2789-S5.zip';
	unzip ni.2789-S5.zip;
	mv SupplementaryData_TutorialPackage/BTM_for_GSEA_20131008.gmt ./;
	rm -r SupplementaryData_TutorialPackage ni.2789-S5.zip;
	echo 'Blood transcription modules gene sets downloaded';
fi

## Download, filter and reformat reactome data set
if [ $SETS = 'reactome' ] || [ $SETS = 'all' ] ; then
	Rscript $SRCDIR/r/init-reactome-gmt.r;
	echo 'Reactome gene sets downloaded';
fi

## Download Kegg genes, map to pathways and format for GMT
if [ $SETS = 'kegg' ] || [ $SETS = 'all' ] ; then
	read -p 'The KEGG API that is used as part of this program is provided for academic use by academic users belonging to academic institutions or those with a commercial license. Are you an academic user or have a commercial license? [y/n]:  ' YN
	if [ $YN = 'y' ] ; then
		Rscript $SRCDIR/r/init-kegg-gmt.r;
		echo 'KEGG gene sets downloaded';
	else 
		exit 0;
	fi
fi