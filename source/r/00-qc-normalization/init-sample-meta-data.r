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
# Program:  init-sample-meta-data.r
# Version:  RSEQREP 2.1.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose: 	merge sample metadata, rseqc, and featureCounts results
# Input:    <pre_dir>/sample_metadata.csv
#			<rseqc_dir>/bam_<qc|rc|jc|gc>_parsed.tab
#			data/fragment_count_matrix.tab.gz
# Output:  	data/sample_metadata_rseqc.csv
#############################################################################################################

## get command line argument
args = commandArgs(trailingOnly = TRUE)
if(length(args)>0) {
	dta.dir =  args[1];
	rseqc.dir = args[2];
	pre.dir = args[3];
	setwd(dta.dir);
}


## read base information including sample information
dta = read.csv(paste(pre.dir,'/sample_metadata.csv',sep=''),header=T,stringsAsFactors=F,sep=',');

####################################################################################
## MERGE RSEQC GC INFORMATION
####################################################################################
dta.rseqc.gc=read.table(paste(dta.dir,'/bam_gc_parsed.tab',sep=''),
		header=T,stringsAsFactors=F,sep='\t');
colnames(dta.rseqc.gc)=c('samid','rseqc.gc.min','rseqc.gc.q1','rseqc.gc.median','rseqc.gc.mean','rseqc.gc.q3','rseqc.gc.max');
dta.rseqc.gc=dta.rseqc.gc[,c(-2,-6)];
dta =merge(x=dta,y=dta.rseqc.gc,by.x='samid',by.y='samid',sort=F)

####################################################################################
## MERGE RSEQC QC INFORMATION
####################################################################################
dta.rseqc.qc= read.table(paste(dta.dir,'/bam_qc_parsed.tab',sep=''),
		header=T,stringsAsFactors=F,sep='\t');
rseqc.qc.names= paste('rseqc.qc.',names(dta.rseqc.qc)[-1],sep='');
colnames(dta.rseqc.qc)=c('samid',rseqc.qc.names);
dta.rseqc.qc$rseqc.qc.total = dta.rseqc.qc$rseqc.qc.unique + dta.rseqc.qc$rseqc.qc.non_unique
dta=merge(x=dta,y=dta.rseqc.qc,by.x='samid',by.y='samid',sort=F)

####################################################################################
## MERGE RSEQC JC INFORMATION
####################################################################################
dta.rseqc.jc= read.table(paste(dta.dir,'/bam_jc_parsed.tab',sep=''),
		header=T,stringsAsFactors=F,sep='\t');
rseqc.jc.names= paste('rseqc.jc.',names(dta.rseqc.jc)[-1],sep='');
colnames(dta.rseqc.jc)=c('samid',rseqc.jc.names);
dta=merge(x=dta,y=dta.rseqc.jc,by.x='samid',by.y='samid',sort=F)

####################################################################################
## MERGE RSEQC READ DISTRIBUTION INFORMATION
####################################################################################
if (file.exists(paste(dta.dir,'/bam_rc_parsed.tab',sep=''))) {
	dta.rseqc.rc= read.table(paste(dta.dir,'/bam_rc_parsed.tab',sep=''),
		header=T,stringsAsFactors=F,sep='\t');
	rseqc.rc.names= paste('rseqc.rc.',names(dta.rseqc.rc)[-1],sep='');
	colnames(dta.rseqc.rc)=c('samid',rseqc.rc.names);
	dta=merge(x=dta,y=dta.rseqc.rc,by.x='samid',by.y='samid',sort=F)
}

####################################################################################
## MERGE FEATURE COUNTS INFORMATION
####################################################################################
in.file =paste('fragment_count_matrix.tab.gz',sep='');
dta.fcounts= read.table(in.file,header=T,stringsAsFactors=F,sep='\t',row.names=1);

dta.fcounts=dta.fcounts[,match(dta$samid,colnames(dta.fcounts))];
dta$feature.counts.total = apply(dta.fcounts,2,sum);

## order the data by subid, specimen type, treatment, and time
dta=dta[order(dta$subid, dta$spct, dta$trt, dta$time),]
write.csv(dta,file='sample_metadata_rseqc.csv',row.names=F,quote=T,eol="\n")
