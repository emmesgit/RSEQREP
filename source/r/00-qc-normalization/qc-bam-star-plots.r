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
# Program:  qc-bam-star-plots.r
# Version:  RSEQREP 2.0.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Plot starplots to evaluate treatment and time effects for each specimen type
# Input:    N/A
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

bam.qc.vars = c(
		'rseqc.qc.total',
		'rseqc.qc.unmapped',
		'rseqc.qc.unique',
		'rseqc.qc.plus_strand',
		'rseqc.qc.minus_strand',
		'rseqc.gc.median',
		'feature.counts.total');

bam.qc.labs = c(
		'#Total Mapped Reads',
		'#Unmapped Reads',
		'#Uniquely Mapped Reads',
		'#Uniquely Mapped + Strand',
		'#Uniquely Mapped - Strand',
		'Median GC Content',
		'Counted Fragments');

## loop through specimen types
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];
	spcLabl = spcLabls[s];
	
	## limit metadata to specimen type
	mta.spc=mta[which(mta$spct==spcFlag),]

	## Define subjects+treatment and timepoints
	subids.trt = unique(mta.spc[,c('subid','trt')])
	time   = unique(mta.spc[,'time'])
	
	## generate factor of timepoints
	mta.spc$timen = as.numeric(factor(mta.spc$time))
		
	## set x (time) and y (subid+trt) coordinates
	locations=matrix(ncol=2,nrow=nrow(mta.spc))
	for(r in 1:nrow(mta.spc)) {
		locations[r,1]=mta.spc[r,'timen'];
		locations[r,2]=which(subids.trt[,1]==mta.spc[r,'subid'] & subids.trt[,2]==mta.spc[r,'trt']);
	}
			
	## prepare plot
	plot(NA,NA,ylim=c(0,nrow(subids.trt)+2),xlim=c(0,length(time)+1),main='',yaxt="n",axes = F,ylab='',xlab='');
		
	## plot stars
	stars(mta.spc[,bam.qc.vars],key.loc = c((length(time)+1)/2,nrow(subids.trt)+1.5), 
			labels=NA,cex=0.5,
			flip.labels = F,col.stars=mta.spc$trtc,frame.plot=F,
			locations=as.matrix(locations),	draw.segments = F, len=0.52,
			key.labels = bam.qc.labs,
			axes=F,add=T);
	axis( 1, at = 1:length(time),   labels =time,cex.axis=0.8,line=0,tick=T);
	axis( 2, at = 1:nrow(subids.trt), labels =subids.trt[,1],cex.axis=0.8,line=0,tick=T,las=2,cex.axis=0.6);
	mtext( side = 2, text='Subject', line = 3.4, cex=0.8)
	mtext( side = 1, text='Study Visit Day', line = 2.2, cex=0.8)
	mtext( side = 3, text=paste(spcLabl,'Reference Alignment Statistics',sep='\n'), cex=1.5 ,line=1.4); ## add subtitle
	legend('topleft',legend=unique(mta.spc$trtl),fill=unique(mta.spc$trtc),cex=0.7,border = FALSE,bty = "n");
	box();
}