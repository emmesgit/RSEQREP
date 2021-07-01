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
# Program:  ma-plots.r
# Version:  RSEQREP 2.1.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Plot concentration (average log counts per million) by log fold change. capture 0.1%-99.9%.  
#			If all are printed, extreme outliers pull the attention away from the important part of the plot. 
# Input:    analysis/glm/<scp>_<trt>_<time>_glm.RData
#			analysis/glm/<scp>_<trt>_<time>_glm_sig.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## set directories
in.dir.glm	 = paste(res.dir,'glm',sep='/');

## specify plotting parameters
par(mfrow=c(5,2),mar=c(3,3.05,2,0.2))

## obtain min and max lofFC and average logCPM
xlim=c(0,0)
ylim=c(0,0)

## loop through specimen types;
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];
	
	## loop through treatment groups;
	for (v in 1:length(trtFlags)) {
		trtFlag = trtFlags[v];
		
		## loop through each post_treatment time
		for(t in 1:length(postb.times)) {
			time = postb.times[t];
			
			## load glm R objects 	
			in.file.glm = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm.RData',sep='');
			load(in.file.glm);
			
			## update min and max values -- use 99.9% and 0.1% intervals
			if(quantile(glm$DGELRT$table$logCPM,0.999) > xlim[2]) {xlim[2]=ceiling(quantile(glm$DGELRT$table$logCPM,0.999))}
			if(quantile(glm$DGELRT$table$logCPM,0.001) < xlim[1]) {xlim[1]=floor(quantile(glm$DGELRT$table$logCPM,0.001))}
			if(quantile(glm$DGELRT$table$logFC,0.999) > ylim[2]) {ylim[2]=ceiling(quantile(glm$DGELRT$table$logFC,0.999))}
			if(quantile(glm$DGELRT$table$logFC,0.001) < ylim[1]) {ylim[1]=floor(quantile(glm$DGELRT$table$logFC,0.001))}
		}
	}
}

## loop through specimen types;
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];
	spcLabl = spcLabls[s];
	
	## loop through treatment groups;
	for (v in 1:length(trtFlags)) {
		trtFlag = trtFlags[v];
		trtLabl = trtLabls[v];
		
		## loop through each post_treatment time
		for(t in 1:length(postb.times)) {
			time = postb.times[t];
			timel = postb.times.l[t];
			
			## load significant genes
			dge.spc.time.sig = c()
			in.file.spc.time.sig = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz',sep='')
			if (file.exists(in.file.spc.time.sig)) {
				dge.spc.time.sig = read.table(in.file.spc.time.sig,header=T,stringsAsFactors=F);
			}
			
			## load glm R objects 	
			in.file.glm = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm.RData',sep='');
			load(in.file.glm);
			
			plotSmear(	glm$DGELRT,
						de.tags=rownames(dge.spc.time.sig),
						xlim=xlim,ylim=ylim,cex.main=1,cex.axis=0.8,cex=0.3,
						main=paste(spcLabl,', ',trtLabl,'\n(',timel,' vs. Pre-treatment)',sep=''),
						panel.first=abline(h=c(log2(1/as.numeric(glm.sdeg.fold)),log2(as.numeric(glm.sdeg.fold))), col="blue",ylab='',xlab=''));
				mtext(expression('average log'[2]*' cpm'),side=1,line=1.8,cex=0.7)
				mtext(expression('log'[2]*' fold change'),side=2,line=1.9,cex=0.7)
		}
		## add blank plots for unused space.
		if((t %% 5)!=0){rep(plot.new(),(5 - (t %% 5))*2)}
	}
}
