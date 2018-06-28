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
# Program:  tmm-box-plots.r 
# Version:  RSEQREP 1.1.3
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Plot boxplots for  pre and post normalization data for each 
#				specimen type colored by treatment.
# Input:    analysis/lcpm/all_alld_lcpm_<nrm>_normalized_filtered.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

in.dir.lcpm = paste(res.dir,'lcpm',sep='/');

## order metadata by spcimen type, patid, and treatment group
mta = mta[order(mta$spct, mta$subid, mta$trt, mta$time),]

## loop through normalization types
for(n in 1:length(nrmFlags)) {
	nrmFlag = nrmFlags[n];
	nrmLabl= nrmLabls[n];
	
	## read data
	in.file = paste(in.dir.lcpm,'/all_alltp_lcpm_',nrmFlag,'_normalized_filtered.tab.gz',sep='');	
	dta = read.csv(in.file,header=T,stringsAsFactors=F,check.names = FALSE,sep='\t');
	
	## match cpm library with meta data samples 
	dta = dta[,mta$samid];
	
	## define plot limits -- based on boxplot() function in R
	ylim = c(floor(min(apply(dta,2,function(x){max(min(x), summary(x)[2] - 1.5 * (summary(x)[5]-summary(x)[2]))}))),
			ceiling(max(apply(dta,2,function(x){min(max(x), summary(x)[5] + 1.5 * (summary(x)[5]-summary(x)[2]))}))))
	
	## ploting parameters
	par(oma = c(0, 2, 2.5, 0));
	par(mar = c(0.5,2.1,0.5,0));
	par(mfrow=c(5,1));
	
	## loop through time points
	for(t in 1:length(times)) {
		time = times[t];
		time.l = times.l[t];
		
		## set boxplot data, labels, and colors
		box.dta= dta[,mta[mta$time==time,'samid']];
		box.lab= mta[mta$time==time,'subid'];
		box.col= mta[mta$time==time,'spctc'];
		
		## Plot boxplot and add labels + axis
		boxplot(box.dta,labels='',ylim=ylim,outline=F,col=box.col,axes=F,las=2,cex.axis=0.9);
		text(1:length(box.lab),rep(ylim[1]-diff(ylim)*0.05,length(box.lab)),labels=box.lab,xpd=T,cex=0.85,las=1)
		axis(2,ylim=c(min(box.dta,max(box.dta))),cex.axis=0.9);
		mtext(side=3,time.l,cex=0.9,line=-0.7);
		mtext(side=2,expression('log'[2]*' cpm'),cex=0.9,line=2.1);
		
		## for every 5 timepoints, make a legend
		if (t %% 5 == 1) {
			legend((length(box.lab)/2+0.5),ylim[2]*1.4, xjust=0.5,inset=c(0,0),legend=spcLabls,fill=unique(mta$spctc),cex=0.95,horiz=T,xpd=NA)
		}
	}
}