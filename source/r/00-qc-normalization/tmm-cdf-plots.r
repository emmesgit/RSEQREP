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
# To cite this software, please reference doi:10.12688/f1000research.13049.1
#
# Program:  tmm-cdf-plots.r 
# Version:  RSEQREP 1.1.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Plot cpm distributions pre and post normalization data for each 
#				specimen type colored by treatment.
# Input:    analysis/lcpm/all_alld_lcpm_<nrm>_normalized_filtered.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## set directories
in.dir.lcpm = paste(res.dir,'lcpm',sep='/');

## ploting parameters
par(oma = c(0, 0, 3.5, 0));
par(mar = c(3.2,3.5,1,0.5));
par(mfrow=c(min(length(times),5),2));

## loop through time points
for(t in 1:length(times)) {
	time = times[t];
	time.l = times.l[t];
		
	## loop through normalization types
	for(n in 1:length(nrmFlags)) {
		nrmFlag = nrmFlags[n];
		nrmLabl= nrmLabls[n];
			
		in.file = paste(in.dir.lcpm,'/all_alltp_lcpm_',nrmFlag,'_normalized_filtered.tab.gz',sep='');
		dta = read.csv(in.file,header=T,stringsAsFactors=F,check.names = FALSE,sep='\t');
			
		## get range of x values
		x = seq(floor(min(dta)),ceiling(max(dta)),1.1);
		
		## set ecdf data, labels, and colors
		ecdf.dta= dta[,as.character(mta[mta$time==time,'samid'])];
		ecdf.lab= mta[mta$time==time,'subid'];
		ecdf.col= mta[mta$time==time,'spctc'];
		
		for(c in 1:ncol(ecdf.dta)) {
			efunc = ecdf(ecdf.dta[,c]);
			if( c==1) {
				plot(NA,NA,cex=0.3,axes=F,xaxt="n",ylim=c(0,1),xlim=range(x),xlab='',ylab='');
				axis(1,cex.axis=0.9)
				axis(2,cex.axis=0.9,at=seq(0,1,0.2),las=2)
				mtext(side = 1, expression('x=log'[2]*' cpm'), line = 2.2, cex=0.75);
				mtext(side = 2, "F(x)", line = 2.3, cex=0.75)
			}				
				
			lines(x,efunc(x),col=ecdf.col[c],lwd=0.5)
			box();
		}
		
		## add plot title
		mtext(paste(time.l,' (',nrmLabl,')',sep=''),cex=0.9);
		
		## for every 5 timepoints, make a legend
		if (n==1 & t %% 5 == 1) {
			legend(max(x)*1.2,1.65, xjust=0.5,inset=c(0,0),legend=spcLabls,fill=unique(mta$spctc),cex=0.95,horiz=T,xpd=NA)
		}
	}	
}