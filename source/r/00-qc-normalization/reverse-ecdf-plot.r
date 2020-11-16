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
# Program:  reverse-ecdf-plot.r
# Version:  RSEQREP 2.0.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  reverse cumulative distribution function plots summarizing the percentage of genes that 
#			exceeded a certain logCPM cut off for the number of samples specified.  A line will be generated
#			for each specimen type on the plot
# Input:    analysis/lcpm/<spc>_alld_lcpm_tmm_normalized_unfiltered.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')
options(scipen=999)

## set directories
in.dir.lcpm = paste(res.dir,'lcpm',sep='/');

## loop through lab types
for (c in 1:length(spcFlags)) {
	spcFlag = spcFlags[c];
	
	# subset metadta to lab and set post_treatment time
	mta.spc = mta[mta$spct == spcFlag & !mta$samid %in% outliers,]
	pb.time = unique(mta.spc$time[mta.spc$time>0])
	
	## import tmm normalized log counts per million
	in.file =paste(res.dir,'/lcpm/all_alltp_lcpm_tmm_normalized_unfiltered.tab.gz',sep='');
	dta.lab = read.csv(in.file,header=T,stringsAsFactors=F,sep='\t',row.names=1,check.names = FALSE);
	dta.lab = dta.lab[,mta.spc$samid]
	
	## generate a vector of values ranging from min lcpm to max lcpm
	lcpm.range = floor(min(dta.lab)):ceiling(max(dta.lab))
	
	if (c==1) {
		## initalize plot
		plot(1,1,col='white',ylab='genes [%]',xlab=expression('log'[2]*' counts per million cutoff'),
				xaxt = "n",ylim=c(0,100),xlim=c(min(lcpm.range),max(lcpm.range)))
		axis(side=1,labels=seq(min(lcpm.range),max(lcpm.range)),at=seq(min(lcpm.range),max(lcpm.range)))
		## create a rectangle representing minimum and maximum number of genes allowed in final set
		rect.size = c(flt.gene.range[1]/nrow(dta.lab),flt.gene.range[2]/nrow(dta.lab)) *100
		rect(min(lcpm.range),rect.size[1],max(lcpm.range),rect.size[2],col='grey90')
		## create a vertical line at the actual value
		abline(v=flt.lcpm.cut,col='black',lwd=2)
	}
	
	## reverse cummulative distribution function
	if (filterFun=='filterLowExpressedGenesMax') {
		rcdf = ecdf(apply(dta.lab,1,max));
	} else if (filterFun=='filterLowExpressedGenesAvg'){
		rcdf = ecdf(apply(dta.lab,1,mean));
	}
	
	lines(y=100-rcdf(lcpm.range)*100,x=lcpm.range,type='l',col=unique(mta.spc$spctc),lwd=2);
}

## Get # genes per specimen type -- record in legend
legend.text = c()

## loop through specimen types
for(c in 1:length(spcFlags)) {
	spcFlag = spcFlags[c];
	spcLabl = spcLabls[c];
	
	## only pursue 'posttp'
	for(d in length(postb.times.m)) {
		time = postb.times.m[d];
		
		## specify infile name
		in.file= paste(in.dir.lcpm,'/',spcFlag,'_',time,'_analysis_gene_set.tab.gz',sep='')
		
		## if file exists, open it and update matrix
		if (file.exists(in.file)) {
			dta = read.table(paste(in.dir.lcpm,'/',spcFlag,'_',time,'_analysis_gene_set.tab.gz',sep=''));
			legend.text = c(legend.text,paste(spcLabl,' ',nrow(dta),' of ',nrow(dta.lab),' (',round(nrow(dta)/nrow(dta.lab)*100,0),'%) selected',sep=''))
		}
	}
}


## add legend
legend('topright',legend=legend.text,bty='n',col=unique(mta$spctc),lwd=2,cex=0.9)
