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
# Program:  hcl-plots.r
# Version:  RSEQREP 2.1.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  generate hierarchical clustering plots based on per_subject data for outlier analysis
# Input:    analysis/dist/<spc>_alltp_euc_dist_tmm_normalized_<std>.RData
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## set directories
in.dir.dist = paste(res.dir,'dist',sep='/');

## update specimen flag to inclue a variable across all
## types if there is more than 1 specimen type
if (length(spcFlags)>1) {
	spcFlags = c('all',spcFlags)
	spcLabls = c('All Specimen Types',spcLabls)
}

## loop through specimen types
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];
	
	## plotting parameters
	par(mfrow=c(1,2),oma=c(0,0,4,1));
	
	## subset metadata to specimen type
	if (spcFlag=='all') {
		mta.spc = mta
	} else {
		mta.spc = mta[mta$spct==spcFlag,]
		mta.spc = mta.spc[!mta.spc$samid %in% outliers,]
	}
	## loop through standardization types
	for(d in 1:length(stdFlags)) {
		stdFlag = stdFlags[d];
		stdLabl = stdLabls[d];
		
		## import data
		filename = paste(in.dir.dist,'/',spcFlag,'_alltp_euc_dist_tmm_normalized_',stdFlag,'.RData',sep='');
		load(file=filename);
		
		## construct HCL object and print as dendrogram
		hcl.res=hclust(dist.res,method=cluster.method);
		hcl.ids = attr(dist.res, "labels")[hcl.res$order]
		hcl.dnd=as.dendrogram(hcl.res);
		plot(colDendo(hcl.dnd,mta.spc$spctc[match(attr(dist.res, "labels"),mta.spc$samid)][hcl.res$order],hcl.ids),axes=T,horiz=T,xlab='Euclidean Distance',xlim=c(max(hcl.res$height),max(hcl.res$height)*-.1));
		if(d==1) {
			legend(max(hcl.res$height)*.75,1.12*nrow(mta.spc),legend=unique(mta.spc$spctl),fill=unique(mta.spc$spctc),cex=0.8,border = FALSE,xpd=TRUE,horiz=T);
		}
		mtext(stdLabl,cex=0.8);	
		
		## Mark outliers
		outs.idx = which(hcl.ids %in% unlist(outliers))
		plot.edge = par()$usr[2]-(abs(par()$usr[1])+abs(par()$usr[2]))/par()$fin[1]*(par()$mai[4]*1.3)
		points(rep(plot.edge,length(outs.idx)),outs.idx,pch=1.5,col='blue',xpd=T)
	}
}