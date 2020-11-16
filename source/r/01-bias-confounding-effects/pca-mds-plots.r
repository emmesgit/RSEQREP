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
# Program:  pca-mds-plots.r
# Version:  RSEQREP 2.0.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Group level principal component analysis (PCA) and multi dimentional 
#			scaling (MDS) plots based on Euclidean and Spearman distances.
# Input:    analysis/pca/<spc>_alltp_pca_tmm_normalized_std.RData
#			analysis/dist/<spc>_alltp_euc_dist_tmm_normalized_std.RData
#			analysis/dist/<spc>_alltp_spm_dist_tmm_normalized_<std>.RData
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## set directories
in.dir.pca  = paste(res.dir,'pca',sep='/');
in.dir.dist = paste(res.dir,'dist',sep='/');

## update specimen flag to inclue a variable across all
## types if there is more than 1 specimen type
if (length(spcFlags)>1) {
	spcFlags = c('all',spcFlags)
	spcLabls = c('All Specimen Types',spcLabls)
}

## set variables
visLabl = '';
alpha = 0.05;

## loop through specimen types
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];
	spcLabl = spcLabls[s];
	
	## subset metadata to specimen type
	if (spcFlag=='all') {
		mta.spc = mta
	} else {
		mta.spc = mta[mta$spct==spcFlag,]
	}
	
	## plotting parameters
	par(mfrow=c(3,1)); 
	par(mar=c(3.5,3.5,3.5,0.2),oma=c(0,0,1.5,0));
	
	## PCA, Euc Dist MDS, and Spearman MDS
	plotPcaByCell(mta,in.dir.pca,spcFlag,spcLabl,'alltp',visLabl);
	plotMdsByCell(mta,in.dir.dist,type='euc',spcFlag,spcLabl,'alltp',visLabl,text='Euclidean Distance Standardized Variables');
	plotMdsByCell(mta,in.dir.dist,type='spm',spcFlag,spcLabl,'alltp',visLabl,text='1-Spearman Cor. Distance Standardized Variables');
	
	## add legend
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend('top',legend=unique(mta.spc$spctl),fill=unique(mta.spc$spctc),horiz=T,cex=1.1,border = FALSE,xpd=TRUE);
}