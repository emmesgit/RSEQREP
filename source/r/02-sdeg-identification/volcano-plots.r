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
# Program:  volcano-plots.r
# Version:  RSEQREP 1.1.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Generate volcano plots. Color scheme: red=Significant, black=fold change 
#				outside of +_1.5, dark grey=non significant inside of +_1.5 fold change
# Input:    analysis/glm/<scp>_<trt>_<time>_glm_<all|sig>.tab.gz
# Output:  	N/A
#############################################################################################################

# specify source
source('../r/init-analysis.r')
options(scipen=999)

# specify in directory
in.dir.glm = paste(res.dir,'/glm',sep='')

#####################################
#
## Volcano Plots
#
#####################################

## loop through treatment groups;
for (v in 1:length(trtFlags)) {
	trtFlag = trtFlags[v];
	trtLabl = trtLabls[v];
	
	## loop through specimen types
	for(s in 1:length(spcFlags)) {
		spcFlag = spcFlags[s];
		spcLabl = spcLabls[s];
		
		# Create a grid to plot Volcano plots
		par(mgp=c(1.4,0.4,0),mar=c(2.5,2.8,2,1.2),oma=c(0,0,0,0),mfrow=c(5,2))
		
		## loop through post treatment timepoints
		for(t in 1:length(postb.times)) {
			time = postb.times[t];
			timel = postb.times.l[t];
			
			# open Sig file
			in.file.sig  = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz',sep='')
			if (file.exists(in.file.sig)) {
				dta.sig = read.table(in.file.sig,sep='\t',header=T,row.names=1)
			} else {
				dta.sig = matrix()
			}
			
			# open all file
			dta = read.table(paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_all.tab.gz',sep=''),sep='\t',header=T,row.names=1)
				
			# assign colors
			dta$color = "grey50"
			dta$color[abs(dta$logFC)>=log2(glm.sdeg.fold)] = "black"
			dta$color[which(rownames(dta) %in% rownames(dta.sig))] = "red"
			
			# Plot Volcano plot
			plot(dta$logFC,-log10(dta$PValue),main=paste(spcLabl,', ',trtLabl,', ',timel, sep = ''),cex.main=1.1,
					xlim=c(floor(min(dta$logFC)), ceiling(max(dta$logFC))), ylim=c(0, ceiling(max(-log10(dta$PValue)))),xlab=expression('log'[2]*' fold change'), 
					ylab=expression('-log'[10]*'(p value)'), col=dta$color,cex=.4,cex.axis=0.8,cex.lab=1)
			abline(v=log2(glm.sdeg.fold),lty=2,col='blue');
			abline(v=log2(1/glm.sdeg.fold),lty=2,col='blue');
				
		}
	}
}