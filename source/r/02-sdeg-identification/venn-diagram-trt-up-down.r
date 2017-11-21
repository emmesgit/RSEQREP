############################################################################################
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
# Program:  venn-diagram-trt-up-down.r
# Version:  RSEQREP 1.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Generate Venn diagrams showing overlap in treatment groups
# Input:    analysis/glm/<scp>_<trt>_<time>_glm_sig.tab.gz
# Output:  	N/A
############################################################################################

source('../r/init-analysis.r')

## only generate plots if there are between 2 and 5 treatment groups
if (length(trtFlags) > 1 & length(trtFlags) <= 5) {
	
	## set directories
	in.dir.glm = paste(res.dir,'glm',sep='/');
	
	## loop through specimen types
	for(s in 1:length(spcFlags)) {
		spcFlag = spcFlags[s];
		spcLabl = spcLabls[s];
	
		## Graphing parameters
		par(mar=c(0.2, 0.2, 2, 0.2),oma=c(0,0,0,0),mfrow=c(5,2),lheight=0.65);
		plot.num = 0
		
		## get unique list of significant genes across time points
		gen.sig.time = c()
		
		## loop through post treatment timepoints
		for(t in 1:length(postb.times)) {
			time = postb.times[t];
			timel = postb.times.l[t];
			
			## get unique list of significant genes
			gen.sig = c()
			
			## loop through treatment groups
			for (v in 1:length(trtFlags)) {
				trtFlag = trtFlags[v];
				trtLabl = trtLabls[v];
				
				in.file = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz',sep='');
				if (file.exists(in.file)) {
					sdeg.res = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T);
					gen.sig = unique(c(gen.sig,rownames(sdeg.res)));
					gen.sig.time = unique(c(gen.sig.time,rownames(sdeg.res)));
				}
			}
			
			## generate binary matrix
			matrx = matrix(rep(0,(length(trtFlags)*length(gen.sig))),ncol=length(trtFlags))
			rownames(matrx) = gen.sig
			## loop through treatment groups
			for (v in 1:length(trtFlags)) {
				trtFlag = trtFlags[v];
				trtLabl = trtLabls[v];
				
				in.file = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz',sep='');
				if (file.exists(in.file)) {
					sdeg.res = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T);
					matrx[rownames(sdeg.res[sdeg.res$logFC>0,]),v] = 1
					matrx[rownames(sdeg.res[sdeg.res$logFC<0,]),v] = -1
				}
			}
			
			## Print significant genes diagram
			vennDiagramMod(abs(matrx), names=trtLabls,lwd=2,cex=1.4,include='both',counts.col=c('black'),show.include=T)
			title(main=paste(spcLabl,'\n',timel,sep=''),line=-0.5,xpd=T,cex.main=1.5)	
			## Print Up and down regulated diagram
			vennDiagramMod(matrx, names=trtLabls,lwd=2,cex=1.4,include=c("up","down"),counts.col=lfc.cols[c(3,1)],show.include=T)
			title(main=paste(spcLabl,'\n',timel,sep=''),line=-0.5,xpd=T,cex.main=1.5)
			plot.num = plot.num + 2
		}
		
		## all post_treatment time plot
		## generate binary matrix
		matrx = matrix(rep(0,(length(trtFlags)*length(gen.sig.time))),ncol=length(trtFlags))
		rownames(matrx) = gen.sig.time
		
		## If there is more than one timepoint
		if(length(postb.times)>1) {
			## loop through treatment groups
			for (v in 1:length(trtFlags)) {
				trtFlag = trtFlags[v];
				trtLabl = trtLabls[v];
				
				gen.sig = c()
				## loop through post treatment timepoints
				for(t in 1:length(postb.times)) {
					time = postb.times[t];
					timel = postb.times.l[t];
					
					in.file = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz',sep='');
					if (file.exists(in.file)) {
						sdeg.spc.trt.time.sig = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T);
						gen.sig = unique(c(gen.sig,rownames(sdeg.spc.trt.time.sig)))
					}
				}
				matrx[match(gen.sig,rownames(matrx)),v] = 1
			}
			
			## Print significant genes diagram
			vennDiagramMod(matrx, names=trtLabls,lwd=2,cex=1.4,include='both',counts.col=c('black'),show.include=T)
			title(main=paste(spcLabl,'\n',postb.times.l[length(postb.times.l)],sep=''),line=-0.5,xpd=T,cex.main=1.2)
			plot.num = plot.num + 1
		}
		
		## Add blank plots if there are not at least 4.
		if ((plot.num %% 10) == 1) {
			plot.new()
			points(1,0,col='grey99')
		}
	}
}