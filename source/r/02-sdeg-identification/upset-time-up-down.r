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
# Program:  upset-time-up-down.r
# Version:  RSEQREP 1.1.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Generate Upset plots showing overlap in time points
#			Each plot is added to a list and the list is printed at the end of the program
# Input:    analysis/glm/<scp>_<trt>_<time>_glm_sig.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## only generate plots if there is 3 or more timepoints
if (length(postb.times) >= 2) {
	
	## set directories
	in.dir.glm = paste(res.dir,'glm',sep='/');
	
	## loop through specimen types
	for(s in 1:length(spcFlags)) {
		spcFlag = spcFlags[s];
		spcLabl = spcLabls[s];
		
		## initialize list to hold grobs/viewports
		vps = list()
		
		##  get unique list of genes across timepoints
		gen.sig.time = c()
		
		sig.trt = c()
				
		## loop through treatment groups;
		for (v in 1:length(trtFlags)) {
			trtFlag = trtFlags[v];
			trtLabl = trtLabls[v];
			
			##  get unique list of genes per specimen type
			gen.sig = c()
			
			## loop through post treatment timepoints
			for(t in 1:length(postb.times)) {
				time = postb.times[t];
				timel = postb.times.l[t];
				
				in.file = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz',sep='');
				if (file.exists(in.file)) {
					sdeg.res = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T);
					gen.sig = unique(c(gen.sig,rownames(sdeg.res)));
					gen.sig.time = unique(c(gen.sig.time,rownames(sdeg.res)));
				}
			}
			
		 	## generate binary matrix
		 	matrx = as.data.frame(matrix(rep(0,(length(postb.times)*length(gen.sig))),ncol=length(postb.times)),stringsAsFactors=F)
		 	rownames(matrx) = gen.sig
			## loop through post treatment timepoints
			for(t in 1:length(postb.times)) {
				time = postb.times[t];
				timel = postb.times.l[t];
				
				in.file = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz',sep='');
				if (file.exists(in.file)) {
		 			sdeg.res = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T);
		 			matrx[rownames(sdeg.res[sdeg.res$logFC>0,]),t] = 1
		 			matrx[rownames(sdeg.res[sdeg.res$logFC<0,]),t] = -1
				}
		 	}
		 	
		 	## Print significant genes diagram
			colnames(matrx) = paste(postb.times.l[-length(postb.times.l)],' (',colSums(abs(matrx)),')',sep='')
			matrx = matrx[,ncol(matrx):1]
			if (length(which(colSums(abs(matrx))>0))>1) {
				matrx.all = matrx[,colSums(matrx)!=0]
				vps = append(vps,list(grid.grabExpr(c(
						upset(as.data.frame(abs(matrx.all)),mainbar.y.label='SDEG Intersection Size',sets.x.label='SDEG Set Size',mainbar.y.max=ceiling(max(table(apply(abs(matrx),1,paste,collapse=';')))*1.2),
									matrix.color=lfc.cols[2],main.bar.color=lfc.cols[2],sets.bar.color=lfc.cols[2],sets=colnames(matrx.all),keep.order=T,order.by="freq"),
								grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "black")),
								grid.text(paste(trtLabl,', ',spcLabl,'\n(Up- and down-regulated DE genes)',sep=''),x = unit(0.68, "npc"), y = unit(0.96, "npc"),gp=gpar(font=2,cex=1),just='top'),
								grid.edit('arrange',name='arrange2')))))
				sig.trt = c(sig.trt,trtLabl)
				
				## Print Up regulated diagram
				matrx.up = as.data.frame(apply(matrx,2,function(x)(as.numeric(gsub('-1','0',x)))),stringsAsFactors=F)
				matrx.up = matrx.up[!rowSums(matrx.up)==0,]
				colnames(matrx.up) = paste(rev(postb.times.l[-length(postb.times.l)]),' (',colSums(abs(matrx.up)),')',sep='')
				if (length(which(colSums(matrx.up)>0))>1) {
					matrx.up = matrx.up[,colSums(matrx.up)!=0]
					vps = append(vps,list(grid.grabExpr(c(
							upset(matrx.up,mainbar.y.label='SDEG Intersection Size',sets.x.label='SDEG Set Size',mainbar.y.max=ceiling(max(table(apply(abs(matrx.up),1,paste,collapse=';')))*1.2),
										matrix.color=lfc.cols[3],main.bar.color=lfc.cols[3],sets.bar.color=lfc.cols[3],sets=colnames(matrx.up),keep.order=T,order.by="freq"),
									grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "black")),
									grid.text(paste(trtLabl,', ',spcLabl,'\n(Up-regulated DE genes)',sep=''),x = unit(0.68, "npc"), y = unit(0.96, "npc"),gp=gpar(font=2,cex=1),just='top'),
									grid.edit('arrange',name='arrange2')))))
				} else {
					vps = append(vps,list(grid.grabExpr(grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "white")))))
				}
				
				## Print Down regulated diagram
				matrx.down = as.data.frame(apply(matrx,2,function(x)(as.numeric(gsub('^1$','0',x)))),stringsAsFactors=F)
				matrx.down = abs(matrx.down[!rowSums(matrx.down)==0,])
				colnames(matrx.down) = paste(rev(postb.times.l[-length(postb.times.l)]),' (',colSums(abs(matrx.down)),')',sep='')
				if (length(which(colSums(matrx.down)>0))>1) {
					matrx.down = matrx.down[,colSums(matrx.down)!=0]
					vps = append(vps,list(grid.grabExpr(c(
							upset(matrx.down,mainbar.y.label='SDEG Intersection Size',sets.x.label='SDEG Set Size',mainbar.y.max=ceiling(max(table(apply(abs(matrx.down),1,paste,collapse=';')))*1.2),
										matrix.color=lfc.cols[1],main.bar.color=lfc.cols[1],sets.bar.color=lfc.cols[1],sets=colnames(matrx.down),keep.order=T,order.by="freq"),
									grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "black")),
									grid.text(paste(trtLabl,', ',spcLabl,'\n(Down-regulated DE genes)',sep=''),x = unit(0.68, "npc"), y = unit(0.96, "npc"),gp=gpar(font=2,cex=1),just='top'),
									grid.edit('arrange',name='arrange2')))))
				} else {
					vps = append(vps,list(grid.grabExpr(grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "white")))))
				}
			}
		
			## if there is more than one treatment group
			if (length(trtFlags)>1) {
				
				## all post_treatment time plot
				## generate binary matrix
				matrx = as.data.frame(matrix(rep(0,(length(postb.times)*length(gen.sig.time))),ncol=length(postb.times)),stringsAsFactors=F)
				rownames(matrx) = gen.sig.time
				
				## loop through post treatment timepoints
				for(t in 1:length(postb.times)) {
					time = postb.times[t];
					timel = postb.times.l[t];
					
					gen.sig = c()
					
					## loop through treatment groups
					for (v in 1:length(trtFlags)) {
						trtFlag = trtFlags[v];
						trtLabl = trtLabls[v];
							
						in.file = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz',sep='');
						if (file.exists(in.file)) {
							sdeg.spc.trt.time.sig = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T);
							gen.sig = unique(c(gen.sig,rownames(sdeg.spc.trt.time.sig)))
						}
					}
					matrx[match(gen.sig,rownames(matrx)),t] = 1
				}
				
				## Print significant genes diagram
				colnames(matrx) = paste(postb.times.l[-length(postb.times.l)],' (',colSums(matrx),')',sep='')
				matrx = matrx[,colSums(matrx)!=0]
				matrx = matrx[,ncol(matrx):1]
				if (length(which(colSums(matrx)>0))>1) {
					vps = append(vps,list(grid.grabExpr(c(
						upset(as.data.frame(matrx),nsets=ncol(matrx),mainbar.y.label='SDEG Intersection Size',sets.x.label='SDEG Set Size',mainbar.y.max=ceiling(max(table(apply(abs(matrx),1,paste,collapse=';')))*1.2),
									matrix.color=lfc.cols[2],main.bar.color=lfc.cols[2],sets.bar.color=lfc.cols[2],sets=colnames(matrx),keep.order=T,order.by="freq"),
								grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "black")),
								grid.text(paste('All Treatment Groups\n',spcLabl,sep=''),x = unit(0.68, "npc"), y = unit(0.96, "npc"),gp=gpar(font=2,cex=1),just='top'),
								grid.edit('arrange',name='arrange2')))))
					sig.trt = c(sig.trt,'All treatment groups')
				}
			}	
			## plot if there are one or more plots
			if (length(vps)>0) {
				for (i in 1:ceiling(length(vps) / 3)) {
					grid.newpage()
					do.call("grid.arrange", c(vps[((i*3)-2):min((i*3),length(vps))], ncol=3,nrow=1,newpage=F))
					upsetLablsTimeSdeg = c(upsetLablsTimeSdeg,paste(spcLabl,sig.trt,sep=', '))
				}
			}
		}
	}
}