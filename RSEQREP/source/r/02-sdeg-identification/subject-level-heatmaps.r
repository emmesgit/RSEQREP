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
# Program:  subject-level-heatmaps.r
# Version:  RSEQREP 2.2.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Generate subject level heatmaps for each spc and time combination across all trt groups
# Input:    data/annot/filtered_gene_annotations.tab
#			analysis/glm/<scp>_<trt>_<time>_glm_sig.tab.gz
#			analysis/lcpm_fc/<spc>_posttp_lcpm_fold_change_tmm_normalized_filtered.tab.gz
#			analysis/pvclust/<spc>_<trt>_<time>_pvclust_sig.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## set directories
in.dir.glm	= paste(res.dir,'glm',sep='/');
in.dir.lfc	= paste(res.dir,'lcpm_fc',sep='/');
in.dir.pvc  = paste(res.dir,'pvclust',sep='/');

## load gene annotations
gen = read.csv(paste(dta.dir,'/annot/filtered_gene_annotations.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T);

## define colors
sdeg.fold = 1.2
colors = c(seq(lfc.col.range[1],log2(1/sdeg.fold),length=100),seq(log2(1/sdeg.fold)+0.00000001,log2(sdeg.fold)-0.00000001,length=100),seq(log2(sdeg.fold),lfc.col.range[length(lfc.col.range)],length=100))
heatcolor = colorRampPalette(lfc.cols)(n =299)

## define color legend
colLevels = rep(NA,length(lfc.col.range))

for(c in 1:length(lfc.col.range)) {
	pos = floor(mean(which(round(colors ,1)== lfc.col.range[c])))
	pos = min(299,pos);
	colLevels[c]=heatcolor[pos]
}

## loop through specimen types
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];
	spcLabl = spcLabls[s];
	
	## import lfc matrix
	dta.lfc.infile = paste(in.dir.lfc,'/',spcFlag,'_posttp_lcpm_fold_change_tmm_normalized_filtered.tab.gz',sep='');
	dta.lfc = read.csv(dta.lfc.infile,header=T,stringsAsFactors=F,check.names = FALSE,sep='\t');
	
	## loop through post treatment timepoints
	for(t in 1:length(postb.times)) {
		time = postb.times[t];
		timel = postb.times.l[t];
		
		## subset metadata
		mta.time.spc = mta[mta$time == time & mta$spct == spcFlag & !mta$samid %in% outliers,]
		
		## collect significant genes
		gen.sig = c();
		
		## loop through treatment groups;
		for (v in 1:length(trtFlags)) {
			trtFlag = trtFlags[v];
			trtLabl = trtLabls[v];
			
			in.file = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz',sep='');
			if (file.exists(in.file)) {
				sdeg.res = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T);
				gen.sig = unique(c(gen.sig,rownames(sdeg.res)));
			}
		}
	
		## subset significant genes
		dta.lfc.time = dta.lfc[gen.sig,mta.time.spc$samid]
		
		## Incorporate pvclusters association coloring to rows
		cls.labs = rep('white',length(rownames(dta.lfc.time)));
		names(cls.labs) = rownames(dta.lfc.time);
		in.file.cls = paste(in.dir.pvc,'/',spcFlag,'_tp',time,'_pvclust_sig.tab.gz',sep='');
		if(file.exists(in.file.cls)) {
			dta.cls = read.csv(in.file.cls,header=T,sep='\t',stringsAsFactor=F);
			for(lab.x in 1:length(cls.labs)) {
				lab.gene = names(cls.labs[lab.x]);
				
				if(lab.gene %in% dta.cls$gene_id) {
					cls.labs[lab.x] = dta.cls[dta.cls$gene_id==lab.gene,'cluster_id'];
				}
			}
		}
		cls.labs[cls.labs!='white'] = as.character(as.numeric(factor(cls.labs[cls.labs!='white']))+1);
		cls.lab.col = as.character(cls.labs);
	
		## specify column coloring by treatment group
		col.color = mta.time.spc$trtc[match(colnames(dta.lfc.time),mta.time.spc$samid)]
		
		## Update row and column names
		colnames(dta.lfc.time) = mta.time.spc$subid[match(colnames(dta.lfc.time),mta.time.spc$samid)]
		rownames(dta.lfc.time) = make.unique(gen$gene_name_lab[match(rownames(dta.lfc.time),gen$gene_id)])
		
		## ensure data frame has at least 3 rows
		if (nrow(dta.lfc.time)>=2) {
		
			## plotting parameters
			par(mar=c(0,0,0,0))
			par(cex.main=0.8)
			if (nrow(dta.lfc.time) < 20) {
				cex.row = 1
			} else {
				cex.row = max(round((1/(min(550,nrow(dta.lfc.time))/3.5))*10,2),0.11);
			}
		
			## Title
			main=paste('Gene Response Heatmap (',spcLabl,', Day ',time,' vs. ',b.times.l,')',sep='');
			
			## print heatmap
			heatmap.2(as.matrix(dta.lfc.time),main=main,
					xlab='',ylab='',cex.axis=0.8,trace="none",col=heatcolor, distfun = function(x) heatmap.dist.fun(x,by.row=T),key=F,
					hclustfun = hclust2,cexRow=cex.row,keysize = 1,scale = c("none"),
					cexCol=0.6, density.info=c('none'),margins=c(6,13), 
					dendrogram =c('both'), symm=F,symkey=F,symbreaks=T, breaks=colors,
					labRow=rownames(dta.lfc.time),RowSideColors=cls.lab.col,ColSideColors=mta.time.spc$trtc,lhei=c(1,4),lwid=c(1,4));
			
			legend('topleft',legend=unique(mta$trtl),fill=unique(mta$trtc),cex=0.7,title='Treatment',horiz=F)
			legend("bottom",legend=lfc.col.rangel,fil=colLevels,pch=NA,y.intersp=0.5,cex=0.6,title=expression('Log'[2]*' Fold Change From Pre Treatment'),title.adj = c(0.5,-0.6), horiz=T,bty = "n")
			
			## caption add-on for KnitR code
			HeatmapLablSdeg = c(HeatmapLablSdeg,paste(spcLabl,timel,sep=', '))
			}
		}
	}
