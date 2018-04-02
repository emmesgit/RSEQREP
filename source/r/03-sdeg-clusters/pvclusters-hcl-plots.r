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
# Program:  pvclusters-hcl-plots.r
# Version:  RSEQREP 1.1.2
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Plot significant clusters using dendrograms
# Input:    data/annot/filtered_gene_annotations.tab
#			analysis/glm/<spc>_<trt>_<time>_glm_sig.tab
#			analysis/pvclust/<spc>_<trt>_<time>_pvclust_sig.tab.gz
#			analysis/pvclust/<spc>_<trt>_<time>_pvclust.RData
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## set directories
in.dir.pvc = paste(res.dir,'pvclust',sep='/');
in.dir.glm = paste(res.dir,'glm',sep='/');

## load gene annotations
gen = read.csv(paste(dta.dir,'/annot/filtered_gene_annotations.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T);

## loop through specimen types
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];
	spcLabl = spcLabls[s];
		
	## loop through post treatment timepoiunts
	## including all post treatment days
	for(t in 1:length(postb.times.m)) {
		postb.time = postb.times.m[t];
		postb.timel = postb.times.l[t];
		
		
		## GET SIGNIFICANT GENES BY TIMEPOINT
		gen.sig = c();
		for (v in 1:length(trtFlags)) {
			trtFlag = trtFlags[v];
			glm.infile = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_',postb.time,'_glm_sig.tab.gz',sep='');
			if(file.exists(glm.infile)) {
				sdeg.tbl.time = read.table(glm.infile,sep='\t',stringsAsFactors=F,header=T);
				gen.sig = unique(c(gen.sig,rownames(sdeg.tbl.time)));
			}
			if(postb.time =='posttp' & length(postb.times)>1){	
				for(tx in 1:(length(postb.times.m)-1)) {
					timex = postb.times.m[tx];
					glm.infile = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_',timex,'_glm_sig.tab.gz',sep='');
					if(file.exists(glm.infile)) {
						sdeg.cel.time = read.table(glm.infile,sep='\t',stringsAsFactors=F,header=T);
						gen.sig = unique(c(gen.sig,rownames(sdeg.cel.time)));
					}
				}
			}
		}
		
		## if there are any signifcant genes
		if(length(gen.sig) > 0) {	
			
			## specify infiles
			pvc.infile = paste(in.dir.pvc,'/',spcFlag,'_',postb.time,'_pvclust.RData',sep='');	
			pvc.infile.sig = paste(in.dir.pvc,'/',spcFlag,'_',postb.time,'_pvclust_sig.tab.gz',sep='');	
			
			## if there is a significant file, generate dendrogram
			if(file.exists(pvc.infile.sig)) {
				load(file=pvc.infile);
			
				## scale based on gene size
				factor = max(min(280,length(res.pvc$hclust$labels))/40,1);
				
				## plotting parameters
				par(mar=rep(0,4));
				all.cex= 0.31/factor;
				cex.lab= 2*factor;	
				cex.pv=1.3*(factor/2);
				cex.axis=0.8*factor;
				linew = 0.8/factor;
				lwd=0.8/factor
				par(cex.main=0.8*3*factor)
				par(mar=c(80*factor,2*factor+2,2*factor,1))
				par(cex=all.cex);
				
				## Update labels
				labs = gen[match(res.pvc$hclust$labels,gen$gene_id),'gene_name_lab'];
				for(l in 1:length(res.pvc$hclust$labels)){
					if(res.pvc$hclust$labels[l] %in% gen.sig) {
						labs[l]=paste(labs[l],'*',sep='');
					}
				}	
				
				## plot dendrograms
				plot(res.pvc,main=paste('Gene cluster dendrogram with bootstrap probabilities (',spcLabl,', ',postb.timel,')',sep=''),sub='',xlab='',ylab='',
						labels = labs,cex.axis=cex.axis,cex.lab=cex.lab,cex.pv=cex.pv,print.num=F,col.pv=c('black','green','blue'),lwd = linew);
				pvrect(res.pvc,alpha=pvclust.pval,pv="au",type="geq",lwd=lwd,lty=1,border='blue',max.dist=pvclust.max.dist)
				box();
				
				## caption add-on for KnitR code
				DendrLabls = c(DendrLabls,paste(spcLabl,', ',postb.timel,sep=''))
			}
		}
	}
}