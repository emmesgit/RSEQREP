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
# Program:  pvclusters-time-trend.r
# Version:  RSEQREP 1.1.3
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Visualize cluster log fold change from treatment trends
#			for significant clusters
# Input:    data/annot/filtered_gene_annotations.tab
#			analysis/glm/<spc>_<trt>_<time>_glm_sig.tab
#			analysis/pvclust/<spc>_<trt>_<time>_pvclust_sig.tab.gz
#			analysis/pvclust/<spc>_<trt>_<time>_pvclust.RData
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

# only plot if there are more than 1 post treatment time points
if(length(postb.times) > 1) {
	
	## set directories
	in.dir.glm  = paste(res.dir,'glm',sep='/');
	in.dir.lfc  = paste(res.dir,'lcpm_fc',sep='/');
	in.dir.pvc = paste(res.dir,'pvclust',sep='/');
	
	## load gene annotations
	gen = read.csv(paste(dta.dir,'/annot/filtered_gene_annotations.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T);
	
	## loop through spc types
	for(s in 1:length(spcFlags)) {
		spcFlag = spcFlags[s];
		spcLabl = spcLabls[s];
		
		## specify plotting parameters
		par(mfrow=c(4,3));	
		par(mar=c(2.5,2.5,2.5,0.2));
		plot.num = 0
		
		## get log fold change from treatment; align meta data
		in.file = paste(in.dir.lfc,'/',spcFlag,'_posttp_lcpm_fold_change_tmm_normalized_filtered.tab.gz',sep='');
		dta.lfc = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T);
		mta.lfc = mta[match(colnames(dta.lfc),mta$samid),]
		
		## loop through post treatment timepoiunts
		## including all post treatment days
		for(t in length(postb.times.m)) {
			postb.time = postb.times.m[t];
			postb.timel = postb.times.l[t];
				
			## specify infile
			pvc.infile.sig = paste(in.dir.pvc,'/',spcFlag,'_',postb.time,'_pvclust_sig.tab.gz',sep='');	
			
			## if there are signficant clusters
			if(file.exists(pvc.infile.sig)) {
					
				## import significant clusters
				clr.res = read.csv(pvc.infile.sig,header=t,sep='\t',stringsAsFactors=F)
				clusters = tapply(clr.res$gene_id,clr.res$cluster_id,list)
				
				## Loop through each significant cluster
				cls.sig.num=0;
				first.in.line=1;
				for(cx in 1:length(clusters)) {
					cluster = clusters[[cx]];		
					cls.sig.num=cls.sig.num+1;
					
					## set Time Labels
					if(postb.time == 'posttp') {
						time.lab = 'ALLTP';
					} else {
						time.lab = toupper(postb.time);
					}
					
					## Set cluster IDs and gene anmes
					clusterId = paste(toupper(spcFlag),time.lab,sprintf("%03d", cls.sig.num),' (n=',length(cluster),')',sep='');
					geneNames = gen[match(cluster,gen$gene_id),'gene_name_lab'];
												
					## generate trt based results
					res = t(dta.lfc[rownames(dta.lfc)%in%cluster,]);
					res = cbind(res,mta.lfc[match(rownames(res),mta.lfc$samid),c('trt','subid','time','samid')])
					
					## average lfc across genes
					res.avg   = cbind(gen_mean=apply(res[,cluster],1,mean),res[,c('trt','subid','time','samid')])
					y.avg.trt = tapply(res.avg[,'gen_mean'],list(res$trt,res$time),mean);
										
					## lfc limits
					max.mean = max(sapply(1:length(cluster),function(x){tapply(res[,cluster[x]],list(res$trt,res$time),mean)}))
					min.mean = min(sapply(1:length(cluster),function(x){tapply(res[,cluster[x]],list(res$trt,res$time),mean)}))
					span = max.mean-min.mean
					ylim =  c(min.mean-(span*0.1),max.mean+(span*0.2))
					
					## Loop through each gene in the cluster
					for(genx in 1:length(cluster)) {
						
						## individual gene lfc
						y.trt = tapply(res[,cluster[genx]],list(res$trt,res$time),mean);
									
						## for the first gene in a cluster, do
						if(genx==1) {
							plot(x=postb.times,main=clusterId,cex.main=0.8,xlab='',ylab='', y=y.avg.trt[1,],
										xlim=c(min(postb.times)-1,max(postb.times)*1.6),ylim=ylim,pch=20,col='white',yaxt="n",axes=F,cex=0.8,cex.axis=0.5);
							if(first.in.line==1) {
								first.in.line=0;
								mtext( side = 2, text=expression('mean log'[2]*' fold change'), cex=0.49 ,line=1.5); ## add subtitle
							} 
							if(cx %% 3 ==0) {
								first.in.line=1
							}
							mtext( side = 3, text=paste(spcLabl,', ',postb.timel,sep=''), cex=0.49 ,line=1.5); ## add subtitle
							
							trtCols = unique(mta.lfc$trtc[order(mta.lfc$trt)])
							for (l in 1:nrow(y.avg.trt)) {
								## average over genes
								points(x=postb.times,y=y.avg.trt[l,],col=trtCols[l],pch=20,cex=0.8)
								lines(x=postb.times,y=y.avg.trt[l,],col=trtCols[l],lwd=1.2);
								## per gene
								points(x=postb.times,y=y.trt[l,],col=makeTransparent(trtCols[l],alpha=70),pch=20,cex=0.6)
								lines(x=postb.times,y=y.trt[l,],col=makeTransparent(trtCols[l],alpha=70),cex=0.5,lwd=0.8);
							}
							
							axis(1,at=postb.times, labels=postb.times.l[1:length(postb.times)],cex.axis=0.5,las=2);
							axis(2,cex.axis=0.5) 

							box();
							abline(h=0,col='darkgrey',lty=1,lwd=0.5);
							legend('topleft',legend=unique(mta.lfc$trtl),fill=trtCols,cex=0.55,border = FALSE,bty = "n");
							legend('bottomright',legend=geneNames,title='Genes:',cex=0.53,border = FALSE,bty = "n");
						} else {
							for (l in 1:nrow(y.avg.trt)) {
								## per gene
								points(x=postb.times,y=y.trt[l,],col=makeTransparent(trtCols[l],alpha=70),pch=20,cex=0.6)
								lines(x=postb.times,y=y.trt[l,],col=makeTransparent(trtCols[l],alpha=70),cex=0.5,lwd=0.8);
							}
						}	
					}
					plot.num = plot.num + 1
				} 
			}
		}
		## caption add-on for KnitR code
		if (cx<12) {
			add.on = ''
		} else {
			add.on = paste(1:ceiling(cx/12),' of ',ceiling(cx/12),' ',sep='')
		}
		TrendLabls = c(TrendLabls,paste(add.on,'(',spcLabl,')',sep=''))
		
		## Add blank plots if there are not at least 4.
		if ((plot.num %% 12) > 0 & (plot.num %% 12) < 3) {
			for (i in 1:(3-(plot.num %% 12))) {
				plot.new()
				if (i==(3-(plot.num %% 12))) {
					points(1,0,col='grey99')
				}
			}
		}
	}
}