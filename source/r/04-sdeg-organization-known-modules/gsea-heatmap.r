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
# Program:  gsea-heatmap.r
# Version:  RSEQREP 2.1.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Generate heatmaps for the top 50 enriched pathways based on jaccard index sums across conditions 
# 			enriched in at least 2 conditions (time,trt,spc).
# Input:    analysis/gsea/<spc>_<trt>_<time>_<gene-set-type>_gsea_all.tab.gz
#			analysis/glm/<scp>_<trt>_<time>_glm_sig.tab.gz
#			data/gene_sets/all.ensembl.tab.gz
#			data/annot/filtered_gene_annotations.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r');
options(scipen=999)

## specify gene sets file
sets.infile = paste(dta.dir,'/gene_sets/all.ensembl.tab.gz',sep='')

if (file.exists(sets.infile)) {
	
	## set directories
	in.dir.gsea  = paste(res.dir,'gsea',sep='/');
	in.dir.glm  = paste(res.dir,'glm',sep='/');
	
	## import gene set information
	sets = read.table(sets.infile,sep='\t',stringsAsFactors=F,header=F,quote="",fill=T,row.names=NULL);
	colnames(sets) = c('category_type','gene_id','category_name','category_link')
	
	# load gene annotations
	gen = read.csv(paste(dta.dir,'/annot/filtered_gene_annotations.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T);
	gen.sets = sets[sets$gene_id %in% gen$gene_id,]
	
	## initialize counters
	joined=c();
	countsJoined = c();
	
	## loop through specimen types
	for(s in 1:length(spcFlags)) {
		spcFlag = spcFlags[s];
		spcLabl = spcLabls[s];
		
		## loop through treatment groups;
		for (v in 1:length(trtFlags)) {
			trtFlag = trtFlags[v];
			trtLabl = trtLabls[v];
			
			## loop through post treatment timepoints
			for(t in 1:length(postb.times)) {
				time = postb.times[t];
				timel = postb.times.l[t];
				
				## loop through each gene set category
				for(k in 1:length(geneset.types)) {
					geneset.type = geneset.types[k];
					
					## specify infiles
					in.file = paste(in.dir.gsea,'/',spcFlag,'_',trtFlag,'_tp',time,'_',gsub(' ','_',tolower(geneset.type)),'_gsea_all.tab.gz',sep='');
						
					if (file.exists(in.file)) {
						
						## read gene set enrichment results
						res = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T,quote="",fill=T);
						
						counts= data.frame(
								category_name = res$category_name,
								category_name = res$category_name,
								category_type = rep(geneset.type,nrow(res)),
								sample_id = rep(paste(spcLabl,', ',trtLabl,' (',timel,')',sep=''),nrow(res)),
								count = res$sdeg_genes,
								jidx = res$jidx,
								es_or=-1*log(res$pval_or,10),
								es_ur=-1*log(res$pval_ur,10));
						countsJoined=rbind(countsJoined,counts)
						
						res = res[res$fdr_or <= goseq.fdr.cutoff,]
	
						if(nrow(res)>0) {	
							res$category_type = geneset.type;
							res$spc_type = spcFlag;
							res$trt_type = trtFlag;
							res$time = time;
							res$sample_id = paste(spcLabl,', ',trtLabl,' (',timel,')',sep='');
							joined = rbind(joined,res);
						}
					}
				}
			}
		}
	}
	
	## if no significant results, joined is just an empty vector
	joined = as.data.frame(joined)
	if (nrow(joined) > 0) {
	
		## get gene sets that are enriched in at least 2 condidtions (specimen type, treatment type, or time point)
		cat.gt2 = sqldf('select category_type,category_name,category_name,count(distinct spc_type) as spc_count, count(distinct trt_type) as trt_count,
		count(distinct time) as time_count, count(*) as total_count, group_concat(distinct spc_type) as spc_types, group_concat(distinct trt_type) as trt_types, 
		group_concat(distinct time) as time_count from joined group by category_type,category_name having spc_count>1 or trt_count>1 or time_count>1 
		order by total_count desc');
		
		## loop through each gene set category
		for(k in 1:length(geneset.types)) {
			geneset.type = geneset.types[k];
			
			## get subset of data for the selected gene set category
			cat.gt2.sel = cat.gt2[cat.gt2$category_type==geneset.type,];
			
			## there must be at least 3 entries to plot results on heatmap
			if(nrow(cat.gt2.sel)>=3) {
				
				## assign sample names
				samples = as.vector(unique(countsJoined$sample_id));
				
				## intit matrix to hold sdeg counts
				mtx.cnt = matrix(rep(0,length(unique(cat.gt2.sel$category_name))*length(samples)),ncol=length(samples),nrow=length(unique(cat.gt2.sel$category_name)))
				
				## init matrix to hold binary absence/presence information
				mtx.bin = matrix(rep(0,length(unique(cat.gt2.sel$category_name))*length(samples)),ncol=length(samples),nrow=length(unique(cat.gt2.sel$category_name)))
				
				## init matrix to hold
				mtx.not = matrix(rep(NA,length(unique(cat.gt2.sel$category_name))*length(samples)),ncol=length(samples),nrow=length(unique(cat.gt2.sel$category_name)))
		
				categories = sort(unique(cat.gt2.sel$category_name));
				
				## Loop through each sample
				for(n in 1:length(samples)) {
				
					## limit to specific subject
					sam.counts = countsJoined[ countsJoined$sample_id==samples[n] & countsJoined$category_type==geneset.type,];
					sam.joined = joined[joined$sample_id==samples[n],]
					sig.idx = which(categories %in% joined[joined$sample_id==samples[n] & joined$category_type==geneset.type,'category_name'])
					
					## update matrices
					mtx.not[,n] = sam.counts[match(categories,sam.counts$category_name),'count'];		# SDEG Counts
					mtx.cnt[,n] = sam.counts[match(categories,sam.counts$category_name),'es_or'];		# SDEG-geneset enrichment score
					mtx.not[sig.idx,n] = paste('[',mtx.not[sig.idx,n],']',sep='')						# add brackets if Enriched
					mtx.bin[sig.idx,n] = 1																# is enrighed? 1==T
				}
				
				## update category names for printing
				categories = gsub('_',' ',categories);
				categories[nchar(categories)>60] = paste(substring(categories[nchar(categories)>60],1,60),'...');	
		
				## assign row and column names to result matrices
				rownames(mtx.cnt) = categories;
				rownames(mtx.bin) = categories;
				colnames(mtx.cnt) = samples;
				colnames(mtx.bin) = samples;
					
				## Report only top 50 (Jaccard Index) if there are more than 50 entries
				if(nrow(mtx.cnt)>50) {
					idx = head(order(rowSums(mtx.cnt),decreasing=T),50);
					x = mtx.cnt[idx,];
					mtx.bin = mtx.bin[idx,]
					mtx.not = mtx.not[idx,]
					cex.cnt =0.35;
				} else {
					x = mtx.cnt;
					cex.cnt = 0.6;
				}
				
				## adding to ensure max(x) is finite
				if(!is.finite(max(x))){next}
				
				## assign heatmap colors gsea.col
				heatcolor 		=  colorRampPalette(gsea.cols)(n =20);
				heatcolor[1] 	= '#A8A8A8'; # color 0 values grey
				color.breaks	= c(-max(x)/(length(heatcolor)-1),(round(seq(0,max(x),max(x)/(length(heatcolor)-1)),3)));
				color.breaks[2] = 0.001
				range = round(seq(0,max(x),max(x)/length(heatcolor)),3);
				
				## plotting parameters
				cex.row = round(min((1/(min(550,nrow(x))/1))*14,0.6),2);
				par(oma=c(0,1,0,0));
				par(mar=c(0,0,0,0))
				par(cex.main=0.55)
				par(cex=0.6)
				
				## plot
				heatmap.2(x,
						main=paste('Pathway Enrichment Heatmap\n(',geneset.type,')',sep=''),  
						Rowv=FALSE, 
						Colv=FALSE,
						xlab='',
						ylab='',
						dendrogram=c('row'),
						cex.axis=0.8,
						trace="none",
						col=heatcolor, 
						distfun= function(x) dist(x,method='euclidean'), 
						key=T, 
						keysize=0.01,
						hclustfun = function(x) hclust(x,method=cluster.method),
						cexRow=cex.row,
						scale = c("none"),
						cexCol=0.7, 
						cellnote=mtx.not,
						notecol='black',
						notecex=cex.cnt,
						density.info=c('none'),
						margins=c(12,15), 
						symm=F,
						symkey=F,
						symbreaks=F,
						breaks=color.breaks,
						colsep=1:ncol(x),
						rowsep=1:nrow(x),
						sepcolor='white',
						sepwidth=c(0.01,0.01),
						lmat=matrix(c(4,3,2,1),ncol=2,byrow=T),
						lwid = c(1,4),
						lhei = c(1.3,6),
						labRow=rownames(x),
						cex.main=0.8, 
						srtCol=45,
						srtRow=0);
				
				## caption add-on for KnitR code
				HeatmapLablGsea = c(HeatmapLablGsea,geneset.type)
			}
		}
	}
}