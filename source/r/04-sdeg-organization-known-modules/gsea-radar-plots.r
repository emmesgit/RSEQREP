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
# Program:  gsea-radar-plots.r
# Version:  RSEQREP 2.0.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Generate radar plots for enriched gene sets by gene set group.
# Input:    analysis/gsea/<spc>_<trt>_<time>_<gene-set-type>_gsea_all.tab.gz
#			analysis/glm/<scp>_<trt>_<time>_glm_sig.tab.gz
#			data/gene_sets/all.ensembl.tab.gz
#			data/annot/filtered_gene_annotations.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## specify gene sets file
sets.infile = paste(dta.dir,'/gene_sets/all.ensembl.tab.gz',sep='')

if (file.exists(sets.infile)) {
	
	## set directories
	in.dir.gsea  = paste(res.dir,'gsea',sep='/');
	in.dir.glm  = paste(res.dir,'glm',sep='/');
	
	## import gene set information
	sets = read.table(sets.infile,sep='\t',stringsAsFactors=F,header=F,quote="",fill=T,row.names=NULL);
	colnames(sets) = c('category_type','gene_id','category_name','category_link')
	
	## load gene annotations
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
								category_type = rep(geneset.type,nrow(res)),
								category_name = res$category_name,
								sample_id = rep(paste(spcLabl,', ',trtLabl,' (',timel,')',sep=''),nrow(res)),
								count = res$sdeg_genes,
								spc = spcFlag,
								time = time,
								trt = trtFlag,
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
	if (is.data.frame(joined) > 0) {
		
		## get gene sets that are enriched in at least 2 condidtions (specimen type, treatment type, or time point)
		cat.gt2 = sqldf('select category_type,category_name,category_name,count(distinct spc_type) as spc_count, count(distinct trt_type) as trt_count,
						count(distinct time) as time_count, count(*) as total_count, group_concat(distinct spc_type) as spc_types, group_concat(distinct trt_type) as trt_types, 
						group_concat(distinct time) as time_types from joined group by category_type,category_name having spc_count>1 or trt_count>1 or time_count>1 
						order by total_count desc');
		
		
		## plotting parameters
		par(mar=c(3.2,3.2,3.2,3.2),
				oma=c(0,0,0,0),
				mfrow=c(1,2),
				xpd=NA);
		
		## loop through each gene set category
		for(k in 1:length(geneset.types)) {
			geneset.type = geneset.types[k];
			
			## loop through specimen types
			for(s in 1:length(spcFlags)) {
				spcFlag = spcFlags[s];
				spcLabl = spcLabls[s];
			
				## loop through treatment groups;
				for (v in 1:length(trtFlags)) {
					trtFlag = trtFlags[v];
					trtLabl = trtLabls[v];
		
					## create dataset
					sig.2.types.es = countsJoined[countsJoined$spc==spcFlag & countsJoined$trt==trtFlag & countsJoined$category_type %in% geneset.type & countsJoined$category_name %in% cat.gt2$category_name,]
					sig.2.types.es = sig.2.types.es[order(sig.2.types.es$time,sig.2.types.es$category_name),]
					genesets.es = do.call(rbind.data.frame, tapply(sig.2.types.es$es_or,sig.2.types.es$time,list))
					
					## there must be 3 or more points to make a good plot
					if (ncol(genesets.es) > 2) {
						names = gsub('_',' ',unique(sig.2.types.es$category_name))[rev(order(colSums(genesets.es)))][1:min(26,ncol(genesets.es))]
						genesets.es = genesets.es[rev(order(colSums(genesets.es)))][1:min(26,ncol(genesets.es))]
						colnames(genesets.es) = toupper(letters[1:ncol(genesets.es)])
						
						## call webplot data
						day.cols = unique(mta$timec[mta$time %in% postb.times])
						
						radarplot(genesets.es[1,min(26,ncol(genesets.es)):1],main=paste(spcLabl,trtLabl,sep=', '),data.row=rownames(genesets.es)[1],
								col=day.cols[1],lty=1,scale=F,labels=rev(colnames(genesets.es)),
								max=max(genesets.es),increment=1,lwd=2,txt.cex=1.3,cex.main=2)
						
						for (d in 2:(nrow(genesets.es)-1)) {
							radarplot(genesets.es[d,min(26,ncol(genesets.es)):1],data.row=rownames(genesets.es)[d],add=T,col=day.cols[d],lty=1,scale=F,
									max=max(genesets.es),increment=1,lwd=2)
						}
							
						radarplot(genesets.es[nrow(genesets.es),min(26,ncol(genesets.es)):1],data.row=rownames(genesets.es)[nrow(genesets.es)],add=T,col=day.cols[d],lty=1,scale=F,
								max=max(genesets.es),increment=1,lwd=2,last=T)
						
						
						plot.new()
						legend('topleft', lty = c(1, 1, 1, 1), ncol=ceiling(length(postb.times)/2), lwd = 4, col = day.cols, legend=postb.times.l[1:length(postb.times)],title=expression('Enrichment Score: -1 x log'[10]*'(FDR-adjusted p-value)'), bty = "n",cex=1)
						legend('bottomleft',legend=paste(colnames(genesets.es),': ',names,sep=''),bty='n',cex=0.9)
						
						## caption add-on for KnitR code
						RadarLablGsea = c(RadarLablGsea,paste(geneset.type,' (',spcLabl,', ',trtLabl,')',sep=''))
					} 
				}
			}
		}
	}
}
