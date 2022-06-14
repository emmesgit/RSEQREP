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
# Program:  init-gsea-sampling.r
# Version:  RSEQREP 2.3.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Run gene set enrichment analysis using GOseq to adjust for 
#			length bias in calling genes differentially expressed.
#			http://www.bioconductor.org/packages/release/bioc/html/goseq.html
# Input:    data/annot/filtered_gene_annotations.tab.gz
#			data/msigdb/all.ensembl.tab.gz
#			analysis/glm/<spc>_<trt>_<time>_glm_sig.tab.gz
# Output:  	analysis/gsea/<spc>_<trt>_<time>_<gene-set-type>_gsea_all.tab.gz
#############################################################################################################

source('init-analysis.r')

## specify gene sets file
sets.infile = paste(dta.dir,'/gene_sets/all.ensembl.tab.gz',sep='')

if (file.exists(sets.infile)) {
	## set directories
	in.dir.glm  = paste(res.dir,'glm',sep='/');
	out.dir.gsea= paste(res.dir,'gsea',sep='/');
	
	## import gene set information
	sets = read.table(sets.infile,sep='\t',stringsAsFactors=F,header=F,quote="",fill=T,row.names=NULL);
	colnames(sets) = c('category_type','gene_id','category_name','category_link')
	
	## load gene annotations and update gene sets to reflect only genes in annotations
	gen = read.csv(paste(dta.dir,'/annot/filtered_gene_annotations.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T);
	gen.sets = sets[sets$gene_id %in% gen$gene_id,]
	
	# parallel implementation
	registerDoParallel(cores=ncores) # use all avaliable cores (~1.5Gb/thread)
	splits = split(expand.grid(1:length(spcFlags),1:length(trtFlags),1:length(postb.times),1:length(geneset.types)),1:(length(spcFlags)*length(trtFlags)*length(postb.times)*length(geneset.types)))
	foreach (z=1:length(splits)) %dopar% {
		split = unlist(splits[[z]])
		spcFlags = spcFlags[split[1]]
		trtFlags = trtFlags[split[2]]
		postb.times = postb.times[split[3]]
		geneset.types = geneset.types[split[4]]
		
		## loop through specimen types
		for(s in 1:length(spcFlags)) {
			spcFlag = spcFlags[s];
			
			## loop through treatment groups;
			for (v in 1:length(trtFlags)) {
				trtFlag = trtFlags[v];
				
				## loop through post treatment timepoints
				for(t in 1:length(postb.times)) {
					time = postb.times[t];
						
					## specify infile
					in.file.sdeg = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz',sep='');
				
					## ensure infile exists
					if(file.exists(in.file.sdeg)) {
					
						## read gc adjusted glm results
						sdeg.cel.time.sig = read.table(gzfile(in.file.sdeg),sep='\t',stringsAsFactors=F,header=T);
						
						## all significatn genes
						gen.sig = unique(rownames(sdeg.cel.time.sig));
						
						## upregulated genes
						gen.upr = unique(rownames(sdeg.cel.time.sig)[2^sdeg.cel.time.sig$logFC> 1]);
						
						## downregulated genes
						gen.dwr = unique(rownames(sdeg.cel.time.sig)[2^sdeg.cel.time.sig$logFC< 1]);
						
						## 0 = not significant, 1 = significant
						gen.bin = as.integer(gen$gene_id %in% gen.sig)
						names(gen.bin)= gen$gene_id;
						
						## OBTAIN NULL DISTRIBUTION
						pwf = nullp(DEgenes=gen.bin, genome=NA, id=NA, bias.data=gen$max_trans_length, plot.fit=FALSE);
						
						## run the analysis for each gene set category
						for(c in 1:length(geneset.types)) {
							geneset.type = geneset.types[c];
							
							gen.set = gen.sets[gen.sets$category_type==geneset.type,];
							
							set.seed(seed)
							res = goseq(
									pwf,
									genome=NA,
									id=NA, 
									gene2cat=gen.set[,c('gene_id','category_name')],
									method = goseq.method,
									repcnt = goseq.randomizations);
								
							## ADJUST FOR MULTIPLE TESTING
							res$fdr.or = p.adjust(res$over_represented_pvalue,method='BH');
							res$fdr.ur = p.adjust(res$under_represented_pvalue,method='BH');		
							
							## create list of gene sets ordered by results and calculate Jaccard similarity index
							gene.set.list = tapply(gen.set$gene_id,gen.set$category_name,list)[res$category]
							res$jidx = round(unlist(lapply(gene.set.list,jaccardIdx,y=gen.sig)),4)
							
							names(res) = c('category_name','pval_or','pval_ur','sdeg_genes',
									'category_genes','fdr_or','fdr_ur','jidx') 
							
							## add up/down regulated genes
							for(r in 1:nrow(res)) {
								cat.id = res[r,'category_name'];
								res[r,'sdeg_genes_upr']=sum(gen.set[gen.set$category_name == cat.id,'gene_id'] %in% gen.upr);
								res[r,'sdeg_genes_dwr']=sum(gen.set[gen.set$category_name == cat.id,'gene_id'] %in% gen.dwr);
							}
							
							## write pwf results
							if (c==1) {
								out.file.pwf = paste(out.dir.gsea,'/',spcFlag,'_',trtFlag,'_tp',time,'_gsea_pwf.tab',sep='');
								write.table(pwf,file=out.file.pwf,sep='\t',quote=F,row.names=T);
								R.utils::gzip(out.file.pwf,overwrite=TRUE);
							}
							
							## sort table
							res = sqldf("select * from res order by fdr_or asc, jidx desc")
							
							## write all results
							out.file.res = paste(out.dir.gsea,'/',spcFlag,'_',trtFlag,'_tp',time,'_',gsub(' ','_',tolower(geneset.type)),'_gsea_all.tab',sep='');
							write.table(res,file=out.file.res,sep='\t',quote=F,row.names=F);
							R.utils::gzip(out.file.res,overwrite=TRUE);
							
							## write significant results
							res.sig = res[res$fdr_or <= goseq.fdr.cutoff,]
							if (nrow(res.sig)>0) {
								out.file.res.sig = paste(out.dir.gsea,'/',spcFlag,'_',trtFlag,'_tp',time,'_',gsub(' ','_',tolower(geneset.type)),'_gsea_sig.tab',sep='');
								write.table(res.sig,file=out.file.res.sig,sep='\t',quote=F,row.names=F);
								R.utils::gzip(out.file.res.sig,overwrite=TRUE);
							}
						}
					}
				}
			}
		}
	}
}