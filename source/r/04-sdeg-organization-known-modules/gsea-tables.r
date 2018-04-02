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
# Program:  gsea-tables.r
# Version:  RSEQREP 1.1.2
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Tabulate significant GSEA results
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
	in.dir.gsea = paste(res.dir,'gsea',sep='/');
	in.dir.glm = paste(res.dir,'glm',sep='/');
	
	## import gene set information
	sets = read.table(sets.infile,sep='\t',stringsAsFactors=F,header=F,quote="",fill=T,row.names=NULL);
	colnames(sets) = c('category_type','gene_id','category_name','category_link')
	
	## load gene annotations
	gen = read.csv(paste(dta.dir,'/annot/filtered_gene_annotations.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T);
	gen.sets = sets[sets$gene_id %in% gen$gene_id,]
	
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
				postb.time = postb.times[t];
				postb.timel = postb.times.l[t];
				
				## loop through each gene set category
				for(k in 1:length(geneset.types)) {
					geneset.type = geneset.types[k];
					
					## specify infile
					in.file = paste(in.dir.gsea,'/',spcFlag,'_',trtFlag,'_tp',postb.time,'_',gsub(' ','_',tolower(geneset.type)),'_gsea_sig.tab.gz',sep='');
					
					if(file.exists(in.file)) {
						
						## read gene set enrichment results
						res = read.table(in.file,sep='\t',stringsAsFactors=F,header=T,quote="",fill=T);
						res = res[with(res, order(fdr_or, -sdeg_genes)),]
						
						## handle p_values
						res$PVAL[res$pval_or<0.0001]='<0.0001';
						res$PVAL[res$pval_or==1]='>0.9999';
						res$PVAL[res$pval_or>=0.0001 & res$pval_or<1]=as.character(round(res$pval_or[res$pval_or>=0.0001 & res$pval_or<1],4));
						
						## add number and percent of overall SDEG genes and by up/down regulation
						res$SDEG_STR=paste(res$sdeg_genes,' (',round(res$sdeg_genes/res$category_genes*100,1),')',sep='');
						res$SDEG_UPR_STR=paste(res$sdeg_genes_upr,' (',round(res$sdeg_genes_upr/res$category_genes*100,1),')',sep='');
						res$SDEG_DWR_STR=paste(res$sdeg_genes_dwr,' (',round(res$sdeg_genes_dwr/res$category_genes*100,1),')',sep='');
						
						res$category_type = geneset.type;
						res$spc_type = spcLabl;
						res$trt_group = trtLabl;
						res$time = postb.timel;
						res$sample_id = paste(spcLabl,', ',trtLabl,' (',postb.timel,')',sep='');
						
						## extract significant genes and get a total list of filtered genes per time/cell type
						sigfile = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',postb.time,'_glm_sig.tab.gz',sep='')
						if (file.exists(sigfile)) {	
							sig = rownames(read.table(sigfile,header=T,stringsAsFactors=F,sep='\t'))
						}
						
						## order by FDR then by Jaccard index
						res = sqldf("select * from res order by fdr_or asc, jidx desc")
						
						res$FDR[res$fdr_or<0.0001]='<0.0001';
						res$FDR[res$fdr_or==1]='>0.9999';
						res$FDR[res$fdr_or>=0.0001 & res$fdr_or<1]=as.character(round(res$fdr_or[res$fdr_or>=0.0001 & res$fdr_or<1],4));
						
						## update category string
						res$category_name_new=gsub('_',' ',res$category_name)
						res$category_name_new=gsub('#','\\\\#',res$category_name_new)
						res$category_name_new=gsub('&','and',res$category_name_new)
						res[which(nchar(res$category_name_new)>50),'category_name_new'] = paste(substr(res[which(nchar(res$category_name_new)>50),'category_name_new'],1,50),'...',sep='')
						
						## add hyperlinks
						res$category_link = paste('\\parbox[t][][t]{4.7in}{\\','href{',gen.sets$category_link[match(res$category_name,gen.sets$category_name)],'}{',res$category_name_new,'}}',sep='');
						
						
						if (nrow(res)>50) {
							res = res[1:50,];
							labtxt = 'Top 50 results are listed.';
						} else {
							labtxt='';
						}
						
						res = res[,c(
										'category_name_new',
										'category_link',
										'category_genes',
										'SDEG_STR',
										'SDEG_UPR_STR',
										'SDEG_DWR_STR',
										'jidx',
										'PVAL',
										'FDR'
								)]
						
						res.exp.tbl = res[,-2];
						
						colnames(res) = c(
								'Category Name',
								'Gene Set Name',
								'Gene Set Genes \\#',
								'DE Genes N (\\%)',
								'Up-reg. \\newline DE Genes \\newline N (\\%)',
								'Down-reg. DE Genes N (\\%)',
								'Jaccard Index',
								'P-Value',
								'FDR \\newline Adjusted \\newline P-Value');	
						
						
						colnames(res.exp.tbl) = c(
								'Gene Set Name',
								'Gene Set Genes #',
								'Differentially Expressed Genes N (%)',
								'Up-regulated Differentially Expressed Genes N (%)',
								'Down-regulated Differentially Expressed Genes N (%)',
								'Jaccard Index',
								'P-Value',
								'FDR Adjusted P-Value');
											
						write.tbl(res.exp.tbl,rownames=F);
						
						caption.short = paste('Enriched ',geneset.type,' (',spcLabl,', ',trtLabl,', ',postb.timel,')',sep='');
						caption.long  = paste('Enriched ',geneset.type,' (',spcLabl,', ',trtLabl,', ',postb.timel,'). Results are sorted by FDR adjusted p-value and Jaccard similarity index. ',labtxt,sep='');
						
						x = xtable(res[,c(-1)],align=c('c','p{4.2in}','p{0.45in}','p{0.5in}','p{0.55in}','p{0.55in}','p{0.38in}','p{0.38in}','p{0.43in}'),digits=c(0,0,0,0,0,0,3,4,4),
								type='latex',floating='F',label=c(paste('tab:gsea_',gsub(' ','_',tolower(geneset.type)),'_',spcFlag,'_',postb.time,sep='')),
								caption=c(caption.long,caption.short)
						);
						
						print(x,tabular.environment='longtable',
								floating=F,
								include.rownames=F,
								hline.after = c(-1,nrow(x)),
								size= 'scriptsize',
								add.to.row = list(pos = list(0),command = "\\hline \\endhead "),
								sanitize.text.function = identity,
								sanitize.rownames.function = identity)
					}
				}
			}
		}
	}
}