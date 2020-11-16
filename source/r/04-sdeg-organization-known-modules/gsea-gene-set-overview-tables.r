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
# Program:  gsea-gene-set-overview-tables.r
# Version:  RSEQREP 2.0.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Tabulate an overview of the gene sets used in enrichment analysis
# Input:    analysis/gsea/<spc>_<trt>_<time>_<gene-set-type>_gsea_all.tab.gz
#			data/gene_sets/all.ensembl.tab.gz
#			data/annot/filtered_gene_annotations.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r');

## specify gene sets file
sets.infile = paste(dta.dir,'/gene_sets/all.ensembl.tab.gz',sep='')

if (file.exists(sets.infile)) {
	
	## set directories
	in.dir.glm = paste(res.dir,'glm',sep='/');
	
	## import gene set information
	sets = read.table(sets.infile,sep='\t',stringsAsFactors=F,header=F,quote="",fill=T,row.names=NULL);
	colnames(sets) = c('category_type','gene_id','category_name','category_link')
	
	## load gene annotations and update gene sets to reflect only genes in annotations
	gen = read.csv(paste(dta.dir,'/annot/filtered_gene_annotations.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T);
	sets = sets[sets$gene_id %in% gen$gene_id,]
	
	## get gene count per category id
	gen.stat = sqldf('select category_type,category_name,count(distinct gene_id) as gene_count from sets 
					group by category_name order by category_type,category_name');
	
	## gen median gene count across category type
	cat.med = round(tapply(gen.stat$gene_count,factor(gen.stat$category_type),median));
	
	## get number of categories per category type and overall number of distinct genes
	x = sqldf('select category_type,count(distinct category_name) as count,count(distinct gene_id) as genes from sets 
	group by category_type order by category_type');
	
	for(r in 1:nrow(x)) {
		cat.type = x[r,'category_type'];
		cat.gen.med =  cat.med[names(cat.med)==cat.type];
		x[r,'medgenes'] = cat.gen.med;
	}
	
	colnames(x) = c(
			'Category Type',
			'Categories',
			'Distinct \\#Genes In Sets',
			'Median \\#Genes Per Set'
	);	
	
	write.tbl(x,rownames=F);
	
	caption.short = paste('Overview of gene sets used for the enrichment analysis (RNA-Seq).',sep='');
	caption.long  = paste('Overview of gene sets used for the enrichment analysis (RNA-Seq).  Genes within gene sets are fltered to reflect only those that exist in filtered Ensembl version ',ensembl.version,' anotations obtained using biomaRt.',sep='');
	
	x = xtable(x,align=c('c','p{2.5in}','p{0.6in}','p{0.4in}','p{0.4in}'),digits=rep(0,5),
			type='latex',floating='F',label='tab:gsea_gene_sets',
			caption=c(caption.long,caption.short)
	);
	
	print(x,tabular.environment='longtable',
			floating=F,
			include.rownames=F,
			hline.after = c(-1,nrow(x)),
			size= 'footnotesize',
			add.to.row = list(pos = list(0),command = "\\hline \\endhead "),
			sanitize.text.function = identity,
			sanitize.rownames.function = identity)
}