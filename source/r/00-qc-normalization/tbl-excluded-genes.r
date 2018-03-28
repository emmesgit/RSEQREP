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
# Program:  tbl-excluded-genes.r 
# Version:  RSEQREP 1.1.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Tabulate genes that were removed prior to the analysis.
# Input:    data/annot/excluded_gene_counts.tab
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

dta =read.csv(paste(dta.dir,'/annot/excluded_gene_counts.tab.gz',sep=''),stringsAsFactors=F,header=T,sep='\t');

dta.tapply = rbind(
		cbind(names(tapply(dta$gene_id,dta$gene_type,length)),tapply(dta$gene_id,dta$gene_type,length)),
		c('Total',nrow(dta)))


colnames(dta.tapply)=c('Category/Chromosome','#Genes')
dta.tapply[,1] = c('mitochondrial ribosomal RNA (Mt_rRNA)','mitochondrial transfer RNA (Mt_tRNA)','ribosomal RNA (rRNA)','Total')

x = xtable(dta.tapply,align=c('l','l','l'),type='latex',floating='F',label="tab:exgenes",caption='Number of excluded genes by gene type (RNA-Seq).');

write.tbl(x,rownames=F);

## print latex table
print(x,tabular.environment='longtable',
		floating=F,
		include.rownames=F,
		hline.after = c(-1,nrow(x)),
		size= 'small',
		add.to.row = list(pos = list(0),command = "\\hline \\endhead "))

