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
# Program:  tbl-retained-genes.r
# Version:  RSEQREP 2.1.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Print # genes after filtering
# Input:    analysis/lcpm/<spc>_<time>_analysis_gene_set.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## set directories
in.dir.lcpm = paste(res.dir,'lcpm',sep='/');

## create empty matrix for results
res = matrix(ncol=1 ,nrow=length(spcFlags) )

## loop through specimen types
for(c in 1:length(spcFlags)) {
	spcFlag = spcFlags[c];
	
	## only pursue 'posttp'
	for(d in length(postb.times.m)) {
		time = postb.times.m[d];
		
		## specify infile name
		in.file= paste(in.dir.lcpm,'/',spcFlag,'_',time,'_analysis_gene_set.tab.gz',sep='')
		
		## if file exists, open it and update matrix
		if (file.exists(in.file)) {
			dta = read.table(paste(in.dir.lcpm,'/',spcFlag,'_',time,'_analysis_gene_set.tab.gz',sep=''));
			res[c,1] = nrow(dta);
		}
	}
}

## update column and rownames
rownames(res) = spcLabls;
colnames(res) = '#Genes';

## print results in tabular format
write.tbl(res,rownames=T);

## specify caption
if(exists('reportAssayFlag')){
	assay.label = ifelse(reportAssayFlag=='rna_seq','RNA-Seq','RP')
	caption = paste0('Number of genes that passed the low expression cut off (',assay.label,').')
}else{
	caption = 'Number of genes that passed the low expression cut off (RNA-Seq).'
}


## create xtable and print
x = xtable(res,align=c('p{2.2in}','p{0.5in}'),type='latex',floating='F',
		label="tab:retainedgenes",caption=caption);

print(x,tabular.environment='longtable',
		floating=F,
		include.rownames=T,
		hline.after = c(-1,nrow(x)),
		size= 'small',
		add.to.row = list(pos = list(0),command = "\\hline \\endhead "))