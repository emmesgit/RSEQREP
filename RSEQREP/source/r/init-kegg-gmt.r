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
# Program:  init-kegg-gmt.r
# Version:  RSEQREP 2.2.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose: 	generate Kegg pathways based GMT formatted gene sets for enrichment analysis
# Input:    http://rest.kegg.jp/<get|link|list>/
# Output:  	./kegg.gmt
#############################################################################################################

ptw2kegg = data.frame(do.call('rbind', strsplit(system('wget -q -O- http://rest.kegg.jp/link/pathway/hsa',intern=T), '\t', fixed=TRUE)),stringsAsFactors=F)
ptw = data.frame(do.call('rbind', strsplit(system('wget -q -O- http://rest.kegg.jp/list/pathway',intern=T), '\t', fixed=TRUE)),stringsAsFactors=F)

ptw2ens = as.data.frame(matrix(nrow=0,ncol=3),stringsAsFactors=F)
colnames(ptw2ens) = c('pathway_id','pathway_link','ensembl_id')
for (i in 1:nrow(ptw2kegg)) {
	ens.call = paste('wget -q -O- http://rest.kegg.jp/get/',ptw2kegg[i,1],'  | grep -Po "ENSG[0-9]{11}"',sep='')
	ens.ids = system(ens.call,intern=T)
	
	ptw.long = ptw[match(gsub('hsa','map',ptw2kegg[i,2]),ptw[,1]),2]
	ptw.link = paste('http://www.genome.jp/kegg-bin/show_pathway?',ptw2kegg[i,2],sep='')
	
	if (length(ens.ids)>0) {
		for (e in 1:length(ens.ids)) {
			ptw2ens = rbind(ptw2ens,c(ptw.long,ptw.link,ens.ids[e]))
		}
	}
	if (i %% 100 == 0) {
		print(c(i,'Complete of',nrow(ptw2kegg)))
	}
}

## make lists of ensembl genes for each pathway
ptw2ens.list = tapply(ptw2ens$ensembl_id,ptw2ens$pathway_id,list)

## generate a pathway to webpage mapping -- order to list
ptw2ens.web.map = unique(cbind(ptw2ens$pathway_id,ptw2ens$pathway_link))
ptw2ens.web.map = ptw2ens.web.map[match(names(ptw2ens.list),ptw2ens.web.map[,1]),]
		
## for each list entry, print to file
for (i in 1:length(ptw2ens.list)) {
	tmp = matrix(c(names(ptw2ens.list)[i],ptw2ens.web.map[i,2],ptw2ens.list[[i]]),nrow=1)
	write.table(tmp,'kegg.gmt',append=T,quote=F,col.names=F,row.names=F,sep='\t')
}