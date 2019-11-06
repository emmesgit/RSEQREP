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
# Program:  init-reactome-gmt.r
# Version:  RSEQREP 2.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose: 	generate Reactome based GMT formatted gene sets for enrichment analysis
# Input:    http://www.reactome.org/download/current/Ensembl2Reactome.txt
# Output:  	./reactome.gmt
#############################################################################################################

raw.dta = read.csv('http://www.reactome.org/download/current/Ensembl2Reactome.txt',header=F,sep='\t', stringsAsFactors=F)

## apply column names
colnames(raw.dta) = c('ensembl_id','ptw_short','ptw_web','ptw_desc','tmp','species')

## subset to human only reactrome pathways
raw.dta.homo = raw.dta[raw.dta$species=='Homo sapiens',]

## remove white spaces from ptw_desc
raw.dta.homo$ptw_desc = gsub('^ | $','',raw.dta.homo$ptw_desc)

## make lists of ensembl genes for each pathway
homo.ptw.list = tapply(raw.dta.homo$ensembl_id,raw.dta.homo$ptw_desc,list)

## generate a pathway to webpage mapping -- order to list
homo.ptw.web.map = unique(cbind(raw.dta.homo$ptw_desc,raw.dta.homo$ptw_web))
homo.ptw.web.map = homo.ptw.web.map[match(names(homo.ptw.list),homo.ptw.web.map[,1]),]
		
## for each list entry, print to file
for (i in 1:length(homo.ptw.list)) {
	tmp = matrix(c(names(homo.ptw.list)[i],homo.ptw.web.map[i,2],homo.ptw.list[[i]]),nrow=1)
	write.table(tmp,'reactome.gmt',append=T,quote=F,col.names=F,row.names=F,sep='\t')
}