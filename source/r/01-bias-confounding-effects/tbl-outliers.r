############################################################################################
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
# Program:  tbl-outlierss.r
# Version:  RSEQREP 1.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  list outliers
# Input:    N/A
# Output:  	N/A
############################################################################################

source('../r/init-analysis.r')

if (!any(is.na(outliers))) {

	mta.out = mta[match(outliers,mta$samid),c('subid','samid','spctl','timel')]

	mta= mta.out[order(mta.out$spctl,mta.out$time),]
	colnames(mta.out)=c('Subject ID','Sequence Library ID','Specimen Type','Timepoint');

	write.tbl(mta.out,rownames=F);

	print(xtable(mta.out,align=c('c','p{0.5in}','p{1in}','p{0.6in}','p{0.9in}'),type='latex',floating='F',label="tab:outliers",
			caption='Outlying observations (RNA-Seq).'),include.rownames=FALSE,size='footnotesize');	
}