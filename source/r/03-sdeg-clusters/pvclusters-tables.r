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
# To cite this software, please reference doi:10.12688/f1000research.13049.1
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# Program:  pvclusters-tables.r
# Version:  RSEQREP 2.0.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Tabulate significant clusters
# Input:    data/annot/filtered_gene_annotations.tab
#			analysis/glm/<spc>_<trt>_<time>_glm_sig.tab
#			analysis/pvclust/<spc>_<trt>_<time>_pvclust_sig.tab.gz
#			analysis/pvclust/<spc>_<trt>_<time>_pvclust.RData
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## set directories
in.dir.pvc = paste(res.dir,'pvclust',sep='/');
in.dir.glm = paste(res.dir,'glm',sep='/');
in.dir.lcpmfc = paste(res.dir,'lcpm_fc',sep='/');

## load gene annotations
gen = read.csv(paste(dta.dir,'/annot/filtered_gene_annotations.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T);

## loop through specimen types
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];
	spcLabl = spcLabls[s];
	
	## loop through post treatment timepoiunts
	## including all post treatment days
	for(t in 1:length(postb.times.m)) {
		postb.time = postb.times.m[t];
		postb.timel = postb.times.l[t];
		
		## specify infile
		pvc.infile.sig = paste(in.dir.pvc,'/',spcFlag,'_',postb.time,'_pvclust_sig.tab.gz',sep='');	
		
		## if there is a significant file, generate xtable
		if(file.exists(pvc.infile.sig)) {
			clr.res = read.csv(pvc.infile.sig,header=T,sep='\t',stringsAsFactors=F)
			
			## update underscores to dashes:
			clr.res$cluster_id = gsub('_','-',clr.res$cluster_id)
		
			## update column names
			colnames(clr.res) = c('Cluster ID','Cluster Size','Gene ID', 'Gene Name','Gene Description',paste('$Log_{2}$ FC ',toupper(trtFlags)));
			
			## remove trt fields for all timepoints -- set alignemnt strings
			if(postb.time !='postd') {
				align=c('c','p{0.85in}','p{0.25in}','p{0.9in}','p{0.6in}','p{2.7in}',rep('p{0.25in}',length(trtFlags)))
			} else if(postb.time =='postd') {
				clr.res= clr.res[,-c(6,7)];
				align=c('c','p{0.85in}','p{0.25in}','p{0.9in}','p{0.6in}','p{2.7in}')	
			}
			
			## write tabular results
			write.tbl(clr.res,rownames=F);
			
			## specify caption
			caption = paste('Co-expressed gene clusters (',spcLabl,', ',postb.timel,')',sep='');
			
			## Genrate Xtable
			x = xtable(clr.res,align=align,type='latex',floating='F',
					label=c(paste('tab:clusters',spcFlag,postb.time,sep='_')),caption=caption);
			
			## print latex table
			print(x,tabular.environment='longtable',
					floating=F,
					include.rownames=F,
					hline.after = c(-1,nrow(x)),
					size= 'scriptsize',
					add.to.row = list(pos = list(0),command = "\\hline \\endhead "),
					sanitize.colnames.function = identity,
					sanitize.text.function = identity,
					sanitize.rownames.function = identity);
		}
	}
}