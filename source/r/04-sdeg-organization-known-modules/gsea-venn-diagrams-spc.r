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
# Program:  gsea-venn-diagrams-spc.r
# Version:  RSEQREP 1.1.3
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Generate Venn Diagram plots by specimen type for significant GSEA results
# Input:    analysis/gsea/<spc>_<trt>_<time>_<gene-set-type>_gsea_all.tab.gz
#			data/gene_sets/all.ensembl.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## specify gene sets file
sets.infile = paste(dta.dir,'/gene_sets/all.ensembl.tab.gz',sep='')

if (file.exists(sets.infile)) {
	
	## only generate plots if there are between 2 and 5 spcFlags
	if (length(spcFlags) > 1 & length(spcFlags) <= 5) {
	
		## set directories
		in.dir.gsea = paste(res.dir,'gsea',sep='/');
	
		## loop through treatment groups;
		for (v in 1:length(trtFlags)) {
			trtFlag = trtFlags[v];
			trtLabl = trtLabls[v];
			
			## loop through each gene set category
			for(k in 1:length(geneset.types)) {
				geneset.type = geneset.types[k];
			
				## Graphing parameters
				par(mar=c(0.2, 0.2, 2, 0.2),oma=c(0,0,0,0),mfrow=c(5,2),lheight=0.65);
				plot.num = 0
				
				##  get unique list of genes across timepoints
				gsea.sig.time = c()
				
				## loop through each post treatment timepoint
				for(t in 1:length(postb.times)) {
					time = postb.times[t];
					timel = postb.times.l[t];
					
					##  get unique list of gsea
					gsea.sig = c()
					for(s in 1:length(spcFlags)) {
						spcFlag = spcFlags[s];
						spcLabl = spcLabls[s];
						
						in.file = paste(in.dir.gsea,'/',spcFlag,'_',trtFlag,'_tp',time,'_',gsub(' ','_',tolower(geneset.type)),'_gsea_sig.tab.gz',sep='');
						if (file.exists(in.file)) {
							gsea.res = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T,quote="",fill=T);
							gsea.sig = unique(c(gsea.sig,gsea.res$category_name));
							gsea.sig.time = unique(c(gsea.sig.time,gsea.res$category_name));
						}
					}
					
					## gseaerate binary matrix
					matrx = matrix(rep(0,(length(spcFlags)*length(gsea.sig))),ncol=length(spcFlags))
					rownames(matrx) = gsea.sig
					for(s in 1:length(spcFlags)) {
						spcFlag = spcFlags[s];
						spcLabl = spcLabls[s];
						
						in.file = paste(in.dir.gsea,'/',spcFlag,'_',trtFlag,'_tp',time,'_',gsub(' ','_',tolower(geneset.type)),'_gsea_sig.tab.gz',sep='');
						if (file.exists(in.file)) {
							gsea.res = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T,quote="",fill=T);
							matrx[gsea.res[,'category_name'],s] = 1
						}
					}
					
					## Print significant gsea diagram
					vennDiagramMod(abs(matrx), names=spcLabls,lwd=2,cex=1.4,include='both',counts.col=c('black'),show.include=T)
					title(main=paste(trtLabl,', ',timel,'\n',geneset.type,sep=''),line=-0.5,xpd=T,cex.main=1.5)
					plot.num = plot.num + 1
				}
			
				## generate binary matrix
				matrx = matrix(rep(0,(length(spcFlags)*length(gsea.sig.time))),ncol=length(spcFlags))
				rownames(matrx) = gsea.sig.time
				
				## If there is more than one timepoint
				if(length(postb.times)>1) {
					
					## loop through specimen types
					for(s in 1:length(spcFlags)) {
						spcFlag = spcFlags[s];
						spcLabl = spcLabls[s];
						
						gsea.sig = c()
						## loop through post treatment timepoints
						for(t in 1:length(postb.times)) {
							time = postb.times[t];
							
							in.file = paste(in.dir.gsea,'/',spcFlag,'_',trtFlag,'_tp',time,'_',gsub(' ','_',tolower(geneset.type)),'_gsea_sig.tab.gz',sep='');
							if(file.exists(in.file)) {
								gsea.res = read.table(gzfile(in.file),sep='\t',stringsAsFactors=F,header=T,quote="",fill=T);
								gsea.sig = unique(c(gsea.sig,gsea.res[,'category_name']))
							}
						}
						matrx[match(gsea.sig,rownames(matrx)),s] = 1
					}
					
					## Print significant gsea diagram
					vennDiagramMod(matrx, names=spcLabls,lwd=2,cex=1.4,include='both',counts.col=c('black'),show.include=T)
					title(main=paste(trtLabl,', ',postb.times.l[length(postb.times.l)],'\n',geneset.type,sep=''),line=-0.5,xpd=T,cex.main=1.2)
					plot.num = plot.num + 1
				}
				
				## Add blank plots if there are not at least 4.
				if ((plot.num %% 10) == 1) {
					plot.new()
					points(1,0,col='grey99')
				}
			}
		}
	}
}