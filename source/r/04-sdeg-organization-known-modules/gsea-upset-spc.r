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
# Program:  gsea-updset-spc.r
# Version:  RSEQREP 2.3.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Generate Upset plots by specimen type for significant GSEA results
# Input:    analysis/gsea/<spc>_<trt>_<time>_<gene-set-type>_gsea_all.tab.gz
#			data/gene_sets/all.ensembl.tab.gz
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## specify gene sets file
sets.infile = paste(dta.dir,'/gene_sets/all.ensembl.tab.gz',sep='')

if (file.exists(sets.infile)) {
	
	## only generate plots if there is more than 5 spcFlags
	if (length(spcFlags) >= 2) {
		
		## set directories
		in.dir.gsea = paste(res.dir,'gsea',sep='/');
		
		## loop through treatment groups;
		for (v in 1:length(trtFlags)) {
			trtFlag = trtFlags[v];
			trtLabl = trtLabls[v];
			
			## loop through each gene set category
			for(k in 1:length(geneset.types)) {
				geneset.type = geneset.types[k];
				
				##  get unique list of genes across timepoints
				gsea.sig.time = c()
				
				## initialize list to hold grobs/viewports
				vps = list()
				
				## loop through each post treatment timepoint
				for(t in 1:length(postb.times)) {
					time = postb.times[t];
					timel = postb.times.l[t];
					
					##  get unique list of gsea
					gsea.sig = c()
					
					## Loop threough specimen types
					for(s in 1:length(spcFlags)) {
						spcFlag = spcFlags[s];
						spcLabl = spcLabls[s];
						
						in.file = paste(in.dir.gsea,'/',spcFlag,'_',trtFlag,'_tp',time,'_',gsub(' ','_',tolower(geneset.type)),'_gsea_sig.tab.gz',sep='');
						if (file.exists(in.file)) {
							gsea.res = read.table(in.file,sep='\t',stringsAsFactors=F,header=T,quote="",fill=T);
							gsea.sig = unique(c(gsea.sig,gsea.res$category_name));
							gsea.sig.time = unique(c(gsea.sig.time,gsea.res$category_name));
						}
					}
					
					## gseaerate binary matrix
					matrx = matrix(rep(0,(length(spcFlags)*length(gsea.sig))),ncol=length(spcFlags))
					rownames(matrx) = gsea.sig
					## Loop threough specimen types
					for(s in 1:length(spcFlags)) {
						spcFlag = spcFlags[s];
						spcLabl = spcLabls[s];
						
						in.file = paste(in.dir.gsea,'/',spcFlag,'_',trtFlag,'_tp',time,'_',gsub(' ','_',tolower(geneset.type)),'_gsea_sig.tab.gz',sep='');
						if (file.exists(in.file)) {
							gsea.res = read.table(in.file,sep='\t',stringsAsFactors=F,header=T,quote="",fill=T);
							matrx[gsea.res[,'category_name'],s] = 1
						}
					}
					
					## Print significant genes diagram
					colnames(matrx) = paste(spcLabls,' (',colSums(abs(matrx)),')',sep='')
					matrx.all = as.data.frame(matrx)[,ncol(matrx):1]
					if (length(gsea.sig)>0) {
						vps = append(vps,list(grid.grabExpr(c(
								upset.custom(matrx.all,mainbar.y.label='GSEA Intersection Size',sets.x.label='GSEA Set Size',mainbar.y.max=ceiling(max(table(apply(matrx.all,1,paste,collapse=';')))*1.25),
											matrix.color=lfc.cols[2],main.bar.color=lfc.cols[2],sets.bar.color=lfc.cols[2],sets=colnames(matrx.all),keep.order=T,order.by="freq"),
										grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "black")),
										grid.text(paste(trtLabl,'\n',timel,'\n',geneset.type,sep=''),x = unit(0.68, "npc"), y = unit(0.96, "npc"),gp=gpar(font=2,cex=1),just='top'),
										grid.edit('arrange',name='arrange2')))))
					}
				}
				
				## generate binary matrix
				matrx = matrix(rep(0,(length(spcFlags)*length(gsea.sig.time))),ncol=length(spcFlags))
				rownames(matrx) = gsea.sig.time
				
				## Loop threough specimen types
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
					matrx[rownames(matrx) %in% gsea.sig,s] = 1
				}
				
				## Print significant GSEA diagram
				colnames(matrx) = paste(spcLabls,' (',colSums(matrx),')',sep='')
				matrx = as.data.frame(matrx)[,ncol(matrx):1]
				if (length(gsea.sig.time)>0) {
					vps = append(vps,list(grid.grabExpr(c(
						upset.custom(matrx,mainbar.y.label='GSEA Intersection Size',sets.x.label='GSEA Set Size',col=lfc.cols[2],
									keep.order=T,order.by="freq",mainbar.y.max=ceiling(max(table(apply(matrx,1,paste,collapse=';')))*1.25)),
								grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "black")),
								grid.text(paste(trtLabl,'\n',postb.times.l[length(postb.times.l)],'\n',geneset.type,sep=''),x = unit(0.68, "npc"), y = unit(0.96, "npc"),gp=gpar(font=2,cex=1),just='top'),
								grid.edit('arrange',name='arrange2')))))
				
				}
				## plot if there are one or more plots
				if (length(vps)>0) {
					add.on.val = if(length(vps)<10){''}else{paste(1:(ceiling(length(vps) / 9)),'of',(ceiling(length(vps) / 9)))}
					for (i in 1:ceiling(length(vps) / 9)) {
						grid.newpage()
						do.call("grid.arrange", c(vps[((i*9)-8):min((i*9),length(vps))], ncol=3,nrow=3,newpage=F))
						if(exists('reportAssayFlag')){
							assay.label = ifelse(reportAssayFlag=='rna_seq','RNA-Seq','RP')
							upsetLablsSpcGsea = c(upsetLablsSpcGsea,paste0(add.on.val,' (',assay.label,', ',trtLabl,', ',geneset.type,').'))
						}else{
							upsetLablsSpcGsea = c(upsetLablsSpcGsea,paste(add.on.val,' (',trtLabl,', ',geneset.type,').',sep=''))
						}
					}
				}
			}
		}
	}
}