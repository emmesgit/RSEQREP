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
# Program:  init-tmm-normalization-fragments.r
# Version:  RSEQREP 1.1.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Calculate moderated log2 counts per million (before/after TMM normalization)
#				and filter lowly expressed genes per specimen type.
# Input:    data/fragment_count_matrix.tab
#			data/annot/filtered_gene_annotations.tab.gz
# Output:  	analysis/lcpm/<spc>_alltp_lcpm_<norm>_normalized_<unfiltered|filtered>.tab
#############################################################################################################

source('init-analysis.r')

## set directories
out.dir.lcpm = paste(res.dir,'/lcpm',sep='');
out.dir.glm  = paste(res.dir,'/glm',sep='');

## load gene list
gen = read.csv(paste(dta.dir,'/annot/filtered_gene_annotations.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T);

## load the count matrix
dta = read.csv(paste(dta.dir,'/fragment_count_matrix.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T,row.names=1)

## update specimen flag to inclue a variable across all types
spcFlags = c('all',spcFlags)
	
## loop through specimen types
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];

	## remove outliers for per specimen type for normalization;
	if (spcFlag=='all') {
		mta.spc = mta;
	} else {
		mta.spc=mta[!mta$samid %in% outliers & mta$spct==spcFlag,];
	}

	##################################################################################
	##
	## CREATE DATA MATRIX
	##
	##################################################################################

	dta.mtx.spc = data.matrix(dta[, match(mta.spc$samid,colnames(dta))]);
	
	##################################################################################
	##
	## EXECUTE TMM NORMALIZATION 
	## CALCULATE MODERATED LOG2 COUNTS PER MILLION
	## SAVE DGE LIST OBJECT
	##
	##################################################################################
	
	## convert the data into an edgeR object
	dge.spc = DGEList(counts=dta.mtx.spc,lib=colSums(dta.mtx.spc),group=rep(1:ncol(dta.mtx.spc)));
	
	## calculate TMM scaling factors
	dge.spc = calcNormFactors(dge.spc);
	
	## save edgeR DGE R object for differential and exploratory analysis
	filename = paste(out.dir.glm,'/',spcFlag,'_alltp_edge_r_dge.RData',sep='');
	save(dge.spc,file=filename,compress=T);
	
	## save TMM scaling factors
	tmm.factors = data.frame(samid=rownames(dge.spc$samples),total_mapped=dge.spc$samples$lib.size,norm_factors=dge.spc$samples$norm.factors)
	out.file.scal.fct = paste(out.dir.lcpm,'/',tolower(spcFlag),'_alltp_edge_r_scaling_factors.tab',sep='');
	write.table(tmm.factors,out.file.scal.fct,sep='\t',quote=F,row.names=F);
	R.utils::gzip(out.file.scal.fct,overwrite=TRUE)
		
	##################################################################################
	##
	## calculate and save moderated log2 counts per million (before and after TMM normalization);
	##
	##################################################################################
	
	lcpm.spc.spc = filterUnwantedGenes(cpm(dge.spc,normalized.lib.sizes=F,prior.count=prior.count,log=T),gen);
	lcpm.spc.nrm = filterUnwantedGenes(cpm(dge.spc,normalized.lib.sizes=T,prior.count=prior.count,log=T),gen);
	out.file.spc = paste(out.dir.lcpm,'/',tolower(spcFlag),'_alltp_lcpm_not_normalized_unfiltered.tab',sep='');
	out.file.nrm = paste(out.dir.lcpm,'/',tolower(spcFlag),'_alltp_lcpm_tmm_normalized_unfiltered.tab',sep='');
	
	## write results
	write.table(lcpm.spc.spc,out.file.spc,sep='\t',quote=F,row.names=T);
	write.table(lcpm.spc.nrm,out.file.nrm,sep='\t',quote=F,row.names=T);
	
	## gzip
	R.utils::gzip(out.file.spc,overwrite=TRUE)
	R.utils::gzip(out.file.nrm,overwrite=TRUE)
	
	##################################################################################
	##
	## filter genes that are lowly expressed for post_treatment times and save filtered results
	##
	##################################################################################
	
	lcpm.spc.spc.flt = get(filterFun)(spcFlag,mta.spc,lcpm.spc.spc,save=F);
	lcpm.spc.nrm.flt = get(filterFun)(spcFlag,mta.spc,lcpm.spc.nrm,save=T);
	out.file.spc.flt = paste(out.dir.lcpm,'/',tolower(spcFlag),'_alltp_lcpm_not_normalized_filtered.tab',sep='');
	out.file.nrm.flt = paste(out.dir.lcpm,'/',tolower(spcFlag),'_alltp_lcpm_tmm_normalized_filtered.tab',sep=''); 
	
	## write results
	write.table(lcpm.spc.spc.flt, out.file.spc.flt, sep='\t',quote=F,row.names=T);
	write.table(lcpm.spc.nrm.flt, out.file.nrm.flt, sep='\t',quote=F,row.names=T);

	## gzip
	R.utils::gzip(out.file.spc.flt,overwrite=TRUE)
	R.utils::gzip(out.file.nrm.flt,overwrite=TRUE)
}