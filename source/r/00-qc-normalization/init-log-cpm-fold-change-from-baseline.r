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
# Program:  init-log-cpm-fold-change-from-baseline.r 
# Version:  RSEQREP 1.1.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose: 	For each specimen type, calculate log fold change from pre-treatment.  
#				If multiple treatment points exist, utilize the mean.
# Input:    analysis/lcpm/<spc>_alltp_lcpm_tmm_normalized_<unfiltered|filtered>.tab.gz
# Output:  	analysis/lcpm_fc/<spc>_posttp_lcpm_fold_change_tmm_normalized_<unfiltered|filtered>.tab.gz
#############################################################################################################

source('init-analysis.r')

## set directories
in.dir.lcpm = paste(res.dir,'lcpm',sep='/');
out.dir.lcpm.fc = paste(res.dir,'lcpm_fc',sep='/'); 

## update specimen flag to inclue a variable across all
## types if there is more than 1 specimen type
if (length(spcFlags)>1) {
	spcFlags = c('all',spcFlags)
}
## loop through filtered and unfiltered sets
for (f in 1:length(fltFlags)) {
	fltFlag = fltFlags[f]
	
	## loop through specimen types
	for(s in 1:length(spcFlags)) {
		spcFlag = spcFlags[s];
		
		## read TMM normalized filtered log counts per million
		in.file = paste(in.dir.lcpm,'/',spcFlag,'_alltp_lcpm_tmm_normalized_',fltFlag,'.tab.gz',sep='');
		dta.spc = read.csv(in.file,header=T,stringsAsFactors=F,sep='\t');
		
		## align meta data and samples
		if (spcFlag=='all') {
			mta.spc=mta[match(colnames(dta.spc),mta$samid),]
		} else {
			mta.spc=mta[intersect(match(colnames(dta.spc),mta$samid),which(mta$spct==spcFlag)),]
		}
		
		## prepare result data frame to store fold changes -- select psot-treatment samples where 
		## that subject has a treatment sample.
		colNames = mta.spc[mta.spc$time %in% postb.times & mta.spc$subid %in% mta.spc$subid[mta.spc$timeb %in% b.times],'samid'];
		rowNames = rownames(dta.spc);
		dta.spc.fld = data.frame(matrix(rep(NA,length(colNames)*length(rowNames)),ncol=length(colNames)));
		colnames(dta.spc.fld) = colNames;
		rownames(dta.spc.fld) = rowNames;
		
		## calculate fold change within subjects that have treatment time
		subs = unique(mta.spc$subid[mta.spc$time %in% b.times]);
		
		## for each subject obtain treatment and post_treatment
		## log_fold change (difference on log scale)
		for(i in 1:length(subs)) {
			
			## get pre_treatment ids
			ids.sub.btp = mta.spc[mta.spc$time %in% b.times & mta.spc$subid==subs[i] ,'samid'];
			
			## get post_treatment ids
			ids.sub.ptp = mta.spc[mta.spc$time %in% postb.times & mta.spc$subid==subs[i] ,'samid'];
			
			## calculate mean of treatment log counts per million
			## using pre/post treatment ids.  If no pre_treatment exists, 
			## caclulate mean across all pre_treatment samples
			if(length(ids.sub.btp)==1) {
				mean.btp = dta.spc[,as.character(ids.sub.btp)]
			} else if (length(ids.sub.btp)>0){
				mean.btp = apply(dta.spc[,as.character(ids.sub.btp)],1,mean);
			} else {
				all.btp = mta.spc[mta.spc$times <=0 ,'samid']
				mean.btp = apply(dta.spc[,as.character(all.btp)],1,mean)
			}
			## calculate log fold change for post-treatment ids
			fold.ptp = dta.spc[,as.character(ids.sub.ptp)]-mean.btp;
	
			## update result data frame
			dta.spc.fld[,as.character(ids.sub.ptp)] = fold.ptp;
		}
		
		## save results
		out.file = paste(out.dir.lcpm.fc,'/',spcFlag,'_posttp_lcpm_fold_change_tmm_normalized_',fltFlag,'.tab',sep='');
		write.table(dta.spc.fld,out.file, sep='\t',quote=F,row.names=T);
		R.utils::gzip(out.file,overwrite=TRUE)
	}
}