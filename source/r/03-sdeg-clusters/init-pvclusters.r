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
# Program:  init-pvclusters.r
# Version:  RSEQREP 2.1.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Identify significant clusters via boostrap analysis of hierarchical clustering 
#			results. Executes multiscale boostrap resampling (pvclust) for hierarchical 
#			clustering results based on treatment log fold changes.  Please note that if parallel=T is set for pvclust(), 
# 			results may or may not be 100% reproducable, as the parallelization process does not implement
# 			the seed in a reproducable way when threading the processes.
# Input:    analysis/glm/<spc>_<trt>_<time>_glm_sig.tab
#			analysis/lcpm_fc/<spc>_psottp_lcpm_fold_change_tmm_normalized_filtered.tab
# Output:  	analysis/pvclust/<spc>_<trt>_<time>_pvclust.RData
#			analysis/pvclust/<spc>_<trt>_<time>_pvclust_sig.tab.gz
#############################################################################################################

source('init-analysis.r')

## set directories
in.dir.glm  = paste(res.dir,'glm',sep='/');
in.dir.fld  = paste(res.dir,'lcpm_fc',sep='/');
out.dir.pvc = paste(res.dir,'pvclust',sep='/');

## load gene annotations
gen = read.csv(paste(dta.dir,'/annot/filtered_gene_annotations.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T);

# parallel implementation
registerDoParallel(cores=ncores) # (~1Gb mem/thread)
splits = split(expand.grid(1:length(spcFlags),1:length(postb.times.m)),1:(length(spcFlags)*length(postb.times.m)))
foreach (z=1:length(splits)) %dopar% {
	split = unlist(splits[[z]])
	spcFlags = spcFlags[split[1]]
	postb.times.m = postb.times.m[split[2]]
	
	## loop through specimen types
	for(s in 1:length(spcFlags)) {
		spcFlag = spcFlags[s];
		
		## GET treatment LOG FOLD CHANGES
		in.file = paste(in.dir.fld,'/',spcFlag,'_posttp_lcpm_fold_change_tmm_normalized_filtered.tab.gz',sep='');
		dta.lfc = read.csv(in.file,sep='\t',stringsAsFactors=F,header=T,check.names=F);
		
		## align meta data and samples 
		mta.lfc=mta[unique(match(colnames(dta.lfc),mta$samid)),]
		
		## loop through post treatment timepoiunts
		## including all post treatment days
		for(t in 1:length(postb.times.m)) {
			postb.time = postb.times.m[t];
			
			## Get significant genes across all post treatment 
			## timepoints and treatment types
			gen.sig = c();
			for (v in 1:length(trtFlags)) {
				trtFlag = trtFlags[v];
				for(tx in 1:length(postb.times)) {
					timex = postb.times[tx];
					glm.infile = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',timex,'_glm_sig.tab.gz',sep='');
					if(file.exists(glm.infile)) {
						sdeg.cel.time = read.table(glm.infile,sep='\t',stringsAsFactors=F,header=T);
						gen.sig = unique(c(gen.sig,rownames(sdeg.cel.time)));
					}
				}
			}
				
			## ennsure there are significant genes
			if(length(gen.sig) > 1) {
				
				## update metadata according to timepoint flag
				if (postb.time=='posttp') {
					mta.lfc.time =mta.lfc
				} else {
					mta.lfc.time = mta.lfc[mta.lfc$time==gsub('tp','',postb.time),];
				}
				
				## prepare data matrix [libraries x genes]
				x = t(data.matrix(dta.lfc[rownames(dta.lfc) %in% gen.sig,mta.lfc.time$samid]));
				rownames(x) = mta.lfc.time$patid;
				colnames(x)= gen[match(colnames(x),gen$gene_id),'gene_id'];
				
				## determine proportions to use for r
				if(nrow(x)<=3) {
					print(paste('ERROR Pvclusters',spcFlag,'Time point',postb.time,'Not enough samples!! Need 4 or more'))
				} else if (nrow(x)==4) {
					pvclust.prop = c(0.5,0.75,1,1.25,1.5)
				} else if (nrow(x)==5) {
					pvclust.prop = c(0.6,0.8,1,1.2,1.4)
				} else if (nrow(x)>=6 & nrow(x)<10) {
					pvclust.prop = c((nrow(x)-3)/nrow(x),(nrow(x)-2)/nrow(x),(nrow(x)-1)/nrow(x),1,(nrow(x)+1)/nrow(x),(nrow(x)+2)/nrow(x),(nrow(x)+3)/nrow(x))
				} else {
					pvclust.prop = seq(0.6,1.4,by=0.1)			
				}
				
				## remove entries where values have no variation -- this will cause an error with pvclust
				remove.ids = names(which(apply(x,2,function(y){all(y[2:length(y)] %in% y[1])})))
				x = x[,!colnames(x) %in% remove.ids]
				
				## RUN PVCLUST MULTISCALE BOOTSTRAP RESAMPLING
				res.pvc = pvclust(x,method.hclust=cluster.method,method.dist=function(x) pvclust.dist.fun(x,by.row=F),iseed=seed,nboot=pvclust.boot,r=pvclust.prop,parallel=F);
				
				## save pvclust object 
				if (nrow(res.pvc$edges)>1) {
					filename = paste(out.dir.pvc,'/',spcFlag,'_',postb.time,'_pvclust.RData',sep='');
					save(res.pvc,file=filename,compress=T);
				}	
				
				## get significant genes
				gen.sig = c();
				for (v in 1:length(trtFlags)) {
					trtFlag = trtFlags[v];
					glm.infile = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_',postb.time,'_glm_sig.tab.gz',sep='');
					if(file.exists(glm.infile)) {
						sdeg.tbl.time = read.table(glm.infile,sep='\t',stringsAsFactors=F,header=T);
						gen.sig = unique(c(gen.sig,rownames(sdeg.tbl.time)));
					}
					if(postb.time =='posttp' & length(postb.times)>1){	
						for(tx in 1:length(postb.times)) {
							timex = postb.times[tx];
							glm.infile = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',timex,'_glm_sig.tab.gz',sep='');
							if(file.exists(glm.infile)) {
								sdeg.cel.time = read.table(glm.infile,sep='\t',stringsAsFactors=F,header=T);
								gen.sig = unique(c(gen.sig,rownames(sdeg.cel.time)));
							}
						}
					}
				}
				
				## get significant clusters
				clr.res = data.frame(cluster_id=c(NA),cluster_size=c(NA),gene_id=c(NA));
				clusters = pvpick(res.pvc,alpha=pvclust.pval,max.dist=pvclust.max.dist)$clusters;
				cls.sig.num = 0;
				for(cx in 1:length(clusters)) {
					cluster = clusters[[cx]];
					if(length(intersect(cluster,gen.sig))>0) {
						cls.sig.num=cls.sig.num+1;
						if(postb.time == 'posttp') {
							time.lab = 'ALLTP';
						} else {
							time.lab = toupper(postb.time);
						}
						
						clusterId = paste(toupper(spcFlag),time.lab,'_',sprintf("%03d", cls.sig.num),sep='');
						clusterSize = length(cluster);
						for(g in 1:clusterSize) {
							gene = cluster[g];
							clr.res = rbind(clr.res,c(clusterId,clusterSize,gene));
						}
					}
				}
				clr.res=clr.res[!is.na(clr.res$cluster_id),];
				
				## ensure there is at least one cluster
				if(nrow(clr.res)>2) {	
					
					## add gene name and description fields to data
					clr.res$gene_name= gen[match(clr.res$gene_id,gen$gene_id),'gene_name'];
					clr.res$gene_desc= gen[match(clr.res$gene_id,gen$gene_id),'gene_desc'];
					
					## add treatment based log fold changes
					for (v in 1:length(trtFlags)) {
						trtLabl = trtLabls[v]
						subs = mta$samid[mta$time==gsub('tp','',postb.time) & mta$spct==spcFlag & mta$trt==trtFlag & !mta$samid %in% outliers]
						clr.res[,paste('log_fc',v,sep='')] = round(rowMeans(dta.lfc[clr.res$gene_id,match(subs,colnames(dta.lfc))]),2);
					}
					
					## order data
					clr.res = clr.res[order(clr.res$cluster_id,clr.res$gene_name),]
					
					## output significant clusters if file does not exists			
					out.file = paste(out.dir.pvc,'/',spcFlag,'_',postb.time,'_pvclust_sig.tab',sep='');
					write.table(clr.res,out.file,quote=F,row.names=F,sep='\t');
					R.utils::gzip(out.file,overwrite=TRUE)
				}
			}	
		}
	}
}