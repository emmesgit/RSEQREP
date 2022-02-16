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
# Program:  init-analysis.r
# Version:  RSEQREP 2.2.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Initializes important global variabels and functions. This includes knitR setup.
# Input:    data/sample_metadata_rseqc.csv
#			data/analysis_config.csv
#			data/gene_sets/all.ensembl.tab.gz
# Output:  	N/A
#############################################################################################################

##################################################################################
##
## LOAD PACKAGES AND RESET VARIABLES
##
##################################################################################

libs = c('edgeR','biomaRt','goseq','xtable','MASS','gplots',
		 'sqldf','pvclust','car','knitr','vegan','R.utils',
		 'gtools','stringr','plyr','grid','gridExtra','UpSetR','Cairo',
		 'doParallel');
	     
## load libraries using the shared R library
lapply(libs, require, character.only=T)

## reset variables
rm(list=ls(all=TRUE));

##################################################################################
##
## SPECICIFY DIRECTORIES
##
##################################################################################

dir.dta = read.csv('../../dir.csv',header=F,sep=',',stringsAsFactors=F)

bas.dir = dir.dta[2,1]
pre.dir = dir.dta[1,1]
src.dir = paste(getwd(),'/..',sep='')

## specify main directories
res.dir  = paste(bas.dir,'analysis',sep='/');
fig.dir  = paste(bas.dir,'report','figs',sep='/');
tbl.dir  = paste(bas.dir,'report','tables',sep='/');
dta.dir  = paste(bas.dir,'data',sep='/');

##################################################################################
##
## LOAD DATA SETS
##
##################################################################################

## load rna_seq meta data
mta = read.csv(paste(dta.dir,'sample_metadata_rseqc.csv',sep='/'),stringsAsFactors=F,header=T)

##################################################################################
##
## ANALYSIS CONFIGURATION
##
##################################################################################

## Import analysis configuration
analysis.config = read.csv(paste(dta.dir,'/analysis_config.csv',sep=''),header=T,stringsAsFactors=F)

## number of cores to use for threadable processes:
ncores = as.numeric(analysis.config[analysis.config$Name=='ncores','Value']);

## Normalization filtering and 0 value imputation paramters
flt.lcpm.cut = as.numeric(analysis.config[analysis.config$Name=='gene_filter_cutoff','Value']);
flt.gene.range = as.numeric(unlist(str_split(analysis.config[analysis.config$Name=='gene_filter_range','Value'],';')))
prior.count = as.numeric(analysis.config[analysis.config$Name=='prior_count','Value']);
filterFun = analysis.config[analysis.config$Name=='gene_filter_fun','Value'];

## cut offs for the selecting differentially expressed genes
glm.sdeg.fold = as.numeric(analysis.config[analysis.config$Name=='glm_model_fc_cutoff','Value']);
glm.sdeg.qval = as.numeric(analysis.config[analysis.config$Name=='glm_model_fdr_cutoff','Value']);
glm.model.paired = as.numeric(analysis.config[analysis.config$Name=='glm_model_paired','Value']);

## maximum distance allowed to be called cluster;
pvclust.max.dist = as.numeric(analysis.config[analysis.config$Name=='pvclust_max_dist','Value']);
pvclust.pval = as.numeric(analysis.config[analysis.config$Name=='pvclust_pval','Value']);
pvclust.boot = as.numeric(analysis.config[analysis.config$Name=='pvclust_boot','Value']);
pvclust.dist.fun = analysis.config[analysis.config$Name=='pvclust_dist_fun','Value'];
cluster.method = analysis.config[analysis.config$Name=='cluster_method','Value']

## random seed
seed = as.numeric(analysis.config[analysis.config$Name=='seed','Value']);

## goseq parameters
goseq.method = analysis.config[analysis.config$Name=='goseq_method','Value'];
goseq.randomizations = as.numeric(analysis.config[analysis.config$Name=='goseq_randomizations','Value']);
goseq.fdr.cutoff = as.numeric(analysis.config[analysis.config$Name=='goseq_fdr_cutoff','Value']);

## ensembl version and genes to remove for analysis
ensembl.version = analysis.config[analysis.config$Name=='ensembl_version','Value'];
ensembl.gene.types.remove = paste(unlist(str_split(analysis.config[analysis.config$Name=='ensembl_gene_types_remove','Value'],';')),collapse='","')

## Read distribution flag 1==T 0==F
## Cram flag 1==T 0==F
rdFlag = analysis.config[analysis.config$Name=='run_read_dist','Value'] == 1;
cramFlag = analysis.config[analysis.config$Name=='save_cram','Value'] == 1;

## define colors and color ranges
lfc.col.range.raw = as.numeric(unlist(str_split(gsub('\\(','',gsub('\\)','',analysis.config[analysis.config$Name=='heatmap_lfc_col_range','Value'])),';')));
lfc.col.range = seq(lfc.col.range.raw[1],lfc.col.range.raw[2],by=1);
lfc.col.rangel = c(paste('<',lfc.col.range[1],sep=''),lfc.col.range[2:(length(lfc.col.range)-1)],paste('>',lfc.col.range[length(lfc.col.range)],sep=''));
lfc.cols = unlist(str_split(analysis.config[analysis.config$Name=='heatmap_lfc_cols','Value'],';'));
gsea.cols = unlist(str_split(analysis.config[analysis.config$Name=='gsea_cols','Value'],';'));

## distance function for heatmaps (use get later -- for now store as string)
heatmap.dist.fun = analysis.config[analysis.config$Name=='heatmap_dist_fun','Value'];

## report variables
rseq_rep_title_short = gsub(';',',',analysis.config[analysis.config$Name=='rseq_rep_title_short','Value']);
rseq_rep_title_long = gsub(';',',',analysis.config[analysis.config$Name=='rseq_rep_title_long','Value']);
study_design_text = gsub(';',',',analysis.config[analysis.config$Name=='study_design_text','Value']);
experiment_text = gsub(';',',',analysis.config[analysis.config$Name=='experiment_text','Value']);
gene_filter_method =  gsub(';',',',analysis.config[analysis.config$Name=='gene_filter_fun','Value']);

if(gene_filter_method =='filterLowExpressedGenesMax'){
	gene_filter_method='maximum';
} else if(gene_filter_method =='filterLowExpressedGenesAvg'){
	gene_filter_method='average';
}

pvclust_cor_name = if(pvclust.dist.fun=='pearsonucdist'){'uncentered pearson correlation'} else if(pvclust.dist.fun=='pearsondist'){'pearson correlation'}else if (pvclust.dist.fun=='spearmandist'){'spearman correlation'}
ensemb_genes_remove = gsub('_','-',unlist(str_split(analysis.config[analysis.config$Name=='ensembl_gene_types_remove','Value'],';')))
aligner_and_version = paste('\\textit{Hisat2} splice-aware read aligner (Version ',system("hisat2 --version | head -1 | grep -o -P ' [0-9.]+'",intern=T),sep='')


##################################################################################
##
## SPECIFY OUTLIERS (OUTLIERS ARE REMOVED BEFORE NORMALIZAION)
##
##################################################################################

outliers = unlist(str_split(analysis.config[analysis.config$Name=='outliers','Value'],';'));

##################################################################################
##
## ANALYSIS VARIABLES
##
##################################################################################

## normalization
nrmFlags = c('not','tmm');
nrmLabls = c('Before TMM Normalization','After TMM Normalization');

## filtering
fltFlags = c('filtered','unfiltered');
fltLabls = c('Filtered Genes','All Genes');

## standardization
stdFlags = c('org','std');
stdLabls = c('Original Variables','Standardized Variables');

## Specimen type variables
spcFlags = unique(mta$spct);
spcLabls = unique(mta$spctl);

## Vaccine/Treatment type variables
trtFlags = unique(mta$trt);
trtLabls = unique(mta$trtl);

## All time
times = unique(mta$time)
times.l = unique(mta$timel)

## treatment time
b.times = unique(mta[mta$timeb==T,'time']);
b.times.l = unique(mta[mta$timeb==T,'timel']);

## post treatment/treatment
postb.times   = unique(mta[mta$timeb==F,'time'])
postb.timesl   = unique(mta[mta$timeb==F,'timel'])
postb.times.m = c(paste('tp', unique(mta[mta$timeb==F,'time']),sep=''),'posttp');
postb.times.l = c(unique(mta[mta$timeb==F,'timel']),'All post-treatment time points');

## Gene set labels
geneset.file = paste(dta.dir,'/gene_sets/all.ensembl.tab.gz',sep='')
geneset.types = system(paste('zcat ',geneset.file," | awk -F '\t' '{print $1}' | sort -u",sep=''),intern=T); # much faster than opening file in R

## KnitR figure caption labels
times.knitr.rep = ceiling(length(times)/5)
times.knitr.nrm.lab = if(length(times)>5) {paste(1:ceiling(length(times)/5),' of ',ceiling(length(times)/5),' ',sep='')} else {''}
postb.times.knitr.ma.lab = if(length(postb.times)>10) {paste(1:ceiling(length(postb.times)/10),' of ',ceiling(length(postb.times)/10),' ',sep='')} else {''}
postb.times.knitr.venn.glm.lab = if(length(postb.times.m)>5) {paste(1:ceiling(length(postb.times.m)/5),' of ',ceiling(length(postb.times.m)/5),' ',sep='')} else {''}
postb.times.knitr.venn.gsea.lab = if(length(postb.times.m)>10) {paste(1:ceiling(length(postb.times.m)/10),' of ',ceiling(length(postb.times.m)/10),' ',sep='')} else {''}
trt.knitr.venn.glm.lab = if(length(trtFlags)>5) {paste(1:ceiling(length(trtFlags)/5),' of ',ceiling(length(trtFlags)/5),' ',sep='')} else {''}
trt.knitr.venn.gsea.lab = if(length(trtFlags)>10) {paste(1:ceiling(length(trtFlags)/10),' of ',ceiling(length(trtFlags)/10),' ',sep='')} else {''}

## upset plot labels
upsetLablsSpcSdeg = c()
upsetLablsTimeSdeg = c()
upsetLablsTrtSdeg = c()
upsetLablsSpcGsea = c()
upsetLablsTimeGsea = c()
upsetLablsTrtGsea = c()

if(rdFlag==TRUE) {
	bpLabel = c('genome','gene model')
} else {
	bpLabel = c('genome')
}
if (length(spcFlags)>1) {
	spcLablsBias = c('All Specimen Types',spcLabls)
} else {
	spcLablsBias = c(spcLabls)
}
HeatmapLablSdeg = c()
TrendLabls = c()
DendrLabls = c()
HeatmapLablGsea = c()
RadarLablGsea = c()


##################################################################################
##
## CUSTOM FUNCTIONS
##
##################################################################################

resetPar <- function() {
	dev.new()
	op <- par(no.readonly = TRUE)
	dev.off()
	op
}

write.tbl = function(x,rownames=F) {
	
	## remove latex characters/labels from first row of DF
	col.names = colnames(x)
	x = as.data.frame(t(apply(x,1,gsub,pattern='\\\\',replacement='')),stringsAsFactors=F) # this results in loss of colnames sometimes
	colnames(x) = col.names
	
	## remove latex characters/labels from columns
	colnames(x) = gsub('\\\\%','%',colnames(x));
	colnames(x) = gsub('\\\\#','#',colnames(x));
	colnames(x) = gsub('\\$Log_\\{2\\}\\$','Log2',colnames(x));
	colnames(x) = gsub('\\\\newline','',colnames(x));
	colnames(x) = gsub('\\\\','',colnames(x));
	colnames(x) = gsub("\\$",'',colnames(x));
	
	if(rownames==T) {
		## remove latex characters/labels from rows
		rownames(x) = gsub('\\\\%','%',rownames(x));
		rownames(x) = gsub('\\\\#','#',rownames(x));
		rownames(x) = gsub("\\$",'',rownames(x));
	}
	
	# if the tables directory does not exist, create it.
	if (!dir.exists(tbl.dir)) {
		system(paste('mkdir',tbl.dir))
	}
	
	tbl.no = length(list.files(tbl.dir))+1;
	xtable.out.file = paste(tbl.dir,'/table_A',tbl.no,'.csv',sep='');
	if(rownames==F) {
		write.table(x,xtable.out.file,quote=F,row.names=F,sep=',',col.names=T);
	} else if(rownames==T) {
		write.table(x,xtable.out.file,quote=F,row.names=T,sep=',',col.names=NA);
	}
	R.utils::gzip(xtable.out.file,overwrite=TRUE);
}

## function to convert ANOVA summary object to data frame
aov.as.dataframe = function(x) {
	
	old.scip = options()$scipen;
	old.dig  = options()$digits;
	
	if(length(x) == 1) {
		res = as.data.frame(x[[1]])
	} else {
		res = lapply(unlist(x, FALSE), as.data.frame)
	}
	
	## update p_values
	if(length(which(res[,5]>= 0.0001))>0) {
		res[which(res[,5]>= 0.0001),5] = round(res[which(res[,5]>= 0.0001),5],4);
	}
	if(length(which(res[,5]< 0.0001))>0) {
		res[which(res[,5]< 0.0001),5] = '<0.0001';
	}
	
	## round values
	res[,2] = round(res[,2],2);
	res[,3] = round(res[,3],2);
	res[,4] = round(res[,4],2);
	
	options("scipen"=100, "digits"=4);
	res[is.na(res)] ='';
	options("scipen"=old.scip, "digits"=old.dig);
	
	return(res);
}

## calculate spearman correlation distance 
spearmandist = function(x,by.row=T,...) {
	if (by.row==T) {
		x  <- as.matrix(t(x))
	} else {
		x  <- as.matrix(x)
	}
	res = 1-cor(x,method='spearman');
	res <- as.dist(res)
	attr(res, "method") <- "1-Spearman Cor."
	return(res)
};

## calculate pearson correlation distance 
pearsondist = function(x,by.row=T,...) {
	if (by.row==T) {
		x  <- as.matrix(t(x))
	} else {
		x  <- as.matrix(x)
	}
	res = 1-cor(x,method='pearson');
	res <- as.dist(res)
	attr(res, "method") <- "1-Pearson Cor."
	return(res)
};

## calculate uncentered pearson correlation distance 
## between rows of x (pvclust)
pearsonucdist = function(x,by.row=T,...) {
	if(sum(is.na(x)) > 0){
		x <- na.omit(x)
		warning("Rows including NAs were omitted")
	}
	if (by.row==T) {
		x  <- as.matrix(t(x))
	} else {
		x  <- as.matrix(x)
	}
	P  <- crossprod(x)
	qq <- matrix(diag(P),ncol=ncol(P))
	Q  <- sqrt(crossprod(qq))
	res <- as.dist(1 - P/Q)
	attr(res,"method") <- "1-Uncentered Pearson Cor."
	return(res)
};

## define distance matrices
hclust2  		= function(x, method='complete') hclust(x, method=method)
vegdist2 		= function(x, method='euclidean') vegdist(x, method=method);

makeTransparent=function(someColor, alpha=60)
{
	newColor=col2rgb(someColor)
	apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
						blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

## function to color dendogram
colDendo = function(hcd,col,label) {
	mycols = col;
	mylabs = label
	i = 0 
	colLab <- function(n) {
		if(is.leaf(n)) {
			i <<- i + 1
			a <- attributes(n)
			attr(n, "nodePar") =
					c(a$nodePar, list(lab.col = mycols[i],lab.bg='grey50',lab.cex=0.5))
			attr(n, "frame.plot") = TRUE;
			attr(n, "label") = mylabs[i];
		}
		n
	}
	return(dendrapply(hcd, colLab));
} 


## calculate jaccard index
jaccardIdx = function(x,y) {
	n1 = length(x)
	n2 = length(y)
	num.int = length(intersect(x,y))
	num.union = length(union(x,y))
	if (num.int==0) {
		return(0)
	} else {
		return(num.int/num.union)
	}
}

## retain genes in the cpm.matrix that are part of the clean.genes file
filterUnwantedGenes = function(cpm.matrix,clean.genes) {
	return(cpm.matrix[rownames(cpm.matrix) %in% clean.genes$gene_id,]);
}

## save genes which have at least one sample value > CPM cutoff
filterLowExpressedGenesMax = function(flt.spc,flt.mta,flt.cpm,save) {
	
	## get gene IDs that have at least one value larger than the the cut offs
	alltp.gen = sort(rownames(flt.cpm)[which((apply(flt.cpm,1,function(x){max(x)})>flt.lcpm.cut))]);
	
	## save gene id list of genes meeting the CPM cutoff criteria
	if(save==T) {
		out.file.gene.set = paste(out.dir.lcpm,'/',tolower(flt.spc),'_posttp_analysis_gene_set.tab',sep='');
		write.table(alltp.gen,out.file.gene.set,sep='\t',quote=F,row.names=F,col.names = F);
		R.utils::gzip(out.file.gene.set,overwrite=TRUE);
	}
	## return data matrix with genes passing filter criteria
	return(flt.cpm[rownames(flt.cpm) %in% alltp.gen,]);
}

## save genes which the mean value across samples > CPM cutoff
filterLowExpressedGenesAvg = function(flt.spc,flt.mta,flt.cpm,save) {
	
	## get gene IDs that have at least one value larger than the the cut offs
	alltp.gen = sort(rownames(flt.cpm)[which((apply(flt.cpm,1,function(x){mean(x)})>flt.lcpm.cut))]);
	
	## save gene id list of genes meeting the CPM cutoff criteria
	if(save==T) {
		out.file.gene.set = paste(out.dir.lcpm,'/',tolower(flt.spc),'_posttp_analysis_gene_set.tab',sep='');
		write.table(alltp.gen,out.file.gene.set,sep='\t',quote=F,row.names=F,col.names = F);
		R.utils::gzip(out.file.gene.set,overwrite=TRUE);
	}
	## return data matrix with genes passing filter criteria
	return(flt.cpm[rownames(flt.cpm) %in% alltp.gen,]);
}

plotPcaByCell = function(mta,in.dir.pca,spcFlag,spcLabl,visFlag,visLabl) {
	
	## read the data
	filename = paste(in.dir.pca,'/',spcFlag,'_',visFlag,'_pca_tmm_normalized_std.RData',sep='');
	load(file=filename);
	
	## get ids and align metadata;
	ids = rownames(pca.res$x);
	mta.pca=mta[match(ids,mta$samid),]
	
	## get first and second principal component
	c1 = pca.res$x[,1];
	c2 = pca.res$x[,2];
	
	out.ids = which(mta.pca$samid %in% outliers);
	out.x = c1[out.ids];
	out.y = c2[out.ids];
	
	n=length(c1);
	
	xlab = paste('PC1 (',pca.res$pvar[1],'%)',sep='');
	ylab = paste('PC2 (',pca.res$pvar[2],'%)',sep='');	
	
	dataEllipse(c1,c2,grid=T,lwd=1,levels = c(1-alpha),xlab='',ylab='',
			groups=as.factor(mta.pca$spct),col=unique(mta.pca$spctc),
			group.labels=c(),cex.axis=0.9,cex.lab=0.8,cex=0.4,id=list(method="mahal",n=3,cex=0.4,col=1,location="lr",labels=ids),robust=T);
	mtext(xlab,side=1,line=2.2,cex=0.7)
	mtext(ylab,side=2,line=2.2,cex=0.7)
	mtext(paste(spcLabl,'\n(PCA ',stdLabls[2],', n=',n,')',sep=''),cex=0.8,line=0.4);
	points(out.x,out.y,col='blue',lwd=1.5,cex=1.3);
}

plotMdsByCell = function(mta,in.dir.dist,type,spcFlag,spcLabl,visFlag,visLabl,text) {
	
	if(type=='euc') {
		file.name = paste(in.dir.dist,'/',spcFlag,'_',visFlag,'_euc_dist_tmm_normalized_std.RData',sep='');
		
	} else if(type=='spm') {
		file.name = paste(in.dir.dist,'/',spcFlag,'_',visFlag,'_spm_dist_tmm_normalized_std.RData',sep='');
		
	}
	
	load(file=file.name);
	
	ids=attr(dist.res,'labels');
	mta.mds=mta[match(ids,mta$samid),]
	
	mds=isoMDS(dist.res,trace=F);
	stress = round(mds$stress,1);
	
	mds$points = mds$points
	
	out.ids = which(mta.mds$samid %in% outliers);
	out.x = mds$points[out.ids,1];
	out.y = mds$points[out.ids,2];
	
	dataEllipse(mds$points,grid=T,lwd=1,levels = c(),xlab="",ylab="",
			groups=as.factor(mta.mds$spct),col=unique(mta.mds$spctc),id=list(method="mahal",n=3,cex=0.4,col=1,location="lr",labels=ids),
			group.labels=c(),cex.axis=0.9,cex.lab=0.8,cex=0.4);
	mtext("MDS Coordinate 1",side=1,line=2.2,cex=0.7)
	mtext("MDS Coordinate 2",side=2,line=2.2,cex=0.7)
	mtext(paste(spcLabl,'\n(',text,', Stress=',stress,'%, n=',nrow(mds$points),')',sep=''),cex=0.75,line=0.4);
	points(out.x,out.y,col='blue',lwd=1.5,cex=1.3);
}

## call radarplot function modified from
## (http://www.statisticstoproveanything.com/2013/11/spider-web-plots-in-r.html)
radarplot = function(data, data.row = NULL, y.cols = NULL, main = NULL, add = F, 
		col = "red", lty = 1, scale = F,max=5,increment=0.5,labels = NULL,txt.cex=1.2,cex.main=3,grey.cir = NULL,lwd=1.2,last=FALSE) {
	if (!is.matrix(data) & !is.data.frame(data)) 
		stop("Requires matrix or data.frame")
	if (is.null(y.cols)) 
		y.cols = colnames(data)[sapply(data, is.numeric)]
	if (sum(!sapply(data[, y.cols], is.numeric)) > 0) {
		out = paste0("\"", colnames(data)[!sapply(data, is.numeric)], "\"", 
				collapse = ", ")
		stop(paste0("All y.cols must be numeric\n", out, " are not numeric"))
	}
	if (is.null(data.row)) 
		data.row = 1
	if (is.character(data.row)) 
		if (data.row %in% rownames(data)) {
			data.row = which(rownames(data) == data.row)
		} else {
			stop("Invalid value for data.row:\nMust be a valid rownames(data) or row-index value")
		}
	if (is.null(main)) 
		main = rownames(data)[data.row]
	if (scale == T) {
		data = scale(data[, y.cols])
		data = apply(data, 2, function(x) x/max(abs(x)))
	}
	max = (ceiling(max*(1/increment)) * increment)+1
	data = as.data.frame(data)
	n.y = length(y.cols)
	min.rad = 360/n.y
	polar.vals = (90 + seq(0, 360, length.out = n.y + 1)) * pi/180
	
	if (add == F) {
		plot(0, xlim = c((-max*1.1), max*1.3), ylim = c((-max*1.1), max*1.3), type = "n", axes = F, 
				xlab = "", ylab = "")
		mtext(main,line=1.3,cex=cex.main)
		lapply(polar.vals, function(x) lines(c(0, max * cos(x)), c(0, max * sin(x))))
		lapply(1:n.y, function(x) text((max*1.27) * cos(polar.vals[x]), (max*1.27) * sin(polar.vals[x]), 
							labels[x], cex = txt.cex, font=2))
		
		lapply(seq(increment,max,by=increment), function(x) lines(x * cos(seq(0, max * pi, length.out = 100)), 
							x * sin(seq(0, max * pi, length.out = 100)), lwd = 0.5, lty = 2, col = "gray60"))
		lines(cos(seq(0, max * pi, length.out = 100)), sin(seq(0, max * pi, length.out = 100)), 
				lwd = 1.2, col = "gray50")
		symbols(0, 0, circles=1, inches=F,bg="white",add=T)
		if(!is.null(grey.cir)) {
			symbols(0, 0, circles=grey.cir+1, inches=F,bg=makeTransparent("grey30"),add=T)
		}
	}
	
	r = 1 + data[data.row, y.cols]
	xs = r * cos(polar.vals)
	ys = r * sin(polar.vals)
	xs = c(xs, xs[1])
	ys = c(ys, ys[1])
	
	lines(xs, ys, col = col, lwd = lwd, lty = lty)
	
	if (last == T) {
		shadowtext <- function(x, y=NULL, labels, col='black', bg='white', 
				theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
			
			xy <- xy.coords(x,y)
			xo <- r*strwidth('A')
			yo <- r*strheight('A')
			
			## draw background text with small shift in x and y in background colour
			for (i in theta) {
				text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
			}
			## draw actual text in exact xy position in foreground colour
			text(xy$x, xy$y, labels, col=col, ... )
		}
		
		shadowtext(x=rep(-0.24,length(seq(1,max,by=increment)+((max-1)/10))),y=seq(1,max,by=increment)+((max-1)/20),labels=seq(0,max-1,by=increment),cex=0.8)
	}
}


# custom version of upset that loows for the printing of empty sets -- 
# the only real change is that "FindStartEnd(data)" was replaced with "1:ncol(data)"
upset.custom = function (data, nsets = 5, nintersects = 40, sets = NULL, keep.order = T, 
		set.metadata = NULL, intersections = NULL, matrix.color = "gray23", 
		main.bar.color = "gray23", mainbar.y.label = "Intersection Size", 
		mainbar.y.max = NULL, sets.bar.color = "gray23", sets.x.label = "Set Size", 
		point.size = 2.2, line.size = 0.7, mb.ratio = c(0.7, 0.3), 
		expression = NULL, att.pos = NULL, att.color = main.bar.color, 
		order.by = "freq", decreasing = c(T, F), show.numbers = "yes", 
		number.angles = 0, group.by = "degree", cutoff = NULL, queries = NULL, 
		query.legend = "none", shade.color = "gray88", shade.alpha = 0.25, 
		matrix.dot.alpha = 0.5, empty.intersections = NULL, color.pal = 1, 
		boxplot.summary = NULL, attribute.plots = NULL, scale.intersections = "identity", 
		scale.sets = "identity", text.scale = 1, set_size.angles = 0) 
{
	startend <- 1:ncol(data)
	first.col <- startend[1]
	last.col <- startend[2]
	if (color.pal == 1) {
		palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", 
				"#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", 
				"#17BECF")
	} else {
		palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
				"#0072B2", "#D55E00", "#CC79A7")
	}
	if (is.null(intersections) == F) {
		Set_names <- unique((unlist(intersections)))
		Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
		New_data <- UpSetR:::Wanted(data, Sets_to_remove)
		Num_of_set <- UpSetR:::Number_of_sets(Set_names)
		if (keep.order == F) {
			Set_names <- UpSetR:::order_sets(New_data, Set_names)
		}
		All_Freqs <- UpSetR:::specific_intersections(data, first.col, 
				last.col, intersections, order.by, group.by, decreasing, 
				cutoff, main.bar.color, Set_names)
	} else if (is.null(intersections) == T) {
		Set_names <- sets
		if (is.null(Set_names) == T || length(Set_names) == 0) {
			Set_names <- UpSetR:::FindMostFreq(data, first.col, last.col, 
					nsets)
		}
		Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
		New_data <- UpSetR:::Wanted(data, Sets_to_remove)
		Num_of_set <- UpSetR:::Number_of_sets(Set_names)
		if (keep.order == F) {
			Set_names <- UpSetR:::order_sets(New_data, Set_names)
		}
		All_Freqs <- UpSetR:::Counter(New_data, Num_of_set, first.col, 
				Set_names, nintersects, main.bar.color, order.by, 
				group.by, cutoff, empty.intersections, decreasing)
	}
	Matrix_setup <- UpSetR:::Create_matrix(All_Freqs)
	labels <- UpSetR:::Make_labels(Matrix_setup)
	att.x <- c()
	att.y <- c()
	if (is.null(attribute.plots) == F) {
		for (i in seq_along(attribute.plots$plots)) {
			if (length(attribute.plots$plots[[i]]$x) != 0) {
				att.x[i] <- attribute.plots$plots[[i]]$x
			} else if (length(attribute.plots$plots[[i]]$x) == 
					0) {
				att.x[i] <- NA
			}
			if (length(attribute.plots$plots[[i]]$y) != 0) {
				att.y[i] <- attribute.plots$plots[[i]]$y
			} else if (length(attribute.plots$plots[[i]]$y) == 
					0) {
				att.y[i] <- NA
			}
		}
	}
	BoxPlots <- NULL
	if (is.null(boxplot.summary) == F) {
		BoxData <- UpSetR:::IntersectionBoxPlot(All_Freqs, New_data, first.col, 
				Set_names)
		BoxPlots <- list()
		for (i in seq_along(boxplot.summary)) {
			BoxPlots[[i]] <- UpSetR:::BoxPlotsPlot(BoxData, boxplot.summary[i], 
					att.color)
		}
	}
	customAttDat <- NULL
	customQBar <- NULL
	Intersection <- NULL
	Element <- NULL
	legend <- NULL
	EBar_data <- NULL
	if (is.null(queries) == F) {
		custom.queries <- UpSetR:::SeperateQueries(queries, 2, palette)
		customDat <- UpSetR:::customQueries(New_data, custom.queries, 
				Set_names)
		legend <- UpSetR:::GuideGenerator(queries, palette)
		legend <- UpSetR:::Make_legend(legend)
		if (is.null(att.x) == F && is.null(customDat) == F) {
			customAttDat <- UpSetR:::CustomAttData(customDat, Set_names)
		}
		customQBar <- UpSetR:::customQueriesBar(customDat, Set_names, 
				All_Freqs, custom.queries)
	}
	if (is.null(queries) == F) {
		Intersection <- UpSetR:::SeperateQueries(queries, 1, palette)
		Matrix_col <- UpSetR:::intersects(QuerieInterData, Intersection, 
				New_data, first.col, Num_of_set, All_Freqs, expression, 
				Set_names, palette)
		Element <- UpSetR:::SeperateQueries(queries, 1, palette)
		EBar_data <- UpSetR:::ElemBarDat(Element, New_data, first.col, 
				expression, Set_names, palette, All_Freqs)
	} else {
		Matrix_col <- NULL
	}
	Matrix_layout <- UpSetR:::Create_layout(Matrix_setup, matrix.color, 
			Matrix_col, matrix.dot.alpha)
	Set_sizes <- UpSetR:::FindSetFreqs(New_data, first.col, Num_of_set, 
			Set_names, keep.order)
	Bar_Q <- NULL
	if (is.null(queries) == F) {
		Bar_Q <- UpSetR:::intersects(QuerieInterBar, Intersection, New_data, 
				first.col, Num_of_set, All_Freqs, expression, Set_names, 
				palette)
	}
	QInter_att_data <- NULL
	QElem_att_data <- NULL
	if ((is.null(queries) == F) & (is.null(att.x) == F)) {
		QInter_att_data <- UpSetR:::intersects(QuerieInterAtt, Intersection, 
				New_data, first.col, Num_of_set, att.x, att.y, expression, 
				Set_names, palette)
		QElem_att_data <- UpSetR:::elements(QuerieElemAtt, Element, New_data, 
				first.col, expression, Set_names, att.x, att.y, palette)
	}
	AllQueryData <- UpSetR:::combineQueriesData(QInter_att_data, QElem_att_data, 
			customAttDat, att.x, att.y)
	ShadingData <- NULL
	if (is.null(set.metadata) == F) {
		ShadingData <- UpSetR:::get_shade_groups(set.metadata, Set_names, 
				Matrix_layout, shade.alpha)
		output <- UpSetR:::Make_set_metadata_plot(set.metadata, Set_names)
		set.metadata.plots <- output[[1]]
		set.metadata <- output[[2]]
		if (is.null(ShadingData) == FALSE) {
			shade.alpha <- unique(ShadingData$alpha)
		}
	}
	if (is.null(ShadingData) == TRUE) {
		ShadingData <- UpSetR:::MakeShading(Matrix_layout, shade.color)
	}
	Main_bar <- suppressMessages(UpSetR:::Make_main_bar(All_Freqs, Bar_Q, 
					show.numbers, mb.ratio, customQBar, number.angles, EBar_data, 
					mainbar.y.label, mainbar.y.max, scale.intersections, 
					text.scale, attribute.plots))
	Matrix <- UpSetR:::Make_matrix_plot(Matrix_layout, Set_sizes, All_Freqs, 
			point.size, line.size, text.scale, labels, ShadingData, 
			shade.alpha)
	Sizes <- UpSetR:::Make_size_plot(Set_sizes, sets.bar.color, mb.ratio, 
			sets.x.label, scale.sets, text.scale, set_size.angles)
	UpSetR:::Make_base_plot(Main_bar, Matrix, Sizes, labels, mb.ratio, 
			att.x, att.y, New_data, expression, att.pos, first.col, 
			att.color, AllQueryData, attribute.plots, legend, query.legend, 
			BoxPlots, Set_names, set.metadata, set.metadata.plots)
}

##################################################################################
##
## Set Knitr Configuration options
##
##################################################################################

## knitr configuration: https://github.com/lazappi/lazappi/blob/master/R/rmarkdown.R
opts_chunk$set(
		## CODE EVALUATION
		eval           = TRUE,         # Whether to evaluate the chunk
		
		## TEXT RESULTS
		collapse       = TRUE,         # Collapse results to single block
		echo           = FALSE,        # Whether to include source
		error          = FALSE,        # Continue after error
		include        = TRUE,         # Whether to include ouput
		message        = FALSE,        # Whether to include messages
		split          = FALSE,        # Split output into multiple files
		strip.white    = TRUE,         # Remove start/end white lines
		warning        = FALSE,        # Whether to include warnings
		
		## CODE DECORATION
		background     = "#F7F7F7",     # Chunk background in LaTeX
		comment        = "##",         # Prefix before output
		highlight      = TRUE,         # Highlight source
		indent         = "  ",         # Indent for markdown output
		prompt         = FALSE,        # Whether to add prompt character
		size           = "normalsize", # Font size for LaTeX output
		tidy           = FALSE,        # Tidy source using formatR
		
		## CACHE
		autodep        = TRUE,         # Automatically set dependencies
		cache          = TRUE,        # Whether to cache output
		cache.comments = FALSE,        # Comments change cache
		cache.lazy     = TRUE,         # Whether to lazy load cache
		cache.path     = "cache/",     # Prefix used for cache files
		cache.rebuild  = FALSE,        # Rebuild cache
		cache.vars     = NULL,         # Variable to save in cache
		
		## PLOTS
		## dev           =               # Plotting device
		dev.args       = NULL,         # Arguments to device
		dpi            = 300,           # Dots Per Inch
		external       = TRUE,         # Whether to externalise tikz
		fig.align      = "center",    # Output figure alignments
		fig.cap        = NULL,         # Figure caption for LaTeX
		fig.env        = "figure",     # LaTeX environment for figures
		fig.ext        = NULL,         # Extension for figure output
		fig.height     = 7,            # Figure height, in inches
		fig.keep       = "high",       # How plots are kept
		fig.lp         = "fig:",       # Prefix for figure labels
		fig.path       = "figs/",      # Prefix for figure filenames
		fig.pos        = "H",           # Figure position arrangement
		fig.process    = NULL,         # Function to post-process figure file
		fig.retina     = 1,            # Adjustment for retina displays
		fig.scap       = NULL,         # Short caption for LaTeX
		fig.show       = "asis",       # How to show plots
		fig.showtext   = NULL,         # Whether to call show.text before plot
		fig.subcap     = NULL,         # Captions for subfigures
		fig.width      = 7,            # Figure width, in inches
		out.extra      = NULL,         # Extra options for figure output
		out.height     = NULL,         # Figure height in output file
		out.width      = NULL,         # Figure width in output file
		resize.height  = NULL,         # Resize height for tikz output in LaTeX 
		resize.width   = NULL,         # Resize width for tikz output in LaTeX
		sanitise       = FALSE,        # Whether to sanitise tikz
		crop 		   = TRUE,		   # Crop figures.  Crop type specified at end of script
		
		## TABLES
		results        = "markup",     # How to display results	
		table.placement= 'H',		 # Where to place tables 
		
		## ANIMATION
		aniopts        = "controls, loop", # Extra options of animations
		interval       = 1,            # Seconds between animation frames
		
		## CODE CHUNK
		code           = NULL,         # Override code in chunk
		ref.label      = NULL,         # Inherit code from other chunks
		
		## CHILD DOCUMENTS
		child          = NULL,         # Filenames of child documents to input
		
		## LANGUAGE ENGINES
		engine         = "R",          # Language of the chunk
		
		## OPTION TEMPLATES
		opts.label     = NULL,         # The label of options in opts_template
		
		## EXTRACT SOURCE
		purl           = TRUE,         # Whether to include chunk in purl
		
		## OTHER
		R.options      = NULL,          # Local R options
		base.dir	   = paste(res.dir,"/../report",sep='') # Where to put the tex file
)

## PDFcrop
knit_hooks$set(crop=hook_pdfcrop);

## assign functions from name
pvclust.dist.fun = get(pvclust.dist.fun);
heatmap.dist.fun = get(heatmap.dist.fun);

##################################################################################
##
## OTHER MODULES
##
##################################################################################

## load custom pvclust code to allow additional distance matrices
## http://www.is.titech.ac.jp/~shimo/prog/pvclust/
source(paste(src.dir,'/r/functions/pvclust-internal.R',sep=''));
source(paste(src.dir,'/r/functions/pvclust-max-dist.R',sep=''));
source(paste(src.dir,'/r/functions/venn-diagram-mod.R',sep=''));