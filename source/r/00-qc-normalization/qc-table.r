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
# Program:  qc-table.r
# Version:  RSEQREP 2.1.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Evaluate correlation between principal components and technical/biological variables in table form
# Input:    N/A
# Output:  	N/A
#############################################################################################################

source('../r/init-analysis.r')

## update specimen flag to inclue a variable across all
## types if there is more than 1 specimen type
if (length(spcFlags)>1) {
	spcFlags = c('all',spcFlags)
	spcLabls = c('All Specimen Types',spcLabls)
}

## if read distribution is run, plot results, otherwise omit
if (rdFlag==T) {
	varFlags = c('tec','gen');
	varLabls = c('genome','gene model');
} else {
	varFlags = c('tec');
	varLabls = c('genome');
}

## define technical variables
tec.vars = c(
		'rseqc.qc.total',
		'rseqc.qc.unmapped',
		'rseqc.qc.unique',
		'rseqc.qc.unique_perc',
		'feature.counts.total',
		'rseqc.qc.plus_strand',
		'rseqc.qc.minus_strand',
		'rseqc.gc.median',
		'rseqc.gc.mean');
tec.lab = c(
		'Total Mapped Reads [$10^6$]',
		'Unmapped Reads [$10^6$]',
		'Uniquely Mapped Reads [$10^6$]',
		'Uniquely Mapped Reads [\\%]',
		'Counted Fragments [$10^6$]',
		'Uniquely Mapped Reads + Strand [$10^6$]',
		'Uniquely Mapped Reads - Strand [$10^6$]',
		'Median GC [\\%]',
		'Mean GC [\\%]');

## define gene model distribution variables 
gen.vars = c(
		'rseqc.rc.perc_exon_tags',
		'rseqc.rc.perc_intron_tags',
		'rseqc.rc.perc_intergenic_tags',
		'rseqc.rc.cds_exons_tags_kb',
		'rseqc.rc.intron_tags_kb',
		'rseqc.rc.tp_utr_exons_tags_kb',
		'rseqc.rc.fp_utr_exons_tags_kb',
		'rseqc.rc.intergenic_up_tags_kb',
		'rseqc.rc.intergenic_down_tags_kb',
		'rseqc.jc.total_splicing_events',
		'rseqc.jc.total_splicing_junctions');
gen.lab = c(
		"Exon Tags [\\%]",
		"Intron Tags [\\%]",
		"Intergenic Tags [\\%]",
		"CDS Exons Tags Per Kb",
		"Intron Tags Per Kb",
		"3' UTR Tags Per Kb",
		"5' UTR Tags Per Kb",
		"TSS upstream 10Kb Tags Per Kb",
		"TES downstream 10Kb Tags Per Kb",
		'Splicing Events [$10^5$]',
		'Splicing Junctions [$10^5$]');

## Update fields -- convert to more reasonable units
mta$rseqc.qc.unique_perc = round(mta$rseqc.qc.unique/mta$rseqc.qc.total,3)*100;
mta$rseqc.qc.total=mta$rseqc.qc.total/10^6;
mta$rseqc.qc.unmapped=mta$rseqc.qc.unmapped/10^6;
mta$rseqc.qc.unique=mta$rseqc.qc.unique/10^6;
mta$feature.counts.total=mta$feature.counts.total/10^6;
mta$rseqc.qc.plus_strand=round(mta$rseqc.qc.plus_strand/10^6,2);
mta$rseqc.qc.minus_strand=round(mta$rseqc.qc.minus_strand/10^6,2);
mta$rseqc.jc.total_splicing_events=mta$rseqc.jc.total_splicing_events/10^5;
mta$rseqc.jc.total_splicing_junctions=mta$rseqc.jc.total_splicing_junctions/10^5;

for(v in 1:length(varFlags)) {
	
	varFlag = varFlags[v];
	varLabl = varLabls[v];
	
	if(varFlag =='tec') {
		var.sel = tec.vars;
		lab.sel = tec.lab;
	} else if (varFlag == 'gen') {
		var.sel = gen.vars;
		lab.sel = gen.lab;
	}
	
	for (s in 1:length(spcFlags)) {
		spcFlag = spcFlags[s];
		spcLabl = spcLabls[s];
		
		if (spcFlag == 'all') {
			mta.spc = mta
		} else {
			mta.spc = mta[mta$spct==spcFlag,]
		}
	
		cols = var.sel;
		res = as.data.frame(matrix(nrow=length(var.sel),ncol=0),stringsAsFactors=F);
		
		res$min =   round(apply(mta.spc[,cols],2,min),2);
		res$q1=     round(apply(mta.spc[,cols],2,quantile,0.25),2);
		res$median= round(apply(mta.spc[,cols],2,median),2);
		res$mean=   round(apply(mta.spc[,cols],2,mean),2);
		res$q3=     round(apply(mta.spc[,cols],2,quantile,0.75),2);
		res$max =   round(apply(mta.spc[,cols],2,max),2);
		res$sd=     round(apply(mta.spc[,cols],2,sd),2);
		res$mad=    round(apply(mta.spc[,cols],2,mad),2);
		res$n=      nrow(mta.spc)
		
		res = res[,c('min','q1','median','mean','q3','max','sd','mad','n')]
		colnames(res) = c('Min','Q1','Median','Mean','Q3','Max','SD','MAD','N');
		rownames(res) =  lab.sel; 
		
		write.tbl(res,rownames=T);
		
		## specify caption
		if(exists('reportAssayFlag')){
			assay.label = ifelse(reportAssayFlag=='rna_seq','RNA-Seq','RP')
			caption = paste0('Summary human reference ',varLabls[v],' alignment statistics for study samples (',assay.label,', ',spcLabl,')')
		}else{
			caption = paste0('Summary human reference ',varLabls[v],' alignment statistics for study samples (',spcLabl,')')
		}
		
		## print latex table
		x = xtable(res,caption=caption,align=c('p{2.1in}','p{0.35in}','p{0.35in}','p{0.38in}','p{0.38in}','p{0.35in}','p{0.38in}','p{0.3in}','p{0.3in}','p{0.3in}'),digits=2,
				label=paste('tab:stats_',varFlag,'_',spcFlag,sep=''));
		cat('\\*');cat('\\nopagebreak');
		print(x,tabular.environment='longtable',
				floating=F,
				include.rownames=T,
				hline.after = c(-1,nrow(x)),
				size= 'footnotesize',
				add.to.row = list(pos = list(0),command = "\\hline \\endhead "),sanitize.text.function=function(x){x})
	}
}