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
# Program:  qc-boxplots.r 
# Version:  RSEQREP 2.1.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Evaluate correlation between principal components and technical/biological variables
# Input:    N/A
# Output:  	N/A
#############################################################################################################

options(scipen = 50) 
source('../r/init-analysis.r')

## if read distribution is run, plot results, otherwise omit
if (rdFlag==T) {
	varFlags = c('tec','gen');
	varLabls = c('genome','gene model');
} else {
	varFlags = c('tec');
	varLabls = c('genome');
}

cor.method = 'spearman';

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
		expression('Total Mapped Reads [10'^6*']'),
		expression('Unmapped Reads [10'^6*']'),
		expression('Uniquely Mapped Reads [10'^6*']'),
		'Uniquely Mapped Reads [%]',
		expression('Counted Fragments [10'^6*']'),
		expression('Uniquely Mapped Reads + Strand [10'^6*']'),
		expression('Uniquely Mapped Reads - Strand [10'^6*']'),
		'Median %GC',
		'Mean %GC');

## if all results are paired end, add metrics
if (all(!is.na(mta$fastq_file_2))) {
	tec.vars = c(tec.vars,'rseqc.qc.proper_pairs','rseqc.qc.non_proper_pairs')
	tec.lab = c(tec.lab,expression('Proper Pair Mapped Reads [10'^6*']'),expression('Non-Proper Pair Mapped Reads [10'^6*']'))
	mta$rseqc.qc.proper_pairs=mta$rseqc.qc.proper_pairs/10^6;
	mta$rseqc.qc.non_proper_pairs=mta$rseqc.qc.non_proper_pairs/10^6;
}

## define gene model distribution variables 
gen.vars = c('rseqc.rc.perc_exon_tags',
		'rseqc.rc.perc_intron_tags',
		'rseqc.rc.perc_intergenic_tags',
		'rseqc.rc.cds_exons_tags_kb',
		'rseqc.rc.intron_tags_kb',
		'rseqc.rc.tp_utr_exons_tags_kb',
		'rseqc.rc.fp_utr_exons_tags_kb',
		'rseqc.jc.total_splicing_events',
		'rseqc.jc.total_splicing_junctions');
gen.lab = c(
		"Exon Tags [%]",
		"Intron Tags [%]",
		"Intergenic Tags [%]",
		"CDS Exons: Tags Per Kb",
		"Intron Tags Per Kb",
		"3' UTR Tags Per Kb",
		"5' UTR Tags Per Kb",
		expression('Splicing Events [10'^5*']'),
		expression('Splicing Junctions [10'^5*']'))
 
## limit to those actually in the metadata (RC may not have been executed)
gen.lab = gen.lab[gen.vars %in% colnames(mta)]
gen.vars = gen.vars[gen.vars %in% colnames(mta)]

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
	## plotting parameters -- cex controls outlier label text size
	par(mar=c(6,4,4,1),mfrow=c(3,3),cex=0.35,cex.axis=2)
	if(varFlag=='tec') {
		for(t in 1:length(tec.vars)) {
			tec.var = tec.vars[t];
			tec.label = tec.lab[t];
			ylim=range(mta[,tec.var]);
			Boxplot(mta[,tec.var]~mta$spctl,labels=mta$samid,ylim=ylim,pch=20,las=2,main='',ylab='',xlab='',col=unique(mta$spctc));
			mtext(tec.label,line=1,side=3,cex=0.7)
		}
	} else if (varFlag=='gen') {
		for(t in 1:length(gen.vars)) {
			gen.var = gen.vars[t];
			gen.label = gen.lab[t];
			ylim=range(mta[,gen.var]);
			boxplot(mta[,gen.var]~mta$spctl,labels=mta$samid,ylim=ylim,pch=20,las=2,main='',ylab='',xlab='',col=unique(mta$spctc));
			mtext(gen.label,line=1,side=3,cex=0.7)
		}
	}
}