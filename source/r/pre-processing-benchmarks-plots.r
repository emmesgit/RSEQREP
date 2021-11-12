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
# Program:  preprocessing-benchmarks-plot.r
# Version:  RSEQREP 2.2.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  generate plots on peak CPU time, Wall cloc time, RAM, and VM usage for each process
# Input:    <pre-process-dir>/preprocess.db
# Output:   N/A
#############################################################################################################

source('../r/init-analysis.r')
options(scipen=999)

## Import data
sqlite.db.file = paste0(src.dir,'/../',grep('*.db',system(paste0('ls ', src.dir,'/..'),intern=T),value=T))

## if only a directory is specified quit without error.
if(unlist(str_split(sqlite.db.file,''))[length(unlist(str_split(sqlite.db.file,'')))] %in% c('/','.')) {
	quit(save='no',status=0);
}

## if the file exists and is not empty, plot:
if (file.exists(sqlite.db.file) & unlist(str_split(system(paste('ls -l',sqlite.db.file),intern=T),' '))[5]!=0) {
	
	dta = dbGetQuery(dbConnect(dbDriver("SQLite"), dbname = sqlite.db.file), 
			paste("select * from benchmark left join file on benchmark.process_id=file.process_id where benchmark.return_code = 0"))
	
	## convert file size from bytes to Gib, and ram from Kib to Gib
	dta$file_bytes = dta$file_bytes/(1024^3)
	dta$virtual = dta$virtual/(1024^2)
	dta$resident_set_size = dta$resident_set_size/(1024^2)
	
	## standard variables	
	benchmark.types = c('cpu_time','wc_time','virtual','resident_set_size','file_bytes')
	benchmark.labs = c('CPU Time [s]','Wall Clock Time [s]','Peak Virtual Memory [Gib]','Peak Resident Set Size Memory [Gib]','Result File Size [Gib]')
	
	## define processes
	process.types = c('encrypted fastq','fastq','STAR mapped and unmapped reads bam file',
			'mapped and unmapped CRAM file','FastQC','reseqc bam read distribution','rseqc bam statistics',
			'rseqc bam junctions','rseqc bam gc%','feature counts')
	process.labs = c('Fastq S3 Download (AWS CLI)','Fastq Decryption (OpenSSL)','Reference Alignments (STAR)',
			'BAM Compression (CRAM)','Read QC (FASTQC)','RSeQC Read Distribution','RSeQC BAM Statistics',
			'RSeQC Junction Annotation','RSeQC Read GC Content','Expression Quantification (Subread)')
	
	if (rdFlag==F) {
		process.types = process.types[which(!process.types=='reseqc bam read distribution')]
		process.labs = process.labs[which(!process.labs=='RSeQC Read Distribution')]
	}
	
	if (length(which(dta$file_type=='encrypted fastq'))==0) {
		process.types = process.types[which(!process.types=='encrypted fastq')]
		process.labs = process.labs[which(!process.labs=='Fastq S3 Download (AWS CLI)')]
	}
	
	if (length(which(dta$file_type=='fastq'))==0) {
		process.types = process.types[which(!process.types=='fastq')]
		process.labs = process.labs[which(!process.labs=='Fastq Decryption (OpenSSL)')]
	}
	
	if (cramFlag==F) {
		process.types = process.types[which(!process.types=='mapped and unmapped CRAM file')]
		process.labs = process.labs[which(!process.labs=='BAM Compression (CRAM)')]
	}
	
	if (length(which(dta$file_type=='FastQC'))==0) {
		process.types = process.types[which(!process.types=='FastQC')]
		process.labs = process.labs[which(!process.labs=='Read QC (FASTQC)')]
	}
	
	## hisat or STAR?
	if(length(grep('HISAT2',dta$file_type))>0) {
		process.types[which(process.types=='STAR mapped and unmapped reads bam file')] = 'HISAT2 mapped and unmapped reads BAM file'
		process.labs[which(process.labs=='Reference Alignments (STAR)')] = 'Reference Alignments (HISAT2)'
	}
	
	## define process colors
	colors.fun = colorRampPalette(c("red",'grey70',"blue"))
	process.colors = colors.fun(length(process.types))
	
	## define plotting area
	par(mar=c(2.5,4,0.5,0.5),mfrow=c(3,2),oma=c(0,0,0,0))
	
	## generate plots:
	for (i in 1:length(benchmark.types)) {
		benchmark.type = benchmark.types[i];
		benchmark.lab  = benchmark.labs[i];
		
		## describe data -- mean and standard deviation
		dta.benchmark.median = tapply(dta[,benchmark.type],dta[,'file_type'],median)[process.types]
		dta.benchmark.max = tapply(dta[,benchmark.type],dta[,'file_type'],max)[process.types]
		dta.benchmark.min = tapply(dta[,benchmark.type],dta[,'file_type'],min)[process.types]
		barCenters = seq(0.75,((length(process.types)*1.2)-0.5),length.out=length(process.types))	
		
		## Plot
		barplot(height=dta.benchmark.median,names.arg=1:length(dta.benchmark.median),ylim=c(0,max(dta.benchmark.max,na.rm=T)*1.2),main='',col=process.colors,ylab=benchmark.lab)
		
		## add min/max segment
		segments(barCenters,dta.benchmark.min,barCenters,dta.benchmark.max,lwd = 2,xpd=T)
		
	}
	plot.new()
	legend('center',col=process.colors,legend=paste(1:length(dta.benchmark.median),process.labs,sep=': '),cex=1.2,pch=15,pt.cex=1.5,bty='n')
}