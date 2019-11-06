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
# Program:  parse-rnaseq-configuration.r
# Version:  RSEQREP 2.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Run sanity checks on configuration, create GMT formatted gene sets, download genome/transcriptome
#				parse workflow and analysis configurations, create work space directories.
# Input:    <COMMAND LINE INPUTS>
# Output:   data/gene_sets/all.ensembl.tab.gz
#			data/analysis_config.csv
#			data/sample_metadata.csv
#############################################################################################################

libs = c('openxlsx','stringr','biomaRt','yaml');
invisible(suppressPackageStartupMessages(lapply(libs, require, character.only=T)))


## Get command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args)>0) {
	infile       = args[1]
	source.dir   = args[2]
}

############################
#
# Import data and perform 
#       Sanity Checks
#
############################

## import metadata
metadata = read.xlsx(infile,sheet=1)

## remove unused fields
metadata = metadata[c(3:nrow(metadata)),]

if (nrow(metadata) == 0) stop('No samples were found\n')

## Import Workflow Data
workflow.config = read.xlsx(infile,sheet=2)[,c('Name','Value')] # key value fields only
workflow.config[is.na(workflow.config)] = '' # replace missing values with empty string

## Import analysis Data
analysis.config = read.xlsx(infile,sheet=3)[,c('Name','Value')] # key value fields only

## define variables
genome.file = workflow.config$Value[workflow.config$Name=='genome_file']
gtf.file = workflow.config$Value[workflow.config$Name=='gtf_file']
threads = workflow.config$Value[workflow.config$Name=='num_threads']
ensembl.version = workflow.config$Value[workflow.config$Name=='ensembl_version']
base.dir = dirname(analysis.config$Value[analysis.config$Name=='report_dir'])
gmt.files = unlist(str_split(workflow.config$Value[workflow.config$Name=='gmt_entrez_files'],';'))
aws.s3.fq.files = c(metadata$fastq_file_1[grep('^s3://',metadata$fastq_file_1)],metadata$fastq_file_2[grep('^s3://',metadata$fastq_file_2)])
rseqc.progs = paste(workflow.config$Value[workflow.config$Name=='rseqc_dir'],
		c('bam_stat.py','junction_annotation.py','read_distribution.py','read_GC.py'),sep='/')
aws.prog = workflow.config$Value[workflow.config$Name=='aws_prog']


###############
# Sanity Checks
###############
cat('Performing sanity checks on configuration options\n')

## ensure all programs exist
progs = c(workflow.config$Value[match(c('aws_prog','openssl_prog','samtools_prog',
		'fcts_prog','fastqc_prog'),workflow.config$Name)],rseqc.progs)
for (i in 1:length(progs)) {
	if(!file.exists(progs[i])) {
		cat(paste('Program does not exist at this location!',progs[i],'\n'))
		q(status=25,save='no')
	}
}

## ensure only one mapping program is valid and exists
progs = c(workflow.config$Value[match(c('star_prog','hisat_prog'),workflow.config$Name)])
if(!(file.exists(progs[1])) & !(file.exists(progs[2]))) {
	cat(paste('Please specify a valid mapping software program location!',progs[1],progs[2],'\n'))
	q(status=25,save='no')
} else if (file.exists(progs[1]) & file.exists(progs[2])) {
	cat(paste('Please only specify ONE mapping software program location!',progs[1],progs[2],'\n'))
	q(status=25,save='no')
} 

## ensure the base result/preprocessing directory exists
if (!dir.exists(base.dir)) {
	cat(paste('Base directory does not exist!',base.dir,'\n'))
	q(status=25,save='no')
}

## ensure colors are defined correctly
areColors = function(x) {sapply(x, function(X) {tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)})}
col.test = areColors(c(metadata$timec,metadata$spctc,metadata$trtc))
if(any(col.test)==F) {
	cat(paste('The following colors are not recognized in the sample metadata:',names(col.test[col.test==F]),'\n'))
	q(status=25,save='no')
}

## ensure ensembl version is valid
ensembl.test = tryCatch(is(useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=ensembl.version),'Mart'), error = function(e) FALSE)
if(ensembl.test==F) {
	cat(paste('The specified Ensembl version is invalid:',ensembl.version,'\n'))
	q(status=25,save='no')
}

############################
#
# Parse Analysis Configuration
#
############################

cat('Parsing analysis configuration\n')

results.dir = analysis.config$Value[analysis.config$Name=='report_dir']
pre.dir = workflow.config$Value[workflow.config$Name=='pre_dir']

## Create preprocessing hierarchy
if(!dir.exists(pre.dir)) {
	system(paste('mkdir',pre.dir))
}

if(!dir.exists(pre.dir)) {
	system(paste('mkdir',pre.dir))
}

if(!dir.exists(paste(pre.dir,'genome',sep='/'))) {
	system(paste('mkdir',paste(pre.dir,'genome',sep='/')))
}
if(!dir.exists(paste(pre.dir,'annot',sep='/'))) {
	system(paste('mkdir',paste(pre.dir,'annot',sep='/')))
}


## Create results / data hierarchy
if(!dir.exists(results.dir)) {
	system(paste('mkdir',results.dir))
}
if(!dir.exists(paste(results.dir,'data',sep='/'))) {
	system(paste('mkdir',paste(results.dir,'data',sep='/')))
}

dta.dir = paste(results.dir,'data',sep='/')
if(!dir.exists(dta.dir)) {
	system(paste('mkdir',dta.dir))
}

if(!dir.exists(paste(dta.dir,'annot',sep='/'))) {
	system(paste('mkdir',paste(dta.dir,'annot',sep='/')))
}
if(!dir.exists(paste(dta.dir,'gene_sets',sep='/'))) {
	system(paste('mkdir',paste(dta.dir,'gene_sets',sep='/')))
}

if(!dir.exists(paste(results.dir,'analysis',sep='/'))) {
	system(paste('mkdir',paste(results.dir,'analysis',sep='/')))
}

analysis.dir = paste(results.dir,'analysis',sep='/')
if(!dir.exists(analysis.dir)) {
	system(paste('mkdir',analysis.dir))
}

if(!dir.exists(paste(analysis.dir,'lcpm',sep='/'))) {
	system(paste('mkdir',paste(analysis.dir,'lcpm',sep='/')))
}
if(!dir.exists(paste(analysis.dir,'lcpm_fc',sep='/'))) {
	system(paste('mkdir',paste(analysis.dir,'lcpm_fc',sep='/')))
}
if(!dir.exists(paste(analysis.dir,'dist',sep='/'))) {
	system(paste('mkdir',paste(analysis.dir,'dist',sep='/')))
}
if(!dir.exists(paste(analysis.dir,'pca',sep='/'))) {
	system(paste('mkdir',paste(analysis.dir,'pca',sep='/')))
}
if(!dir.exists(paste(analysis.dir,'glm',sep='/'))) {
	system(paste('mkdir',paste(analysis.dir,'glm',sep='/')))
}
if(!dir.exists(paste(analysis.dir,'gsea',sep='/'))) {
	system(paste('mkdir',paste(analysis.dir,'gsea',sep='/')))
}
if(!dir.exists(paste(analysis.dir,'pvclust',sep='/'))) {
	system(paste('mkdir',paste(analysis.dir,'pvclust',sep='/')))
}

## Create report directory 
if(!dir.exists(paste(results.dir,'report',sep='/'))) {
	system(paste('mkdir',paste(results.dir,'report',sep='/')))
}

## add ensembl version to analysis variables
ensembl.version = workflow.config$Value[workflow.config$Name=='ensembl_version']
analysis.config = rbind(analysis.config,c('ensembl_version',ensembl.version))
run.read.dist = workflow.config$Value[workflow.config$Name=='run_read_dist']
analysis.config = rbind(analysis.config,c('run_read_dist',run.read.dist))
star.prog = workflow.config$Value[workflow.config$Name=='star_prog']
analysis.config = rbind(analysis.config,c('star_prog',star.prog))
cram.flag = workflow.config$Value[workflow.config$Name=='save_cram']
analysis.config = rbind(analysis.config,c('save_cram',cram.flag))

## Write analysis configuration 
write.table(analysis.config,paste(dta.dir,'analysis_config.csv',sep='/'),sep=',',quote=T,row.names=F)

############################
#
# Parse Metadata
#
############################

## write CSV formatted metadata file
cat('Parsing sample metadata\n')
write.table(metadata,paste(pre.dir,'sample_metadata.csv',sep='/'),na='',sep=',',quote=F,row.names=F)
write.table(metadata,paste(dta.dir,'sample_metadata.csv',sep='/'),na='',sep=',',quote=F,row.names=F)


############################
#
# Parse Workflow Configuration
# and create indexes, BED files, 
# and gene sets files
#
############################

cat('Parsing workflow configuration\n')

## define variables
gene.set.labels = unlist(str_split(workflow.config$Value[workflow.config$Name=='gmt_entrez_files_labels'],';'))
ensembl.version = workflow.config$Value[workflow.config$Name=='ensembl_version']

## download GTF from ensembl FTP
gtf.file = paste(pre.dir,'/annot/Homo_sapiens.ensembl.version',ensembl.version,'.chr.gtf',sep='')
if(!file.exists(gtf.file)) {
	cat('Downloading GTF annotations from Ensembl FTP site\n')
	system(paste(source.dir,'/shell/download-ensembl-gtf.sh ',ensembl.version,' ',dirname(gtf.file),sep=''))
}

	
## Create BED file if not exists
bed.file = gsub('gtf$','bed',gtf.file)
if(!file.exists(bed.file)) {
	cat('Creating BED File from GTF annotations\n')
	system(paste('perl ',paste(source.dir,'/perl/gtf2bed.pl',sep='/'),gtf.file,'>',bed.file))
}

## download from ensembl FTP
genome.file = paste(pre.dir,'/genome/Homo_sapiens.ensembl.version',ensembl.version,'.genome.fa',sep='')
if(!file.exists(genome.file)) {
	cat('Downloading Genome annotations from Ensembl FTP site\n')
	system(paste(source.dir,'/shell/download-ensembl-genome.sh ',ensembl.version,' ',dirname(genome.file),sep=''))
}



## Create Gene Set tab delimited file
gene.sets.file = paste(dta.dir,'gene_sets/all.ensembl.tab',sep='/')
if(!file.exists(paste(gene.sets.file,'.gz',sep='')) & length(gmt.files)>0 & gmt.files[1]!='') {
	cat('Preparing GSEA Gene Sets\n')
	
	## import ensembl annotations
	ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version= ensembl.version)
	att =  c('ensembl_gene_id','entrezgene','hgnc_symbol')
	ens.ano = getBM(attributes=att,mart=ensembl)
	colnames(ens.ano) = c('gene_id','entrez_id','gene_name')
	
	## write mapping table
	outfile = paste(dta.dir,'annot/ensembl_entrez_genename.tab',sep='/')
	write.table(ens.ano,outfile,quote=F,row.names=F,sep='\t')
	R.utils::gzip(outfile,overwrite=TRUE);
	
	## Get labels
	gmt.labs = unlist(str_split(workflow.config$Value[workflow.config$Name=='gmt_entrez_files_labels'],';'))
	
	## initialize data set
	dta = matrix(ncol=4,nrow=0)
	
	## for each gene set file
	for (i in 1:length(gmt.files)) {
		
		## Migrate files
		new.gmt.file = paste(dta.dir,'gene_sets',basename(gmt.files[i]),sep='/')
		system(paste('cp',gmt.files[i],new.gmt.file))
		
		## extract GMT information
		gmt = scan(new.gmt.file, what="", sep="\n")
		
		dta.gmt = matrix(ncol=4,nrow=0)
		
		## for each gene set
		for (j in 1:length(gmt)) {
			
			## split GMT elements
			gene.set.str = unlist(str_split(gmt[[j]],'\t'))
			
			## determine what format the gene IDs are in
			if (j==1) {
				## for each entry, determine mappings to Ensembl IDs, gene symbols, or entrez IDs.
				## The entry with the most mapped annotations will be used.
				type = colnames(ens.ano)[which.max(sapply(1:length(ens.ano),function(x) {length(which(ens.ano[,x] %in% gene.set.str[3:length(gene.set.str)]))}))]
			}
			
			## get ensembl IDs
			ens.ids = ens.ano$gene_id[ens.ano[,type] %in% gene.set.str[3:length(gene.set.str)]]
			
			if (length(ens.ids) > 0) {
				## generate table: 1) gene set category label 2) gene ids 3) gene set name 4) link
				gene.set.tbl = cbind(gmt.labs[[i]],ens.ids,gene.set.str[1],gene.set.str[2])
			
				## add to category type -- multiple levels of rbind steps speed up the process dramatically.
				dta.gmt = rbind(dta.gmt,gene.set.tbl)
			}
		}
		## merge category data -- multiple levels of rbind steps speed up the process dramatically.
		dta = rbind(dta,dta.gmt)
	}
	
	## remove duplicate entries
	dta.u = unique(na.omit(dta))
	colnames(dta.u) = c('','','','')
	
	## write results
	write.table(dta.u,gene.sets.file,sep='\t',quote=F,row.names=F,col.names=F)
	
	## Gzip all files in the gene_sets directory
	system(paste('gzip',paste(dta.dir,'gene_sets/*',sep='/')))
}

## remove/update/add values to workflow configuration
workflow.config = workflow.config[!workflow.config$Name %in% c('gmt_entrez_files','gmt_entrez_files_labels'),]
workflow.config = rbind(workflow.config,c('genome_file',genome.file))
workflow.config = rbind(workflow.config,c('gtf_file',gtf.file))
workflow.config = rbind(workflow.config,c('bed_file',bed.file))
workflow.config = rbind(workflow.config,c('maxthreads','4'))


## Convert workflow config to YAML
workflow.config.yaml = list()
for (x in 1:nrow(workflow.config)) {
	workflow.config.yaml[[workflow.config$Name[x]]] = workflow.config$Value[x]
} 

## Add sample IDs, paths to fastqs
workflow.config.yaml[['samid']]  = metadata$samid	
workflow.config.yaml[['fastq1']] = metadata$fastq_file_1

if (all(is.na(metadata$fastq_file_2))) {
	workflow.config.yaml[['fastq2']] = ''
	cat('Working in single-end mode\n')
} else {
	workflow.config.yaml[['fastq2']] = metadata$fastq_file_2
	cat('Working in paired-end mode\n')
}

## Parse adapter sequences 
if (all(is.na(metadata$threep_adapter_seq))) {
	cat("Skipping 3' adapter trimming... some or all adapters omitted.\n")
} 
if (all(is.na(metadata$fivep_adapter_seq))) {
	cat("Skipping 5' adapter trimming... some or all adapters omitted.\n")
} 

## Check to make sure an appropriate quality trimming value was specified 
if (is.na(workflow.config.yaml[['quality_trim']])) {
	stop("No quality trimming value specified")
} else if (!as.numeric(workflow.config.yaml[['quality_trim']]) %in% 0:40) {
	stop("Quality trimming value specified is outside allowed range: ", workflow.config.yaml[['quality_trim']])
}

## Check to make sure Bowtie2 remapping flag was set 
if (is.na(workflow.config.yaml[['ion_torrent']])) {
	stop("Ion Torrent remapping flag not set")
} else if (!workflow.config.yaml[['ion_torrent']] %in% 0:1) {
	stop("Ion Torrent remapping flag incorrectly set: ", workflow.config.yaml[['ion_torrent']])
}



metadata$threep_adapter_seq[is.na(metadata$threep_adapter_seq)] = 'NA'
metadata$fivep_adapter_seq[is.na(metadata$fivep_adapter_seq)]   = 'NA'
workflow.config.yaml[['tp_adapter_seq']]  = metadata$threep_adapter_seq	
workflow.config.yaml[['fp_adapter_seq']]  = metadata$fivep_adapter_seq	

## Write out preprocessing config in YAML
sink(paste(pre.dir,'preprocess_config.yaml',sep='/'))
cat(as.yaml(workflow.config.yaml))
sink()


## Write temporary file with paths to directories
## workflow directory,analysis directory,workflow configuration, and sample_metadata csv
tmp.dta = c(pre.dir, analysis.config[analysis.config$Name=='report_dir','Value'], paste(pre.dir,'preprocess_config.yaml',sep='/'),paste(pre.dir,'sample_metadata.csv',sep='/'))
write.table(tmp.dta,paste(source.dir,'/../dir.csv',sep=''),sep=',',quote=F,row.names=F,col.names=F)