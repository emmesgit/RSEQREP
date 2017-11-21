############################################################################################
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
# Program:  init-export-annotations.r
# Version:  RSEQREP 1.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose: 	Export Ensembl annotation data
# Input:    N/A
# Output:  	analysis/annot/filtered_gene_annotations.tab
#			analysis/annot/excluded_gene_counts.tab
############################################################################################

source('init-analysis.r')

## import ensembl annotations
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=ensembl.version)
att =  c('ensembl_gene_id','start_position','end_position','strand','gene_biotype','percentage_gc_content','transcript_count','description','chromosome_name','source','entrezgene','hgnc_id','hgnc_symbol','transcript_length','ensembl_transcript_id','ensembl_exon_id')
if (ensembl.version > 87) {
	att[att=='percentage_gc_content'] = 'percentage_gene_gc_content'
}
dta = getBM(attributes=att,mart=ensembl,uniqueRows=F)

## First filter to count the number of exons, then count number of transcripts, and the maximum length of transcripts per gene.
## Finally, filter for entries that exist on a chromosome of intrest
colnames(dta) = c('gene_id','gene_start','gene_end','strand','gene_type','gc_content','num_trans','gene_desc','chr','source','entrez_id','hugo_gene_name','gene_name','trans_length','transcript_id','exon_id')
dtaFilt = sqldf(paste('select gene_id, gene_start, gene_end, strand, gene_type, gc_content, num_trans, gene_desc, chr, source, entrez_id, hugo_gene_name, gene_name, trans_length, transcript_id, count(exon_id) as num_exons from dta where (gene_type not in ("',ensembl.gene.types.remove,'") and chr in ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT")) group by transcript_id',sep=''))
dta.filtered.chr = sqldf('select gene_id, gene_start, gene_end, strand, gene_type, gc_content, num_trans, gene_desc, chr, source, entrez_id, hugo_gene_name, gene_name, max(trans_length) as max_trans_length, count(distinct(transcript_id)) as num_trans, max(num_exons) as max_num_exons from dtaFilt group by gene_id')
dta.filtered.chr$gene_name_lab = dta.filtered.chr$gene_name
dta.filtered.chr$gene_name_lab[which(dta.filtered.chr$gene_name=='')] = dta.filtered.chr$gene_id[which(dta.filtered.chr$gene_name=='')]

## write resulting table
outfile = paste(dta.dir,'/annot/filtered_gene_annotations.tab',sep='')
write.table(dta.filtered.chr,outfile,quote=F,row.names=F,sep='\t')
R.utils::gzip(outfile,overwrite=TRUE);

## grab excluded genes
dta.excluded = dta[!dta$gene_id %in% dta.filtered.chr$gene_id,]
dta.excluded.chr = dta.excluded[dta.excluded$chr %in% c(1:22,'X','Y','MT'),]
dta.excluded.chr.unique = dta.excluded.chr[which(!duplicated(dta.excluded.chr$gene_id)),]

## write excluded genes
outfile2 = paste(dta.dir,'/annot/excluded_gene_counts.tab',sep='')
write.table(dta.excluded.chr.unique,outfile2,quote=F,row.names=F,sep='\t')
R.utils::gzip(outfile2,overwrite=TRUE);