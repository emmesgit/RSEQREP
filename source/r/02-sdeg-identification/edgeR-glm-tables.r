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
# Program:  edgeR-glm-tables.r
# Version:  RSEQREP 1.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Tabulate differential gene results
# Input:    analysis/glm/<scp>_<trt>_<time>_glm_sig.tab.gz
# Output:  	N/A
############################################################################################

source('../r/init-analysis.r')

## dont print in scientific notation -- use 6 digits following decomal
options("scipen"=100, "digits"=6)

## set directories
in.dir.glm = paste(res.dir,'glm',sep='/');

## load gene annotations
gen = read.csv(paste(dta.dir,'/annot/filtered_gene_annotations.tab.gz',sep=''),sep='\t',stringsAsFactors=F,header=T);

## loop through specimen types
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];
	spcLabl = spcLabls[s];
		
	## loop through treatment groups;
	for (v in 1:length(trtFlags)) {
		trtFlag = trtFlags[v];
		trtLabl = trtLabls[v];
		
		## loop through timepoints
		for(t in 1:length(postb.times)) {
			time = postb.times[t];
			timel = postb.times.l[t];
			
			## specify infile
			in.file = paste(in.dir.glm,'/',spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz',sep='');
			
			## ensure there are significant results
			if (file.exists(in.file)) {
				tbl = read.csv(in.file,sep='\t',stringsAsFactors=F,header=T);			
			
				## add fields to data from annotations
				tbl$geneName = gen[match(rownames(tbl),gen$gene_id),'gene_name'];
				tbl$geneName = gsub('_','\\',tbl$geneName);
				tbl$desc = gen[match(rownames(tbl),gen$gene_id),'gene_desc'];
				tbl$desc = gsub('_','\\',tbl$desc);
				tbl$gene_type = gen[match(rownames(tbl),gen$gene_id),'gene_type'];
				tbl$gene_type = gsub('_',' ',tbl$gene_type);
				tbl$source = gsub('_',' ',gen[match(rownames(tbl),gen$gene_id),'source']);
				tbl$num_trans = gen[match(rownames(tbl),gen$gene_id),'num_trans'];
				tbl$num_exons = gen[match(rownames(tbl),gen$gene_id),'max_num_exons'];
				tbl$gene_id = rownames(tbl)
				
				## sort by absolute log fold change
				if (nrow(tbl)>1) {
					tbl$abs_fc = abs(tbl$logFC);
					tbl = tbl[order(tbl$abs_fc,decreasing=T),]			
				}

				## handle p_value/FDR presentation
				tbl$PString = as.character(tbl$PValue);
				tbl$PString[tbl$PValue<0.0001]='<0.0001';
				tbl$PString[tbl$PValue==1]='>0.9999';
				tbl$PString[tbl$PValue>=0.0001 & tbl$PValue<1]=as.character(round(tbl$PValue[tbl$PValue>=0.0001 & tbl$PValue<1],4));
				tbl$FDRString[tbl$FDR<0.0001]='<0.0001';
				tbl$FDRString[tbl$FDR==1]='>0.9999';
				tbl$FDRString[tbl$FDR>=0.0001 & tbl$FDR<1]=as.character(round(tbl$FDR[tbl$FDR>=0.0001 & tbl$FDR<1],4));
				
				## select data and update column names
				tbl = tbl[,c('gene_id','geneName','desc','gene_type','logFC','logCPM','LR','PValue','FDR')]
				colnames(tbl) = c('Ensembl Gene ID', 'Ensembl Gene Name','Ensembl Gene Description','Gene Type',paste('$Log_{2}$ Fold Change (',timel,' vs. pre-treatment)',sep=''),'Average $Log_{2}$ CPM','Likelihood Ratio Test Statistic','P-Value','FDR Adjusted P-Value');
				
				## print results in tabular format
				write.tbl(tbl,rownames=F);
				
				## set captions
				caption.short = paste('Genes differentially expressed at ',timel,' compared to pre-treatment (',spcLabl,', ',trtLabl,').',sep='');
				caption.long  = paste('Genes differentially expressed at ',timel,' compared to pre-treatment  (',spcLabl,', ',trtLabl,'). Sorted by descending absolute $log_{2}$ fold change (Day ',time,' vs. pre-treatment). Gene model summaries and annotations are based on Ensembl Version ',ensembl.version,'.',sep=''); 
				
				## add link to ensembl ID
				tbl[,'Ensembl Gene ID'] = paste('\\parbox{0.9in}{\\','href{http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=',
						tbl[,'Ensembl Gene ID'],'}{',tbl[,'Ensembl Gene ID'],'}}',sep='');
				
				## generate xtable and print
				x = xtable(tbl,align=c('c','p{0.9in}','p{0.6in}','p{2.35in}','p{0.5in}','p{0.75in}','p{0.45in}','p{0.45in}','p{0.45in}','p{0.45in}'),type='latex',floating='F',label=c(paste('tab:sdeg',spcFlag,'d',time,sep='')),
						caption=c(caption.long,caption.short));
				
				print(x,tabular.environment='longtable',
						floating=F,
						include.rownames=F,
						hline.after = c(-1,nrow(x)),
						size= 'scriptsize',
						add.to.row = list(pos = list(0),command = "\\hline \\endhead "),
						sanitize.colnames.function = identity,
						sanitize.text.function = identity,
						sanitize.rownames.function = identity)
			}
		}
	}
}