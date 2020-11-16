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
# Program:  session-info.r
# Version:  RSEQREP 2.0.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Tabulate packages and versions.  print to 3 tier latex table.
# Input:    N/A
# Output:   N/A
#############################################################################################################

source('../r/init-analysis.r');

## get session info object
session = sessionInfo()
loaded.names = names(session$loadedOnly)
extended.names = names(session$otherPkgs)

## extract information to DF
tbl = matrix(ncol=2,nrow=0)
colnames(tbl) = c('Package Name','Version')
for (i in 1:length(loaded.names)) {
	tbl= rbind(tbl,c(session$loadedOnly[[i]]$Package,session$loadedOnly[[i]]$Version))
}
for (i in 1:length(extended.names)) {
	tbl= rbind(tbl,c(session$otherPkgs[[i]]$Package,session$otherPkgs[[i]]$Version))
}

## push into 3 tiers
idx = ceiling(nrow(tbl)/3)
tbl3wide = cbind(tbl[1:idx,],tbl[(idx+1):(2*idx),],rbind(tbl[(2*idx+1):nrow(tbl),],matrix(ncol=2,nrow=(idx*3)-nrow(tbl))))

write.tbl(tbl,rownames=F);

caption.short = 'List of R packages and versions used for the analyses presented in this report.'
caption.long = gsub('_','-',paste(caption.short,' ',session$R.version$version.string," '",session$R.version$nickname,sep=''))

x = xtable(tbl3wide,align=c('c',rep(c('|p{1.0in}','p{0.5in}|'),3)),
		type='latex',floating='F',label='tab:session_info',
		caption=c(caption.long,caption.short)
);

print(x,tabular.environment='longtable',
		floating=F,
		include.rownames=F,
		size= 'scriptsize',
		sanitize.text.function = identity,
		sanitize.rownames.function = identity)