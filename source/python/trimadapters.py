#############################################################################################################
# RSEQREP: RNA-Seq Reports, an open-source cloud-enabled framework for reproducible
# RNA-Seq data processing, analysis, and result reporting
# 
# https://github.com/emmesgit/RSEQREP
#
# Copyright (C) 2019 The Emmes Corporation 
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
# Program:  trimadapters.py
# Version:  RSEQREP 2.2.0
# Author:   William F Hooper, Travis L. Jensen, Johannes B. Goll
# Purpose:  Cutadapt wrapper 
# Input:    N/A
# Output:   N/A
#############################################################################################################

import sys
import os
import subprocess
import utils
import argparse


## Get the correct adapter sequences for this sample
def getAdapters(sam, fp_adapters, tp_adapters, samid_all):
	fp_adapter = [fp_adapters[i[0]] for i in enumerate(samid_all) if samid_all[i[0]] == sam][0]
	tp_adapter = [tp_adapters[i[0]] for i in enumerate(samid_all) if samid_all[i[0]] == sam][0]
	return({'tp':tp_adapter, 'fp':fp_adapter})


def trimAdapters(cutadapt, infile, outfile, adapt3='NA', adapt5='NA'):
	
	## Determine how many adapters there are, and the ended-ness of the data
	is_paired_end  = False if len(infile) == 1 else True
	is_three_prime = False if adapt3 == 'NA' else True
	is_five_prime  = False if adapt5 == 'NA' else True
	
	## Check to make sure at least one adapter is present
	if (not is_three_prime and not is_five_prime):
		raise RuntimeError('Missing adapters to use for trimming')
	
	## If both adapters are present and this is single ended data, use linked trimming
	if (is_three_prime and is_five_prime and not is_paired_end):
		trim_string = '-a '+adapt5+'...'+adapt3+' -o '+outfile[0]+' '+infile[0]
		
	## If only single ended and 5' adapter, use anchored 5' trimming
	elif (is_five_prime and not is_three_prime and not is_paired_end):
		trim_string = '-g ^'+adapt5+' -o '+outfile[0]+' '+infile[0]
		
	
	## If only single ended and 3' adapter, use default 3' trimming
	elif (is_three_prime and not is_five_prime and not is_paired_end):
		trim_string = '-a '+adapt3+' -o '+outfile[0]+' '+infile[0]
		
	
	## If paired end, make sure that both adapters are present
	elif (is_three_prime and is_five_prime and is_paired_end):
		trim_string = '-a '+adapt3 + ' -A '+adapt5+' -o '+outfile[0]+' -p '+outfile[1]+' '+infile[0]+' '+infile[1]
		
	trim_cmd = cutadapt+' '+trim_string+' -m 15 -e 0.1 -O 5'
	bench_obj = utils.logging_call(trim_cmd, shell=True)
	
	return(bench_obj)