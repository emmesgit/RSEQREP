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
# Program:  putfile.py
# Version:  RSEQREP 2.0.1
# Author:   William F Hooper, Travis L. Jensen, Johannes B. Goll
# Purpose:  Interface for uploading files to the cloud
# Input:    N/A
# Output:   N/A
#############################################################################################################

import subprocess
import logging
import utils



## Upload to an Amazon S3 bucket, preserving folder structure
def s3Upload(file, destination, prog='aws'):
	if destination[-1] != '/':
		destination += '/'
	
	if not isinstance(file, list):
		file = file.split(' ')
	
	for f in file:
		cmd = '%s s3 cp %s %s' % (prog, f, destination+f)
		
		## Attempt to upload
		try:
			utils.logging_call(cmd, shell=True)
		except subprocess.CalledProcessError:
			logging.error('S3 upload failed. See above for more details.')
			exit(1)
			
	return(True)



## Use different functions depending on cloud provider (For now it's just AWS)
def upload(file, destination, cloud, prog):
	upload_func = {'aws': s3Upload}
	upload_func[cloud](file=file, destination=destination, prog=prog)