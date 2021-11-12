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
# Program:  utils.py
# Version:  RSEQREP 2.2.0
# Author:   William F Hooper, Travis L. Jensen, Johannes B. Goll
# Purpose:  Miscellaneous utility functions
# Input:    N/A
# Output:   N/A
#############################################################################################################

import datetime
import shutil
import snakemake
import subprocess
import logging
import os
import glob
import hashlib
import psutil
import time
from itertools import chain


## Simple benchmark class
## Essentially a stripped down version of the one provided by snakemake
class Benchmark:
    def __init__(self):
        self.running_time = 0
        self.cpu_seconds = 0
        self.max_rss = 0
        self.max_vms = 0 
        self.return_code = None
        
        
    ## Take measurements
    def update(self, pid):
    
        # Memory measurements
        rss, vms = 0, 0
        
        # Iterate over process and all children
        try:
            main = psutil.Process(pid)
            
            ## Sum virtual memory, resident set size over all threads
            ## For a given update, sum CPU times across threads and replace recorded time with sum
            for proc in chain((main,), main.children(recursive=True)):
                self.cpu_seconds = sum([x for x in proc.cpu_times()])
                meminfo = proc.memory_full_info()
                rss += meminfo.rss
                vms += meminfo.vms
                
            rss /= 1024 * 1024
            vms /= 1024 * 1024
            
            # Update benchmark record's RSS and VMS
            self.max_rss = max(self.max_rss, rss)
            self.max_vms = max(self.max_vms, vms)
            
        except psutil.Error as e:
            return
            
       

## Dummy log handler function -- snakemake output gets doubled if updating the logging config 
## without redirecting log output
def logHandler(x):
    return



## Run a process and redirect process output to log
def logging_call(popenargs, **kwargs):
    
    ## Log command
    logging.debug(popenargs)
    
    ## Start process, get pid
    start_time = time.time()
    process = subprocess.Popen(popenargs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **kwargs)
    pid = process.pid
    
    ## Check/capture process output    
    def check_io():
        output = process.stdout.readline().decode()
        if output:
            logging.info(output)
    
    ## Init benchmark object
    bench_obj = Benchmark()

    ## keep checking stdout/stderr until the child exits
    while process.poll() is None:
        check_io()
        bench_obj.update(pid)
        time.sleep(0.5)
    
    ## Compute wall clock time
    bench_obj.running_time = time.time() - start_time 
    
    ## Check return code
    bench_obj.return_code = process.returncode
    if bench_obj.return_code != 0:
        raise subprocess.CalledProcessError(returncode=bench_obj.return_code, cmd=popenargs)

    ## Return benchmarking object
    return(bench_obj)



## Insert object into its own list if it isn't already of class list
def toList(x):
    if isinstance(x, str):
        x = [x]
    elif not isinstance(x, list):
        x = list(x)
    return(x) 
    
    

## cat together files
def cat(in_files, out_file):
    f = ' '.join(in_files)
    cmd = 'cat '+f+' > '+out_file
    subprocess.run(cmd, shell=True, check=True)
    
    
    
## Decrypt file(s), return path to result
def decryptFile(in_file, openssl, password, hash):
    res = []
    for f in in_file:
        cmd = openssl+' aes-256-cbc -md '+hash+' -d -pass pass:'+password+' < '+f+' > ' + os.getcwd() + '/' +os.path.basename(f)[:-4]
        #subprocess.run(cmd, shell=True, check=True)
        logging_call(cmd, shell=True)
        res.append(os.getcwd() + '/' +os.path.basename(f)[:-4])
    return(res)



## Encrypt file(s), return path to result
def encryptFile(file, openssl, password, hash):
    res = []
    
    for f in file:
        cmd = openssl+' aes-256-cbc -md '+hash+' -pass pass:'+password+' < '+f+' > ' + f + '.enc'
        #subprocess.run(cmd, shell=True, check=True)
        logging_call(cmd, shell=True)
        res.append(f + '.enc')
    return(res)
    


## Compute an sha256sum/md5sum for a file 
## Source: https://docs.python.org/3/library/hashlib.html
def sha256(fname):
    hash_sha256 = hashlib.sha256()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(2 ** 20), b""):
            hash_sha256.update(chunk)
    return hash_sha256.hexdigest()
def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(2 ** 20), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()



## Catch last return code, write to log, pass back to snakemake so 
## that it knows to quit if a non-zero exit status was returned
def returnCode(log, sample='', process='', outfile=''):
    
    ## Build log string
    t = datetime.datetime.now().strftime('%m/%d/%Y %I:%M:%S %p')
    log_str = 'echo "[rseqrep] ' + t + '- INFO - Return Code for ' + sample + ' ' + process + ': $? ">> ' + log 
    return('; e=$?;  ' + log_str + '; exit $e')




## Timestamp and compress the log upon exit
def archiveLog(log):
    t = datetime.datetime.now().strftime('%Y-%m-%d.h%H-m%M-s%S')
    new_log = log+'.'+t
    try:
        shutil.copy2(log, new_log)
        snakemake.shell('gzip -9 '+new_log)
    except (FileNotFoundError, PermissionError) as e:
        pass



## Shorthand for touch - takes list<str> or str
def touch(file):
    if type(file) is list:
        [snakemake.shell('touch ' + x) for x in file]
    else:
        snakemake.shell('touch {file}')
    


## Write process info so we can correctly merge benchmarks later
def writeBenchmarkInfo(procid, samid, bench_file, return_code):
    with open(bench_file, 'w') as f:
        f.write('procid\tsamid\tbench_file\treturn_code\n')
        f.write('%s\t%s\t%s\t%s\n' % (procid, samid, bench_file[:-5], return_code))
    


## Copy length counts to the data directory
def addResultsToDataDir(datadir, resdir):
    destdir = datadir + '/' + resdir 
    if (not os.path.isdir(destdir)):
        os.mkdir(destdir)
    [shutil.copy2(x, destdir) for x in glob.glob(resdir+'/*_counts.RData')]
    
    
    
## Determine the outputs from HISAT2 -- changes depending on whether we're remapping or not
def hisatOutput(ends, iontorrent):
    res = {'aln' : snakemake.io.temp('tmp/{sample}.sam')}
            
    if (iontorrent):
        res['un_aln'] = snakemake.io.temp(expand('tmp/{{sample}}_iontorrent_{pe}.fastq.gz', pe=ends))
        
    return(res)
    


## Determine bam output based on whether we're remapping/compressing
def sortBamOutput(iontorrent, cram, temp):
    if iontorrent:
        return(snakemake.io.temp('tmp/{sample}_sorted.bam'))
    else:
        if (cram == 1 or temp):
            return(snakemake.io.temp('bam/{sample}.bam'))
        else: 
            return('bam/{sample}.bam')



## Merge RSeQC results, return paths to merged results
def mergeRSEQC(srcdir, run_read_dist):
    retvals = ['sample_metadata.csv','rseqc/bam_qc_parsed.tab', 'rseqc/bam_gc_parsed.tab', 'rseqc/bam_jc_parsed.tab', 'feature_counts/fragment_count_matrix.tab.gz','feature_counts/gene_lengths.tab']
    
    ## QC
    snakemake.shell("find rseqc | grep -P 'bam_qc.txt$' > bam_qc_outfiles.txt")
    snakemake.shell("perl {srcdir}/perl/parse-rseqc-bam-qc-results.pl bam_qc_outfiles.txt > rseqc/bam_qc_parsed.tab")
    snakemake.shell("rm bam_qc_outfiles.txt")

    ## GC
    snakemake.shell("find rseqc | grep -P 'bam_gc.txt$' > bam_gc_outfiles.txt")
    snakemake.shell("perl {srcdir}/perl/parse-rseqc-bam-gc-results.pl bam_gc_outfiles.txt > rseqc/bam_gc_parsed.tab")
    snakemake.shell("rm bam_gc_outfiles.txt")

    ## JC
    snakemake.shell("find rseqc | grep -P 'bam_jc.txt$' > bam_jc_outfiles.txt")
    snakemake.shell("perl {srcdir}/perl/parse-rseqc-bam-jc-results.pl bam_jc_outfiles.txt > rseqc/bam_jc_parsed.tab")
    snakemake.shell("rm bam_jc_outfiles.txt")
    
    ## RC
    if (run_read_dist) == 1:
        snakemake.shell("find rseqc | grep -P 'bam_rc.txt$' > bam_rc_outfiles.txt")
        snakemake.shell("perl {srcdir}/perl/parse-rseqc-read-distribution-results.pl bam_rc_outfiles.txt > rseqc/bam_rc_parsed.tab")
        snakemake.shell("rm bam_rc_outfiles.txt")
        retvals.append('rseqc/bam_rc_parsed.tab')

    ## featureCounts
    snakemake.shell("find feature_counts | grep -P '_count.tab$' > feature_counts_outfiles.txt")
    snakemake.shell("perl {srcdir}/perl/parse-read-count-matrix-subread-1.4.6.pl feature_counts_outfiles.txt > feature_counts/fragment_count_matrix.tab")
    snakemake.shell("gzip -f feature_counts/fragment_count_matrix.tab")
    
    ## Export gene lengths from a single featurecounts file 
    snakemake.shell("head -1 feature_counts_outfiles.txt | xargs tail -n +2 | awk '{{print $1,$6}}' > feature_counts/gene_lengths.tab")
    snakemake.shell("rm feature_counts_outfiles.txt")
    
    return(retvals)
