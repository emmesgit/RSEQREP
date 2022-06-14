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
# Program:  stage2/Snakefile.sh
# Version:  RSEQREP 2.3.0
# Author:   William F. Hooper, Travis L. Jensen, Johannes B. Goll
# Purpose:  Run second stage of preprocessing: run HISAT2 alignment, perform optional remap with Bowtie2,
#            generate QC metrics
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Import modules
import shutil
import logging as _logging

##############################
#        CONFIGURATION         #
##############################

## Specify YAML config file location
configfile: "preprocess_config.yaml"

## Program locations
RSEQC         = config["rseqc_dir"]
SAMTOOLS      = config["samtools_prog"]
HISAT2        = config["hisat_prog"]
HISAT2BUILD   = HISAT2 + '-build'
HISAT2SPLICE  = HISAT2 + '_extract_splice_sites.py'
BOWTIE2       = config["bowtie_prog"]
BOWTIE2BUILD  = BOWTIE2+'-build'
FEATURECOUNTS = config["fcts_prog"]
AWS              = config["aws_prog"]
OPENSSL      = config["openssl_prog"]
FASTQC       = config["fastqc_prog"]
CUTADAPT     = config["cutadapt_prog"]
FASTQDUMP    = config["fastqdump_prog"]
TRIMMOMATIC  = config["trimmomatic_prog"]

## Directories
SOURCEDIR     = config["srcdir"]
DATADIR       = config["datadir"]

## Use the source dir to import helper modules
sys.path.append(SOURCEDIR+'/python')
import trimadapters
import getfile  
import putfile  
import utils 
import sqlite
import time 

## Adapter sequences
FP_ADAPTERS   = [x.strip() for x in utils.toList(config["fp_adapter_seq"])]
TP_ADAPTERS   = [x.strip() for x in utils.toList(config["tp_adapter_seq"])]

## Ensembl version number
ENSEMBL = config["ensembl_version"]

## Logging config
LOG_LEVEL = config["log_level"]
LOG_FILE  = config["log_file"]
LOG_DB       = config["log_file"] + '.db'

## Save intermediate files?
REMOVEINTFILES = not config["saveintlocalfiles"]

## Reference datasets
INDEXSEQ = 'genome/Homo_sapiens.ensembl.version'+str(ENSEMBL)+'.genome.fa'
ANNOTATIONS_BED = 'annot/Homo_sapiens.ensembl.version'+str(ENSEMBL)+'.chr.bed'
ANNOTATIONS_GTF = 'annot/Homo_sapiens.ensembl.version'+str(ENSEMBL)+'.chr.gtf'

## Remap reads? 
IONTORRENT   = int(config["ion_torrent"]) == 1

## Generate CRAM files?
CRAM         = int(config["save_cram"]) 
CRAMFLAG     = ['bam','cram'][CRAM]

## run/not run certain steps
RUN_READ_DIST = int(config["run_read_dist"])
RUN_FASTQC = int(config["run_fastqc"])

## List of samples to process
SAMID        = utils.toList(config["samid"]) 

## List of input files
FASTQ_1      = utils.toList(config["fastq1"])
FASTQ_2      = utils.toList(config["fastq2"])

## Some steps allow paired end samples to be run separately
if (FASTQ_2 != ['']):
    ENDS  = ['1','2']
else:
    ENDS  = ['1']

## strandedness of experiment
STRANDED = int(config["stranded"]) 

## Determine whether adapters should be trimmed or not
TRIM_FP = sum([x == 'NA'  for x in FP_ADAPTERS]) == 0
TRIM_TP = sum([x == 'NA'  for x in TP_ADAPTERS]) == 0
TRIM_ADAPTERS_OUTPUT = '.fastq.gz' if (TRIM_FP or TRIM_TP) else '.skipped'

## Configure uploading 
ARCHIVE = config["archive_bucket"]
DOARCHIVE = len(ARCHIVE) != 0

## Configure encryption
DOENCRYPT    = int(config["encrypt_local"]) == 1
DECRYPT_PASS = config["decrypt_pass"]
ENCRYPT_PASS = config["decrypt_pass"]

## Configure quality trimming level
QUAL_CUTOFF = config["quality_trim"]

## number of cores to use
NCORES  = int(config["ncores"])

## hash for encryption/decryption & Checksums
HASH = config["hash"]

## Set up logging
_logging.basicConfig(level=LOG_LEVEL, 
                    format='[rseqrep] %(asctime)s - %(levelname)s - %(message)s', 
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    handlers=[_logging.FileHandler(LOG_FILE),
                             _logging.StreamHandler(sys.stdout)])


## Define final output
OUTPUT = [expand('rseqc/{sample}_bam_qc.txt', sample=SAMID),
    expand('rseqc/{sample}_bam_gc.txt', sample=SAMID),
    expand('rseqc/{sample}_bam_jc.txt', sample=SAMID),
    expand('feature_counts/{sample}_count.tab', sample=SAMID)]
if CRAM == 1 and REMOVEINTFILES:
    OUTPUT.append(expand('progress/{sample}_cram.done', sample=SAMID))
if not REMOVEINTFILES:
    OUTPUT.append(expand(CRAMFLAG+'/{sample}.'+CRAMFLAG, sample=SAMID))
if RUN_READ_DIST == 1:
    OUTPUT.append(expand('rseqc/{sample}_bam_rc.txt', sample=SAMID))
if RUN_FASTQC == 1:
    OUTPUT.append(expand('fastqc/{sample}_fastqc.tar.gz', sample=SAMID))


# RULE DEFINITIONS #
onstart:
    ## Initialize the SQL metadata database, add samples
    sqlite.initSqliteDb(db=LOG_DB)
    [sqlite.addSample(db=LOG_DB, sample_name=x) for x in SAMID]
onsuccess:
    ## merge featurecounts, rseqc results
    merged_results = utils.mergeRSEQC(SOURCEDIR, RUN_READ_DIST)
    
    ## Copy to datadir
    [shutil.copy2(x, DATADIR) for x in merged_results]
    
    ## Encrypt and/or upload merged results
    if DOENCRYPT: merged_results = [utils.encryptFile(file=x, openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH) for x in merged_results]
    if DOARCHIVE: [putfile.upload(file=x, destination=ARCHIVE, cloud='aws', prog=AWS) for x in merged_results]



## Define final output
rule all:
    input:
        OUTPUT
    
    
## Set up directory structure
## Ignores non-zero exit status returned when any directories already exist
rule directory_setup: 
    output: 
        'progress/dirs.done'
    threads:1
    run:
        cmd = "mkdir cram annot genome index input bam feature_counts rseqc fastqc tmp cutadapt rqual_filter progress benchmark -p 2> /dev/null"
        shell(cmd+utils.returnCode(process='Directory setup', log=LOG_FILE), shell=True)
        utils.touch(output)

    
    
## Build HISAT2 index
rule build_hisat_index:
    input:
        rules.directory_setup.output
    output:
        'progress/hisat2_index_built.done'
    benchmark:
        'benchmark/build_index.tab'
    threads: max(1,NCORES)
    priority: 1000
    params:
        benchmark='benchmark/build_index.tab.info'
    run:
        ## Run command
        cmd = '%s -p %s %s index/hisat2_genome_index && %s %s > annot/splicesites.txt' % (HISAT2BUILD, threads, INDEXSEQ, HISAT2SPLICE, ANNOTATIONS_GTF)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='HISAT2 Build', log=LOG_FILE), shell=True)    
        end_time = time.time()
        utils.touch(output)
        
        ## Add process and output file to output (just use the first sample in the sample list)
        samid = sqlite.getSamid(db=LOG_DB, sample_name=SAMID[0])
        procid = sqlite.addProcess(db=LOG_DB, samid=sqlite.getSamid(db=LOG_DB, sample_name=SAMID[0]), cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=SAMID[0], file_type='build_hisat_index', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)



## Build Bowtie2 index
rule build_bowtie_index:
    input:
        rules.directory_setup.output
    output:
        'progress/bowtie2_index_built.done'
    benchmark:
        'benchmark/bowtie2_index.tab'
    threads: max(1,NCORES)
    params:
        benchmark='benchmark/bowtie2_index.tab.info'
    run:
        ## Run command
        cmd = '%s -f %s index/bowtie2_genome_index --threads %s' % (BOWTIE2BUILD, INDEXSEQ, threads)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='Bowtie2 Build', log=LOG_FILE), shell=True)
        end_time = time.time()
        utils.touch(output)
        
        ## Add process and output file to output (just use the first sample in the sample list)
        samid = sqlite.getSamid(db=LOG_DB, sample_name=SAMID[0])
        procid = sqlite.addProcess(db=LOG_DB, samid=sqlite.getSamid(db=LOG_DB, sample_name=SAMID[0]), cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=SAMID[0], file_type='build_bowtie_index', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)


## Index reference genome with faidx (only necessary if CRAM compression is enabled)
rule faidx: 
    output:
        'progress/faidx.done'
    benchmark:
        'benchmark/faidx.tab'
    params:
        benchmark='benchmark/faidx.tab.info'
    threads: 1
    run:
        ## Run command
        cmd = '%s faidx %s' % (SAMTOOLS, INDEXSEQ)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='Bowtie2 Build', log=LOG_FILE), shell=True)
        end_time = time.time()
        utils.touch(output)
        
        ## Add process and output file to output (just use the first sample in the sample list)
        samid = sqlite.getSamid(db=LOG_DB, sample_name=SAMID[0])
        procid = sqlite.addProcess(db=LOG_DB, samid=sqlite.getSamid(db=LOG_DB, sample_name=SAMID[0]), cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=SAMID[0], file_type='faidx', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
    
    
    
## Get input
rule getfile:
    input:
        rules.directory_setup.output,
        rules.build_hisat_index.output
    output:
        temp(expand('input/{{sample}}_{pe}.fastq.gz',  pe=ENDS))
    params:
        sample='{sample}'
    priority: 1
    threads: max(1,min(8,NCORES))
    run:
        ## Get files to pull
        if ENDS == ['1']:
            in_file = [FASTQ_1[i[0]] for i in enumerate(SAMID) if SAMID[i[0]] == params.sample]
        else:
            in_file_1 = [FASTQ_1[i[0]] for i in enumerate(SAMID) if SAMID[i[0]] == params.sample][0]
            in_file_2 = [FASTQ_2[i[0]] for i in enumerate(SAMID) if SAMID[i[0]] == params.sample][0]
            in_file = [in_file_1, in_file_2] 
        
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        
        ## Run command
        for f in range(len(in_file)):
            start_time = time.time()
            return_code = getfile.getFile(in_file=in_file[f], out_file=utils.toList(output)[f], aws_prog=AWS, fastq_dump_prog=FASTQDUMP, openssl_prog=OPENSSL, pw=DECRYPT_PASS, hash=HASH)
            end_time = time.time()
            
            ## Add process and output file to output
            procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd='getfile.getFile()', return_code=return_code, start_time=start_time, end_time=end_time)
            sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output)[f], sample_name=params.sample, file_type='getfile', hash=HASH)
            


## Trim adapters with cutadapt
rule trimadapters:
    input:
        fa=expand('input/{{sample}}_{pe}.fastq.gz',  pe=ENDS)
    output:
        temp([x + TRIM_ADAPTERS_OUTPUT for x in expand('cutadapt/{{sample}}_{pe}', pe=ENDS)])
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_trim_adapters.tab.info'
    threads: max(1,min(8,NCORES))
    priority: 2
    benchmark:
        'benchmark/{sample}_trim_adapters.tab'
    run: 
        if TRIM_FP or TRIM_TP:
            adapters = trimadapters.getAdapters(sam=params.sample, fp_adapters=FP_ADAPTERS, tp_adapters=TP_ADAPTERS, samid_all=SAMID)
            
            ## Run command
            start_time = time.time()
            bench_obj = trimadapters.trimAdapters(infile=input.fa, outfile=output, adapt5=adapters['fp'], adapt3=adapters['tp'], cutadapt=CUTADAPT, threads=threads)
            end_time = time.time()
            
            ## Add process and output file to output
            samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
            procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd='trimadapters.trimAdapters()', return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
            sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='trimadapters', hash=HASH)
            sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
        else:
            utils.touch(output)
    
    
    
## Optionally quality-trim reads via Trimmomatic
rule qualityfilter:
    input:
        rules.trimadapters.output if TRIM_FP or TRIM_TP else expand('input/{{sample}}_{pe}.fastq.gz',  pe=ENDS)
    output:
    	temp(expand('rqual_filter/{{sample}}_{pe}{paired}_qual.fastq.gz', pe=ENDS, paired=['P','U'])) if len(ENDS)==2 else temp(expand('rqual_filter/{{sample}}_{pe}_qual.fastq.gz', pe=ENDS)) 
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_quality_filter.tab.info'
    threads: max(1,min(8,NCORES))
    priority: 4
    benchmark:
        'benchmark/{sample}_quality_filter.tab'
    run:
        ## Build command (if quality cutoff is 0 just copy files to their new locations)
        if (int(QUAL_CUTOFF) > 0):
            paired_end_flag = 'PE' if len(ENDS) == 2 else 'SE'
            cmd_input = ' '.join(utils.toList(input))
            cmd_output = ' '.join(utils.toList(output))
            cmd = 'java -jar %s %s -threads %s %s %s SLIDINGWINDOW:4:%s' % (TRIMMOMATIC, paired_end_flag, threads, cmd_input, cmd_output, QUAL_CUTOFF)                            
        else:
            ## Join together all copy commands
            cmd = ['cp %s %s' % (x,y) for x,y in zip(utils.toList(input), utils.toList(output))]
            cmd = '; '.join(cmd)
            
        ## Run command
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(sample='{params.sample}',process='Trimmomatic', log=LOG_FILE, outfile=utils.toList(output)[0]), shell=True)
        end_time = time.time()
        
        ## Add process and output file to output
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='qualityfilter', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)



## Map reads to the reference genome using HISAT2 
rule run_hisat:
    input:
        tch=rules.build_hisat_index.output,
        fa=rules.qualityfilter.output
    output:
        **utils.hisatOutput(ends=ENDS, iontorrent=IONTORRENT)
    benchmark:
        'benchmark/{sample}_run_hisat.tab'
    threads: max(1,min(8,NCORES))
    priority: 5
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_run_hisat.tab.info'
    run:
        ## Construct input based on read ended-ness
        if (len(ENDS) == 1):
            in_fa_str = "-U "+input.fa[0]
            
            if (IONTORRENT):
                unmapped_str = "--un-gz tmp/"+params.sample+"_iontorrent_1.fastq.gz"
            else:
                unmapped_str = ''
            
        elif (len(ENDS) == 2):
            in_fa_str = '-1 '+input.fa[0]+' -2 '+input.fa[2]
            
            if (IONTORRENT):
                unmapped_str = "--un-conc-gz tmp/"+params.sample+"_iontorrent_%.fastq.gz"
            else:
                unmapped_str = ''
        
        ## Run command
        cmd = '%s -p %s -x index/hisat2_genome_index %s %s --known-splicesite-infile annot/splicesites.txt -S %s' % (HISAT2, threads, in_fa_str, unmapped_str, output.aln)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='HISAT2 Align', sample='{params.sample}', log=LOG_FILE), shell=True)
        end_time = time.time()
        
        ## Add process and output file to output
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='hisat2 sam', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)        

rule run_bowtie:
    input:
        rules.build_bowtie_index.output,
        fa=expand('tmp/{{sample}}_iontorrent_{pe}.fastq.gz', pe=ENDS)
    output:
        temp('tmp/{sample}_iontorrent.bam')
    benchmark:
        'benchmark/{sample}_run_bowtie.tab'
    threads: max(1,min(8,NCORES))
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_run_bowtie.tab.info'
    run:
        ## Construct input based on read ended-ness
        if (len(ENDS) == 1):
            in_fa_str = "-U "+input.fa[0]
        elif (len(ENDS) == 2):
            in_fa_str = '-1 '+input.fa[0]+' -2 '+input.fa[1]
        
        ## Run command
        cmd = '%s --threads %s -x index/bowtie2_genome_index %s --sensitive-local | %s view -h -u -b | %s sort -@ %s - > %s' % (BOWTIE2, threads, in_fa_str, SAMTOOLS, SAMTOOLS, threads, output)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='Bowtie2 Remapping', sample='{params.sample}', log=LOG_FILE), shell=True)
        end_time = time.time()
        
        ## Add process and output file to output 
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='sorted bam', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)        



## Convert SAM to BAM
rule sam_to_bam:
    input:
        rules.run_hisat.output.aln
    output:
        temp('tmp/{sample}.bam')
    benchmark:
        'benchmark/{sample}_sam_to_bam.tab'
    priority: 5
    threads: max(1,min(8,NCORES))
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_sam_to_bam.tab.info'
    run:
        ## Run command
        cmd = '%s view -@ %s -Sbh %s > %s' % (SAMTOOLS, threads, input, output)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='SAM to BAM', sample='{params.sample}', log=LOG_FILE, outfile=utils.toList(output)[0]), shell=True)
        end_time = time.time()
        
        ## Add process and output file to output 
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='HISAT2 mapped and unmapped reads BAM file', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
        


## Sort BAM
rule sort_bam:
    input:
        bam=rules.sam_to_bam.output,
    output:
        utils.sortBamOutput(iontorrent=IONTORRENT, cram=CRAM, temp=REMOVEINTFILES)
    benchmark:
        'benchmark/{sample}_sort_bam.tab'
    priority: 5
    threads: max(1,min(8,NCORES))
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_sort_bam.tab.info'
    run:
        ## Run command
        cmd = '%s sort -@ %s %s > %s' % (SAMTOOLS, threads, input.bam, output)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='Sort & Index BAM', sample='{params.sample}', log=LOG_FILE, outfile=utils.toList(output)[0]), shell=True)
        end_time = time.time()
        
        ## Add process and output file to output 
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='sorted bam', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
        

## Merge sorted alignments if we re-mapped with bowtie
rule merge_bam:
    input:
        bowtie=rules.run_bowtie.output,
        hisat=rules.sort_bam.output
    output:
        temp('bam/{sample}.bam') if (CRAM==1) else 'bam/{sample}.bam'
    benchmark:
        'benchmark/{sample}_merge_bam.tab'
    threads: 1
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_merge_bam.tab.info'
    run:
        ## Run command
        cmd = '%s merge %s %s %s' % (SAMTOOLS, output, input.hisat, input.bowtie)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='Merge BAM', sample='{params.sample}', log=LOG_FILE, outfile=utils.toList(output)[0]), shell=True)
        end_time = time.time()
        
        ## Add process and output file to output 
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='merged sorted bam', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
        
        
## Index BAM
rule index_bam:
    input:
        bam='bam/{sample}.bam'
    output:
        temp('bam/{sample}.bam.bai') if REMOVEINTFILES else 'bam/{sample}.bam.bai'
    benchmark:
        'benchmark/{sample}_index_bam.tab'
    priority: 5
    threads: max(1,min(8,NCORES))
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_index_bam.tab.info'
    run:
        ## Run command
        cmd = '%s index -@ %s %s' % (SAMTOOLS, threads, input.bam)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='Index BAM', sample='{params.sample}', log=LOG_FILE, outfile=utils.toList(output)[0]), shell=True)
        end_time = time.time()
        
        ## Add process and output file to output 
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='index bam', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)


## Convert bam to cram
rule bam_to_cram:
    input:
        bam='bam/{sample}.bam',
        faidx=rules.faidx.output
    output:
        cram=temp('cram/{sample}.cram') if REMOVEINTFILES else 'cram/{sample}.cram',
        prog=touch('progress/{sample}_cram.done')
    benchmark: 
        'benchmark/{sample}_bam_to_cram.tab'
    threads: max(1,min(8,NCORES))
    priority: 5
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_bam_to_cram.tab.info'
    run:
        ## Run command
        cmd = '%s view -C -@ %s -T %s -o %s %s' % (SAMTOOLS, threads, INDEXSEQ, output.cram, input.bam)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='CRAM Compression', sample='{params.sample}', log=LOG_FILE, outfile=utils.toList(output.cram)[0]), shell=True)
        end_time = time.time()
        
        ## Add process and output file to output 
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output.cram), sample_name=params.sample, file_type='cram', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
        
        ## Encrypt and/or upload if necessary
        if DOENCRYPT: output.cram = utils.encryptFile(file=utils.toList(output.cram), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output.cram), destination=ARCHIVE, cloud='aws', prog=AWS)
        if DOARCHIVE: putfile.upload(file=utils.toList(INDEXSEQ), destination=ARCHIVE, cloud='aws', prog=AWS, rm=False)    

## Run FASTQC
rule fastqc: 
    input:
        'bam/{sample}.bam'
    output:
        temp('fastqc/{sample}_fastqc.tar.gz') if REMOVEINTFILES else 'fastqc/{sample}_fastqc.tar.gz'
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_fastqc.tab.info'
    threads: 1
    priority: 5
    benchmark:
        'benchmark/{sample}_fastqc.tab'
    run:
        ## Determine intermediate files
        fq_base = 'fastqc/%s_fastqc' % (params.sample)
        fq_zip  = fq_base + '.zip'
        fq_html = fq_base + '.html'
        
        ## Run command
        cmd = '%s %s -q -o fastqc' % (FASTQC, input)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(sample='{params.sample}',process='FastQC', log=LOG_FILE),shell=True)
        end_time = time.time()
        
        ## Add process metadata to DB
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
        
        ## unzip, remove zipped results, HTML duplicate, and tarball results
        cmd = 'unzip %s -d %s && tar -zvcf %s %s && rm -r %s %s %s' % (fq_zip, fq_base, output, fq_base, fq_zip, fq_html, fq_base)
        utils.logging_call(cmd, shell=True)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='FastQC', hash=HASH)
        
        ## Encrypt and/or upload if necessary
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)



## Run RSEQC bam_stat.py
rule bam_qc:
    input:
        bam='bam/{sample}.bam',
        idx=rules.index_bam.output
    output:
        'rseqc/{sample}_bam_qc.txt'
    benchmark:
        'benchmark/{sample}_bam_qc.tab'
    threads: 1
    priority: 5
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_bam_qc.tab.info'
    run:
        ## Run command
        cmd = '%s/bam_stat.py -i %s > %s' % (RSEQC, input.bam, output)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='BAM QC', sample='{params.sample}', log=LOG_FILE), shell=True)    
        end_time = time.time()
        
        ## Add process and output file to output 
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='rseqc bam statistics', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
        
        ## Encrypt and/or upload if necessary
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)        
    
            
            
## Run RSEQC read_gc.py 
rule bam_gc:
    input:
        bam='bam/{sample}.bam',
        idx=rules.index_bam.output
    output:
        r=temp('rseqc/{sample}.GC_plot.r') if REMOVEINTFILES else 'rseqc/{sample}.GC_plot.r',
        txt='rseqc/{sample}_bam_gc.txt'
    benchmark:
        'benchmark/{sample}_bam_gc.tab'
    threads: 1
    priority: 5
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_bam_gc.tab.info'
    run:
        ## Run RSEQC read_GC
        cmd = '%s/read_GC.py -i %s -o rseqc/%s' % (RSEQC, input.bam, params.sample)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='BAM GC', sample='{params.sample}', log=LOG_FILE), shell=True)    
        end_time = time.time()
        
        ## Add some extra code to output R script, run it
        cmd = """echo "out=as.vector(summary(gc));dta = data.frame('%s',out[1],out[2],out[3],out[4],out[5],out[6]);write.table(dta,file='%s',sep='\t',row.names=F,col.names=F,quote=F);" >> %s""" % (params.sample, output.txt, output.r)
        utils.logging_call(cmd, shell=True)
        utils.logging_call('Rscript --vanilla --quiet '+output.r+utils.returnCode(process='BAM GC Rscript', sample='{params.sample}', log=LOG_FILE), shell=True)
        
        ## Add process and output file to output 
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='rseqc bam gc%', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
        
        ## Encrypt and/or upload if necessary
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        
        ## Remove intermediate read_GC files after s3 upload
        if REMOVEINTFILES:
            cmd = 'rm rseqc/%s.GC.xls rseqc/%s.GC_plot.pdf' % (params.sample, params.sample)
            utils.logging_call(cmd, shell=True)

## Run RSEQC junction_annotation.py 
rule bam_jc:
    input:
        bam='bam/{sample}.bam',
        idx=rules.index_bam.output
    output:
        'rseqc/{sample}_bam_jc.txt'
    benchmark:
        'benchmark/{sample}_bam_jc.tab'
    threads: 1
    priority: 5
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_bam_jc.tab.info'
    run:
        ## Run command
        cmd = '%s/junction_annotation.py -i %s -o rseqc/%s -r %s 2> %s' % (RSEQC, input.bam, params.sample, ANNOTATIONS_BED, output)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='BAM JC', sample='{params.sample}', log=LOG_FILE),shell=True)    
        end_time = time.time()
        
        ## Add process and output file to output 
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='rseqc bam junctions', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
        
        ## Encrypt and/or upload if necessary
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)

## Run RSEQC junction_annotation.py 
rule read_distribution:
    input:
        bam='bam/{sample}.bam',
        idx=rules.index_bam.output
    output:
        'rseqc/{sample}_bam_rc.txt'
    benchmark:
        'benchmark/{sample}_bam_rc.tab'
    threads: 1
    priority: 5
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_bam_rc.tab.info'
    run:
        ## Run command
        cmd = '%s/read_distribution.py -i %s -r %s > %s' % (RSEQC, input.bam, ANNOTATIONS_BED, output)
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='BAM Read Distribution', sample='{params.sample}', log=LOG_FILE),shell=True)    
        end_time = time.time()
        
        ## Add process and output file to output 
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='reseqc bam read distribution', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
        
        ## Encrypt and/or upload if necessary
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        
## Run featureCounts
rule feature_counts:
    input:
        'bam/{sample}.bam'
    output:
        'feature_counts/{sample}_count.tab'
    benchmark:
        'benchmark/{sample}_feature_counts.tab'
    threads: max(1,min(8,NCORES))
    priority: 5
    params:
        sample='{sample}',
        benchmark='benchmark/{sample}_feature_counts.tab.info'
    run:
        ## Run command - if paired add args
        if len(ENDS)>1:
            cmd = '%s -B -p -C -T %s -s %s -a %s -o %s %s' % (FEATURECOUNTS, threads, STRANDED, ANNOTATIONS_GTF, output, input)
        else:
            cmd = '%s -T %s -s %s -a %s -o %s %s' % (FEATURECOUNTS, threads, STRANDED, ANNOTATIONS_GTF, output, input)
            
        start_time = time.time()
        bench_obj = utils.logging_call(cmd+utils.returnCode(process='FeatureCounts', sample='{params.sample}', log=LOG_FILE),shell=True)
        end_time = time.time()
        
        ## Add process and output file to output 
        samid = sqlite.getSamid(db=LOG_DB, sample_name=params.sample)
        procid = sqlite.addProcess(db=LOG_DB, samid=samid, cmd=cmd, return_code=bench_obj.return_code, start_time=start_time, end_time=end_time)
        sqlite.addFile(db=LOG_DB, samid=samid, procid=procid, file=utils.toList(output), sample_name=params.sample, file_type='feature counts', hash=HASH)
        sqlite.addBenchmark(db=LOG_DB, procid=procid, samid=samid, bench_obj=bench_obj)
        
        ## Encrypt and/or upload if necessary
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        
