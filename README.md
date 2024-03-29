# RSEQREP

*RNA-Seq Reports, an open-source cloud-enabled framework for reproducible RNA-Seq data processing, analysis, and result reporting*

## INSTALLATION
 
### Option 1 (local Ubuntu server): 

* install LTS Ubuntu desktop version 18.04.02 on your machine (http://releases.ubuntu.com/18.04.02/, install from bootable USB)
* copy/clone the RSEQREP github source code to your Ubuntu machine (git clone https://github.com/emmesgit/RSEQREP.git). Ensure read write and execute permissions are set for RSEQREP (chmod -R u+rwx RSEQREP).
* execute our installation shell script to install the software on your own Ubuntu machine (sh RSEQREP/ubuntu/install-software-v2.0.0.sh)

### Option 2 (RSEQREP AWS AMI):

* initialize RSEQREP AMI (https://aws.amazon.com, AMI ID: RPREP RSEQREP (Ribosome Profiling and RNA-Seq Reports) v2.0 (ami-0dfef5a604d41a286)).  To do this, create a AWS account (https://aws.amazon.com), log into the AWS console, and navigate to the EC2 resources.  Next select AMIs in the navigation pane and select public images.  Finally search for RSEQREP, find "RPREP RSEQREP (Ribosome Profiling and RNA-Seq Reports) v2.0" and launch the AMI (ensure you are in the US EAST (N. Virginia) region).
* using an ssh command line connection or X2GO GUI supported connection client (https://wiki.x2go.org/doku.php/download:start) (XFCE session type) log into the RSEQREP AMI using username: repuser and password: repuser2019. The IP address of your instance can be found on the AWS console EC2->Instances->"IPv4 Public IP".  Please note this IP may change every time the machine is started/stopped.
* copy/clone the RSEQREP github source code to your Ubuntu machine (git clone https://github.com/emmesgit/RSEQREP.git).  Ensure read write and execute permissions are set for RSEQREP (chmod -R u+rwx RSEQREP).

Additional information on AWS and AMI configuration can be found in RSEQREP/aws/aws_instructions.docx

### Option 3 (RSEQREP Docker Image):

* On an Ubuntu 18.04.02 machine with Docker software installed (https://docs.docker.com/install), pull RSEQREP image from the Docker repository (docker pull emmesdock/rseqrep) [https://hub.docker.com/r/emmesdock/rseqrep/]. Run the docker image as a container in interactive mode (docker run --name rseqrep1 -i -t emmesdock/rseqrep /bin/bash).
* copy/clone the RSEQREP github source code to the container (git clone https://github.com/emmesgit/RSEQREP.git).  Ensure read write and execute permissions are set for RSEQREP (chmod -R u+rwx RSEQREP).

## EXECUTION

* Download GMT formatted gene sets for pathway enrichment (if no gene sets are specified, pathway enrichment is not performed).  We provide a custom script (https://github.com/emmesgit/RSEQREP/blob/master/source/shell/download-gene-sets.sh) to download Blood Transcription Modules (option btm), Reactome pathways (reactome option), and KEGG pathways (kegg option). Note for the KEGG option you need to either be an academic user or have a commercial KEGG license (http://www.kegg.jp/kegg/rest).  Gene sets can also be downloaded from MSigDB after email registration (http://software.broadinstitute.org/gsea/msigdb).  Specify the locations and labels of the gene sets in the gmt_entrez_files and gmt_entrez_files_labels fields of the worklfow_config tab of the configuration file (RSEQREP/config/config.xlsx).
* fill out the configuration file (RSEQREP/config/config.xlsx).  RSEQREP/case-study/config-henn.xlsx is an example of a complete configuration file.
* execute start-to-end analysis (sh RSEQREP/run-all.sh) or for a particular component (sh RSEQREP/run-pre-processing.sh; sh RSEQREP/run-analysis.sh; sh RSEQREP/run-report.sh). Specify at minimum the number of threads (-@ <threads>) and the full path to the configuration file (-c <config-file-path>) after the shell script name (ex: sh RSEQREP/run-all.sh -@ 32 -c /home/repuser/config.xlsx). 
 
## TROUBLESHOOTING

Upon inspecting the initial run of the report, you may find that the configuration option that you initially chose does not fit your data.  For example, you inspect the reverse cumulative distribution function plot comparing log count per million cutoffs with the number of retained genes.  You find that the cutoff you had originally selected resulted in too few genes for the analysis.  To remedy this, you would update the configuration file to reflect a more appropriate log counts per million cutoff.  You then determine the steps of the analysis that will be affected by this configuration change.  In this case, the analysis and report steps are affected.  Re-execute the analysis and report steps (sh RSEQREP/run-analysis.sh; sh RSEQREP/run-report.sh).  The configuration file is re-parsed each time a RSEQREP/run-* script is executed. Please see https://github.com/emmesgit/RSEQREP/wiki for additional information.
 
## BUGS/ISSUES

If you identify bugs/issues, please submit them via the GitHub Issue Tracker 
https://github.com/emmesgit/RSEQREP/issues
 
## ADDITIONAL RESOURCES

For a detailed explanation of the software and its capabilities please navigate to our F1000 research publication:
https://f1000research.com/articles/6-2162/v2

## RELEASE NOTES 

### RSEQREP RNA-Seq Reports - Version 2.3.0

#### New Features:

* Utilizing multicore functionality of cutadapt within the trimadapters rule.  
* Added a --dryrun arg which can passed to the pre-processing execution to print the list of jobs to be completed without actually executing the workflow.

### RSEQREP RNA-Seq Reports - Version 2.2.1

#### Bug Fixes:

* typo in the invocation of rseqrep python script caused failure because the number of cores not properly set -- "cores" used rather than intended "ncores". This is now repaired.  
* The paired end functionality of quality trimming with trimmomatic using v39 is not working properly.  This is now repaired.
* Updated the rule qualityfilter in the snakemake worflow to work correctly with both v38 and v39 trimmomatic.

### RSEQREP RNA-Seq Reports - Version 2.2.0

#### Bug Fixes:

* "ncores" specified in the workflow configuration was not correctly reflecting the number of cores actually used by snakemake on the backend.  The number of cores now must be specifed through the command line (--threads / -@). "ncores" is no longer supported in the configuration file. if used, the program will write an error to alert the user.
* If the remove intermediate files flag was specified in earlier versions, the reference genome files would be removed after the first sample, making them unavaliable for the remaining samples.  This is now repaired.
* Between v2.0.1 and v2.1.0, the paired end functionality of the reature_counts rule had lost paired end functionality.  This functionality is now restored.
* The strandedness option provided in the configuration file was never properly linked to the feature_counts rule in snakemake.  This functionality is now working properly.
* added "import trimadapters" into snakemake.sh script to properly execute trimming.
* added indexing for sorted bam files to remove non-impacting error message from RSEQC.

### RSEQREP RNA-Seq Reports - Version 2.1.2

#### Bug Fixes:

* OpenSSL encryption hash and checksum hash can now be specified as md5 or sha256 in the configuration xlsx file.  previously, md5 was hard-coded.

### RSEQREP RNA-Seq Reports - Version 2.1.0

#### New Features:

* Pathway enrichment, unsupervised gene clustering, and DE gene identification are now using threading via foreach/doparallel R packages.
* An addtional flag has been added to reduce the local footrint.  By default, this is active (intermediate files are removed).
* Added flag to run/not run RSEQC read distribution and FASTQC preprocessing steps.
* Added "ncores" variable for preprocessing steps (Samtools, Hisat2, Feature Counts) to specify the number of cores to thread across for pre-processing steps. Max value should not exceed avaliable cores.
* Added "ncores" variable for analysis steps (DE Gene Identification, Gene Clustering, Pathway Enrichment) to specify the number of cores to thread across for analysis steps. Max value should not exceed avaliable cores.

#### Bug Fixes:

* A previous hard coded trimmomatic flag "--discard-untrimmed" has been removed. With this flag, reads without the adapter would be thrown out.  This is undesired default behavior. 
* Repaired issue with spcLabl being used instead of the correct trtLabl in upset-trt-up-down.r.
* updated configuration file variable "glm\_model\_lfc\_cutoff" to "glm\_model\_fc\_cutoff" to correctly represent how the variable is processed in the workflow. The variable has always been processed on the original scale and reflects the fold change.

### RSEQREP RNA-Seq Reports - Version 2.0.1

#### Bug Fixes:
* Feature Counts has been updated to correctly handle paired end reads. Pritor to this fix, featureCounts only runs in single end mode, even if 2 fasta files are provided in the sample metadata file.

### RSEQREP RNA-Seq Reports - Version 2.0.0

#### New Features:

* New Amazon web services machine image based on ubuntu 18.04.02 LTS (<INSTANCE_ID>).  The updated software versions and licenses can be found in SOFTWARE.xlsx.  For simple user access to this machine: username: repuser password: repuser2019
* The pre-processing step has been re-implemented in Snakemake with a command-line interface. 
* Optional adapter trimming (via Cutadapt) and quality trimming (via Trimmomatic) have been added to the pre-processing step. 
* Optional read re-mapping (via Bowtie2) has been added to the pre-processing step. 
* Pre-processing files can now be encrypted and/or archived to AWS S3 cloud storage.
* Please note RSEQREP v2.0.0 is not compatible with the v1.0 AMI/Docker images and RSEQREP v1.X.X is not compatible with the v2.0 AMI/Docker image.  Please use the newer v2.0 publicly available AMI ID: RPREP RSEQREP (Ribosome Profiling and RNA-Seq Reports) v2.0 (ami-0dfef5a604d41a286)

#### Bug Fixes:

* For differing number of sets in the UpSet plots (both DE and GSEA) -- these were not plotting as expected.  UpSet plots now use upset.custom function to plot all sets with at least one significant gene/gene set.

#### Known Bugs:

* CPU time is incorrectly reported as 0 by pre-processing steps

### RSEQREP RNA-Seq Reports - Version 1.1.3

#### New Features:

* Multiplexed fastq files can be merged by specifying a semicolon separated list in the configuration file fastq fields. This will work for all download/decryption methods (Amazon S3, local, and SRA Archive).

#### Bug Fixes:

* Previously unsorted pathway enrichment flat files are now sorted by FDR adjusted P-value (ascending order) followed by Jaccard similarity index (descending order).
* Using an empty sting as input for GMT file locations/names in the configuration file for the pathway enrichment analysis would break the program.  An empty string for GMT file locations can now be supplied. In that case no pathway enrichment analysis will take place.
* SOFTWARE.xlsx link repaired in program headers.

### RSEQREP RNA-Seq Reports - Version 1.1.2

* Updated README file to include information about how to download external gene set information for pathway enrichment analysis and how to include this information in the RSEQREP configuration file

### RSEQREP RNA-Seq Reports - Version 1.1.1

#### Bug Fixes:

* Previously unsorted feature counts matrix is now sorted by ensembl gene ID.
* The code to generate Upset plots summarizing DE genes (reference code)  was not correctly reading input matrices
* A missing "}" was added to RSEQREP/source/r/02-sdeg-identification/upset-trt-up-down.r.
* Benchmark barplots in the report (RSEQREP/source/r/preprocessing-benchmarks-plot.r) now utilize only SQLite database entries with a return code of 0.

### RSEQREP RNA-Seq Reports - Version 1.1.0

#### New Features:

* The run-*.sh scripts can now be executed from anywhere on the file system.
* The download gene sets shell script (RSEQREP/source/shell/download-gene-sets.sh) was updated to support command line arguments for selectively download Blood Transcription Modules (option btm), Reactome pathways (reactome option), and KEGG pathways (kegg option). Note for the KEGG option you need to either be an academic user or have a commercial KEGG license (http://www.kegg.jp/kegg/rest).

#### Bug Fixes:

* SRA file download was switched from FTP URLs to use fastq-dump, a program that is as part of the SRA toolkit.
* Quotation of fields (double quotes) in .csv files for pre-processing steps generated when running RSEQREP/source/r/parse-rnaseq-configuration.r was removed.
* Result and preprocessing directories will only be created if they do not already exist.
* UpSet plots were not showing data if there was an equal number of up and down-regulated DE genes for a particular specimen type, treatment group, or timepoint.
* Two hidden line breaks were repaired in the Henn et.al. configuration file under case-study/config-henn.xlsx.

### RSEQREP RNA-Seq Reports - Version 1.0.0

* Initial release.

 
