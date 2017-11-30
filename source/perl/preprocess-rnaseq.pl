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
# https://github.com/emmesgit/RSEQREP/SOFTWARE.xlsx  
# 
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# To cite this software, please reference doi:10.12688/f1000research.10464.1
#
# Program:  preprocess-rnaseq.pl 
# Version:  RSEQREP 1.0.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:	Clinical RNA-Seq data preprocessing pipeline
# Input:	1) Workflow configuration file
#			2) Sample metadata file
# Output: 	cram/<sample_id>.cram
#			fastqc/<sample_id>_fastqc.tar.gz
#			rseqc/<sample_id>_bam_<qc|rc|jc|gc>.txt
#			feature_counts/
#				<sample_id>_count.tab
#				<sample_id>_count.tab.summary
#			star/<sample_id>.star.*
#			done/*.done
#############################################################################################################

## Specify packages
use strict;
use warnings;
use DBI;	
use File::Basename;
use POSIX qw(:sys_wait_h);
$|++;

## input file locations/names from command line
my $worklowConfigFile		= $ARGV[0];
my $sampleMetadataFile		= $ARGV[1];

## parse workflow configuration file
## initialize container for configuration variables
my %config;
open(my $cfg, '<', $worklowConfigFile) or die "Cannot open workflow configuration file: $!";
while (my $csvkvpair = <$cfg>) {
	# remove trailing whitespace
	$csvkvpair =~ s/\s*\z//;
	# split key and value -- store in an array
	my @kv = split(/,/, $csvkvpair);
	# store key/value
	my $k = $kv[0];
	my $v = $kv[1];
	# hash referenced variable and value
	$config{$k} = $v;
}
	
## Make main, tmp, and out directories if they do not exist
my $tmpDir = "$config{'pre_dir'}/tmp";
my $outDir = "$config{'pre_dir'}";
unless(-d $config{'pre_dir'}) {`mkdir $config{'pre_dir'}`;}
unless(-d $outDir) {`mkdir $outDir`;}
unless(-d "$outDir/cram") {`mkdir $outDir/cram $outDir/star $outDir/rseqc $outDir/fastqc $outDir/feature_counts $outDir/done`;}
unless(-d $tmpDir) {`mkdir $tmpDir`;}

## open log file -- Always append data to existing log
my $logFile = "$outDir/preprocess.log";
unless (-e "$logFile") {`touch $logFile`}
open(LOG, ">>",$logFile);

## Initalize database connection -- Always append data to database
my $dbh = DBI->connect("dbi:SQLite:dbname=$outDir/preprocess.db","","");
initDb();
	
## iterate through sample metadata file
open(my $mta, '<', $sampleMetadataFile) or die "Cannot open sample metadata file: $!";
while (my $entry = <$mta>) {
	
	# dont read header
	unless ($entry=~/subid,samid/) {
		
		# remove newline
		chomp($entry);
		
		# assign variables to entries
		my ($subid,$samid,$tp,$tpl,$tpc,$tpb,$spc,$spcl,$spcc,$trt,$trtcl,$trtc,$fq1,$fq2)= split(/,/, $entry);
	
		## define new file names;
		my $bamEncFile		= "$tmpDir/$samid.bam.enc";
		my $bamFile			= "$tmpDir/$samid.bam";
		my $fq1File			= "$tmpDir/$samid"."_1.fastq.gz";
		my $fq2File 		= "$tmpDir/$samid"."_2.fastq.gz";
		my $cramFile		= "$outDir/cram/$samid.cram";
		my $starFile		= "$outDir/star/$samid.star.";
		my $spliceFile		= "$outDir/annot/splicesites.txt";
		my $bamQcFilePrfx	= "$outDir/rseqc/$samid";
		my $bamQcFile		= "$outDir/rseqc/$samid"."_bam_qc.txt";
		my $bamJcQcFile		= "$outDir/rseqc/$samid"."_bam_jc.txt";
		my $bamRdQcFile		= "$outDir/rseqc/$samid"."_bam_rc.txt";
		my $bamGcQcFile		= "$outDir/rseqc/$samid"."_bam_gc.txt";
		my $bamFastQcDir	= "$outDir/fastqc/$samid"."_fastqc";
		my $cntFile 		= "$outDir/feature_counts/$samid"."_count.tab";
		my $doneDir 		= "$outDir/done/$samid";
		my $sampleId = 		addSample($dbh,$samid);
		
		# If not specified, replace with empty string
		if ($fq2 eq '') {
			$fq2File = '';
		}
		
		# if the file has "s3://" at the beginning of the file name, download it from the s3 bucket.
		# we assume paired reads will be tored in the same place -- Both S3 or both local
		
		if($fq1=~/^s3:\/\//) {
			my $down_done_1 = "$doneDir"."_download_fq1.done";
			unless(-e $down_done_1) {
				my $return_down1 = downloadFileFromS3($dbh,$sampleId,$fq1,"$fq1File.enc",$tmpDir,$config{'aws_prog'});
				if ($return_down1 == 0) {
					`touch $down_done_1`;
				}
			}
			my $down_done_2 = "$doneDir"."_download_fq2.done";
			unless(-e $down_done_2 || $fq2 eq '') {
				my $return_down2 = downloadFileFromS3($dbh,$sampleId,$fq2,"$fq2File.enc",$tmpDir,$config{'aws_prog'});
				if ($return_down2 == 0) {
					`touch $down_done_2`;
				}
			}
		} 
		
		# if the file ends in sra, we will attempt to download it from the SRA FTP site
		if ($fq1=~/sra$/) {
			my $sra_done_1 = "$doneDir"."_download_sra.done";
			unless(-e $sra_done_1) {
				my $return_down1 = downloadFileFromSra($dbh,$sampleId,$fq1,$fq1File,$tmpDir,$config{'fastqdump_prog'});
				if ($return_down1 == 0) {
					`touch $sra_done_1`;
				}
			}
			my $sra_done_2 = "$doneDir"."_download_sra.done";
			unless(-e $sra_done_2 || $fq2 eq '') {
				my $return_down2 = downloadFileFromSra($dbh,$sampleId,$fq2,$fq2File,$tmpDir,$config{'fastqdump_prog'});
				if ($return_down2 == 0) {
					`touch $sra_done_2`;
				}
			}
		}
		
		# if the file is encrypted -- ending in "enc" or "encr" -- decrypt it.
		if ($fq1=~/enc[r]?$/) {
			my $decrypt_done_1 = "$doneDir"."_decrypt_fq1.done";
			unless (-e $decrypt_done_1) {
				my $return_decrypt1 = decryptFile($dbh,$sampleId,"$fq1File.enc",$fq1File,0,$tmpDir,$config{'openssl_prog'},$config{'decrypt_pass'});
				if ($return_decrypt1 == 0) {
					`touch $decrypt_done_1`;
				}
			}
			my $decrypt_done_2 = "$doneDir"."_decrypt_fq2.done";
			unless(-e $decrypt_done_2 || $fq2 eq '') {
				my $return_decrypt2 = decryptFile($dbh,$sampleId,"$fq2File.enc",$fq2File,0,$tmpDir,$config{'openssl_prog'},$config{'decrypt_pass'});
				if ($return_decrypt2 == 0) {
					`touch $decrypt_done_2`;
				}
			}	
		} else {
			`cp $fq1 $fq1File`;
			unless ($fq2 eq '') {
				`cp $fq2 $fq2File`;
			}
		}
		
		
		
		# perform STAR read assembly (indexing performed if needed)
		my $star_done = "$doneDir"."_star.done";
		if(defined($config{'star_prog'})) {
			unless (-e $star_done) {
				my $star_index_done = "$outDir/done/star_index.done";
				unless (-e $star_index_done) {
					my $return_star_index = starIndex($dbh,$config{'num_threads'},$tmpDir,$sampleId,$config{'star_prog'},dirname($config{'genome_file'}),$config{'genome_file'},$config{'gtf_file'});	
					if ($return_star_index == 0) {
						`touch $star_index_done`;
					}
				}
				my $return_star = star($dbh,$fq1File,$fq2File,$config{'num_threads'},$starFile,$bamFile,$tmpDir,$sampleId,$config{'star_prog'},dirname($config{'genome_file'}));
				if ($return_star == 0) {
					`touch $star_done`;
				}
			}
		}
		
		# perform HISAT read assembly (indexing performed if needed)
		my $hisat_done = "$doneDir"."_hisat.done";
		if(defined($config{'hisat_prog'})) {
			unless (-e $hisat_done) {
				my $hisat_index_done = "$outDir/done/hisat_index.done";
				unless (-e $hisat_index_done) {
					my $return_hisat_index = hisatIndex($dbh,$config{'num_threads'},$tmpDir,$sampleId,$config{'hisat_prog'},dirname($config{'genome_file'}),$config{'genome_file'},$config{'gtf_file'},$spliceFile);	
					if ($return_hisat_index == 0) {
						`touch $hisat_index_done`;
					}
				}
				my $return_histat = hisat($dbh,$fq1File,$fq2File,$config{'num_threads'},$bamFile,$tmpDir,$sampleId,$config{'hisat_prog'},dirname($config{'genome_file'}),$spliceFile,$config{'samtools_prog'});
				if ($return_histat == 0) {
					`touch $hisat_done`;
				}
			}
		}
		
		# create the cram file
		unless (-e $cramFile) {
			if ($config{'save_cram'} == 1) {
				compress_bam($dbh,$config{'genome_file'},$config{'num_threads'},$tmpDir,$sampleId,$cramFile,$config{'samtools_prog'},$bamFile);
			}
		}
		
		# Clean up fastq files
		if (-e "$fq1File.enc") {print LOG "rm $fq1File.enc\n";`rm $fq1File.enc`;}
		if (-e "$fq2File.enc") {print LOG "rm $$fq2File.enc\n";`rm $fq2File.enc`;}
		if (-e $fq1File) {print LOG "rm $fq1File\n";`rm $fq1File`;}
		if (-e $fq2File) {print LOG "rm $fq2File\n";`rm $fq2File`;}
		
		## if a BAM file does not exist, make it from the CRAM file
		## This is useful if there was an error in the run, and the process is re-run.
		unless (-e $bamFile) {
			`$config{'samtools_prog'} view -@ $config{'num_threads'} --output-fmt bam --reference $config{'genome_file'} $cramFile > $bamFile`;
			addFile($dbh,$sampleId,"","$bamFile",basename($bamFile),md5sum("$bamFile"),bytes("$bamFile"),'Mapped Reads Bam file')
		}
		
		## execute FastQC on bam file
		my $fqc_done = "$doneDir"."_fastqc.done";
		unless (-e $fqc_done) {
			my $return_fqc = executeFastQC($dbh,"$outDir/fastqc",$sampleId,$bamFile,$tmpDir,$bamFastQcDir,$config{'fastqc_prog'});
			if ($return_fqc == 0) {
				`touch $fqc_done`;
			}
		}
		
		## execute Read distribution
		my $rd_done = "$doneDir"."_read_dist.done";
		if ($config{'run_read_dist'} == 1) {
			unless (-e $rd_done) {
				my $return_rd = executeBamReadDistribution($dbh,$sampleId,$bamFile,$config{'bed_file'},$bamRdQcFile,$tmpDir,"$config{'rseqc_dir'}/read_distribution.py");
				if ($return_rd == 0) {
					`touch $rd_done`;
				}
			}	
		}
		
		## execute BAM STAT QC
		my $st_done = "$doneDir"."_stat_qc.done";
		unless (-e $st_done) {
			my $return_st = executeBamStatQc($dbh,$sampleId,$bamFile,$bamQcFile,$tmpDir,"$config{'rseqc_dir'}/bam_stat.py");
			if ($return_st == 0) {
				`touch $st_done`;
			}
		}
		
		## execute BAM JUNCTIONS QC
		my $ju_done = "$doneDir"."_junct_qc.done";
		unless (-e $ju_done) {
			my $return_ju = executeBamJunctionQc($dbh,$sampleId,$bamFile,$config{'bed_file'},$bamQcFilePrfx,$bamJcQcFile,$tmpDir,"$config{'rseqc_dir'}/junction_annotation.py");
			if ($return_ju == 0) {
				`touch $ju_done`;
			}
		}
		
		## execute BAM GC CONTENT
		my $gc_done = "$doneDir"."_read_qc.done";
		unless (-e $gc_done) {
			my $return_gc = executeBamReadGc($dbh,$samid,$bamFile,$sampleId, $bamQcFilePrfx, $bamGcQcFile,$tmpDir,"$config{'rseqc_dir'}/read_GC.py");
			if ($return_gc == 0) {
				`touch $gc_done`;
			}
		}
		
		## execute Feature Counts
		my $fcts_done = "$doneDir"."_feat_cts.done";
		unless (-e $fcts_done) {
			my $return_fcts = executeFeatureCountsBam($dbh,$samid,$bamFile,$sampleId,$cntFile,$config{'num_threads'},$tmpDir,$config{'fcts_prog'},$config{'gtf_file'},$fq2File,$config{'stranded'});
			if ($return_fcts == 0) {
				`touch $fcts_done`;
			}
		}
		
		## remove non archivable files
		if (-e $bamFile) {print LOG "rm $bamFile\n";`rm $bamFile`;}
		if (-e $bamEncFile) {print LOG "rm $bamEncFile\n";`rm $bamEncFile`;}
	}
}

## download file from S3
sub downloadFileFromS3{
	my($dbh,$sampleId,$key,$encFile,$tmpDir,$aws_cli_prog)=@_;	
	my $timeBefore = time();
	my $cmd = "$aws_cli_prog s3 cp $key $encFile";
	print LOG "$cmd\n";
	my ($bench_ref) = benchmark($cmd,$tmpDir);
	print LOG "S3 download output:$$bench_ref{output}:\n";
	print LOG "S3 download return code:$$bench_ref{return_code}:\n";
	my $downloaded_md5 = md5sum($encFile);
	my $processId = addProcess($dbh,$sampleId, $cmd,$$bench_ref{return_code},$timeBefore,time());
	addFile($dbh,$sampleId,$processId,$encFile,basename($encFile),$downloaded_md5,bytes($encFile),'encrypted fastq');
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};
}

sub downloadFileFromSra{
	my($dbh,$sampleId,$sra,$fqFile,$tmpDir,$sraProg)=@_;	
	my $timeBefore = time();
	(my $sraFile = $fqFile) =~ s/fastq.gz/sra/g;
	my $cmd = "wget $sra -O $sraFile && $sraProg --gzip -O $tmpDir $sraFile && rm $sraFile";
	print LOG "$cmd\n";
	my ($bench_ref) = benchmark($cmd,$tmpDir);
	print LOG "SRA download output:$$bench_ref{output}:\n";
	print LOG "SRA download return code:$$bench_ref{return_code}:\n";
	my $downloaded_md5 = md5sum($fqFile);
	my $processId = addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());
	addFile($dbh,$sampleId,$processId,$fqFile,basename($fqFile),$downloaded_md5,bytes($fqFile),'fastq');
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};
}

## get bam file from encrypted bam file;
sub decryptFile{			
	my($dbh,$sampleId,$encFile,$bamFile,$numRetry,$tmpDir,$openssl_prog,$password)=@_;	
	my $timeBefore= time();
	my $cmd = "$openssl_prog aes-256-cbc -d -pass pass:$password < $encFile > $bamFile";
	print LOG "$cmd\n";
	my ($bench_ref) = benchmark($cmd,$tmpDir);
	if($numRetry==0) {
		print LOG "openssl output:$$bench_ref{output}:\n";
		print LOG "openssl return code:$$bench_ref{return_code}:\n";
	} else {
		print LOG "Retry $numRetry: openssl output:$$bench_ref{output}:\n";
		print LOG "Retry $numRetry: openssl return code:$$bench_ref{return_code}:\n";
	}	
	
	## restart the decryption if it failed with code 256
	if($$bench_ref{return_code}==256 && $numRetry <=3) {
		$numRetry=$numRetry+1;
		decryptFile($dbh,$sampleId,$encFile,$bamFile,$numRetry);
	}
	
	my $processId = addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());	
	addFile($dbh,$sampleId,$processId,$bamFile,basename($bamFile),md5sum($bamFile),bytes($bamFile),'fastq');
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};	
}

## execute fastQC
sub executeFastQC {
	my($dbh,$outDir,$sampleId,$bamFile,$tmpDir,$bamFastQcDir,$fastqc_prog)=@_;		
	my $timeBefore = time();
	my $cmd = "$fastqc_prog -q -o $outDir $bamFile";
	print LOG "$cmd\n";
	my ($bench_ref) = benchmark($cmd,$tmpDir);
	print LOG "FastQC output:$$bench_ref{output}:\n";
	my $processId=addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());	
	# unzip, remove zipped results, HTML duplicate, and tarball results
	my $mvFastq = "unzip $bamFastQcDir.zip -d $bamFastQcDir && tar -zcvf $bamFastQcDir.tar.gz $bamFastQcDir && rm -r $bamFastQcDir.zip $bamFastQcDir.html $bamFastQcDir";
	print LOG "$mvFastq\n";
	`$mvFastq`;
	my $bamFastQcFile = "$bamFastQcDir.tar.gz";
	addFile($dbh,$sampleId,$processId,$bamFastQcFile,basename($bamFastQcFile),md5sum($bamFastQcFile),bytes($bamFastQcFile),'FastQC');	
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};
}

## execute bam_stat.py [bam statistics];
sub executeBamStatQc {
	my($dbh,$sampleId,$bamFile,$bamQcFile,$tmpDir,$bam_stat_prog)=@_;	
	my $timeBefore = time();
	my $cmd= "$bam_stat_prog -i $bamFile > $bamQcFile";
	print LOG "$cmd\n";
	`echo "$cmd" > $bamQcFile.sh`;
	my $command = "sh " . $bamQcFile . ".sh";
	my ($bench_ref) = benchmark($command,$tmpDir);
	print LOG "Mapped Reads Stat QC output:$$bench_ref{output}:\n";
	my $processId=addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());
	addFile($dbh,$sampleId,$processId,$bamQcFile,basename($bamQcFile),md5sum($bamQcFile),bytes($bamQcFile),'rseqc bam statistics');	
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};
}

## execute junction_annotation.py [junction statistics];
sub executeBamJunctionQc {
	my($dbh,$sampleId,$bamFile,$annotFile,$bamQcFilePrfx,$bamJcQcFile,$tmpDir,$jun_anot_prog)=@_;	
	my $timeBefore = time();
	my $cmd = "$jun_anot_prog -i $bamFile -o $bamQcFilePrfx -r $annotFile 2>$bamJcQcFile";
	print LOG "$cmd\n";
	`echo "$cmd" > $bamJcQcFile.sh`;
	my $command = "sh " . $bamJcQcFile . ".sh";
	my ($bench_ref) = benchmark($command,$tmpDir);
	print LOG "Junction QC return code:$$bench_ref{return_code}:\n";	
	my $processId=addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());	
	addFile($dbh,$sampleId,$processId,$bamJcQcFile,basename($bamJcQcFile),md5sum($bamJcQcFile),bytes($bamJcQcFile),'reseqc bam junctions');
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};	
}

## execute read_distribution.py [5', 3' UTR, exons, introns, intergenic];
sub executeBamReadDistribution {
	my($dbh,$sampleId,$bamFile,$annotFile,$bamRdQcFile,$tmpDir,$read_dist_prog)=@_;		
	my $timeBefore = time();
	my $cmd = "$read_dist_prog -i $bamFile -r $annotFile > $bamRdQcFile";
	print LOG "$cmd\n";
	`echo "$cmd" > $bamRdQcFile.sh`;
	my $command = "sh " . $bamRdQcFile . ".sh";
	my ($bench_ref) = benchmark($command,$tmpDir);
	print LOG "Read Distribution output:$$bench_ref{output}:\n";
	my $processId=addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());	
	addFile($dbh,$sampleId,$processId,$bamRdQcFile,basename($bamRdQcFile),md5sum($bamRdQcFile),bytes($bamRdQcFile),'reseqc bam read distribution');	
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};
}

## execute read_GC.py [GC statistics];
sub executeBamReadGc {
	my($dbh,$samid,$bamFile,$sampleId, $bamQcFilePrfx, $bamGcQcFile,$tmpDir,$read_gc_prog)=@_;		
	my $timeBefore = time();
	my $cmd = "$read_gc_prog -i $bamFile -o $bamQcFilePrfx";
	print LOG "$cmd\n";
	my ($bench_ref) = benchmark($cmd,$tmpDir);
	print LOG "read_GC.py output:$$bench_ref{output}:\n";
	print LOG "read_GC.py return $$bench_ref{return_code}:\n";
	my $processId = addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	
	## extract summary information ["sample_id"     "Min."     "1st Qu."     "Median"     "Mean"     "3rd Qu."     "Max."] from R file
	my $rCmd1 = "echo \"out=as.vector(summary(gc));dta = data.frame('$samid',out[1],out[2],out[3],out[4],out[5],out[6]);write.table(dta,file='$bamGcQcFile',sep='\t',row.names=F,col.names=F,quote=F);\" >> $bamQcFilePrfx.GC_plot.r";
	print LOG "$rCmd1\n";		
	`$rCmd1`;		
	my $rCmd2 = "/usr/local/bin/Rscript --quiet --vanilla $bamQcFilePrfx.GC_plot.r";
	print LOG "$rCmd2\n";		
	my $rCmdOut2=qx($rCmd2 2>&1);
	$rCmdOut2 =~ s/^\s+|\s+$//g ; 
	my $rReturnCode2=$?;
	print LOG "R output:$rCmdOut2:\n";
	print LOG "R return code:$rReturnCode2:\n";			
	addFile($dbh,$sampleId,$processId,$bamGcQcFile,basename($bamGcQcFile),md5sum($bamGcQcFile), bytes($bamGcQcFile),'rseqc bam gc%');
	return $$bench_ref{return_code}			
}

## execute subreadfeaturecounts v1.4.6
sub executeFeatureCountsBam{
	my($dbh,$samid,$bamFile,$sampleId,$cntFile,$threads,$tmpDir,$feat_cts_prog,$annotGtfFile,$fastq2,$stranded)=@_;		
	my $timeBefore = time();
	my $feat_cts_addon = '';
	# add on for paired end reads -- count fragments rather than reads
	if ($fastq2 ne '') {
		$feat_cts_addon = '-B -p -C';
	}
	my $cmd ="$feat_cts_prog $feat_cts_addon -T $threads -s $stranded -a $annotGtfFile -o $cntFile $bamFile";
	print LOG "$cmd\n";
	my ($bench_ref) = benchmark($cmd,$tmpDir);
	print LOG "Subread Feature Counts output:$$bench_ref{output}:\n";
	print LOG "Subread Feature Counts return code:$$bench_ref{return_code}:\n";
	my $processId = addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());	
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	addFile($dbh,$sampleId,$processId,$cntFile,basename($cntFile),md5sum($cntFile),bytes($cntFile),'feature counts');
	return $$bench_ref{return_code};		
}

## execute STAR Read mapping
sub star {
	my($dbh,$fastq1,$fastq2,$threads,$starFile,$bam_file,$tmpDir,$sampleId,$star_prog,$genomeDir)=@_;	
	my $timeBefore = time();
	my $cmd = "$star_prog --runThreadN $threads --genomeDir $genomeDir --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --readFilesIn $fastq1 $fastq2 --outFileNamePrefix $starFile --outTmpDir $tmpDir/star/";
	print LOG "$cmd\n";
	my ($bench_ref) = benchmark($cmd,$tmpDir);
	print LOG "star return code:$$bench_ref{return_code}:\n";
	my $processId = addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());
	my $mvBamCall = "mv $starFile" . "Aligned.sortedByCoord.out.bam $bam_file";
	print LOG "$mvBamCall\n";
	`$mvBamCall`; # rename the output
	addFile($dbh,$sampleId,$processId,"$bam_file",basename($bam_file),md5sum("$bam_file"),bytes("$bam_file"),'STAR mapped and unmapped reads bam file');
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};
}
## execute HISAT2 Read mapping
sub hisat {
	my($dbh,$fastq1,$fastq2,$threads,$bam_file,$tmpDir,$sampleId,$hisat_prog,$genomeDir,$spliceFile,$samtools_prog)=@_;	
	my $timeBefore = time();
	my $hisat_reads = "-U $fastq1";
	if($fastq2 ne '') { # for paired reads
		$hisat_reads = "-1 $fastq1 -2 $fastq2"
	}
	(my $sam_file = $bam_file) =~ s/bam$/sam/g;
	my $cmd = "$hisat_prog -p $threads -x $genomeDir/hisat2_genome_index $hisat_reads --known-splicesite-infile $spliceFile -S $sam_file";
	print LOG "$cmd\n";
	my ($bench_ref) = benchmark($cmd,$tmpDir);
	print LOG "Hisat2 return code:$$bench_ref{return_code}:\n";
	# convert sam to bam and remove sam file
	my $mvBamCall = "$samtools_prog view -@ $threads -b -o $bam_file $sam_file && rm $sam_file";
	print LOG "$mvBamCall\n";
	`$mvBamCall`; # convert sam to bam
	my $processId = addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());
	addFile($dbh,$sampleId,$processId,"$bam_file",basename($bam_file),md5sum("$bam_file"),bytes("$bam_file"),'HISAT2 mapped and unmapped reads BAM file');
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};
}
					
sub starIndex {
	my($dbh,$threads,$tmpDir,$sampleId,$star_prog,$genomeDir,$geneomeFile,$annotGtfFile)=@_;	
	my $timeBefore = time();
	my $cmd = "$star_prog --runThreadN $threads --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $geneomeFile --sjdbGTFfile $annotGtfFile";
	print LOG "$cmd\n";
	my ($bench_ref) = benchmark($cmd,$tmpDir);
	print LOG "star index return code:$$bench_ref{return_code}:\n";
	my $processId = addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};
}

sub hisatIndex {
	my($dbh,$threads,$tmpDir,$sampleId,$hisat_prog,$genomeDir,$genomeFile,$annotGtfFile,$spliceFile)=@_;	
	my $timeBefore = time();
	# run both indexing and known splice junctions
	my $cmd = "$hisat_prog-build -p $threads $genomeFile $genomeDir/hisat2_genome_index && $hisat_prog" . "_extract_splice_sites.py $annotGtfFile > $spliceFile";
	print LOG "$cmd\n";
	my ($bench_ref) = benchmark($cmd,$tmpDir);
	print LOG "hisat index return code:$$bench_ref{return_code}:\n";
	my $processId = addProcess($dbh,$sampleId,$cmd,$$bench_ref{return_code},$timeBefore,time());
	addBenchmark($dbh,$processId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};
}

## Compress bam
sub compress_bam {
	my($dbh,$genome,$threads,$tmpDir,$sampleId,$cramFile,$samtools_prog,$bamFile)=@_;	
	my $timeBefore = time();
	my $Cmd = "$samtools_prog view -@ $threads --output-fmt cram --reference $genome -o $cramFile $bamFile";
	print LOG "$Cmd\n";
	my ($bench_ref) = benchmark($Cmd,$tmpDir);
	print LOG "Samtools CRAM return code:$$bench_ref{return_code}:\n";
	my $ProcessId = addProcess($dbh,$sampleId,$Cmd,$$bench_ref{return_code},$timeBefore,time());
	addFile($dbh,$sampleId,$ProcessId,"$cramFile",basename($cramFile),md5sum("$cramFile"),bytes("$cramFile"),'mapped and unmapped CRAM file');
	addBenchmark($dbh,$ProcessId,$sampleId,$bench_ref);
	return $$bench_ref{return_code};
}

## get md5sum per file
sub md5sum {
	my $samid = shift;
	my $md5Cmd= "md5sum $samid";
	my $md5sum=qx($md5Cmd 2>&1);
	($md5sum, undef) = split('\s',$md5sum);
	$md5sum =~ s/^\s+|\s+$//g ; 
	return $md5sum;
}

## get bytes per file
sub bytes {
	my $samid = shift;
	my $bytesCmd= "du -b $samid";
	my $bytes=qx($bytesCmd 2>&1);
	($bytes, undef) = split('\s',$bytes);
	$bytes =~ s/^\s+|\s+$//g ; 
	return $bytes;
}

###############################################################
##
## SQLite database interactions
##
###############################################################

## initializat the database
sub initDb {
	$dbh->do("create table if not exists process(process_id integer primary key autoincrement, process_str text, return_code, start_time integer,
	end_time integer,wc_time integer, sample_id integer,created timestamp default (datetime('now','localtime')))");
	
	$dbh->do("create table if not exists file(file_id integer primary key autoincrement, file_name text, file_md5 text, file_bytes integer, file_type text, 
	process_id integer, sample_id integer, file_path text, created timestamp default (datetime('now','localtime')))");
	
	$dbh->do("create table if not exists sample(sample_id integer primary key autoincrement, sample_name text unique,
	created timestamp default (datetime('now','localtime')))");
	
	$dbh->do("create table if not exists benchmark(benchmark_id integer primary key autoincrement, process_id integer, sample_id integer, virtual integer, 
	resident_set_size integer, cpu_time integer, wc_time integer, return_code integer, created timestamp default (datetime('now','localtime')))")
}

## add a process to the metadata database
sub addProcess {
	my ($dbhr, $sampleId, $cmd,$returnCode,$start,$end)=@_;
	$dbhr->do("insert into process (sample_id,process_str,return_code,start_time,end_time, wc_time) 
	values($sampleId,'$cmd',$returnCode,$start,$end,$end-$start)");
	
	my $processId=undef;

 	my $sth = $dbhr->prepare('select max(process_id) from process');
	 	$sth->execute();
	 	while (my @row = $sth->fetchrow_array) {
	 	$processId = $row[0];
 	}
	return $processId;
}

## add a sample to the metadata database
sub addSample {
	my ($dbhr, $sampleName)=@_;
	
	$dbhr->do("insert or ignore into sample (sample_name) 
	values('$sampleName')");
	
	my $sampleId=undef;

 	my $sth = $dbhr->prepare("select sample_id from sample where sample_name='$sampleName'");
	 	$sth->execute();
	 	while (my @row = $sth->fetchrow_array) {
	 	$sampleId = $row[0];
 	}
	return $sampleId;
}

## add new files and file metadata
sub addFile {
	my ($dbhr, $sampleId,$processId,$filePath,$samid,$fileMd5,$fileBytes,$fileType)=@_;
	$dbhr->do("insert into file (sample_id,process_id,file_path,file_name,file_md5,file_bytes,file_type) 
	values($sampleId,$processId,'$filePath','$samid','$fileMd5','$fileBytes','$fileType')");
}
## add benchmark data to SQLite database
sub addBenchmark {
	my ($dbhr,$processId,$sampleId,$info)=@_;
	$dbhr->do("insert into benchmark (process_id,sample_id,virtual,resident_set_size,cpu_time,wc_time,return_code) 
	values($processId,$sampleId,$$info{max_virtual},$$info{max_rss},$$info{cpu_time},$$info{wc_time},$$info{return_code})");
}

###############################################################
##
## Benchmarking
##
###############################################################

## Main process for benchmarking. Give this subroutine a linux command string and it runs the command, and returns:
## max Virtual memory, max RSS memory, CUP_time, wall clock time, retrun code, and any output from the command pointed at STDOUT/STDERR
sub benchmark{
	my ($command, $tmpDir) = @_;
	# define variables
	my $pid = fork;
	my $isfirstloop = 1;
	my ($cpu_start_str, $timestart) = "";
	my $max_rss_mem = 0;
	my $max_virt_mem = 0;
	
	# if a child process exists, fork
	if ( !$pid ) {
		exec("$command 2> $tmpDir/bench_out.txt");
	} else {
		
		# While the process exists, monitor its CPU and memory
		while ( !waitpid( $pid, WNOHANG )) {
			
			if ($isfirstloop == 1) {
				$cpu_start_str = `cat /proc/$$/stat`;
				$timestart = time();
				$isfirstloop = 0;
			}
			
			# Grab sum of all running threads under $pid and (sum the RSS, then sum the Virtual) every second
			my $rss_mem = `pstree -p -A $pid | grep -o -P "[0-9]{1,5}" | xargs ps -o rss |  awk '{ sum+=\$1} END {print sum}'`;
			my $virt_mem = `pstree -p -A $pid | grep -o -P "[0-9]{1,5}" | xargs ps -o vsz |  awk '{ sum+=\$1} END {print sum}'`;
			if ($rss_mem > $max_rss_mem) {$max_rss_mem = $rss_mem;}
			if ($virt_mem > $max_virt_mem) {$max_virt_mem = $virt_mem;}
			sleep 1;
		}
	}
	my $return_code=$?;
    my $cpu_end_str = `cat /proc/$$/stat`;
    my $timeend = time();
    # grab final stats for CPU time, User time, and wall clock time
    # each read is in ticks: 100 ticks = 1 second in unix
    # [14] = perl program cpu time / Kernel time, [16] = child process CPU time.
    # ([14] + [16]) / 100 = total CPU time for all perl program and all children/threads
    my @cpustart = split(" ", $cpu_start_str);
    my @cpuend = split(" ", $cpu_end_str);
    my $cputime = (($cpuend[14]-$cpustart[14]) + ($cpuend[16]-$cpustart[16])) / 100;
    my $wall_clock_time = ($timeend - $timestart);
    my $out = `cat $tmpDir/bench_out.txt`;
    chomp $out;
    `rm -r $tmpDir/bench_out.txt`;
	# return a referenced hash
    my %rethash = (
    "max_rss"     => "$max_rss_mem",
    "max_virtual" => "$max_virt_mem",
    "cpu_time"    => "$cputime",
    "wc_time"     => "$wall_clock_time",
	"return_code" => "$return_code",
	"output"      => "$out"
	);
	return(\%rethash);
}
