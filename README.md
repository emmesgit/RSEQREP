# RSEQREP

RNA-Seq Reports, an open-source cloud-enabled framework for reproducible RNA-Seq data processing, analysis, and result reporting


INSTALLATION
 
Option 1 (local Ubuntu server): 
a)       install LTS Ubuntu desktop version 16.04.03 on your machine (http://releases.ubuntu.com/16.04.3/, install from bootable USB)
b)      copy/clone the RSEQREP github source code to your Ubuntu machine (git clone https://github.com/emmesgit/RSEQREP.git)
c)       execute our installation shell script to install the software on your own Ubuntu machine (sh RSEQREP/ubuntu/install-software.sh)
Option 2 (RSEQREP AWS AMI):
a)       initialize RSEQPRE AMI ((https://aws.amazon.com, AMI ID: RSEQREP (RNA-Seq Reports) v1.0 (ami-e1708b9b))).  To do this, log into the AWS console, navigate to the EC2 resources.  Next select AMIs in the navigation pane and select public images.  Finally search for RSEQREP and launch the AMI.
b)      copy/clone the RSEQREP github source code to your Ubuntu machine (git clone https://github.com/emmesgit/RSEQREP.git)
EXECUTION
a)       fill out the configuration file (RSEQREP/config/config.xlsx).  RSEQREP/case-study/config-henn.xlsx is an example of a complete configuration file. 
b)      execute start-to-end analysis (sh RSEQREP/run-all.sh) or for a particular component (sh RSEQREP/run-pre-processing.sh; sh RSEQREP/run-analysis.sh; sh RSEQREP/run-report.sh).
 
TROUBLESHOOTING
Upon inspecting the initial run of the report, you may find that the configuration option that you initially chose does not fit your data.  For example, you inspect the reverse cumulative distribution function plot comparing log count per million cutoffs with the number of retained genes.  You find that the cutoff you had originally selected resulted in too few genes for the analysis.  To remedy this, you would update the configuration file to reflect a more appropriate log counts per million cutoff.  You then determine the steps of the analysis that will be affected by this configuration change.  In this case, the analysis and report steps are affected.  Re-execute the analysis and report steps (sh RSEQREP/run-analysis.sh; sh RSEQREP/run-report.sh).  The configuration file is re-parsed each time a RSEQREP/run-* script is executed.
 
BUGS/ISSUES
If you identify bugs/issues, please submit them via the GitHub Issue Tracker 
https://github.com/emmesgit/RSEQREP/issues
 
ADDITIONAL RESOURCES
For a detailed explanation of the software and its capabilities please navigate to our F1000 research publication:
https://f1000research.com/articles/6-2162/v1
 

