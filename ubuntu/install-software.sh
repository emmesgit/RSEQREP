#!/bin/bash
##############################################################################################################################
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
# Program:  install-software.pl 
# Version:  RSEQREP 1.0.0
# Author:   Kevin Conway and Leigh Villarroel
# Purpose:  install software onto 16.0.4 Ubuntu desktop computer required for RSEQREP usage.
# Input:    N/A
# Output:   N/A
##############################################################################################################################

##Add Swap for extra memory if needed 
sudo /bin/dd if=/dev/zero of=/var/swap.1 bs=1M count=1024
sudo /sbin/mkswap /var/swap.1
sudo chmod 600 /var/swap.1
sudo /sbin/swapon /var/swap.1


##Initial setup
## emmes_mnt.sh should be run first 
##create /usr/local/bin/emmes dir if it doesn't already exist
mydir=/usr/local/emmes
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	sudo mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi

##create/cd into emmes installation dir
mydir=/home/ubuntu/emmes_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

##JAVA
myjava=/usr/bin/java
if [ ! -f "$myjava" ]
then
	sudo apt -y install default-jre
fi
#sudo apt-get -y install openjdk-6-jdk 
echo java is INSTALLED

##GCC and other Configure/ Make Components
mygcc=/usr/bin/gcc
if [ ! -f "$mygcc" ] 
then
	sudo apt-get -y install build-essential
fi	
sudo apt-get -y install libncurses5-dev #Needed for SAMtools install
sudo apt-get -y install zlib1g-dev #Needed for SAMtools install
#sudo apt-get -y install libcurl4-nss-dev #Needed for SAMtools install either
#sudo apt-get -y install libnss3-dev #Needed for SAMtools install
sudo apt-get -y install libcurl4-openssl-dev #Needed for SAMtools install
sudo apt-get -y install libssl-dev #Needed for SAMtools install
sudo apt-get -y install autoconf #Needed for htop install
sudo apt-get -y install libncursesw5-dev #Needed for htop install
sudo apt-get -y install gfortran #Needed for R-Base install
sudo apt-get -y install aptitude
sudo aptitude -y install libreadline-dev
sudo apt-get -y install libxml2-dev #Needed for BiocLite('biomaRt') to run successfully
sudo apt-get -y install libcairo2-dev #Needed for Cairo package in R
sudo apt-get -y install libxt-dev #Needed for Cairo package in R
sudo apt-get -y install libgmp3-dev #Needed for iterpc package in R
sudo apt-get -y install liblzma-doc liblzma-dev xz-utils #for samtools/htslib-1.5
#sudo apt-get -y install libbz2-dev
sudo apt-get -y install libtbb2 #for bowtie2
sudo apt-get -y install openjdk-8-jdk #for R 3.4.1
sudo apt-get -y install mesa-common-dev libglu1-mesa-dev  #R package rgl
sudo apt-get -y install evince #PDF viewer



	
echo gcc and other necessary compiler tools INSTALLED

##PYTHON
mypython27=/usr/bin/python2.7
mypython=/usr/bin/python
if [ ! -f "$mypython" ] && [ ! -f "$mypython27" ]
then
	sudo apt-get -y install python2.7 python-dev libpython2.7-dev 
	sudo cp /usr/bin/python2.7 /usr/bin/python
elif [ ! -f "$mypython" ] && [ -f "$mypython27" ]
then
	sudo apt-get -y install python-dev libpython2.7-dev 
	sudo cp /usr/bin/python2.7 /usr/bin/python
fi
echo "Installing required python modules."
sudo apt-get -y install python-dev libpython2.7-dev #Gives the Python.h header files for rseqc install
sudo apt-get -y install libbz2-1.0 libbz2-dev libbz2-ocaml libbz2-ocaml-dev #for samtools/htslib-1.5

#rm -fr /usr/lib/python2.7/dist-packages/distribute*
echo python INSTALLED

##PIP
mypip=/usr/bin/pip
if [ ! -f "$mypip" ]
then
	sudo apt-get -y install python-pip
	sudo pip install --upgrade pip
fi
echo pip INSTALLED

##UNZIP
myfile=/usr/bin/unzip
if [ ! -f "$myfile" ]
then
	echo installing unzip
	sudo apt-get -y install unzip
    if [ ! -f "$myfile" ]
    then
		echo could not install unzip
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"
echo unzip INSTALLED

##GIT
mygit=/usr/bin/git
if [ ! -f "$mygit" ]
then
	sudo apt-get -y install git
fi
echo git INSTALLED

###USER UTILITIES
##Firefox
myfirefox=/usr/bin/firefox
if [ ! -f "$myfirefox" ]
then
	sudo apt-get -y install firefox
fi
echo firefox INSTALLED

##Libre office
mylibre=/usr/bin/libreoffice
if [ ! -f "$mylibre" ]
then
	sudo apt-get -y install libreoffice
fi
echo libreoffice INSTALLED

##X2GO?


###UTILITIES
##Sqlite
#wget https://sqlite.org/2017/sqlite-tools-linux-x86-3200100.zip
# create/cd into sqlite installation dir
mydir=/home/ubuntu/emmes_install/sqlite_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# download sqlite binary tarball
myfile=sqlite-autoconf-3200100.tar.gz #sqlite-tools-linux-x86-3200100.zip 
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget https://sqlite.org/2017/$myfile
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# extract tophat binary tarball
mydir=/home/ubuntu/emmes_install/sqlite_install/sqlite-autoconf-3200100
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	tar xvfz $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot extract "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
 fi
#configure
cd $mydir 

myfile=/usr/local/bin/sqlite3
if [ ! -f "$myfile" ]
then
	echo building sqlite
	./configure 
	make
	sudo make install
	if [ ! -f "$myfile" ]
	then
		echo could not build sqlite
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo sqlite3 INSTALLED

##SRATOOLKIT
#wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-ubuntu64.tar.gz
# create/cd into sqlite installation dir
mydir=/home/ubuntu/emmes_install/sratoolkit_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# download sratoolkit binary tarball
myfile=sratoolkit.2.8.2-1-ubuntu64.tar.gz 
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-ubuntu64.tar.gz
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# extract tophat binary tarball
mydir=/home/ubuntu/emmes_install/sratoolkit_install/sratoolkit.2.8.2-1-ubuntu64
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	tar xvfz $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot extract "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
 fi

if [ ! -f "/usr/bin/sra-sort.2.8.2" ]
then
	cd $mydir/bin/
	echo copying sratoolkit files to /usr/bin
	sudo cp -a * /usr/bin/
	echo sratoolkit is INSTALLED
fi

##TOPHAT
# create/cd into tophat installation dir
mydir=/home/ubuntu/emmes_install/tophat_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# download tophat binary tarball
myfile=tophat-2.1.1.Linux_x86_64.tar.gz
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget http://ccb.jhu.edu/software/tophat/downloads/$myfile
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# extract tophat binary tarball
mydir=/home/ubuntu/emmes_install/tophat_install/tophat-2.1.1.Linux_x86_64
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	tar xvzf $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot extract "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
 fi

# copy new tophat executables from tophat installation dir to /usr/bin
# note this is a little ugly because on very first pass /usr/bin/tophat does not exist
if [ -f "/usr/bin/tophat" ]
then
	/usr/bin/tophat --version | grep 2\.1\.1 || (
		# archive original tophat executables in /usr/bin
		echo tophat version is not 2.1.1, presumed old or not available
		echo ls-ing archived files in /usr/bin starting with tophat
		ls -alh /usr/bin/tophat*
		myfile=/usr/bin/orig_tophat
		if [ ! -f "$myfile" ]
		then
			echo archiving original tophat executables in /usr/bin
			for myitem in bam2fastx bam_merge bed_to_juncs contig_to_chr_coords fix_map_ordering gtf_juncs gtf_to_fasta juncs_db long_spanning_reads map2gtf prep_reads sam_juncs samtools_0.1.18 segment_juncs sra_to_solid tophat tophat2 tophat-fusion-post tophat_reports
			do
				sudo mv /usr/bin/$myitem /usr/bin/orig_$myitem
			done
			if [ ! -f "$myfile" ]
			then
				echo could not archive original tophat executables
				echo exiting with error code 1 ...
				exit 1
			fi
		fi
	)

	/usr/bin/tophat --version | grep 2\.1\.1 && (
	   echo tophat version is INSTALLED and OK at 2.1.1
	)
else 
	##tophat not there, move the files to /usr/bin
	myfile=/home/ubuntu/emmes_install/tophat_install/tophat-2.1.1.Linux_x86_64/tophat
	if [ -f "$myfile" ]
	then
		echo copying new tophat executables into /usr/bin
		for myotheritem in bam2fastx bam_merge bed_to_juncs contig_to_chr_coords fix_map_ordering gtf_juncs gtf_to_fasta juncs_db long_spanning_reads map2gtf prep_reads sam_juncs samtools_0.1.18 segment_juncs sra_to_solid tophat tophat2 tophat-fusion-post tophat_reports
		do 
			sudo cp -a /home/ubuntu/emmes_install/tophat_install/tophat-2.1.1.Linux_x86_64/$myotheritem /usr/bin/
		done
		if [ ! -f /usr/bin/tophat ]
		then
			echo could not copy new tophat executables to /usr/bin
			echo exiting with error code 1 ...
			exit 1
		fi
	fi
fi ##end tophat exists check
echo tophat version is INSTALLED

##STAR
# create/cd into STAR installation dir
mydir=/home/ubuntu/emmes_install/star_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# get STAR via git
# get the latest version (2.5.2a) 
# wget https://github.com/alexdobin/STAR/archive/2.5.2a.tar.gz
# download STAR binary tarball
myfile=2.5.2a.tar.gz
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget wget https://github.com/alexdobin/STAR/archive/$myfile
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# extract STAR binary tarball
mydir=/home/ubuntu/emmes_install/star_install/STAR-2.5.2a
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	tar xvzf $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot extract "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi

#jump to bin file area
cd "$mydir/bin/Linux_x86_64_static"

# install STAR binary
mybinfile=/usr/bin/STAR
if [ ! -f $mybinfile ]
then
	sudo cp -a STA* /usr/bin #copies all STAR files in the above location to bin
	if [ ! -f "$mybinfile" ]
	then
		echo could not install "$mybinfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
#Verify
$mybinfile --version | grep 2\.5\.2a && (
   echo STAR version in /usr/bin is OK and INSTALLED at 2.5.2a
)

$mybinfile --version | grep 2\.5\.2a || (
   echo STAR version in /usr/bin is not 2.5.2a
   echo installing executable obtained via git
   sudo cp -a $myfile $mybinfile
   $mybinfile --version | grep 2\.5\.2a || (
         echo could not install STAR version 2.5.2a
         echo exiting with error code 1 ...
         exit 1
   )
)
#sudo chmod ugo+x /usr/bin/STAR
echo STAR version is INSTALLED


##SAMtools
# create/cd into SAMtools installation dir
mydir=/home/ubuntu/emmes_install/samtools_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# download SAMtools source tarball (1.5)
myfile=samtools-1.5.tar.bz2
if [ ! -f "$myfile" ]
   then
      echo wget-ing file "$myfile"
      wget https://github.com/samtools/samtools/releases/download/1.5/$myfile
      if [ ! -f "$myfile" ]
         then
            echo could not wget file "$mydir"
            echo exiting with error code 1 ...
            exit 1
      fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# extract SAMtools source tarball
mydir=/home/ubuntu/emmes_install/samtools_install/samtools-1.5
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	tar xvjf $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot extract "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi

# build SAMtools in "in-place mode"
cd $mydir
myfile=samtools
if [ ! -f "$myfile" ]
   then
      echo building SAMtools "in-place"
      ./configure --enable-plugins --enable-libcurl --with-plugin-path=$PWD/htslib-1.5
      sudo make all plugins-htslib
      if [ ! -f "$myfile" ]
         then
            echo could not build SAMtools
            echo exiting with error code 1 ...
            exit 1
      fi
	##run the install using sudo apt install samtools
	#echo running apt to install SAMtools
	
	##move the samtools executable to /usr/bin
	
	myfile=/usr/bin/samtools
	myfile2=/usr/bin/orig_samtools
	
	if [ ! -f "$myfile" ] && [ ! -f "$myfile2" ]
	then
		#copy over samtools executable to usr/bin
		sudo cp -a samtools /usr/bin
		echo samtools version is at 1.5 and INSTALLED
		
	elif [ -f "$myfile" ] && [ ! -f "$myfile2" ] 
	then
		##make copies and copy over executables from current install directory
		echo archiving original samtools executables in /usr/bin	
		sudo mv -v "/usr/bin/samtools" "/usr/bin/orig_samtools"
		#move the new samtools version	
		echo copying samtools executables to /usr/bin
		sudo cp -a samtools /usr/bin
		echo samtools version is at 1.5 and INSTALLED
		
	elif [ ! -f "$myfile" ] && [ -f "$myfile2" ] 
	then
		sudo cp -a samtools /usr/bin
		echo samtools version is at 1.5 and INSTALLED
		
	else
		## Both files are there, remove the original, move the old to the orignal then copy over the executable
		sudo rm -f "/usr/bin/orig_samtools"
		##make copies and copy over executables from current install directory
		echo archiving original samtools executables in /usr/bin	
		sudo mv -v "/usr/bin/samtools" "/usr/bin/orig_samtools"
		#move the new samtools version	
		echo copying samtools executables to /usr/bin
		sudo cp -a samtools /usr/bin
		echo samtools version is at 1.5 and INSTALLED
	fi

fi
echo samtools is INSTALLED

##BCFtools
# create/cd into BCFtools installation dir
mydir=/home/ubuntu/emmes_install/bcftools_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# download BCFtools source tarball
myfile=bcftools-1.5.tar.bz2
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget https://github.com/samtools/bcftools/releases/download/1.5/$myfile
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# extract BCFtools source tarball
mydir=/home/ubuntu/emmes_install/bcftools_install/bcftools-1.5
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	tar xvjf $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot extract "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi

# build BCFtools in "in-place mode"
cd $mydir
myfile=bcftools
if [ ! -f "$myfile" ]
then
	echo building BCFtools "in-place"
	./configure --enable-plugins --enable-libcurl --with-plugin-path=$mydir/htslib-1.5 ##changed plugin-path to point to the install area
	make all plugins-htslib
	if [ ! -f "$myfile" ]
	then
		echo could not BCF
		echo exiting with error code 1 ...
		exit 1
	fi

	##move the bcftools executable to /usr/bin ##EDIT
	
	myfile=/usr/bin/bcftools
	myfile2=/usr/bin/orig_bcftools
	
	if [ ! -f "$myfile" ] && [ ! -f "$myfile2" ]
	then
		#copy over bcftools executable to usr/bin
		sudo cp -a bcftools /usr/bin
		echo bcftools version is at 1.5 and INSTALLED
		
	elif [ -f "$myfile" ] && [ ! -f "$myfile2" ] 
	then
		##make copies and copy over executables from current install directory
		echo archiving original bcftools executables in /usr/bin	
		sudo mv -v "/usr/bin/bcftools" "/usr/bin/orig_bcftools"
		#move the new bcftools version	
		echo copying bcftools executables to /usr/bin
		sudo cp -a bcftools /usr/bin
		echo bcftools version is at 1.5 and INSTALLED
		
	elif [ ! -f "$myfile" ] && [ -f "$myfile2" ] 
	then
		sudo cp -a bcftools /usr/bin
		echo bcftools version is at 1.5 and INSTALLED
		
	else
		## Both files are there, remove the original, move the old to the orignal then copy over the executable
		sudo rm -f "/usr/bin/orig_bcftools"
		##make copies and copy over executables from current install directory
		echo archiving original bcftools executables in /usr/bin	
		sudo mv -v "/usr/bin/bcftools" "/usr/bin/orig_bcftools"
		#move the new bcftools version	
		echo copying bcftools executables to /usr/bin
		sudo cp -a bcftools /usr/bin
		echo bcftools version is at 1.5 and INSTALLED
	fi
	
fi
echo BCFtools is INSTALLED

##SUBREAD
# create/cd into subread installation dir
mydir=/home/ubuntu/emmes_install/subread_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# download subread source tarball
myfile=subread-1.5.3-Linux-x86_64.tar.gz
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget https://downloads.sourceforge.net/project/subread/subread-1.5.3/$myfile
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# extract subread tarball
mydir=/home/ubuntu/emmes_install/subread_install/subread-1.5.3-Linux-x86_64/
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	tar xvzf $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot extract "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi

# copy new subread executables from subread installation dir to /usr/bin
# note this is a little ugly because on the very first pass /usr/bin/subread-align does not exist
# note also the fancy footwork with redirects. Apparently subread-align -v
# returns a non-zero error code which is filtered out in this way
# (ultimately by the grep, I guess)

/usr/bin/subread-align -v 2>&1 | grep 1\.5\.3 && (
   echo subread-align version is INSTALLED OK at 1.5.3
)

/usr/bin/subread-align -v 2>&1 | grep 1\.5\.3 || (
	echo subread-align version is not 1.5.3, presumed old
	echo ls-ing archived files in /usr/bin starting with subread
	ls -alh /usr/bin/subread*
	myfile=/home/ubuntu/emmes_install/subread_install/subread-1.5.3-Linux-x86_64/bin/subread-align
	if [ -f "$myfile" ]
	then
		echo copying new subread executables into /usr/bin
		for myotheritem in exactSNP subindel subread-align featureCounts subjunc subread-buildindex
		do 
			#sudo mv /usr/bin/$myotheritem 
			sudo cp -a /home/ubuntu/emmes_install/subread_install/subread-1.5.3-Linux-x86_64/bin/$myotheritem /usr/bin/
		done
		if [ ! -f /usr/bin/subread-align ]
		then
			echo could not copy new subread executables to /usr/bin
			echo exiting with error code 1 ...
			exit 1
		fi
	fi
	echo ls-ing files in /usr/bin/ starting with subread again
	ls -alh /usr/bin/subread*
)
echo Subread-align is INSTALLED

##RSEQC
# create/cd into rseqc installation dir
mydir=/home/ubuntu/emmes_install/rseqc_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
		if [ ! -d "$mydir" ]
		then
			echo cannot create directory "$mydir"
			echo exiting with error code 1 ...
			exit 1
		fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# download rseqc source tarball
# http://sourceforge.net/projects/rseqc/files/RSeQC-2.6.4.tar.gz/download
# hacked download location to a mirror, otherwise hangs or fails
myfile=RSeQC-2.6.4.tar.gz
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	#     used direct link because user friendly link failed
	wget https://sourceforge.net/projects/rseqc/files/RSeQC-2.6.4.tar.gz/download?use_mirror=pilotfiber#
	#     make the dir name less cumbersome and more as expected
	echo "if hanging will not get to here"
	pwd
	mv download?use_mirror=pilotfiber RSeQC-2.6.4.tar.gz
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# extract rseqc source tarball and cd into installation dir
mydir=/home/ubuntu/emmes_install/rseqc_install/RSeQC-2.6.4
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	tar xvzf $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot extract "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# build and install rseqc according to docs/INSTALL file
# expecting to produce file /usr/local/bin/read_GC.py
myfile=/usr/local/bin/read_GC.py
if [ ! -f "$myfile" ]
then
	echo running default setup.py
	#mv setup.py setup.py.orig
	#cp -a ../../../emmes_rseqc.setup.py setup.py 
	sudo python setup.py install 
	if [ ! -f "$myfile" ]
	then
		echo could not generate file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"
echo "save copy of setup.py"
cp -a setup.py ../../../setup-RSeQC-2.6.4.py
echo RSeQC is INSTALLED

## may need to remove old distribute package..
#sudo rm -fr /usr/local/lib/python2.7/dist-packages/distribute*



# install FastQC
# create/cd into fastqc installation dir
mydir=/home/ubuntu/emmes_install/fastqc_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# download fastqc binary tarball
myfile=fastqc_v0.11.5.zip
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/$myfile
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# unzip fastqc binary tarball
# note no version info in directory name, so might quit too soon
mydir=/home/ubuntu/emmes_install/fastqc_install/FastQC
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	unzip $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot unzip "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi

# copy new fastqc directory from fastqc installation dir to /usr/local/emmes
mydir=/usr/local/emmes/FastQC
if [ ! -d "$mydir" ]
then
	echo copying installation FastQC dir and contents to "$mydir"
	sudo cp -rp  /home/ubuntu/emmes_install/fastqc_install/FastQC /usr/local/emmes/
	sudo chmod ugo+x /usr/local/emmes/FastQC/fastqc
	sudo chown -R root.root /usr/local/emmes/FastQC
	sudo mv /usr/local/bin/fastqc /usr/local/bin/fastqc_orig
	sudo ln -s /usr/local/emmes/FastQC/fastqc /usr/local/bin/fastqc
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
else
	echo $mydir already exists, nothing to do
fi
 
##if we got to here, fastqc is either absent or old
#myfile=/usr/local/emmes/FastQC/fastqc
#if [ -f "$myfile" ]
#then
#      echo fastqc is already installed, checking version
#      $myfile --version | grep 0\.11\.5 && (
#         echo fastqc version is 0.11.5, nothing to do
#         #exit 1
#      )
#      $myfile --version | grep 0\.11\.5 || (
#         echo fastqc is installed but version is not 0.11.5, presumed old
#         #exit 1
#      )
#fi
echo fastqc version is 0.11.5 and INSTALLED
 
 
 
##LATEX FONTS
## Check for /usr/share/fonts/type1/texlive-fonts-recommended/ and /usr/share/doc/texlive-fonts-extra/
myfonts1=/usr/share/fonts/type1/texlive-fonts-recommended/
myfonts2=/usr/share/doc/texlive-fonts-extra/
if [ ! -d "$myfonts1" ] || [ ! -d "$myfonts2" ]
then
	sudo apt-get -y install texlive-fonts-recommended && sudo apt-get -y install texlive-fonts-extra
fi
echo Latex fonts are INSTALLED

##Seqtk
# create/cd into seqtk installation dir
mydir=/home/ubuntu/emmes_install/seqtk_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# get seqtk via git
# git clone https://github.com/lh3/seqtk.git

myfile=seqtk/seqtk.c
if [ ! -f "$myfile" ]
then
	echo gitting seqtk from github repository
	git clone https://github.com/lh3/seqtk.git
	if [ ! -f "$myfile" ]
	then
		echo could not git file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

cat $myfile | grep 1\.2-r94 && (
   echo seqtk version obtained via git is OK at 1.2-r94
)

cat $myfile | grep 1\.2-r94 || (
   echo seqtk version obtained via git is not 1.2-r94
   echo exiting with error code 1 ...
   exit 1
)

# build seqtk binary as needed
mydir=/home/ubuntu/emmes_install/seqtk_install/seqtk
cd $mydir
myfile=seqtk
if [ ! -f "$myfile" ]
   then
      echo building seqtk
      make
      if [ ! -f "$myfile" ]
         then
            echo could not build Seqtk
            echo exiting with error code 1 ...
            exit 1
      fi
fi

# install seqtk binary
mybinfile=/usr/bin/seqtk

if [ ! -f $mybinfile ]
   then
   sudo cp -a $myfile $mybinfile
   if [ ! -f "$mybinfile" ]
      then
         echo could not install "$mybinfile"
         echo exiting with error code 1 ...
         exit 1
   fi
fi

$mybinfile 2>&1 | grep 1\.2-r94 && (
   echo seqtk version in /usr/bin is OK at 1.2-r94
)

$mybinfile 2>&1 | grep 1\.2-r94 || (
   echo seqtk version in /usr/bin is not 1.2-r94
   echo installing executable obtained via git
   sudo cp -a $myfile $mybinfile
   $mybinfile  2>&1 | grep 1\.2-r94 || (
         echo could not install seqtk version 1.2-r94
         echo exiting with error code 1 ...
         exit 1
   )
)
echo Seqtk is INSTALLED

# install amazon command line interface software (AWSCLI)
myfile=/usr/bin/aws
if [ -f "$myfile" ]
then
	echo awscli is already installed
else
	#try apt then try pip
	#sudo apt -y install awscli
	
	/usr/bin/python --version | grep 2\.7 && (
	 echo python is installed but major version is not 2.7, presumed old
	 echo resolve before continuing with installation of aws cli
	 exit 1
	)
	/usr/bin/python --version | grep 2\.7 || (
	 echo python is installed and is at major version 2.7
	 /usr/bin/pip --help | grep "Search PyPI" || (
		echo pip is not installed
		echo resolve before continuing with installation of aws cli
		exit 1
	 )
	 /usr/bin/pip --help | grep "Search PyPI" && (
		echo pip is installed, proceed with installation of aws cli
		sudo pip install awscli==1.10.48
	# make aws help generally available
		sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg/EGG-INFO/top_level.txt
	# resolve permissions issue encountered by Travis
		sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg/EGG-INFO/requires.txt
		if [ ! -f "$myfile" ]
		   then
			  echo could not install awscli
			  exit 1
		fi
	 )
	)

fi
echo awscli is INSTALLED

##install fastx_toolkit

# create/cd into fastx_toolkit installation dir
mydir=/home/ubuntu/emmes_install/fastx_toolkit_install
if [ ! -d "$mydir" ]
   then
      echo creating directory "$mydir"
      mkdir $mydir
      if [ ! -d "$mydir" ]
         then
            echo cannot create directory "$mydir"
            echo exiting with error code 1 ...
            exit 1
      fi
 fi
 echo cd-ing into directory "$mydir"
 cd "$mydir"

##Firstly need to get gtextutils from http://hannonlab.cshl.edu/fastx_toolkit/libgtextutils-0.6.1.tar.bz2 before the fastx is installed
myfile=libgtextutils-0.6.1.tar.bz2
if [ ! -f "$myfile" ] #may not work multiple times, need to find a better way to check for this file installed 
then
	echo wget-ing file "$myfile"
	wget http://hannonlab.cshl.edu/fastx_toolkit/libgtextutils-0.6.1.tar.bz2
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
	#quick approach to extract 
	tar -xvjf $myfile
	cd libgtextutils-0.6.1
	./configure
	make
	sudo make install
	cd ..
fi

##Now fastx_toolkit
# getting via git https://github.com/agordon/fastx_toolkit.git  doesnt bring the configure file so getting it via wget
myfile=fastx_toolkit-0.0.14.tar.bz2
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# extract fastx_toolkit source tarball
mydir=/home/ubuntu/emmes_install/fastx_toolkit_install/fastx_toolkit-0.0.14
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	tar -xvjf $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot extract "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi

# build fastx_toolkit binary as needed
cd $mydir
myfile=/usr/local/bin/fastx_uncollapser
if [ ! -f "$myfile" ]
then
	echo building fastx_toolkit
	./configure 
	make
	sudo make install
	if [ ! -f "$myfile" ]
	then
		echo could not build fastx_toolkit
		echo exiting with error code 1 ...
		exit 1
	fi
fi

##Adding the new shared library to the PATH for it to work
sudo ldconfig
#export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH 

mybinfile=/usr/local/bin/fastx_uncollapser
if [ -f "$mybinfile" ]
   then
      echo fastx_toolkit is already installed, checking version
      $mybinfile -h 2>&1 | grep 0\.0\.14 && (
         echo fastx_toolkit version is 0.0.14, nothing to do
         exit 1
      )
      $mybinfile -h 2>&1 | grep 0\.0\.14 || (
         echo fastx_toolkit is installed but version is not 0\.0\.14, presumed old
         exit 1
      )
fi

echo Fastx_toolkit is INSTALLED

##HISAT 2.1.0
#ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/

# create/cd into Hisat installation dir
mydir=/home/ubuntu/emmes_install/hisat_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# download hisat binary tarball
myfile=hisat2-2.1.0-Linux_x86_64.zip
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/$myfile
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# unzip hisat binary tarball
mydir=/home/ubuntu/emmes_install/hisat_install/hisat2-2.1.0
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	unzip $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot unzip "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi

#hisat version check

hisatversioncheck=`/usr/bin/hisat2-align-s --version | grep -o 2\.1\.0 `
myfile=/usr/bin/orig_hisat2
myfile2=/usr/bin/hisat2
if [ ! -f "$myfile2" ] ##doesnt exist
then
	##copy over executables from current install directory
	echo copying new hisat2 executables into /usr/bin
	sudo cp -a /home/ubuntu/emmes_install/hisat_install/hisat2-2.1.0/hisat* /usr/bin/
		if [ ! -f /usr/bin/hisat2 ]
		then
			echo could not copy new hisat2 executables to /usr/bin
			echo exiting with error code 1 ...
			exit 1
		fi
elif [ -f "$myfile2" ] && [ "$hisatversioncheck" != "2.1.0" ] ##exists but not the right version
then
	##make copies and copy over executables from current install directory
	echo hisat2 version is not 2.1.0, presumed old
	echo archiving original hisat2 executables in /usr/bin
	for file in $(ls /usr/bin/hisat2*)
	do
		sudo mv -v "$file" "/usr/bin/orig_${file##*/}"
	done
	if [ ! -f "$myfile" ]
	then
		echo could not archive original hisat2 executables
		echo exiting with error code 1 ...
		exit 1
	fi
else
	## hisat2 should be present and current
	echo ls-ing files in /usr/bin/ containing hisat2
    ls -alh /usr/bin/*hisat2*
	echo hisat2 version is OK at 2.1.0	and INSTALLED
fi
echo HISAT INSTALLED

#Htop
# create/cd into htop installation dir
mydir=/home/ubuntu/emmes_install/htop_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

# get htop via git
myfile=htop/configure.ac
if [ ! -f "$myfile" ]
then
	echo gitting htop from github repository
	git clone https://github.com/hishamhm/htop.git
	if [ ! -f "$myfile" ]
	then
		echo could not git file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

cat $myfile | grep 2\.0\.2 && (
   echo htop version obtained via git is OK at 2.0.2
)

cat $myfile | grep 2\.0\.2 || (
   echo htop version obtained via git is not 2.0.2
   echo exiting with error code 1 ...
   exit 1
)

# build htop binary as needed
mydir=/home/ubuntu/emmes_install/htop_install/htop
cd $mydir
myfile=htop
if [ ! -f "$myfile" ]
   then
      echo building htop
      ./autogen.sh && ./configure && make
      if [ ! -f "$myfile" ]
         then
            echo could not build Htop
            echo exiting with error code 1 ...
            exit 1
      fi
fi

# install htop binary
mybinfile=/usr/bin/htop

if [ ! -f $mybinfile ]
   then
   sudo cp -a $myfile $mybinfile
   if [ ! -f "$mybinfile" ]
      then
         echo could not install "$mybinfile"
         echo exiting with error code 1 ...
         exit 1
   fi
fi

$mybinfile --help 2>&1 | grep 2\.0\.2 && (
   echo htop version in /usr/bin is OK at 2.0.2
)

$mybinfile --help 2>&1 | grep 2\.0\.2 || (
   echo htop version in /usr/bin is not 2.0.2
   echo installing executable obtained via git
   sudo cp -a $myfile $mybinfile
   $mybinfile --help 2>&1 | grep 2\.0\.2 || (
         echo could not install htop version 2.0.2
         echo exiting with error code 1 ...
         exit 1
   )
)
echo htop is INSTALLED

##INSTALL R (r-base)

# may need to add a couple of public keys before updating the repos
# but neither of these worked, nor with keyserver.ubuntu.com
# sudo apt-key adv --keyserver nebc.nerc.ac.uk --recv-keys 679917B3C53271E8
# sudo apt-key adv --keyserver research.cs.wisc.edu --recv-keys 973FC7D2670079F6

# ran these but most recent available version is 3.0.1
# sudo apt-get update
# sudo apt-get install r-base

# if r-base 3.0.1 is installed, remove via apt-get
#    note that removing r-base alone did not remove /usr/bin/R

myfile=/usr/bin/R
if [ -f "$myfile" ]
   then
      /usr/bin/R --version | grep 3\.4\.1 && (
         echo R is installed and is at version 3.4.1, so do no more
      )
      /usr/bin/R --version | grep 3\.4\.1 || (
	  # extra work to get the right version
         echo R is installed but version is not 3.4.1, presumed old
        # /usr/bin/R --version | grep 3\.2\.5 || (
        #    echo R is installed but version is not at 3.2.5, so do not remove it via apt-get
        # )
        # /usr/bin/R --version | grep 3\.0\.1 && (
            echo R version is presumed old and installed via apt-get
            echo so remove it via apt-get
            sudo apt-get remove r-base
            sudo apt-get remove r-base-core
            if [ -f "$myfile" ]
               then
                  /usr/bin/R --version | grep 3\.0\.1 && (
                     echo R version is 3.0.1, but could not remove via apt-get
                     exit 1
                  )
            fi
        # )
      )
fi

# at this point, R should be either missing or at version 3.4.1 if present
# source tarball is available for 3.4.1
# there is no separate tarball for r-base-core
# May require openjdk09
#apt-get install openjdk-9-jdk
#rm -rf /usr/lib/jvm/default-java
#ln -s /usr/lib/jvm/java-9-openjdk-amd64/ /usr/lib/jvm/default-java
#
#
echo at this point, R should be either missing, unable to be removed via apt-get or at version 3.4.1 if present

myfile=/usr/bin/R
if [ ! -f "$myfile" ]
then
echo creating/cd-ing into r-base installation dir
mydir=/home/ubuntu/emmes_install/r-base_install
if [ ! -d "$mydir" ]
 then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	   then
		  echo cannot create directory "$mydir"
		  echo exiting with error code 1 ...
		  exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"

#     download r-base source tarball @ https://cran.r-project.org/src/base/R-3/R-3.4.1.tar.gz
myfile2=R-3.4.1.tar.gz
if [ ! -f "$myfile2" ]
then
	echo wget-ing file "$myfile2"
	wget https://cran.r-project.org/src/base/R-3/R-3.4.1.tar.gz
	if [ ! -f "$myfile2" ]
	   then
		  echo could not wget file "$myfile2"
		  echo exiting with error code 1 ...
		  exit 1
	fi
fi
echo ls-ing file "$myfile2"
ls -alh "$myfile2"
echo extracting r-base source tarball and cd into installation dir
mydir=/home/ubuntu/emmes_install/r-base_install/R-3.4.1
if [ ! -d "$mydir" ]
 then
	echo unzipping file "$myfile2"
	tar xvzf $myfile2
	if [ ! -d "$mydir" ]
	   then
		  echo cannot extract "$myfile2" to create "$mydir"
		  echo exiting with error code 1 ...
		  exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"
#     build and install r-base according to docs/INSTALL file
#     note this is a little ugly because on the very first pass /usr/bin/subread-align does not exist
#     expecting to produce file /usr/local/bin/R
./configure --with-x=no
make
sudo make install
#     make install-info
#     make install-pdf
sudo ln -s /usr/local/bin/R /usr/bin/R
ls -al /usr/local/bin/R /usr/bin/R
if [ ! -f "$myfile" ]
 then
	echo /usr/bin/R still does not exist, installation deemed failed
	exit 1
 else
	/usr/bin/R --version | grep 3\.4\.1 || (
	echo R is not version 3.4.1, installation deemed failed
	)
	/usr/bin/R --version | grep 3\.4\.1 && (
	echo R has been installed from source and is at version 3.4.1
	)
fi
fi
echo R is INSTALLED


##Install Devtools (for specific verion of R packages)
myfile=/usr/local/lib/R/library/devtools/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package devtools is already installed, so do nothing
else
	echo installing package devtools in R
	echo install.packages\(\"devtools\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_devtools.R
	cat install_devtools.R
	sudo R CMD BATCH install_devtools.R
	cat install_devtools.Rout
	rm --interactive=never install_devtools.R install_devtools.Rout
fi
echo R package devtools INSTALLED

##Install RcppEigen (for requirements for car)
#  install.packages("RcppEigen", repos = "http://cran.us.r-project.org")



##install biocLite.R in R (3.5)
##the same two commands can test for a previous biocLite installation
##and install it if not already installed
##so just grep the output for "Using Bioconductor 3.5" to verify a working installation

echo testing for bioconductor version 3.5 in R, and installing if necessary
echo source\(\"http://bioconductor.org/biocLite.R\"\) > install_biocLite.R
echo biocLite\(\) >> install_biocLite.R
cat install_biocLite.R
sudo R CMD BATCH install_biocLite.R
cat install_biocLite.Rout
grep "Using Bioconductor 3.5" install_biocLite.Rout || (
   echo biocLite version 3.5 not detected, installation deemed failed
   exit 1
)
grep "Using Bioconductor 3.5" install_biocLite.Rout && (
   echo biocLite version 3.5 detected
)
pwd
ls -al install_biocLite.R install_biocLite.Rout
rm --interactive=never install_biocLite.R install_biocLite.Rout

##install some biocLite packages (in R) if they are not already installed
# echo libs = c('edgeR','biomaRt','goseq','xtable','MASS','gplots','sqldf','pvclust','car','knitr');
# install Gviz# install impute# install affy

cd /home/ubuntu/emmes_install

echo checking for prior installation of R packages edgeR, goseq, biomaRt, erccdashboard, DESeq2, PROPER, Gviz, impute, affy and limma
echo libs \= c\(\'edgeR\'\,\'goseq\'\,\'biomaRt\'\,\'erccdashboard\'\,\'DESeq2\'\,\'PROPER\'\,\'Gviz\'\,\'impute\'\,\'affy\'\,\'limma\'\)\; > check_R_packages.R
echo lapply\(libs\, require\, character\.only\=T\) >> check_R_packages.R
cat  check_R_packages.R
sudo R CMD BATCH check_R_packages.R
cat check_R_packages.Rout
grep FALSE  check_R_packages.Rout || (
   echo R packages edgeR, goseq, biomaRt, erccdashboard, DESeq2, PROPER, Gviz, impute, affy and limma already installed, so do nothing
  rm --interactive=never check_R_packages.Rout
)
grep FALSE  check_R_packages.Rout && (
   echo installing biocLite packages edgeR, goseq, biomaRt, erccdashboard, DESeq2, PROPER, Gviz, impute, affy and limma
   echo source\(\"http://bioconductor.org/biocLite.R\"\) > install_biocLite_packages.R
   echo biocLite\(\'edgeR\'\) >> install_biocLite_packages.R
   echo biocLite\(\'goseq\'\) >> install_biocLite_packages.R
   echo biocLite\(\'biomaRt\'\) >> install_biocLite_packages.R
   echo biocLite\(\'erccdashboard\'\) >> install_biocLite_packages.R
   echo biocLite\(\'DESeq2\'\) >> install_biocLite_packages.R
   echo biocLite\(\'PROPER\'\) >> install_biocLite_packages.R
   echo biocLite\(\'Gviz\'\) >> install_biocLite_packages.R
   echo biocLite\(\'impute\'\) >> install_biocLite_packages.R
   echo biocLite\(\'affy\'\) >> install_biocLite_packages.R
   echo biocLite\(\'limma\'\) >> install_biocLite_packages.R
   cat install_biocLite_packages.R
   sudo R CMD BATCH install_biocLite_packages.R
   cat install_biocLite_packages.Rout
   egrep "DONE \(edgeR\)" install_biocLite_packages.Rout || (
      echo biocLite package edgeR not installed, installation deemed failed
      exit 1
   )
   egrep "DONE \(edgeR\)" install_biocLite_packages.Rout && (
      echo biocLite package edgeR INSTALLED
   )
   egrep "DONE \(goseq\)" install_biocLite_packages.Rout || (
      echo biocLite package goseq not installed, installation deemed failed
      exit 1
   )
   egrep "DONE \(goseq\)" install_biocLite_packages.Rout && (
      echo biocLite package goseq INSTALLED
   )
   egrep "DONE \(biomaRt\)" install_biocLite_packages.Rout || (
      echo biocLite package biomaRt not installed, installation deemed failed
      exit 1
   )
   egrep "DONE \(biomaRt\)" install_biocLite_packages.Rout && (
      echo biocLite package biomaRt INSTALLED
   )
   egrep "DONE \(erccdashboard\)" install_biocLite_packages.Rout || (
      echo biocLite package erccdashboard not installed, installation deemed failed
      exit 1
   )
   egrep "DONE \(erccdashboard\)" install_biocLite_packages.Rout && (
      echo biocLite package erccdashboard INSTALLED
   )
   egrep "DONE \(DESeq2\)" install_biocLite_packages.Rout || (
      echo biocLite package DESeq2 not installed, installation deemed failed
      exit 1
   )
   egrep "DONE \(DESeq2\)" install_biocLite_packages.Rout && (
      echo biocLite package DESeq2 INSTALLED
   )
   egrep "DONE \(PROPER\)" install_biocLite_packages.Rout || (
      echo biocLite package PROPER not installed, installation deemed failed
      exit 1
   )
   egrep "DONE \(PROPER\)" install_biocLite_packages.Rout && (
      echo biocLite package PROPER INSTALLED
   )
   egrep "DONE \(Gviz\)" install_biocLite_packages.Rout || (
      echo biocLite package Gviz not installed, installation deemed failed
      exit 1
   )
   egrep "DONE \(Gviz\)" install_biocLite_packages.Rout && (
      echo biocLite package INSTALLED
   )
   egrep "DONE \(impute\)" install_biocLite_packages.Rout || (
      echo biocLite package impute not installed, installation deemed failed
      exit 1
   )
   egrep "DONE \(impute\)" install_biocLite_packages.Rout && (
      echo biocLite package impute INSTALLED
   )
   egrep "DONE \(affy\)" install_biocLite_packages.Rout || (
      echo biocLite package affy not installed, installation deemed failed
      exit 1
   )
   egrep "DONE \(affy\)" install_biocLite_packages.Rout && (
      echo biocLite package affy INSTALLED
   )
   egrep "DONE \(limma\)" install_biocLite_packages.Rout || (
      echo biocLite package limma not installed, installation deemed failed
      exit 1
   )
   egrep "DONE \(limma\)" install_biocLite_packages.Rout && (
      echo biocLite package limma INSTALLED
   )
   rm --interactive=never install_biocLite_packages.R install_biocLite_packages.Rout
)


# install some more packages (in R) if they are not already installed the R install.packages way or using devtools

# install pvclust in R
myfile=/usr/local/lib/R/library/pvclust/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package pvclust is already installed, so do nothing
else
	echo installing package pvclust in R
	echo install.packages\(\"pvclust\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_pvclust.R
	cat install_pvclust.R
	sudo R CMD BATCH install_pvclust.R
	cat install_pvclust.Rout
	rm --interactive=never install_pvclust.R install_pvclust.Rout
fi
echo R package pvclust INSTALLED

# install xtable in R
myfile=/usr/local/lib/R/library/xtable/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package xtable is already installed, so do nothing
else
	echo installing package xtable in R
	echo install.packages\(\"xtable\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_xtable.R
	cat install_xtable.R
	sudo R CMD BATCH install_xtable.R
	cat install_xtable.Rout
	rm --interactive=never install_xtable.R install_xtable.Rout
fi
echo R package xtable INSTALLED

# install MASS in R
myfile=/usr/local/lib/R/library/MASS/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package MASS is already installed, so do nothing
else
	echo installing package MASS in R
	echo install.packages\(\"MASS\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_MASS.R
	cat install_MASS.R
	sudo R CMD BATCH install_MASS.R
	cat install_MASS.Rout
	rm --interactive=never install_MASS.R install_MASS.Rout
fi
echo R package MASS INSTALLED

# install gplots in R
myfile=/usr/local/lib/R/library/gplots/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package gplots is already installed, so do nothing
else
	echo installing package gplots in R
	echo install.packages\(\"gplots\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_gplots.R
	cat install_gplots.R
	sudo R CMD BATCH install_gplots.R
	cat install_gplots.Rout
	rm --interactive=never install_gplots.R install_gplots.Rout
fi
echo R package gplots INSTALLED

# install sqldf in R
myfile=/usr/local/lib/R/library/sqldf/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package sqldf is already installed, so do nothing
else
	echo installing package sqldf in R
	echo install.packages\(\"sqldf\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_sqldf.R
	cat install_sqldf.R
	sudo R CMD BATCH install_sqldf.R
	cat install_sqldf.Rout
	rm --interactive=never install_sqldf.R install_sqldf.Rout
fi
echo R package sqldf INSTALLED

# install knitr in R
myfile=/usr/local/lib/R/library/knitr/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package knitr is already installed, so do nothing
else
	echo installing package knitr in R
	echo install.packages\(\"knitr\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_knitr.R
	cat install_knitr.R
	sudo R CMD BATCH install_knitr.R
	cat install_knitr.Rout
	rm --interactive=never install_knitr.R install_knitr.Rout
fi
echo R package knitr INSTALLED

# install car in R
myfile=/usr/local/lib/R/library/car/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package car is already installed, so do nothing
else
	echo installing package car in R
	#install.packages("RcppEigen", repos = "http://cran.us.r-project.org")
	echo install.packages\(\"RcppEigen\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_car.R
	echo install.packages\(\"car\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_car.R
	cat install_car.R
	sudo R CMD BATCH install_car.R
	cat install_car.Rout
	rm --interactive=never install_car.R install_car.Rout
fi
echo R package car INSTALLED

# install R.utils in R
myfile=/usr/local/lib/R/library/R.utils/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package R.utils is already installed, so do nothing
else
	echo installing package R.utils in R
	echo install.packages\(\"R.utils\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_R.utils.R
	cat install_R.utils.R
	sudo R CMD BATCH install_R.utils.R
	cat install_R.utils.Rout
	rm --interactive=never install_R.utils.R install_R.utils.Rout
fi
echo R package R.utils INSTALLED

# install vegan in R
myfile=/usr/local/lib/R/library/vegan/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package vegan is already installed, so do nothing
else
	echo installing package vegan in R
	echo install.packages\(\"vegan\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_vegan.R
	cat install_vegan.R
	sudo R CMD BATCH install_vegan.R
	cat install_vegan.Rout
	rm --interactive=never install_vegan.R install_vegan.Rout
fi
echo R package vegan INSTALLED

# install openxlsx in R
myfile=/usr/local/lib/R/library/openxlsx/libs/openxlsx.so
if [ -f "$myfile" ]
then
	echo R package openxlsx is already installed, so do nothing
else
	echo installing legacy archive package openxlsx in R
	echo first downloading tarball and extracting it into an installation directory
	mkdir /home/ubuntu/emmes_install/openxlsx_install
	cd /home/ubuntu/emmes_install/openxlsx_install
	wget https://cran.r-project.org/src/contrib/openxlsx_4.0.17.tar.gz
	#     tar -zxvf openxlsx_3.0.0.tar.gz #no extraction needed
	#     create R script to install package
	echo pkgFile = \"/home/ubuntu/emmes_install/openxlsx_install/openxlsx_4.0.17.tar.gz\" > install_openxlsx.R
	#     load package to the shared R library
	echo .libPaths\(\"/usr/local/lib/R/library\"\) >> install_openxlsx.R
	#     Install package . Don.t specify repository . build from archive
	echo install.packages\(pkgs\=pkgFile, type\=\"source\", repos\=NULL\) >> install_openxlsx.R
	#     print R script
	cat install_openxlsx.R
	#     Run R script
	sudo R CMD BATCH install_openxlsx.R
	cat install_openxlsx.Rout
	#     clean up
	rm --interactive=never openxlsx_4.0.17.tar.gz install_openxlsx.R
fi
echo R package openxlsx INSTALLED

# install Cairo for R
myfile=/usr/local/lib/R/library/Cairo/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package Cairo is already installed, so do nothing
else
	echo installing package Cairo in R using devtools
	echo require\(devtools\) > install_Cairo.R
	echo install_version\(\"Cairo\"\, version \= \"1\.5\-9\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_Cairo.R
	cat install_Cairo.R
	sudo R CMD BATCH install_Cairo.R
	cat install_Cairo.Rout
	rm --interactive=never install_Cairo.R install_Cairo.Rout
fi
echo R package Cairo INSTALLED

# install seqinr for R
myfile=/usr/local/lib/R/library/seqinr/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package seqinr is already installed, so do nothing
else
	echo installing package seqinr in R using devtools
	echo require\(devtools\) > install_seqinr.R
	echo install_version\(\"seqinr\"\, version \= \"3\.3\-3\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_seqinr.R
	cat install_seqinr.R
	sudo R CMD BATCH install_seqinr.R
	cat install_seqinr.Rout
	rm --interactive=never install_seqinr.R install_seqinr.Rout
fi
echo R package seqinr INSTALLED

# install iterpc for R
myfile=/usr/local/lib/R/library/iterpc/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package iterpc is already installed, so do nothing
else
	echo installing package iterpc in R using devtools
	echo require\(devtools\) > install_iterpc.R
	echo install_version\(\"iterpc\"\, version \= \"0\.3\.0\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_iterpc.R
	cat install_iterpc.R
	sudo R CMD BATCH install_iterpc.R
	cat install_iterpc.Rout
	rm --interactive=never install_iterpc.R install_iterpc.Rout
fi
echo R package iterpc INSTALLED

# install gridExtra for R
myfile=/usr/local/lib/R/library/gridExtra/DESCRIPTION ##change
if [ -f "$myfile" ]
then
	echo R package gridExtra is already installed, so do nothing
else
	echo installing package gridExtra in R using devtools
	echo require\(devtools\) > install_gridExtra.R
	echo install_version\(\"gridExtra\"\, version \= \"2\.2\.1\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_gridExtra.R
	cat install_gridExtra.R
	sudo R CMD BATCH install_gridExtra.R
	cat install_gridExtra.Rout
	rm --interactive=never install_gridExtra.R install_gridExtra.Rout
fi
echo R package gridExtra INSTALLED

# install UpSetR for R
myfile=/usr/local/lib/R/library/UpSetR/DESCRIPTION ##change
if [ -f "$myfile" ]
then
	echo R package UpSetR is already installed, so do nothing
else
	echo installing package UpSetR in R using devtools
	echo require\(devtools\) > install_UpSetR.R
	echo install_version\(\"UpSetR\"\, version \= \"1\.3\.3\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_UpSetR.R
	cat install_UpSetR.R
	sudo R CMD BATCH install_UpSetR.R
	cat install_UpSetR.Rout
	rm --interactive=never install_UpSetR.R install_UpSetR.Rout
fi
echo R package UpSetR INSTALLED

# install mixOmics for R
myfile=/usr/local/lib/R/library/mixOmics/DESCRIPTION ##change
if [ -f "$myfile" ]
then
	echo R package mixOmics is already installed, so do nothing
else
	echo installing package mixOmics in R using devtools
	echo require\(devtools\) > install_mixOmics.R
	echo install_version\(\"mixOmics\"\, version \= \"6\.1\.3\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_mixOmics.R
	cat install_mixOmics.R
	sudo R CMD BATCH install_mixOmics.R
	cat install_mixOmics.Rout
	rm --interactive=never install_mixOmics.R install_mixOmics.Rout
fi
echo R package mixOmics INSTALLED

# install glmnet in R
myfile=/usr/local/lib/R/library/glmnet/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package glmnet is already installed, so do nothing
else
	echo installing package glmnet in R
	echo install.packages\(\"glmnet\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_glmnet.R
	cat install_glmnet.R
	sudo R CMD BATCH install_glmnet.R
	cat install_glmnet.Rout
	rm --interactive=never install_glmnet.R install_glmnet.Rout
fi
echo R package glmnet INSTALLED

# install networkx
##pip install networkx==1.11
##Need to add a part for adding downloaded file to the emmes_install aread
 
myfile=/usr/local/lib/python2.7/dist-packages/networkx
if [ -d "$myfile" ]
then
	echo networkx is already installed
else
	sudo pip install networkx==1.11
	#     make aws help generally available
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg/EGG-INFO/top_level.txt
	#     resolve permissions issue encountered by Travis
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg/EGG-INFO/requires.txt
	if [ ! -d "$myfile" ]
	then
	   echo could not install networkx
	   exit 1
	fi
fi
echo networkx INSTALLED

# install pysam 0.12.0.1
##sudo pip install pysam
myfile=/usr/local/lib/python2.7/dist-packages/pysam
if [ -d "$myfile" ]
then
	echo pysam is already installed
else
	sudo pip install pysam
	#     make aws help generally available
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg/EGG-INFO/top_level.txt
	#     resolve permissions issue encountered by Travis
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg/EGG-INFO/requires.txt
	if [ ! -d "$myfile" ]
	then
	   echo could not install pysam
	   exit 1
	fi
fi
echo pysam INSTALLED

# install cython
##sudo pip install cython 0.27
myfile=/usr/local/lib/python2.7/dist-packages/cython
if [ -d "$myfile" ]
then
	echo cython is already installed
else
	sudo pip install cython
	#     make aws help generally available
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg/EGG-INFO/top_level.txt
	#     resolve permissions issue encountered by Travis
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg/EGG-INFO/requires.txt
	if [ ! -d "$myfile" ]
	then
	   echo could not install cython
	   exit 1
	fi
fi
echo cython INSTALLED

##CLEANUP

#Remove swap
sudo /sbin/swapoff /var/swap.1
sudo rm -f /var/swap.1


echo Script Completed.
##end of script 
exit 1