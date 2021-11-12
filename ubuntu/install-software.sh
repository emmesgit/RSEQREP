#!/usr/bin/env bash
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
# Program:  install-software.pl 
# Version:  RSEQREP 2.1.2
# Author:   Kevin Conway and Leigh Villarroel
# Purpose:  install software onto 18.0.4 Ubuntu desktop computer required for RSEQREP usage.
# Input:    N/A
# Output:   N/A
##############################################################################################################################

##Add Swap for extra memory if needed 
sudo /bin/dd if=/dev/zero of=/var/swap.1 bs=1M count=1024
sudo /sbin/mkswap /var/swap.1
sudo chmod 600 /var/swap.1
sudo /sbin/swapon /var/swap.1


##Initial setup

sudo apt-get update

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
echo installing other required packages for software

sudo apt-get -y install aptitude
sudo aptitude -y install libreadline-dev	
sudo apt-get -y install libncurses5-dev #Needed for htop/SAMtools install
sudo apt-get -y install libncursesw5-dev #Needed for htop install
sudo apt-get -y install zlib1g-dev #Needed for SAMtools install
sudo apt-get -y install libcurl4-openssl-dev #Needed for SAMtools install
sudo apt-get -y install libssl-dev #Needed for SAMtools install
sudo apt-get -y install autoconf #Needed for htop install
sudo apt-get -y install automake #Needed for htop install
sudo apt-get -y install libtool #Needed for htop install
sudo apt-get -y install gfortran #Needed for R-Base install
sudo apt-get -y install libxml2-dev #Needed for BiocLite('biomaRt') to run successfully
sudo apt-get -y install libcairo2-dev #Needed for Cairo package in R
sudo apt-get -y install libxt-dev #Needed for Cairo package in R
sudo apt-get -y install libgmp3-dev #Needed for iterpc package in R
sudo apt-get -y install liblzma-doc liblzma-dev xz-utils #for samtools/htslib-1.5
sudo apt-get -y install libbz2-dev #SAMTools
sudo apt-get -y install libtbb2 #for bowtie2
sudo apt-get -y install openjdk-8-jdk #for R 3.6.1
sudo apt-get -y install default-jdk #for R packages
sudo apt-get -y install mesa-common-dev libglu1-mesa-dev  #R package rgl
#sudo apt-get -y install python-setuptools
sudo apt-get -y install python3-dev python3-setuptools #needed for psutil Utility 
sudo apt-get -y install  libmariadb-client-lgpl-dev #for RMySQL package which is needed for goseq package
sudo apt-get -y install gcc-6 g++-6 g++-6-multilib gfortran-6 #for gcc-6 which could be needed for fastx compiling
sudo apt-get -y install libffi-dev #python3
sudo apt-get -y install tcl-dev tk-dev #Needed for sm R package
	
echo gcc and other necessary compiler tools INSTALLED

##EVINCE
sudo apt-get -y install evince

##IMAGEMAGICK
sudo apt-get -y install imagemagick

##OPENSSL
# create/cd into openssl installation dir
mydir=/home/ubuntu/emmes_install/openssl_install
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

# download openssl binary tarball
myfile=openssl-1.1.1b.tar.gz
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget --no-check-certificate https://www.openssl.org/source/openssl-1.1.1b.tar.gz
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# unzip openssl binary tarball and build
mydir=/home/ubuntu/emmes_install/openssl_install/openssl-1.1.1b
if [ ! -d "$mydir" ]
then
	echo unzipping file "$myfile"
	tar xvfz $myfile
	if [ ! -d "$mydir" ]
	then
		echo cannot unzip "$myfile" to create "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
	
	cd $mydir
	#build openssl
	./config  #./config -Wl,--enable-new-dtags,-rpath,'$(LIBRPATH)'
	sudo make
	#make test
	sudo make install #in /usr/local/ssl/bin/openssl?
fi

## if necessary, archive original openssl executables in /usr/bin
## copy new openssl executables from openssl installation dir to /usr/bin
opensslversioncheck=`/usr/bin/openssl version | grep -o 1\.1\.1b `
myfile=/usr/bin/orig_openssl
myfile2=/usr/bin/openssl
if [ -f "$myfile2" ] && [ "$opensslversioncheck" != "1.1.1b" ] ##exists but not the right version
then
	##make copies and copy over executables from current install directory
	echo openssl version is not 1.1.1b, presumed old
	echo archiving original openssl executables in /usr/bin
	
	sudo mv -v "/usr/bin/openssl" "/usr/bin/orig_openssl"
	
	if [ ! -f "$myfile" ]
	then
		echo could not archive original openssl executables
		echo exiting with error code 1 ...
		exit 1
	else
		#move the new openssl version	
		echo moving openssl executables to /usr/bin
		sudo mv /usr/local/bin/openssl /usr/bin
		echo openssl version is at 1.1.1b and INSTALLED
	fi
else
	## openssl should be present and current
	echo ls-ing files in /usr/bin/
    ls -alh /usr/bin/openssl
	echo openssl version is OK at 1.1.1b and INSTALLED
fi



#PYTHON 2
mypython27=/usr/bin/python2.7
mypython=/usr/bin/python
if [ ! -f "$mypython" ] && [ ! -f "$mypython27" ]
then
	sudo apt-get -y install python2.7 python-dev libpython2.7-dev 
	sudo cp /usr/bin/python2.7 /usr/bin/python2
elif [ ! -f "$mypython" ] && [ -f "$mypython27" ]
then
	sudo apt-get -y install python-dev libpython2.7-dev 
	sudo cp /usr/bin/python2.7 /usr/bin/python2
fi
echo "Installing required python modules."
sudo apt-get -y install python-dev libpython2.7-dev #Gives the Python.h header files for rseqc install
sudo apt-get -y install libbz2-1.0 libbz2-dev libbz2-ocaml libbz2-ocaml-dev #for samtools/htslib-1.5
#Make python bin python2

sudo mv /usr/bin/python /usr/bin/python2

echo python INSTALLED

##PYTHON3 
#Python3.7.3
mypython3=/usr/bin/python3
if [ ! -f "$mypython3" ]
then
	sudo apt-get update
	sudo apt-get -y install python3.7-dev
	sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.6 1
	sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.7 2
	#Point python3.7 to python3
	sudo update-alternatives --config python3
fi

##RSeQC Needs Python3 so fluip the current python with the latest version

sudo cp -a /usr/bin/python3 /usr/bin/python

#Needs another update before pip can be installed
sudo apt-add-repository universe
sudo apt-get update


##PIP
mypip=/usr/bin/pip
if [ ! -f "$mypip" ]
then
	sudo apt-get -y install python-pip 
	sudo pip install --upgrade pip
fi
echo pip INSTALLED

##PIP3
mypip=/usr/bin/pip3
if [ ! -f "$mypip" ]
then
	sudo apt-get -y install python3-pip
fi
echo pip3 INSTALLED

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

##X2GO--> Set up AFTER the script is run


###UTILITIES

# install pysam 0.12.0.1
##sudo pip install pysam
myfile=/usr/local/lib/python3.6/dist-packages/pysam
if [ -d "$myfile" ]
then
	echo pysam is already installed
else
	sudo pip install pysam
	#     make aws help generally available
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.4-py2.7-linux-x86_64.egg/EGG-INFO/top_level.txt
	#     resolve permissions issue encountered by Travis
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.4-py2.7-linux-x86_64.egg/EGG-INFO/requires.txt
	if [ ! -d "$myfile" ]
	then
	   echo could not install pysam
	   exit 1
	fi
fi
echo pysam INSTALLED

# install cython
##sudo pip install cython 0.27
myfile=/usr/local/lib/python3.6/dist-packages/Cython
if [ -d "$myfile" ]
then
	echo cython is already installed
else
	sudo pip install cython
	#     make aws help generally available
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.4-py2.7-linux-x86_64.egg/EGG-INFO/top_level.txt
	#     resolve permissions issue encountered by Travis
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.4-py2.7-linux-x86_64.egg/EGG-INFO/requires.txt
	if [ ! -d "$myfile" ]
	then
	   echo could not install cython
	   exit 1
	fi
fi
echo cython INSTALLED

# install amazon command line interface software (AWSCLI)
myfile=/usr/local/bin/aws
if [ -f "$myfile" ]
then
	echo awscli is already installed
else
	#try apt then try pip
	#sudo apt -y install awscli
		echo pip is installed, proceed with installation of aws cli
		sudo pip install awscli==1.16.152
	# make aws help generally available
		sudo chmod o+r /usr/local/lib/python3.6/dist-packages/
	# resolve permissions issue encountered by Travis
		sudo chmod o+r /usr/local/lib/python3.6/dist-packages/
		if [ ! -f "$myfile" ]
		   then
			  echo could not install awscli
			  exit 1
		fi
	

fi
echo awscli is INSTALLED

##SNAKEMAKE (NEEDS PYTHON3)
myfile=/usr/local/lib/python3.6/dist-packages/snakemake
if [ -d "$myfile" ]
then
	echo snakemake is already installed
else
	#Needs datrie which seems to be having issues  https://bitbucket.org/snakemake/snakemake/issues/934/installation-failed-in-python-37
	sudo pip3 install git+https://github.com/pytries/datrie.git
	sudo pip3 install snakemake==5.4.5
	#     make aws help generally available
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.4-py2.7-linux-x86_64.egg/EGG-INFO/top_level.txt
	#     resolve permissions issue encountered by Travis
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.4-py2.7-linux-x86_64.egg/EGG-INFO/requires.txt
	if [ ! -d "$myfile" ]
	then
	   echo could not install snakemake
	   exit 1
	fi
fi
echo snakemake INSTALLED

##PSUTIL (NEEDS PYTHON3)
myfile=/usr/local/lib/python3.6/dist-packages/psutil #Trying a wildcard instead of the fixed python version to handle complexity
if [ -d "$myfile" ]
then
	echo psutil is already installed
else
	sudo pip3 install psutil==5.6.2
	#     make aws help generally available
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.4-py2.7-linux-x86_64.egg/EGG-INFO/top_level.txt
	#     resolve permissions issue encountered by Travis
	#      sudo chmod o+r /usr/local/lib/python2.7/dist-packages/RSeQC-2.6.4-py2.7-linux-x86_64.egg/EGG-INFO/requires.txt
	if [ ! -d "$myfile" ]
	then
	   echo could not install psutil
	   exit 1
	fi
fi
echo psutil INSTALLED

##INSTALL cutadapt (2.3)
myfile=/usr/local/lib/python3.6/dist-packages/cutadapt
if [ -d "$myfile" ]
then
	echo cutadapt is already installed
else
	sudo python3 -m pip install cutadapt==2.3
	if [ ! -d "$myfile" ]
	then
	   echo could not install cutadapt
	   exit 1
	fi
fi
echo cutadapt INSTALLED


##Sqlite
#wget sqlite-autoconf-3280100.tar.gz
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
myfile=sqlite-autoconf-3280000.tar.gz 
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget --no-check-certificate https://sqlite.org/2019/$myfile
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
mydir=/home/ubuntu/emmes_install/sqlite_install/sqlite-autoconf-3280000
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
#wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz
# create/cd into sratoolkit installation dir
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
myfile=sratoolkit.2.9.6-ubuntu64.tar.gz
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget --no-check-certificate https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/$myfile
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
mydir=/home/ubuntu/emmes_install/sratoolkit_install/sratoolkit.2.9.6-ubuntu64
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

if [ ! -f "/usr/bin/sra-sort.2.9.6" ]
then
	cd $mydir/bin/
	echo copying sratoolkit files to /usr/bin
	sudo cp -a * /usr/bin/
	echo sratoolkit is INSTALLED
fi


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
sudo apt-get -y install libcurl4-openssl-dev #In case it got downgraded..

# download SAMtools source tarball (1.9)
## https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2

myfile=samtools-1.9.tar.bz2
if [ ! -f "$myfile" ]
   then
      echo wget-ing file "$myfile"
      wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.9/$myfile
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
mydir=/home/ubuntu/emmes_install/samtools_install/samtools-1.9
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
      ./configure --enable-plugins --enable-libcurl --with-plugin-path=$PWD/htslib-1.9
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
		echo samtools version is at 1.9 and INSTALLED
		
	elif [ -f "$myfile" ] && [ ! -f "$myfile2" ] 
	then
		##make copies and copy over executables from current install directory
		echo archiving original samtools executables in /usr/bin	
		sudo mv -v "/usr/bin/samtools" "/usr/bin/orig_samtools"
		#move the new samtools version	
		echo copying samtools executables to /usr/bin
		sudo cp -a samtools /usr/bin
		echo samtools version is at 1.9 and INSTALLED
		
	elif [ ! -f "$myfile" ] && [ -f "$myfile2" ] 
	then
		sudo cp -a samtools /usr/bin
		echo samtools version is at 1.9 and INSTALLED
		
	else
		## Both files are there, remove the original, move the old to the orignal then copy over the executable
		sudo rm -f "/usr/bin/orig_samtools"
		##make copies and copy over executables from current install directory
		echo archiving original samtools executables in /usr/bin	
		sudo mv -v "/usr/bin/samtools" "/usr/bin/orig_samtools"
		#move the new samtools version	
		echo copying samtools executables to /usr/bin
		sudo cp -a samtools /usr/bin
		echo samtools version is at 1.9 and INSTALLED
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
myfile=bcftools-1.9.tar.bz2

if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.9/$myfile
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
mydir=/home/ubuntu/emmes_install/bcftools_install/bcftools-1.9
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
	./configure --enable-plugins --enable-libcurl --with-plugin-path=$mydir/htslib-1.9 ##changed plugin-path to point to the install area
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
		echo bcftools version is at 1.9 and INSTALLED
		
	elif [ -f "$myfile" ] && [ ! -f "$myfile2" ] 
	then
		##make copies and copy over executables from current install directory
		echo archiving original bcftools executables in /usr/bin	
		sudo mv -v "/usr/bin/bcftools" "/usr/bin/orig_bcftools"
		#move the new bcftools version	
		echo copying bcftools executables to /usr/bin
		sudo cp -a bcftools /usr/bin
		echo bcftools version is at 1.9 and INSTALLED
		
	elif [ ! -f "$myfile" ] && [ -f "$myfile2" ] 
	then
		sudo cp -a bcftools /usr/bin
		echo bcftools version is at 1.9 and INSTALLED
		
	else
		## Both files are there, remove the original, move the old to the orignal then copy over the executable
		sudo rm -f "/usr/bin/orig_bcftools"
		##make copies and copy over executables from current install directory
		echo archiving original bcftools executables in /usr/bin	
		sudo mv -v "/usr/bin/bcftools" "/usr/bin/orig_bcftools"
		#move the new bcftools version	
		echo copying bcftools executables to /usr/bin
		sudo cp -a bcftools /usr/bin
		echo bcftools version is at 1.9 and INSTALLED
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
myfile=subread-1.6.4-Linux-x86_64.tar.gz
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget --no-check-certificate https://downloads.sourceforge.net/project/subread/subread-1.6.4/$myfile
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
mydir=/home/ubuntu/emmes_install/subread_install/subread-1.6.4-Linux-x86_64/
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

/usr/bin/subread-align -v 2>&1 | grep 1\.6\.4 && (
   echo subread-align version is INSTALLED OK at 1.6.4
)

/usr/bin/subread-align -v 2>&1 | grep 1\.6\.4 || (
	echo subread-align version is not 1.6.4, presumed old
	echo ls-ing archived files in /usr/bin starting with subread
	ls -alh /usr/bin/subread*
	myfile=/home/ubuntu/emmes_install/subread_install/subread-1.6.4-Linux-x86_64/bin/subread-align
	if [ -f "$myfile" ]
	then
		echo copying new subread executables into /usr/bin
		for myotheritem in exactSNP subindel subread-align featureCounts subjunc subread-buildindex
		do 
			#sudo mv /usr/bin/$myotheritem 
			sudo cp -a /home/ubuntu/emmes_install/subread_install/subread-1.6.4-Linux-x86_64/bin/$myotheritem /usr/bin/
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

##PyBigWig needed for RSEQC
sudo pip install pyBigWig

##numpy needed for RSEQC
sudo pip install numpy

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
# http://sourceforge.net/projects/rseqc/files/RSeQC-3.0.0.tar.gz/download
# hacked download location to a mirror, otherwise hangs or fails
myfile=RSeQC-3.0.0.tar.gz
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	#     used direct link because user friendly link failed
	wget --no-check-certificate https://sourceforge.net/projects/rseqc/files/RSeQC-3.0.0.tar.gz/download
	#     make the dir name less cumbersome and more as expected
	echo "if hanging will not get to here"
	pwd
	mv download RSeQC-3.0.0.tar.gz
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
mydir=/home/ubuntu/emmes_install/rseqc_install/RSeQC-3.0.0
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
cp -a setup.py ../../../setup-RSeQC-3.0.0.py
echo RSeQC is INSTALLED


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
myfile=fastqc_v0.11.8.zip
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget --no-check-certificate http://www.bioinformatics.babraham.ac.uk/projects/fastqc/$myfile
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
		#line below may be unnecessary
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
echo fastqc version is 0.11.8 and INSTALLED
 
 
 
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

## check it is currently 1.3-r106

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

cat $myfile | grep 1\.3-r106 && (
   echo seqtk version obtained via git is OK at 1.3-r106
)

cat $myfile | grep 1\.3-r101 || (
   echo seqtk version obtained via git is not 1.3-r106
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

$mybinfile 2>&1 | grep 1\.3-r106 && (
   echo seqtk version in /usr/bin is OK at 1.3-r106
)

$mybinfile 2>&1 | grep 1\.3-r106 || (
   echo seqtk version in /usr/bin is not 1.3-r106
   echo installing executable obtained via git
   sudo cp -a $myfile $mybinfile
   $mybinfile  2>&1 | grep 1\.3-r106 || (
         echo could not install seqtk version 1.3-r106
         echo exiting with error code 1 ...
         exit 1
   )
)
echo Seqtk is INSTALLED


##INSTALL fastx_toolkit
# used sudo apt-get install fastx-toolkit instead of building from source due to an issue with gcc version

# create/cd into fastx_toolkit installation dir
# mydir=/home/ubuntu/emmes_install/fastx_toolkit_install
# if [ ! -d "$mydir" ]
   # then
      # echo creating directory "$mydir"
      # mkdir $mydir
      # if [ ! -d "$mydir" ]
         # then
            # echo cannot create directory "$mydir"
            # echo exiting with error code 1 ...
            # exit 1
      # fi
 # fi
 # echo cd-ing into directory "$mydir"
 # cd "$mydir"

# ##Firstly need to get gtextutils from http://hannonlab.cshl.edu/fastx_toolkit/libgtextutils-0.6.1.tar.bz2 before the fastx is installed
# myfile=libgtextutils-0.6.1.tar.bz2
# if [ ! -f "$myfile" ] #may not work multiple times, need to find a better way to check for this file installed 
# then
	# echo wget-ing file "$myfile"
	# wget --no-check-certificate http://hannonlab.cshl.edu/fastx_toolkit/libgtextutils-0.6.1.tar.bz2
	# if [ ! -f "$myfile" ]
	# then
		# echo could not wget file "$myfile"
		# echo exiting with error code 1 ...
		# exit 1
	# fi
	# #quick approach to extract 
	# tar -xvjf $myfile
	# cd libgtextutils-0.6.1
	# ./configure
	# make
	# sudo make install
	
# #May need to set enviornment variables GTEXTUTILS_CFLAGS and GTEXTUTILS_LIBS

	
	# cd ..
# fi

# ##Now fastx_toolkit
# # getting via git https://github.com/agordon/fastx_toolkit.git  doesnt bring the configure file so getting it via wget
# myfile=fastx_toolkit-0.0.14.tar.bz2
# if [ ! -f "$myfile" ]
# then
	# echo wget-ing file "$myfile"
	# wget --no-check-certificate https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2
	# if [ ! -f "$myfile" ]
	# then
		# echo could not wget file "$myfile"
		# echo exiting with error code 1 ...
		# exit 1
	# fi
# fi
# echo ls-ing file "$myfile"
# ls -alh "$myfile"

# # extract fastx_toolkit source tarball
# mydir=/home/ubuntu/emmes_install/fastx_toolkit_install/fastx_toolkit-0.0.14
# if [ ! -d "$mydir" ]
# then
	# echo unzipping file "$myfile"
	# tar -xvjf $myfile
	# if [ ! -d "$mydir" ]
	# then
		# echo cannot extract "$myfile" to create "$mydir"
		# echo exiting with error code 1 ...
		# exit 1
	# fi
# fi

# # build fastx_toolkit binary as needed
# cd $mydir
# myfile=/usr/local/bin/fastx_uncollapser
# if [ ! -f "$myfile" ]
# then
	# echo building fastx_toolkit
	# ./configure 
	# make
	# sudo make install
	# if [ ! -f "$myfile" ]
	# then
		# echo could not build fastx_toolkit
		# echo exiting with error code 1 ...
		# exit 1
	# fi
# fi

sudo apt-get -y install fastx-toolkit

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
	wget --no-check-certificate ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/$myfile
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

cat $myfile | grep 2\.2\.0 && (
   echo htop version obtained via git is OK at 2.2.0
)

cat $myfile | grep 2\.2\.0 || (
   echo htop version obtained via git is not 2.2.0
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

## install Trimmomatic 
#http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip

# create/cd into trimmomatic installation dir
mydir=/home/ubuntu/emmes_install/trimmomatic_install
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

# download trimmomatic zip
myfile=Trimmomatic-0.38.zip
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget --no-check-certificate http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/$myfile
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# unzip trimmomatic
mydir=/home/ubuntu/emmes_install/trimmomatic_install/Trimmomatic-0.38
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

##Create trimmomatic folder on /usr/local/bin
localdir=/usr/local/bin/trimmomatic
myfile2=/usr/local/bin/trimmomatic/trimmomatic-0.38.jar
if [ ! -f "$myfile2" ] ##doesnt exist
then
	##copy over executables from current install directory
	echo copying new trimmomatic executables into /usr/local/bin
	sudo mkdir $localdir
	sudo cp -a /home/ubuntu/emmes_install/trimmomatic_install/Trimmomatic-0.38/* $localdir
		if [ ! -f $myfile2 ]
		then
			echo could not copy new trimmomatic executables to /usr/bin
			echo exiting with error code 1 ...
			exit 1
		fi
else
	## trimmomatic should be present and current
	echo ls-ing files in /usr/local/bin/ containing trimmomatic
    ls -alh /usr/local/bin/trimmomatic/*
	echo trimmomatic is version 0.38 and INSTALLED
fi
echo Trimmomatic INSTALLED

##BOWTIE2
# create/cd into bowtie installation dir
mydir=/home/ubuntu/emmes_install/bowtie2_install
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

# download bowtie2 binary tarball
myfile=bowtie2-2.3.5-linux-x86_64.zip # replaced with download to improve readability
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget --no-check-certificate https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5/$myfile
	
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# unzip bowtie2 binary tarball
mydir=/home/ubuntu/emmes_install/bowtie2_install/bowtie2-2.3.5-linux-x86_64
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

## if necessary, archive original bowtie2 executables in /usr/bin
## copy new bowtie2 executables from bowtie2 installation dir to /usr/bin
bowtieversioncheck=`/usr/bin/bowtie2 --version | grep -o 2\.3\.5 `
myfile=/usr/bin/orig_bowtie2
myfile2=/usr/bin/bowtie2
if [ ! -f "$myfile2" ] ##doesnt exist
then
	##copy over executables from current install directory
	echo copying new bowtie2 executables into /usr/bin
	sudo cp -a /home/ubuntu/emmes_install/bowtie2_install/bowtie2-2.3.5-linux-x86_64/bowtie2* /usr/bin/
		if [ ! -f /usr/bin/bowtie2 ]
		then
			echo could not copy new bowtie2 executables to /usr/bin
			echo exiting with error code 1 ...
			exit 1
		fi
elif [ -f "$myfile2" ] && [ "$bowtieversioncheck" != "2.3.5" ] ##exists but not the right version #failed out here when installed
then
	##make copies and copy over executables from current install directory
	echo bowtie2 version is not 2.3.5, presumed old
	echo archiving original bowtie2 executables in /usr/bin
	for file in $(ls /usr/bin/bowtie2*)
	do
		sudo mv -v "$file" "/usr/bin/orig_${file##*/}"
	done
	if [ ! -f "$myfile" ]
	then
		echo could not archive original bowtie2 executables
		echo exiting with error code 1 ...
		exit 1
	fi
else
	## bowtie2 should be present and current
	echo ls-ing files in /usr/bin/ containing bowtie2
    ls -alh /usr/bin/*bowtie2*
	echo bowtie2 version is OK at 2.3.5 and INSTALLED
fi


##ENSEMBL-GIT-TOOLS
# create/cd into seqtk installation dir
mydir=/home/ubuntu/emmes_install/ensembl-git-tools_install
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

#git clone https://github.com/Ensembl/ensembl-git-tools.git
#export PATH=$PWD/ensembl-git-tools/bin:$PATH

myfile=ensembl-git-tools/bin/git-ensembl
if [ ! -f "$myfile" ]
then
	echo gitting ensembl-git-tools from github repository
	git clone https://github.com/Ensembl/ensembl-git-tools.git
	if [ ! -f "$myfile" ]
	then
		echo could not git file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

# Symbolic link to /usr/local/bin
mydir=/home/ubuntu/emmes_install/ensembl-git-tools_install/ensembl-git-tools/bin
sudo cp -a ${mydir}/* /usr/local/bin

# check binary
mybinfile=/usr/local/bin/git-ensembl 
if [ ! -f "$mybinfile" ]
   then
      sudo cp -a ${mydir}/* /usr/local/bin
      if [ ! -f "$myfile" ]
         then
            echo could not copy ensembl-git-tools to bin
            echo exiting with error code 1 ...
            exit 1
      fi
else
	echo ensembl-git-tools is INSTALLED
fi

##INSTALL Ensembl Perl API
#Needs
#wget ftp://ftp.ensembl.org/pub/ensembl-api.tar.gz
#wget https://cpan.metacpan.org/authors/id/C/CJ/CJFIELDS/BioPerl-1.6.924.tar.gz 
# create/cd into htop installation dir
mydir=/home/ubuntu/emmes_install/ensembl_install
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
# get and extract files
#wget https://cpan.metacpan.org/authors/id/C/CJ/CJFIELDS/BioPerl-1.6.924.tar.gz 

myfile=BioPerl-1.6.924.tar.gz 
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget --no-check-certificate https://cpan.metacpan.org/authors/id/C/CJ/CJFIELDS/$myfile
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"
echo unzipping file "$myfile"
tar xvzf $myfile

#wget ftp://ftp.ensembl.org/pub/ensembl-api.tar.gz
myfile=ensembl-api.tar.gz
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget --no-check-certificate ftp://ftp.ensembl.org/pub/$myfile
	if [ ! -f "$myfile" ]
	then
		echo could not wget file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"
echo unzipping file "$myfile"
tar xvzf $myfile

#Check 
#mydir=/home/ubuntu/emmes_install/ensembl_install/ensembl

# Manage the ensembl files now
#jump to local file directory where files will be stored
cd "/usr/local/emmes"

# Copy all directories extracted to /usr/local/emmes/ and add stuff
echo copying new ensembl directories into /usr/local/emmes
for myotheritem in BioPerl-1.6.924 ensembl ensembl-compara ensembl-funcgen ensembl-io ensembl-tools ensembl-variation
do 
	sudo cp -a /home/ubuntu/emmes_install/ensembl_install/$myotheritem /usr/local/emmes
	
	#add file locations to PERL5LIB using PWD variable
	#PERL5LIB=${PERL5LIB}:${PWD}/$myotheritem	
done

#export variable
PERL5LIB=${PERL5LIB}:/usr/local/emmes/bioperl-1.6.924
PERL5LIB=${PERL5LIB}:/usr/local/emmes/ensembl/modules
PERL5LIB=${PERL5LIB}:/usr/local/emmes/ensembl-compara/modules
PERL5LIB=${PERL5LIB}:/usr/local/emmes/ensembl-variation/modules
PERL5LIB=${PERL5LIB}:/usr/local/emmes/ensembl-funcgen/modules
PERL5LIB=${PERL5LIB}:/usr/local/emmes/ensembl-io/modules
PERL5LIB=${PERL5LIB}:/usr/local/emmes/ensembl-tools/modules
export PERL5LIB #need to make this persistent

#sudo chmod ugo+x /usr/bin/STAR
echo ensembl version is INSTALLED

##MEME SUITE
#http://meme-suite.org/meme-software/5.0.5/meme-5.0.5.tar.gz

mydir=/home/ubuntu/emmes_install/meme_install
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

#download meme source tarball @ http://meme-suite.org/meme-software/5.0.5/meme-5.0.5.tar.gz
myfile=meme-5.0.5.tar.gz
if [ ! -f "$myfile" ]
then
	echo wget-ing file "$myfile"
	wget --no-check-certificate http://meme-suite.org/meme-software/5.0.5/meme-5.0.5.tar.gz
	if [ ! -f "$myfile" ]
	   then
		  echo could not wget file "$myfile"
		  echo exiting with error code 1 ...
		  exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"
echo extracting meme source tarball and cd into installation dir
mydir=/home/ubuntu/emmes_install/meme_install/meme-5.0.5
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

./configure --prefix=/home/ubuntu/emmes_install/meme_install/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
make
make test
sudo make install

#Copy files to user

sudo cp -a /home/ubuntu/emmes_install/meme_install/meme/bin/* /usr/local/bin/ 


#Still need a check for Meme
ls -al /usr/local/bin/
#if [ ! -f "$myfile" ]
 #then
#	echo /usr/bin/R still does not exist, installation deemed failed
#	exit 1
 #else
#	/usr/bin/R --version | grep 3\.4\.2 || (
#	echo R is not version 3.4.2, installation deemed failed
#	)
#	/usr/bin/R --version | grep 3\.4\.2 && (
#	echo R has been installed from source and is at version 3.4.2
#	)
#fi

echo MEME is INSTALLED


##Pigz
# create/cd into Pigz installation dir
mydir=/home/ubuntu/emmes_install/pigz_install
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

# get pigz via git
# git clone https://github.com/madler/pigz.git

myfile=pigz/pigz.c
if [ ! -f "$myfile" ]
then
	echo gitting pigz from github repository
	git clone https://github.com/madler/pigz.git
	if [ ! -f "$myfile" ]
	then
		echo could not git file "$myfile"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo ls-ing file "$myfile"
ls -alh "$myfile"

cat $myfile | grep 2\.4 && (
   echo pigz version obtained via git is OK at 2.4
)

cat $myfile | grep 2\.4 || (
   echo pigz version obtained via git is not 2.4
   echo exiting with error code 1 ...
   exit 1
)

# build pigz binary as needed
mydir=/home/ubuntu/emmes_install/pigz_install/pigz
cd $mydir
myfile=pigz
if [ ! -f "$myfile" ]
   then
      echo building pigz
      make
      if [ ! -f "$myfile" ]
         then
            echo could not build pigz
            echo exiting with error code 1 ...
            exit 1
      fi
fi

# install pigz binary
mybinfile=/usr/local/bin/pigz

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

$mybinfile 2>&1 | grep 2\.4 && (
   echo pigz version in /usr/bin is OK at 2.4
)

$mybinfile 2>&1 | grep 2\.4 || (
   echo pigz version in /usr/bin is not 2.4
   echo installing executable obtained via git
   sudo cp -a $myfile $mybinfile
   $mybinfile  2>&1 | grep 2\.4 || (
         echo could not install pigz version 2.4
         echo exiting with error code 1 ...
         exit 1
   )
)
echo pigz is INSTALLED

# install pigz binary
myfile=unpigz
mybinfile=/usr/local/bin/unpigz

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

$mybinfile 2>&1 | grep 2\.4 && (
   echo unpigz version in /usr/bin is OK at 2.4
)

$mybinfile 2>&1 | grep 2\.4 || (
   echo unpigz version in /usr/bin is not 2.4
   echo installing executable obtained via git
   sudo cp -a $myfile $mybinfile
   $mybinfile  2>&1 | grep 2\.4 || (
         echo could not install unpigz version 2.4
         echo exiting with error code 1 ...
         exit 1
   )
)
echo unpigz is INSTALLED


##INSTALL R (r-base)

# may need to add a couple of public keys before updating the repos
# but neither of these worked, nor with keyserver.ubuntu.com
# sudo apt-key adv --keyserver nebc.nerc.ac.uk --recv-keys 679917B3C53271E8
# sudo apt-key adv --keyserver research.cs.wisc.edu --recv-keys 973FC7D2670079F6

# ran these but most recent available version is 3.6.0
# sudo apt-get update
# sudo apt-get install r-base

# if r-base 3.6.0 is installed, remove via apt-get
#    note that removing r-base alone did not remove /usr/bin/R

myfile=/usr/bin/R
if [ -f "$myfile" ]
   then
      /usr/bin/R --version | grep 3\.6\.0 && (
         echo R is installed and is at version 3.6.0, so do no more
      )
      /usr/bin/R --version | grep 3\.6\.0 || (
	  # extra work to get the right version
         echo R is installed but version is not 3.6.0, presumed old
		 echo R version is presumed old and installed via apt-get
		 echo so remove it via apt-get
		 sudo apt-get remove r-base
		 sudo apt-get remove r-base-core
            if [ -f "$myfile" ]
               then
                  /usr/bin/R --version | grep 3\.4\.2 && (
                     echo R version is 3.4.2, but could not remove via apt-get
                     exit 1
                  )
            fi
        # )
      )
fi

# at this point, R should be either missing or at version 3.4.2 if present
# source tarball is available for 3.4.2
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

#download r-base source tarball @ https://cran.r-project.org/src/base/R-3/R-3.6.0.tar.gz
myfile2=R-3.6.0.tar.gz
if [ ! -f "$myfile2" ]
then
	echo wget-ing file "$myfile2"
	wget --no-check-certificate https://cran.r-project.org/src/base/R-3/R-3.6.0.tar.gz
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
mydir=/home/ubuntu/emmes_install/r-base_install/R-3.6.0
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
./configure --with-x=no --with-tcl-config=/usr/lib/tclConfig.sh --with-tk-config=/usr/lib/tkConfig.sh
make
sudo make install

sudo ln -s /usr/local/bin/R /usr/bin/R
ls -al /usr/local/bin/R /usr/bin/R
if [ ! -f "$myfile" ]
 then
	echo /usr/bin/R still does not exist, installation deemed failed
	exit 1
 else
	/usr/bin/R --version | grep 3\.6\.0 || (
	echo R is not version 3.6.0, installation deemed failed
	)
	/usr/bin/R --version | grep 3\.6\.0 && (
	echo R has been installed from source and is at version 3.6.0
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




##install biocLite.R in R (3.9)
##the same two commands can test for a previous biocLite installation
##and install it if not already installed
##so just grep the output for "Using Bioconductor 3.9" to verify a working installation
## FOR 3.9 do
#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

	
#install.packages("BiocManager")
#BiocManager::install()


##FOR PACKAGES
#BiocManager::install(c("mixOmics", "msa", "BiocGenerics", "Biobase", "affyPLM", "qvalue", "genefilter", "sva", "edgeR", "goseq", "biomaRt", "erccdashboard", "DESeq2", "PROPER", "Gviz", "impute", "affy", "limma"))


# echo testing for bioconductor version 3.9 in R, and installing if necessary
# echo source\(\"http://bioconductor.org/biocLite.R\"\) > install_biocLite.R
# echo biocLite\(\) >> install_biocLite.R
# cat install_biocLite.R
# sudo R CMD BATCH install_biocLite.R
# cat install_biocLite.Rout
# grep "Using Bioconductor 3.9" install_biocLite.Rout || (
   # echo biocLite version 3.9 not detected, installation deemed failed
   # exit 1
# )
# grep "Using Bioconductor 3.9" install_biocLite.Rout && (
   # echo biocLite version 3.9 detected
# )
# pwd
# ls -al install_biocLite.R install_biocLite.Rout
# rm --interactive=never install_biocLite.R install_biocLite.Rout


cd /home/ubuntu/emmes_install


#if (!requireNamespace("BiocManager", quietly = TRUE))
# echo install.packages\(\"BiocManager\"\) >> install_biocLite.R
# echo BiocManager::install\(\) >> install_biocLite.R


echo checking for prior installation of R packages mixOmics, msa, BiocGenerics, Biobase, affyPLM, qvalue, genefilter, sva, edgeR, goseq, biomaRt, erccdashboard, DESeq2, PROPER, Gviz, impute, affy and limma

BiocManager::install\(c\(\"mixOmics\"\, \"msa\"\, \"BiocGenerics\"\, \"Biobase\"\, \"affyPLM\"\, \"qvalue\"\, \"genefilter\"\, \"sva\"\, \"edgeR\"\, \"goseq\"\, \"biomaRt\"\, \"erccdashboard\"\, \"DESeq2\"\, \"PROPER\"\, \"Gviz\"\, \"impute\"\, \"affy\"\, \"limma\"\)\) >> install_biocLite.R

sudo R CMD BATCH install_biocLite.R

# echo libs \= c\(\'mixOmics\'\,\'msa\'\,\'BiocGenerics\'\,\'Biobase\'\,\'affyPLM\'\,\'qvalue\'\,\'genefilter\'\,\'sva\'\,\'edgeR\'\,\'goseq\'\,\'biomaRt\'\,\'erccdashboard\'\,\'DESeq2\'\,\'PROPER\'\,\'Gviz\'\,\'impute\'\,\'affy\'\,\'limma\'\)\; > check_R_packages.R
# echo lapply\(libs\, require\, character\.only\=T\) >> check_R_packages.R
# cat  check_R_packages.R
# sudo R CMD BATCH check_R_packages.R
# cat check_R_packages.Rout
# grep FALSE  check_R_packages.Rout || (
   # echo R packages mixOmics, msa, BiocGenerics, Biobase, affyPLM, qvalue, genefilter, sva, edgeR, goseq, biomaRt, erccdashboard, DESeq2, PROPER, Gviz, impute, affy and limma already installed, so do nothing
  # rm --interactive=never check_R_packages.Rout
# )
# grep FALSE  check_R_packages.Rout && (
   # echo installing biocLite packages mixOmics, msa, BiocGenerics, Biobase, affyPLM, qvalue, genefilter, sva, edgeR, goseq, biomaRt, erccdashboard, DESeq2, PROPER, Gviz, impute, affy and limma
   # echo source\(\"http://bioconductor.org/biocLite.R\"\) > install_biocLite_packages.R
   # echo biocLite\(\'mixOmics\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'msa\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'BiocGenerics\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'Biobase\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'affyPLM\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'qvalue\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'genefilter\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'sva\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'edgeR\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'goseq\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'biomaRt\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'erccdashboard\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'DESeq2\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'PROPER\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'Gviz\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'impute\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'affy\'\) >> install_biocLite_packages.R
   # echo biocLite\(\'limma\'\) >> install_biocLite_packages.R
   # cat install_biocLite_packages.R
   # sudo R CMD BATCH install_biocLite_packages.R
   # cat install_biocLite_packages.Rout
   # egrep "DONE \(mixOmics\)" install_biocLite_packages.Rout || (
      # echo biocLite package mixOmics not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(mixOmics\)" install_biocLite_packages.Rout && (
      # echo biocLite package mixOmics INSTALLED
   # )
   # egrep "DONE \(msa\)" install_biocLite_packages.Rout || (
      # echo biocLite package msa not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(msa\)" install_biocLite_packages.Rout && (
      # echo biocLite package msa INSTALLED
   # )
   # egrep "DONE \(BiocGenerics\)" install_biocLite_packages.Rout || (
      # echo biocLite package BiocGenerics not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(BiocGenerics\)" install_biocLite_packages.Rout && (
      # echo biocLite package BiocGenerics INSTALLED
   # )
   # egrep "DONE \(Biobase\)" install_biocLite_packages.Rout || (
      # echo biocLite package Biobase not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(Biobase\)" install_biocLite_packages.Rout && (
      # echo biocLite package Biobase INSTALLED
   # )
   # egrep "DONE \(affyPLM\)" install_biocLite_packages.Rout || (
      # echo biocLite package affyPLM not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(affyPLM\)" install_biocLite_packages.Rout && (
      # echo biocLite package affyPLM INSTALLED
   # )
   # egrep "DONE \(qvalue\)" install_biocLite_packages.Rout || (
      # echo biocLite package qvalue not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(qvalue\)" install_biocLite_packages.Rout && (
      # echo biocLite package qvalue INSTALLED
   # )
   # egrep "DONE \(genefilter\)" install_biocLite_packages.Rout || (
      # echo biocLite package genefilter not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(genefilter\)" install_biocLite_packages.Rout && (
      # echo biocLite package genefilter INSTALLED
   # )
   # egrep "DONE \(sva\)" install_biocLite_packages.Rout || (
      # echo biocLite package sva not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(sva\)" install_biocLite_packages.Rout && (
      # echo biocLite package sva INSTALLED
   # )
   # egrep "DONE \(edgeR\)" install_biocLite_packages.Rout || (
      # echo biocLite package edgeR not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(edgeR\)" install_biocLite_packages.Rout && (
      # echo biocLite package edgeR INSTALLED
   # )
   # egrep "DONE \(goseq\)" install_biocLite_packages.Rout || (
      # echo biocLite package goseq not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(goseq\)" install_biocLite_packages.Rout && (
      # echo biocLite package goseq INSTALLED
   # )
   # egrep "DONE \(biomaRt\)" install_biocLite_packages.Rout || (
      # echo biocLite package biomaRt not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(biomaRt\)" install_biocLite_packages.Rout && (
      # echo biocLite package biomaRt INSTALLED
   # )
   # egrep "DONE \(erccdashboard\)" install_biocLite_packages.Rout || (
      # echo biocLite package erccdashboard not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(erccdashboard\)" install_biocLite_packages.Rout && (
      # echo biocLite package erccdashboard INSTALLED
   # )
   # egrep "DONE \(DESeq2\)" install_biocLite_packages.Rout || (
      # echo biocLite package DESeq2 not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(DESeq2\)" install_biocLite_packages.Rout && (
      # echo biocLite package DESeq2 INSTALLED
   # )
   # egrep "DONE \(PROPER\)" install_biocLite_packages.Rout || (
      # echo biocLite package PROPER not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(PROPER\)" install_biocLite_packages.Rout && (
      # echo biocLite package PROPER INSTALLED
   # )
   # egrep "DONE \(Gviz\)" install_biocLite_packages.Rout || (
      # echo biocLite package Gviz not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(Gviz\)" install_biocLite_packages.Rout && (
      # echo biocLite package INSTALLED
   # )
   # egrep "DONE \(impute\)" install_biocLite_packages.Rout || (
      # echo biocLite package impute not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(impute\)" install_biocLite_packages.Rout && (
      # echo biocLite package impute INSTALLED
   # )
   # egrep "DONE \(affy\)" install_biocLite_packages.Rout || (
      # echo biocLite package affy not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(affy\)" install_biocLite_packages.Rout && (
      # echo biocLite package affy INSTALLED
   # )
   # egrep "DONE \(limma\)" install_biocLite_packages.Rout || (
      # echo biocLite package limma not installed, installation deemed failed
      # exit 1
   # )
   # egrep "DONE \(limma\)" install_biocLite_packages.Rout && (
      # echo biocLite package limma INSTALLED
   # )
   # rm --interactive=never install_biocLite_packages.R install_biocLite_packages.Rout
)


# install some more packages (in R) if they are not already installed the R install.packages way or using devtools
# Script loop to check and install pacakges
# library(package)
# if error, install.packages()
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
	wget --no-check-certificate https://cran.r-project.org/src/contrib/openxlsx_4.1.0.tar.gz
	#     tar -zxvf openxlsx_3.0.0.tar.gz #no extraction needed
	#     create R script to install package
	echo pkgFile = \"/home/ubuntu/emmes_install/openxlsx_install/openxlsx_4.1.0.tar.gz\" > install_openxlsx.R
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
	rm --interactive=never openxlsx_4.1.0.tar.gz install_openxlsx.R
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
	echo install_version\(\"Cairo\"\, version \= \"1\.5\-10\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_Cairo.R
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
	echo install_version\(\"seqinr\"\, version \= \"3\.4\-5\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_seqinr.R
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
	echo install_version\(\"iterpc\"\, version \= \"0\.4\.1\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_iterpc.R
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
	echo install_version\(\"gridExtra\"\, version \= \"2\.3\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_gridExtra.R
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


# install iGraph for R
myfile=/usr/local/lib/R/library/iGraph/DESCRIPTION ##change
if [ -f "$myfile" ]
then
	echo R package iGraph is already installed, so do nothing
else
	echo installing package iGraph in R using devtools
	echo require\(devtools\) > install_iGraph.R
	echo install_version\(\"igraph\"\, version \= \"1\.2\.4\.1\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_iGraph.R
	cat install_iGraph.R
	sudo R CMD BATCH install_iGraph.R
	cat install_iGraph.Rout
	rm --interactive=never install_iGraph.R install_iGraph.Rout
fi
echo R package iGraph INSTALLED

# install glmnet in R
myfile=/usr/local/lib/R/library/glmnet/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package glmnet is already installed, so do nothing
else
	echo installing package glmnet in R using devtools
	echo require\(devtools\) > install_glmnet.R
	echo install_version\(\"glmnet\"\, version \= \"2\.0\-16\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_glmnet.R
	cat install_glmnet.R
	sudo R CMD BATCH install_glmnet.R
	cat install_glmnet.Rout
	rm --interactive=never install_glmnet.R install_glmnet.Rout
fi
echo R package glmnet INSTALLED

# install fpca for R (needs xvfb to bypass NO DISPLAY error)
myfile=/usr/local/lib/R/library/fpca/DESCRIPTION ##change
if [ -f "$myfile" ]
then
	echo R package fpca is already installed, so do nothing
else
	echo installing package fpca in R using devtools
	echo require\(devtools\) > install_fpca.R
	echo install_version\(\"fpca\"\, version \= \"0\.2\-1\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_fpca.R
	cat install_fpca.R
	sudo R CMD BATCH install_fpca.R
	cat install_fpca.Rout
	rm --interactive=never install_fpca.R install_fpca.Rout
fi
echo R package fpca INSTALLED

# install smacof for R
myfile=/usr/local/lib/R/library/smacof/DESCRIPTION ##change
if [ -f "$myfile" ]
then
	echo R package smacof is already installed, so do nothing
else
	echo installing package smacof in R using devtools
	echo require\(devtools\) > install_smacof.R
	echo install_version\(\"smacof\"\, version \= \"1\.10\-8\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_smacof.R
	cat install_smacof.R
	sudo R CMD BATCH install_smacof.R
	cat install_smacof.Rout
	rm --interactive=never install_smacof.R install_smacof.Rout
fi
echo R package smacof INSTALLED

# install MBA in R
myfile=/usr/local/lib/R/library/MBA/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package MBA is already installed, so do nothing
else
	echo installing package MBA in R
	echo install.packages\(\"MBA\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_MBA.R
	cat install_MBA.R
	sudo R CMD BATCH install_MBA.R
	cat install_MBA.Rout
	rm --interactive=never install_MBA.R install_MBA.Rout
fi
echo R package MBA INSTALLED

# install reshape2 in R
myfile=/usr/local/lib/R/library/reshape2/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package reshape2 is already installed, so do nothing
else
	echo installing package reshape2 in R
	echo install.packages\(\"reshape2\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_reshape2.R
	cat install_reshape2.R
	sudo R CMD BATCH install_reshape2.R
	cat install_reshape2.Rout
	rm --interactive=never install_reshape2.R install_reshape2.Rout
fi
echo R package reshape2 INSTALLED

# install ggplot2 in R
myfile=/usr/local/lib/R/library/ggplot2/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package ggplot2 is already installed, so do nothing
else
	echo installing package ggplot2 in R
	echo install.packages\(\"ggplot2\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_ggplot2.R
	cat install_ggplot2.R
	sudo R CMD BATCH install_ggplot2.R
	cat install_ggplot2.Rout
	rm --interactive=never install_ggplot2.R install_ggplot2.Rout
fi
echo R package ggplot2 INSTALLED

# install gtools in R
myfile=/usr/local/lib/R/library/gtools/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package gtools is already installed, so do nothing
else
	echo installing package gtools in R
	echo install.packages\(\"gtools\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_gtools.R
	cat install_gtools.R
	sudo R CMD BATCH install_gtools.R
	cat install_gtools.Rout
	rm --interactive=never install_gtools.R install_gtools.Rout
fi
echo R package gtools INSTALLED

# install stringr in R
myfile=/usr/local/lib/R/library/stringr/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package stringr is already installed, so do nothing
else
	echo installing package stringr in R
	echo install.packages\(\"stringr\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_stringr.R
	cat install_stringr.R
	sudo R CMD BATCH install_stringr.R
	cat install_stringr.Rout
	rm --interactive=never install_stringr.R install_stringr.Rout
fi
echo R package stringr INSTALLED

# install dplyr for R
myfile=/usr/local/lib/R/library/dplyr/DESCRIPTION ##change
if [ -f "$myfile" ]
then
	echo R package dplyr is already installed, so do nothing
else
	echo installing package dplyr in R using devtools
	echo require\(devtools\) > install_dplyr.R
	echo install_version\(\"dplyr\"\, version \= \"0\.8\.0\.1\"\, repos \= \"http\://cran.us.r\-project.org\"\, quiet \= F\) >> install_dplyr.R
	cat install_dplyr.R
	sudo R CMD BATCH install_dplyr.R
	cat install_dplyr.Rout
	rm --interactive=never install_dplyr.R install_dplyr.Rout
fi
echo R package dplyr INSTALLED

# install plot3D in R
myfile=/usr/local/lib/R/library/plot3D/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package plot3D is already installed, so do nothing
else
	echo installing package plot3D in R
	echo install.packages\(\"plot3D\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_plot3D.R
	cat install_plot3D.R
	sudo R CMD BATCH install_plot3D.R
	cat install_plot3D.Rout
	rm --interactive=never install_plot3D.R install_plot3D.Rout
fi
echo R package plot3D INSTALLED

# install doParallel in R
myfile=/usr/local/lib/R/library/doParallel/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package doParallel is already installed, so do nothing
else
	echo installing package doParallel in R
	echo install.packages\(\"doParallel\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_doParallel.R
	cat install_doParallel.R
	sudo R CMD BATCH install_doParallel.R
	cat install_doParallel.Rout
	rm --interactive=never install_doParallel.R install_doParallel.Rout
fi
echo R package doParallel INSTALLED

# install foreach in R
myfile=/usr/local/lib/R/library/foreach/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package foreach is already installed, so do nothing
else
	echo installing package foreach in R
	echo install.packages\(\"foreach\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_foreach.R
	cat install_foreach.R
	sudo R CMD BATCH install_foreach.R
	cat install_foreach.Rout
	rm --interactive=never install_foreach.R install_foreach.Rout
fi
echo R package foreach INSTALLED

# install fdapace in R
myfile=/usr/local/lib/R/library/fdapace/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package fdapace is already installed, so do nothing
else
	echo installing package fdapace in R
	echo install.packages\(\"fdapace\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_fdapace.R
	cat install_fdapace.R
	sudo R CMD BATCH install_fdapace.R
	cat install_fdapace.Rout
	rm --interactive=never install_fdapace.R install_fdapace.Rout
fi
echo R package fdapace INSTALLED

# install ape in R
myfile=/usr/local/lib/R/library/ape/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package ape is already installed, so do nothing
else
	echo installing package ape in R
	echo install.packages\(\"ape\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_ape.R
	cat install_ape.R
	sudo R CMD BATCH install_ape.R
	cat install_ape.Rout
	rm --interactive=never install_ape.R install_ape.Rout
fi
echo R package ape INSTALLED

# install png in R
myfile=/usr/local/lib/R/library/png/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package png is already installed, so do nothing
else
	echo installing package png in R
	echo install.packages\(\"png\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_png.R
	cat install_png.R
	sudo R CMD BATCH install_png.R
	cat install_png.Rout
	rm --interactive=never install_png.R install_png.Rout
fi
echo R package png INSTALLED

#install pdfrop and other utilities
sudo apt-get install texlive-extra-utils
sudo apt-get install poppler-utils

##CLEANUP

#Apt get
sudo apt-get clean
sudo apt-get autoclean
sudo apt autoremove


#Remove swap
sudo /sbin/swapoff /var/swap.1
sudo rm -f /var/swap.1

## may need to remove old distribute package..
#sudo rm -fr /usr/local/lib/python2.7/dist-packages/distribute*

#Delete emmes_install to gain back space on server?


echo Script Completed.
##end of script 
exit 0
