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
# Program:  sqlite.py
# Version:  RSEQREP 2.0.0
# Author:   William F Hooper, Travis L. Jensen, Johannes B. Goll
# Purpose:  sqlite3 database interface
# Input:    N/A
# Output:   N/A
#############################################################################################################

import sqlite3
import utils
import os
import logging

## Initialize a SQLite3 database if it doesn't already exist 
def initSqliteDb(db):
	
	## Open the connection, get the cursor
	conn = sqlite3.connect(db)
	c = conn.cursor()
	
	## Create process table
	c.execute('''CREATE TABLE IF NOT EXISTS process(process_id integer primary key autoincrement, process_str text, return_code, start_time integer,
	end_time integer,wc_time integer, sample_id integer,created timestamp default (datetime('now','localtime')))''')
	
	## Create file table
	c.execute('''CREATE TABLE IF NOT EXISTS file(file_id integer primary key autoincrement, file_name text, file_md5 text, file_bytes integer, file_type text, 
	process_id integer, sample_id integer, file_path text, created timestamp default (datetime('now','localtime')))''')
	
	## Create sample table
	c.execute('''CREATE TABLE IF NOT EXISTS sample(sample_id integer primary key autoincrement, sample_name text unique,
	created timestamp default (datetime('now','localtime')))''');
	
	## Create benchmark table
	c.execute('''CREATE TABLE IF NOT EXISTS benchmark(benchmark_id integer primary key autoincrement, process_id integer, sample_id integer, virtual integer, 
	resident_set_size integer, cpu_time integer, wc_time integer, return_code integer, created timestamp default (datetime('now','localtime')))''')
	
	## Commit and close the connection
	conn.commit()
	conn.close()

	

## Add a sample to the sample table
def addSample(db, sample_name):
	
	## Open the connection, get the cursor
	conn = sqlite3.connect(db, timeout=1000)
	c = conn.cursor()
	
	## Try to insert into table
	c.execute("INSERT OR IGNORE INTO sample (sample_name) VALUES (?)", (sample_name, ))
	
	## Commit and close the connection
	conn.commit()
	conn.close()
	


## Pull the database sample ID for a given sample name
def getSamid(db, sample_name):

	## Open the connection, get the cursor
	conn = sqlite3.connect(db, timeout=1000)
	c = conn.cursor()
	
	## Pull samid
	c.execute("SELECT sample_id FROM sample WHERE sample_name=(?)", (sample_name, ))
	samid = c.fetchone()[0]
	
	## Close the connection & return
	conn.close()
	return(samid)
	


## Add a process to the process table
def addProcess(db, samid, cmd, return_code, start_time, end_time):
	
	## Compute wall clock time
	wc = end_time - start_time
	
	## Open the connection, get the cursor
	conn = sqlite3.connect(db, timeout=1000)
	c = conn.cursor()

	## Add process to table
	c.execute('''INSERT INTO process (sample_id, process_str, return_code, start_time, end_time, wc_time) VALUES (?, ?, ?, ?, ?, ?)''', (samid, cmd, return_code, start_time, end_time, wc))

	## Return process ID
	c.execute("SELECT MAX(process_id) from PROCESS")
	procid = c.fetchone()[0]
	
	## Commit, close the connection, return process id
	conn.commit()
	conn.close()
	return(procid)



## Add file(s) to the file table
def addFile(db, samid, procid, file, sample_name, file_type):
	
	## Make sure file is a list
	file = utils.toList(file) 
	
	for f in file:
	
		## Compute md5sum, file size in bytes
		file_md5 = utils.md5(f)
		file_size = os.path.getsize(f)
		
		## Open the connection, get the cursor
		conn = sqlite3.connect(db, timeout=1000)
		c = conn.cursor()
			
		## Add file to table
		c.execute('''INSERT INTO file (sample_id, process_id, file_path, file_name, file_md5, file_bytes, file_type) VALUES (?, ?, ?, ?, ?, ?, ?)''', (samid, procid, f, sample_name, file_md5, file_size, file_type))
		
		## Commit and close the connection
		conn.commit()
		conn.close()

	
	
	

## Parse a single process benchmark into a dictionary
def parseBenchmark(bench_file):
	with open(bench_file) as f:
		keys = f.readline().split()
		values = f.readline().split()
	return(dict(zip(keys, values)))
	


## Add a process benchmark to the benchmark table
def addBenchmark(db, procid, samid, bench_obj):
	
	## Open the connection, get the cursor
	conn = sqlite3.connect(db, timeout=1000)
	c = conn.cursor()

	## Add file to table
	c.execute('''INSERT INTO benchmark (process_id, sample_id, virtual, resident_set_size, cpu_time, wc_time, return_code) VALUES (?, ?, ?, ?, ?, ?, ?)''', (procid, samid, bench_obj.max_vms, bench_obj.max_rss, bench_obj.cpu_seconds, bench_obj.running_time, bench_obj.return_code))

	## Commit and close the connection
	conn.commit()
	conn.close()



## Add all flat-file benchmarks 
def addFolderBenchmarks(db, dir):
	for file in os.listdir(dir):
		if file.endswith('.info'):
			bench_info = parseBenchmark(dir+'/'+file)
			addBenchmark(db=db, procid=bench_info['procid'], samid=bench_info['samid'], bench_file=bench_info['bench_file'], return_code=bench_info['return_code'])
			os.remove(dir+'/'+file)
