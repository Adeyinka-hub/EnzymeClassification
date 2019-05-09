#!/usr/bin/python3.6
# Remmy Chen 04/13/2019

import requests 
import os
import numpy as np
import pandas as pd

RAW_DATA_DIR = 'raw_data'
PROCESSED_DATA_DIR = 'processed_data'


"""
Creates a directory in the current working directory 
if it doesn't exist.

Argument:
	-	dir_name : string that is the name of the 
		directory that the client wants to create.

Return:
	-	entry_filepath : string that is the 
		absolute path to entry_data directory.
"""
def create_entry_dir(dir_name):
	current_dir = os.getcwd()
	entry_filepath = os.path.join(current_dir, dir_name)
	if not os.path.exists(entry_filepath):
		os.makedirs(entry_filepath)
	return entry_filepath

"""
Checks if a file name exists in a directory.

Arguments:
	- 	file_name : string that is the file name. 
		Assumes that the string will be a file name.
	- 	entry_filepath : string that is the absolute 
		path of directory to search in

Return:
	-	EXISTS : bool that is True if file exists and 
		False otherwise.
"""
def check_file_exists(file_name, entry_filepath):
	EXISTS = True
	files = list(os.listdir(entry_filepath)) # lists files AND dirs
	if file_name in files:
		return EXISTS
	return not EXISTS

"""
Writes entry text to file, creating the file if
it doesn't exist.

Arguments:
	- 	file_name : string that is the file name.
	- 	text : string that is the content to write.
"""
def write_to_file(file_name, text):
	with open(file_name, 'w+') as f:
		f.write('%s' % text)
		f.close()

"""
Checks if requested data is stored locally. If not, performs a 
GET request to the KEGG database and stores the list information
in a file. Returns the file path and name. 

Argument:
	-	database : string that is assumed to be an 
		abbreviation of a KEGG database.

Return:
	- 	FILEPATH_AND_NAME : string that is the absolute 
		file path and name that stores the KEGG data for 
		the requested information.

http://rest.kegg.jp/list/ko
"""
def list_KEGG_data(database):
	database = database.lower()
	URL_LIST = 'http://rest.kegg.jp/list/' + database
	FILE_NAME = database + '_list.txt'
	ENTRY_FILEPATH = create_entry_dir(RAW_DATA_DIR)
	FILEPATH_AND_NAME = ENTRY_FILEPATH + '/' + FILE_NAME
	if check_file_exists(FILE_NAME, ENTRY_FILEPATH) != True:
		r = requests.get(url = URL_LIST)
		if r.status_code != 200:
			print("Bad request")
			return ''
		write_to_file(FILEPATH_AND_NAME, r.text)
	return FILEPATH_AND_NAME	

"""
Parses the list file created in list_KEGG_data().

Argument:
	- 	FILEPATH_AND_NAME : string that is the absolute 
		file path and name that stores the KEGG data for 
		the requested information.

Return:
	- 	ko_num_list : a list of KO numbers in the file.
	- 	ko_ec_list : a list of EC numbers associated with KO
	    numbers.
	- 	ko_ec_bin_list: a list of ints where 1 means KO number
	    has an associated EC number and 0 otherwise.
"""
def parse_KO_list(FILEPATH_AND_NAME):
	ko_num_list = []
	ko_ec_list = []
	ko_ec_bin_list = []
	with open(FILEPATH_AND_NAME, 'r') as f:
		for line in f:
			ko_parts = line.split()
			ko_num_list.append(ko_parts[0])
			if '[' in line:
				parts = line.split('[')
				ko_ec_list.append(parts[1][:-2])
				ko_ec_bin_list.append(1)
			else:
				ko_ec_list.append("")
				ko_ec_bin_list.append(0)
	return ko_num_list, ko_ec_list, ko_ec_bin_list

"""
Checks if requested entry is stored locally. If not, performs 
a GET request to the KEGG database and stores the KEGG entry 
in a file. Returns the file path and name. 

Arguments:
	- 	database : string that is assumed to be an 
		abbreviation of a KEGG database.
	-	entry : string that is assumed to be a 
		compound, enzyme, RClass, etc that exists
		in the KEGG database. Assumes alphabetical
		letters to be in upper case.

Return:
	-	FILEPATH_AND_NAME : string that is the absolute 
		file path and name that stores the KEGG data for 
		the requested entry.

ko http://rest.kegg.jp/get/ko:K00001
hsa http://rest.kegg.jp/get/hsa:7108

Note:
	- 	Seems like for organism, only hsa is currently supported.
"""
def get_KEGG_data(database, entry):
	database = database.lower()
	FILE_NAME = database + '_' + entry + '.txt'
	URL_GET = 'http://rest.kegg.jp/get/' + database + ':' + entry
	ENTRY_FILEPATH = create_entry_dir(RAW_DATA_DIR)
	FILEPATH_AND_NAME = ENTRY_FILEPATH + '/' + FILE_NAME	
	if check_file_exists(FILE_NAME, ENTRY_FILEPATH) != True:
		#print("making a GET request to : %s" % URL_GET)
		r = requests.get(url = URL_GET)
		#print("GET request response code : %s %s" % (r.status_code, r.reason))
		if r.status_code != 200:
			print("Bad request")
			return ''
		write_to_file(FILEPATH_AND_NAME, r.text)
		#print("GET request response saved to file : %s" % FILE_NAME)
	#else:
		#print("request %s already stored in : %s" % (entry, FILEPATH_AND_NAME))
	return FILEPATH_AND_NAME


"""
Searches file of KO number for the genes section, looks
for associated gene of specified organism, and retrieves
amino acid sequence of associated gene.

Arguments:
	- 	organism : string that is assumed to be an 
		abbreviation of an organism in the KEGG database.
	-	ko_num : string that is assumed to be a 
		KO number that exists in the KEGG database. 
		Assumes alphabetical letters to be in upper case.

Return:
	-	ko_hsa : 1 if found gene of specified organism
		associated with a specified KO number.
	-	aaseqs : list of strings representing amino acid
		sequences of a specified organism associated with
		a specified KO number.
"""
def get_ko_data(organism, ko_num):
	organism = organism.upper()
	ENTRY_FILEPATH = create_entry_dir(RAW_DATA_DIR)
	FILE_NAME = 'ko_' + ko_num + '.txt'
	FILEPATH_AND_NAME = ENTRY_FILEPATH + '/' + FILE_NAME	
	ko_hsa = 0
	aaseqs = []
	with open(FILEPATH_AND_NAME, 'r') as f:
		for line in f:
			if line[:3] == "///":
				continue
			if line[:12] != "            ":
				keyword = line[:12]
			data = line[12:].strip()
			if keyword == "GENES       ":
				words = data.split()
				if organism.upper() in words[0]:
					ko_hsa = 1
					for word in words[1:]: # can have more than just one gene
						gene = word.split('(')[0]
						FILEPATH_AND_NAME = get_KEGG_data(organism, gene)
						aaseqs.append(get_aa_seq(FILEPATH_AND_NAME))
	return ko_hsa, aaseqs

"""
Parses the gene file created in get_KEGG_data().

Argument:
	- 	FILEPATH_AND_NAME : string that is the absolute 
		file path and name that stores the KEGG data for 
		the requested information.

Return:
	- 	aaseq : string that is the AA sequence found in file.
"""
def get_aa_seq(FILEPATH_AND_NAME):	
	aaseq = ""
	with open(FILEPATH_AND_NAME, 'r') as f:
		for line in f:
			if line[:3] == "///":
				continue
			if line[:12] != "            ":
				keyword = line[:12]
			data = line[12:].strip()
			if keyword == "AASEQ       ":
				aaseq += data
	return aaseq.lstrip('0123456789') # remove the line indicating AA sequence length

"""
Downloads all KO numbers in KEGG database, saves KO numbers and associated
AA sequences in humans (HSA) into ./processed_data/x.csv, saves binary value
indicating presence or absence of EC number into ./processed_data/y.csv.
"""
def run():
	FILEPATH_AND_NAME = list_KEGG_data('ko')
	ko_num_list, ko_ec_list, ko_ec_bin_list = parse_KO_list(FILEPATH_AND_NAME)
	ko_hsa_list = []
	aa_hsa_list = []
	ec_bin_hsa_list = []
	for ko, ec_bin in zip(ko_num_list, ko_ec_bin_list):
		ko_parts = ko.split(':')
		ko_db = ko_parts[0]
		ko_num = ko_parts[1]
		ko_hsa, aaseqs = get_ko_data('HSA', ko_num)
		if ko_hsa == 1: # if hsa found
			for aaseq in aaseqs:
				if aaseq != "": # sometimes gene page only has nucleotide sequence and not amino acid sequence
					ko_hsa_list.append(ko_num)
					aa_hsa_list.append(aaseq)
					ec_bin_hsa_list.append(ec_bin)
	ENTRY_FILEPATH = create_entry_dir(PROCESSED_DATA_DIR)
	df = pd.DataFrame({"KO_num" : np.array(ko_hsa_list), "AA_seq" : np.array(aa_hsa_list)})
	df.to_csv(ENTRY_FILEPATH + "/x.csv", index=False)
	df = pd.DataFrame({"has_EC_num" : np.array(ec_bin_hsa_list)})
	df.to_csv(ENTRY_FILEPATH + "/y.csv", index=False)

if __name__ == "__main__":
	run()



