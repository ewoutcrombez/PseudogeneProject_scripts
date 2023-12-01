# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 13:43:00 2020
@author: Ewout
"""
import sys
import re

def converttoDnaSeg(iadhore_output=sys.argv[1], species=sys.argv[2], chromosome=sys.argv[3]):
	fh = open("%s/genome_lists/%s"%(iadhore_output, chromosome), "r")
	string = 'name\tstart\tend\tstrand\tcol\tfill\n'
	coordinate1 = 0
	coordinate2 = 150
	for line in fh:
		if line == '\n':
			next
		else:
			if re.search("Pseudogene", line):
				col = 'blue'
				fill = 'blue'
			else:
				col = 'red'
				fill = 'red'
			string += line[:-2] + '\t' + str(coordinate1) + '\t' + str(coordinate2) + '\t' + line[-2] + '\t' + col + '\t' + fill + '\n'
			coordinate1 += 400
			coordinate2 += 400
	fhh = open("%s_%s.tab.txt"%(chromosome, species), "w")
	fhh.write(string)

converttoDnaSeg()
