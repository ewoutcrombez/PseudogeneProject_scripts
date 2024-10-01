# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 13:43:00 2020
@author: Ewout
"""
import sys

def converttoDnaSeg(iadhore_input=sys.argv[1], species=sys.argv[2], chromosome=sys.argv[3], output=sys.argv[4]):
	fh = open("%s/genelists/%s/%s"%(iadhore_input, species, chromosome), "r")
	string = 'name\tstart\tend\tstrand\tcol\tfill\n'
	coordinate1 = 0
	coordinate2 = 150
	col = 'red'
	fill = 'red'
	for line in fh:
		if line == '\n':
			next
		else:
			string += line[:-2] + '\t' + str(coordinate1) + '\t' + str(coordinate2) + '\t' + line[-2] + '\t' + col + '\t' + fill + '\n'
			coordinate1 += 400
			coordinate2 += 400
	fhh = open("%s/%s_%s.tab.txt"%(output, chromosome, species), "w")
	fhh.write(string)
	fhh.close()
	fh.close()

converttoDnaSeg()
