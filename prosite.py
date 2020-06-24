#!/usr/bin/env python
# -- coding: utf-8 --

"""

@author: Marina Ortega Angulo
"""


import re

from Bio.ExPASy import Prosite,Prodoc
from Bio import SeqIO


def Parsing_prosite(handle="prosite.dat"):

	with open("prosite.dat", "r") as handle, \
			open("prosite_parsed.tsv", "w") as parsed_out:

		parsed_out.write("NAME\tACCESSION\tDESCRIPTION\tPATTERN\n")
		records = Prosite.parse(handle)

		#Finding each element


		for record in records:
			parsed_out.write("%s\t%s\t%s\t%s\n" % (
				record.name,
				record.accession,
				record.description,
				record.pattern)
							)

		handle.close()
		records.close()
		parsed_out.close()


def convertTOregex(parsed_out="prosite_parsed.tsv"):


	#Replacing the elements that can not be read by regex format.


	toregex = parsed_out.replace(".", "").replace("x", ".").replace("-", "")\
		.replace("{", "[^").replace("}", "]").replace("(", "{").replace(")", "}")\
		.replace("<", "^").replace(">", "$")
	return (toregex)

def Domains(usr_query, parsed_out="prosite_parsed.tsv"):

	for record in SeqIO.parse(usr_query, "fasta"):

		with open(record.id + "_premuscle.fasta", "r") as inputf, \
				open("prosite_parsed.tsv","r") as parsed_out,\
				open(record.id + "_domains.txt", "w") as final:

			splited_i = inputf.read().split("\n")
			splited_p = parsed_out.read().split("\n")

			# Read the file 2 by 2 lines


			for i in range(len(splited_i) // 2):
				SEQ = splited_i[2 * i + 1]
				ID = splited_i[2 * i]
				subID = str("Patterns in sequence: %s\n\n" % (ID[1:]))
				final.write(subID)

				#Parse the patterns without the header


				for pattern in splited_p[1:]:

					try:
						pattern = pattern.split("\t")

						#converting the pattron column to a regex lenguage


						patron_re = convertTOregex(pattern[3])

						if pattern[3] != "" and re.search(patron_re,SEQ):
							final.write("\tName: %s\n\tAccession: %s\n\tDescription: %s\n\t"
										"Pattern: %s\n\n" % (
							pattern[0], pattern[1], pattern[2], pattern[3]))

					except:

						#Not counting the lines that are in blank


						pattern = "NA"

			final.close()

	print("Domains ✓")