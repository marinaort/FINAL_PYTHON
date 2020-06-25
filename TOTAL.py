#!/usr/bin/env python
# -- coding: utf-8 --

"""

@author: Marina Ortega Angulo
"""


import os
import shutil
import sys

from Bio import SeqIO

import funcionblast
import ARBOL
import prosite

def Help():

    class color:
        PURPLE = '\033[95m'
        DARKCYAN = '\033[36m'
        GREEN = '\033[92m'
        YELLOW = '\033[93m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
        END = '\033[0m'

    print("----------------------------")
    print(color.BOLD+"\tSCRIPT USAGE"+color.END)
    print("----------------------------\n\n"
          "What is this script for?\n\n"
          "With this script you can parse genbanks (1 or more) to find the secuence that "
          "we will use as subject for the blastp.After the parse the script will use the "
          "file introduced as query to do the blastp against the genbank_parsed.fasta.")
    print(color.UNDERLINE+"* The query file can have one or more sequence, each "
                          "result file will be named with the name of each query sequence that the user "
                          "have introduced."+color.END)
    print("\n\nThe blastp results will be filtered with the coverage and identity "
          "that the user enters. This filtered file will be converted to a .fasta so we "
          "can use it for the muscle. With muscle we will get a aligned file so we can "
          "make a Neighbor-Joining (N-J) tree.\n\n"
          "The last thing that this script does is to find protein domains in the "
          "Prosite database. For each protein it will show the domain or domains, the "
          "accession, the description and the pattern found.\n\n"
          "What do you need to introduce?\n")
    print("./TOTAL.py"+color.DARKCYAN+"[genbank]"+color.YELLOW+"[query.fasta]"
          +color.GREEN+"[identity]"+color.PURPLE+"[coverage]"+color.END)
    print(color.DARKCYAN+"- Genbank"+color.END+": introduce or a file or a "
                                               "folder with the Genbank files.")
    print(color.YELLOW+"- Query.fasta"+color.END+": with all the querys that "
                                                 "you want to compare with the genbank files.")
    print(color.BOLD+"OPTIONAL:"+color.END)
    print(color.GREEN+"- Identity"+color.END+": it is better to use "
                                             "one greater than 30.0")
    print(color.PURPLE+"- Coverage"+color.END+": it is better to use "
                                              "one greater than 50.0")
    print(color.UNDERLINE+"* If you do not enter any of these two, the script"
                          " will automatically use 30.0 as identity and 50.0 as coverage."+color.END)
    print("\n\n")
    print(color.BOLD+"all the files will be saved in a folder called with the"
                     " name that you want to introduce and _results."+color.END)
    print("\n\nThis files must to be in the same folder as the script "
          "and its modules:\n"
          "- query.fasta\n"
          "- prosite.dat\n"
          "- The folder with the genbank files.\n\n")
    print("\n\nTHANKS FOR USING THIS SCRIPT HELP.\n\n")

#Arguments control


if len(sys.argv) <= 1:

     print("ERROR: wrong arguments, you have to introduce the script name and other"
           " two at least. Please try again (getting out of the script by pressing"
           " [Ctrl+Z]) or press any key:")
     help = input()

     if help == '-h' or '-help':
         Help()
         sys.exit()
     else:
         sys.exit()

else:
    pass

#Create a folder to save all the resulting files


print ("How do you want to name the results folder?")
usr_ip = input()
os.mkdir(str(usr_ip)+"_Results")

#First argument - genbank


genbank = sys.argv[1]


#Function for parsing the genbanks


funcionblast.Parsear(genbank)


#Second argument - query file


usr_query = sys.argv[2]


#Function for doing the blastp


funcionblast.BLASTP(usr_query)


#Function for filtering the blasp result


funcionblast.FILTER(usr_query, 'ID_usr', 'COV_usr')


#Convert a .tsv file to .fasta


ARBOL.TSVtoFASTA(usr_query)


#Function to do the alignment


ARBOL.ALIGN(usr_query, file='record.id+"_muscle.fasta"')


#Function to do the tree


ARBOL.Maketree(usr_query, input_m='record.id+"_align.fasta"')


# Functions for parsing Prosite and finding the domains


prosite.Parsing_prosite(handle="prosite.dat")
prosite.Domains(usr_query, parsed_out="prosite_parsed.tsv")


#Copy files into the results folder.


for record in SeqIO.parse(usr_query, "fasta"):
    shutil.move(record.id+"_bpresults.tsv", str(usr_ip)+"_Results"+"/")
    shutil.move(record.id+"_filtred.tsv", str(usr_ip)+"_Results"+"/")
    shutil.move(record.id+"_premuscle.fasta", str(usr_ip)+"_Results"+"/")
    shutil.move(record.id+"_align.fasta", str(usr_ip)+"_Results"+"/")
    shutil.move(record.id+"_tree.nw", str(usr_ip)+"_Results"+"/")
    shutil.move(record.id + "_domains.txt", str(usr_ip)+"_Results"+"/")

shutil.move("gb_parsed.fasta", str(usr_ip)+"_Results"+"/")
shutil.move("prosite_parsed.tsv", str(usr_ip)+"_Results"+"/")

