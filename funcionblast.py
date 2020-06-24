#!/usr/bin/env python
# -- coding: utf-8 --

"""

@author: Marina Ortega Angulo
"""


import os
import pandas as pd
import sys


from Bio import SeqIO

from subprocess import Popen, PIPE




def Parsear(genbank):
    directorio = os.getcwd()
    try:
        if os.path.isfile("gb_parsed.fasta") == True:
            os.remove("gb_parsed.fasta")
        else:
            pass

        output_f = open("gb_parsed.fasta","a")
        os.chdir(genbank)
        os.listdir()

        for file in os.listdir():
            with open(file, "r") as input_f:

                for seq_record in SeqIO.parse(input_f, "genbank"):

                    for seq_feature in seq_record.features:
                        try:

                            #Finding the elements that we want inside the genbank files


                            if seq_feature.type == 'CDS':
                                output_f.write("> %s@%s\n%s\n" % (
                                    seq_feature.qualifiers['locus_tag'][0],
                                    seq_record.name,
                                    seq_feature.qualifiers['translation'][0])
                                )
                        except:
                            pass

            input_f.close()
        output_f.close()

        print("Parse ✓")

        os.chdir(directorio)
    except:
        print("The input does not have genbank files. Please try again.")
        sys.exit()


def BLASTP(usr_query, subject = 'gb_parsed.fasta'):

    #With the try element I see if the file introduced is fasta or not


    try:
        with open(usr_query,"r") as file:

            for record in SeqIO.parse(file, "fasta"):
                seq = str(">%s\n%s\n" % (record.id, record.seq))
                temp = open(record.id+"temp.fasta", "w")
                temp.write(seq)
                temp.close()

                with open(record.id+"_bpresults.tsv", "w") as blastp_r:

                    blastp = Popen(['blastp', '-query', record.id+"temp.fasta", '-subject',
                                    subject, '-evalue', '0.00001', '-outfmt', "6 qseqid qcovs pident evalue sseqid"
                                                                              " sseq"], stdout=PIPE, stderr=PIPE)
                    header_blastp = str("queryID\tCOVERAGE\tIDENTITY\tEVALUE\t"
                                        "subjectID\tSUBJECTseq\n")
                    result_blastp = blastp.stdout.read().decode("utf-8")
                    blastp_r.write(header_blastp)
                    blastp_r.write(result_blastp)
                    blastp_r.close()

                    os.remove(record.id + "temp.fasta")

            print("Blastp ✓")
    except:
        print("ERROR: the input file is not in fasta format. Try again.")

def FILTER(usr_query, ID_usr, COV_usr):


    #IDENTITY
    try:
        if int(sys.argv[3])>=100:
            print("ERROR: identity has to be in between 0 and 100. Try again.")
            sys.exit()

        elif int(sys.argv[3])>=0 or int(sys.argv[3])<=100:
            ID_usr = sys.argv[3]
            pass
    except:
        ID_usr = 30.0

    #COVERAGE


    try:
        if int(sys.argv[4])>=100:
            print("ERROR: coverage has to be in between 0 and 100. Try again.")
            sys.exit()

        elif int(sys.argv[4])>=0 or int(sys.argv[4])<=100:
            COV_usr = sys.argv[4]
            pass
    except:
        COV_usr = 50.0

    for record in SeqIO.parse(usr_query, "fasta"):

        with open(record.id+"_bpresults.tsv") as tsvfile, \
                open(record.id+"_filtred.tsv", "w") as tsv_filtred:

            tsvreader=pd.read_csv(tsvfile, delimiter='\t')
            trying = tsvreader.loc[(tsvreader['IDENTITY']>=int(ID_usr)) &
                                   (tsvreader['COVERAGE']>=int(COV_usr)), :]
            trying.to_csv(tsv_filtred, sep='\t')
            tsvfile.close()
            tsv_filtred.close()

    print("Filtered ✓")

