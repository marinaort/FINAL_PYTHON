#!/usr/bin/env python
# -- coding: utf-8 --

"""

@author: Marina Ortega Angulo
"""


from Bio import SeqIO

from subprocess import Popen, PIPE


def TSVtoFASTA(usr_query):

    for record in SeqIO.parse(usr_query, "fasta"):
        seq = str(">%s\n%s\n" % (record.id, record.seq))
        test = open(record.id+"_filtred.tsv", "r")
        muscle_out=open(record.id+"_premuscle.fasta","w")
        muscle_out.write(seq)


        with open(record.id+"_premuscle.fasta","a+") as muscle_out:

            for row in test.readlines():
                fields = row.rstrip().split("\t")

                #This if lets me to read the file and printing it without the header.


                if fields[6] != 'SUBJECTseq':
                    subject_ID = fields[5]
                    SUBJECTseq = fields[6]
                    muscle_out.write(">%s\n%s\n" % (subject_ID, SUBJECTseq))

            muscle_out.close()


def ALIGN(usr_query, file='record.id+"_premuscle.fasta"'):

    for record in SeqIO.parse(usr_query, "fasta"):

        with open(record.id+"_premuscle.fasta","r") as file, \
                open(record.id+"_align.fasta", "a") as align:

            process = Popen(['muscle', '-in', record.id+"_premuscle.fasta","-out",
                             record.id+"_align.fasta"], stdout=PIPE, stderr=PIPE)
            OUTPUT = process.stdout.read().decode("utf-8")
            process.stdout.close()
            align.write(OUTPUT)

        align.close()
        file.close()

    print("Alignment ✓")


def Maketree(usr_query, input_m='record.id+"_align.fasta"'):

    for record in SeqIO.parse(usr_query, "fasta"):

        with open(record.id+"_align.fasta", "r") as input_m, \
                open(record.id+"_tree.nw","w") as output_tree:

            tree = Popen(['muscle','-maketree','-in',record.id+"_align.fasta",'-out',
                          record.id+"_tree.nw",'-cluster','neighborjoining'], stderr=PIPE)

    input_m.close()
    output_tree.close()
    print("Tree ✓")