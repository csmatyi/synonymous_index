#!/usr/bin/python

import re, os, sys, getopt, math

def translate(mrna_sequence):
    genetic_code = {"AAA":'K',"AAC":'N',"AAG":'K',"AAT":'N',"ACA":'T',"ACC":'T',"ACG":'T',"ACT":'T',"AGA":'R',
    "AGC":'S',"AGG":'R',"AGT":'S',"ATA":'I',"ATC":'I',"ATG":'M',"ATT":'I',"CAA":'Q',"CAC":'H',
    "CAG":'Q',"CAT":'H',"CCA":'P',"CCC":'P',"CCG":'P',"CCT":'P',"CGA":'R',"CGC":'R',"CGG":'R',
    "CGT":'R',"CTA":'L',"CTC":'L',"CTG":'L',"CTT":'L',"GAA":'E',"GAC":'D',"GAG":'E',"GAT":'D',
    "GCA":'A',"GCC":'A',"GCG":'A',"GCT":'A',"GGA":'G',"GGC":'G',"GGG":'G',"GGT":'G',"GTA":'V',
    "GTC":'V',"GTG":'V',"GTT":'V',"TAA":'*',"TAC":'Y',"TAG":'*',"TAT":'Y',"TCA":'S',"TCC":'S',
    "TCG":'S',"TCT":'S',"TGA":'*',"TGC":'C',"TGG":'W',"TGT":'C',"TTA":'L',"TTC":'F',"TTG":'L',"TTT":'F'}
    mrna_sequence = mrna_sequence.upper().replace('U','T')
    protein_sequence = ""
    for i in range(0,len(mrna_sequence),3):
        if mrna_sequence[i:i+3] in genetic_code.keys():
            protein_sequence = protein_sequence + str(genetic_code[mrna_sequence[i:i+3]])
        else:
            protein_sequence = protein_sequence + "X"
    return protein_sequence

def getsd(si,siavg):
    sd = 0
    for val in si:
        sd += (val - siavg)**2
    sd /= len(si)
    return math.sqrt(sd)

def main(argv):
    codefreq = {'A':4/64.0, 'C':2/64.0, 'D':2/64.0, 'E':2/64.0, 'F':2/64.0, 'G':4/64.0, 'H':2/64.0, 'I':3/64.0, 'K':2/64.0, 'L':6/64.0, 
    'M':1/64.0, 'N':2/64.0, 'P':4/64.0, 'Q':2/64.0, 'R':6/64.0, 'S':6/64.0, 'T':4/64.0, 'V':4/64.0, 'W':1/64.0, 'Y':2/64.0}

    inputfile = ''
    outputfile = ''
    protein = True

    try:
        opts, args = getopt.getopt(argv,"hi:o:n",["ifile=","ofile="])
    except getopt.GetoptError:
        print('Error! Run as: SI2.py -i <inputfile> -o <outputfile> -n <nucleotide>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('SI2.py -i <inputfile> -o <outputfile> -n <nucleotide>')
            sys.exit()
        if opt == '-n':
            protein = False
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    length = 0
    seqid = ''
    SI = []
    nseq = 0
    sequence = ""

    fi = open(inputfile,'r')
    fo = open(outputfile,'w')
    for line in fi:
        lin = line.rstrip('\n')
        line = str(lin)
        if len(line) > 0:
            if line[0] == '>':
                if nseq > 0:
                    # get sequence
                    if protein is False:
                        sequence = translate(sequence)
                    sequence = sequence.upper().replace('*','')
                    length = len(sequence)
                    for aa in sequence:
                        if aa in codefreq.keys():
	                          SI.append(codefreq[aa])
                    
                    # calculate stats
                    si_avg = 0
                    for si in SI:
                        si_avg += si
                    si_total = si_avg
                    si_avg /= length
                    si_sd = getsd(SI,si_avg)
                    fo.write(str(nseq) + '\t' + seqid + '\t' + str(length) + '\t' + str(si_total) + '\t' + str(si_avg) + '\t' + str(si_sd) + '\n')

                    # reset variables                    
                    length = 0
                    sequence = ""
                    SI = []
                nseq += 1
                seqid = line[1:]
            else:
                for letter in line:
                    sequence += letter

    if protein is False:
        sequence = translate(sequence)
    sequence = sequence.upper().replace('*','')
    length = len(sequence)
    for aa in sequence:
        if aa in codefreq.keys():
            SI.append(codefreq[aa])
    
    si_avg = 0
    for si in SI:
        si_avg += si
    si_total = si_avg
    si_avg /= length
    si_sd = getsd(SI,si_avg)
    fo.write(str(nseq) + '\t' + seqid + '\t' + str(length) + '\t' + str(si_total) + '\t' + str(si_avg) + '\t' + str(si_sd) + '\n')

    fi.close()
    fo.close()

if __name__ == "__main__":
    main(sys.argv[1:])
