#! /usr/bin/env python

import argparse
import re
import sys
import gzip

baseRev = {}
baseRev['A'] = 'T'
baseRev['C'] = 'G'
baseRev['G'] = 'C'
baseRev['T'] = 'A'
baseRev['R'] = 'Y'
baseRev['Y'] = 'R'
baseRev['S'] = 'S'
baseRev['W'] = 'W'
baseRev['K'] = 'M'
baseRev['M'] = 'K'
baseRev['B'] = 'V'
baseRev['D'] = 'H'
baseRev['H'] = 'D'
baseRev['V'] = 'B'
baseRev['N'] = 'N'

baseAmb = {}
baseAmb['A'] = 'A'
baseAmb['C'] = 'C'
baseAmb['G'] = 'G'
baseAmb['T'] = 'T'
baseAmb['R'] = '[AG]'
baseAmb['Y'] = '[CT]'
baseAmb['S'] = '[GC]'
baseAmb['W'] = '[AT]'
baseAmb['K'] = '[GT]'
baseAmb['M'] = '[AC]'
baseAmb['B'] = '[CGT]'
baseAmb['D'] = '[AGT]'
baseAmb['H'] = '[ACT]'
baseAmb['V'] = '[ACG]'
baseAmb['N'] = '[ACGTN]'

def getCommandArgs():
    """This funciton processes the sys.argv array to get command line arguments.
    """
    parser = argparse.ArgumentParser(description='Estimate the number of occurences of a given sequence in a fasta file.')
    parser.add_argument('-f','--fasta', metavar='fasta', dest='fasta',  type=str, help='Input fasta file', required=True)
    parser.add_argument('-s','--sequence', metavar='fileOfSequences', dest='seqs',  type=str, help='File with sequences (2 cols, name and sequence)', required=True)
    parser.add_argument('-o','--output', metavar='outFile', dest='out',  type=str, help='Output file', default="")
    parser.add_argument('-b', '--make-bed', metavar="bedFile", dest='bed', type=str, default="", help="Output a bed file with locations of sequence.")
    parser.add_argument('-n', '--no-complement', dest="nocomp", action="store_true")
    args = parser.parse_args()
    return (args)


def revComplement(seq):
    """Get reverse complement of dna sequence.
    """
    revseq = ""
    for base in seq:
        revseq += baseRev[base]
    return(revseq[::-1])

def makeAmbiguous(seq):
    """Add ambiguity in the string using regular expression format.
    """
    newseq = ""
    for base in seq:
        newseq += baseAmb[seq]
    return(newseq)

def processOneChromosome(cname, cseq, seqs, tots, outfile):
    """Process one chromosome before jumping onto the next. 
    Output the details of this one for each of the matching sequences.
    """
    for sname in seqs:
        seq = seqs[sname]
        allfinds = re.findall(seq, cseq, re.I)
        fwdfinds = len(allfinds)
        seq = revComplement(seq)
        allfinds = re.findall(seq, cseq, re.I)
        revfinds = len(allfinds)
        tots[sname][0] += fwdfinds
        tots[sname][1] += revfinds
        if (fwdfinds > 0 or revfinds > 0):
            outfile.write(cname+"\t"+sname+"\t"+str(fwdfinds)+"\t"+str(revfinds)+"\n")

def processOneChromosomeBed(cname, cseq, seqs, bedfile):
    """Process one chromosome, to get all the + and -
    strand matches, and output them as a bed file.
    """
    for sname in seqs:
        seq = seqs[sname]
        alliter = re.finditer(seq, cseq, re.I)
        for match in alliter:
            bedfile.write(cname+"\t"+str(match.start())+"\t"+str(match.end())+"\t"+sname+"_"+cname+"_"+str(match.start())+"\t0\t+\n")
        seq = revComplement(seq)
        alliter = re.finditer(seq, cseq)
        for match in alliter:
            bedfile.write(cname+"\t"+str(match.start())+"\t"+str(match.end())+"\t"+sname+"_"+cname+"_"+str(match.start())+"\t0\t-\n")

        
def readPatterns(sequence):
    """Read in the patterns to be matched.
    """
    seqs = {}
    seqfile = open(sequence)
    for line in seqfile:
        toks = line.strip().split()
        if len(toks) != 2:
            print "Number of tokens on line is not 2. Kill!!!", line
            sys.exit(1)
        seqs[toks[0]] = toks[1]
    seqfile.close()
    return(seqs)

if __name__ == "__main__":
    args = getCommandArgs()
    outfile = sys.stdout
    if args.out != "":
        outfile = open(args.out, "w")
    if args.bed != "":
        bedfile = open(args.bed, "w")
    seqs = readPatterns(args.seqs)
    totalCounts = {}
    for sname in seqs:
        totalCounts[sname] = [0,0]
    if args.fasta[-3:] == ".gz":
        infile = gzip.open(args.fasta)
    else:
        infile = open(args.fasta)
    
    chrname = ""
    chrseq = ""
    ## Read one chr at a time of fasta
    for line in infile:
        line = line.strip()
        if line[0] == ">":
            if chrname != "":
                processOneChromosome(chrname, chrseq, seqs, totalCounts, outfile)
                if args.bed != "":
                    processOneChromosomeBed(chrname, chrseq, seqs, bedfile)
            chrname = line[1:]
            chrseq = ""
        else:
            chrseq += line
    if chrname != "":
        processOneChromosome(chrname, chrseq, seqs, totalCounts, outfile)
        if args.bed != "":
            processOneChromosomeBed(chrname, chrseq, seqs, bedfile)
    infile.close()
    for sname in seqs:
        outfile.write("Total\t"+sname+"\t"+str(totalCounts[sname][0])+"\t"+str(totalCounts[sname][1])+"\n")
    if args.out != "":
        outfile.close()
    if args.bed != "":
        bedfile.close()
