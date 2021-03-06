#! /usr/bin/env python3
import os
import re
import sys
import argparse
import pandas as pd
from Bio import SeqIO
cigarre = re.compile('[0-9]+[ISMD]')

def complement(seq):
    # calculate the complement sequence string
    # input a string
    # return a string
    iupac = {
        "A": "T", "G": "C", "C": "G", "T": "A", "Y": "R", "R": "Y",
        "W": "W", "S": "S", "K": "M", "M": "K", "D": "H", "V": "B",
        "H": "D", "B": "V", "N": "N", "X": "X", "-": "-",
        "a": "t", "g": "c", "c": "g", "t": "a", "y": "r", "r": "y",
        "w": "w", "s": "s", "k": "m", "m": "k", "d": "h", "v": "b",
        "h": "d", "b": "v", "n": "n", "x": "x", "-": "-",
        " ": " ", ".": "."
    }
    return ''.join(iupac[i] for i in seq)


def reverse_complement(seq):
    # calculate the reverse complement sequence string
    # input a string
    # return a string
    return complement(seq)[::-1]

def splitcigar(x):
    return (int(x[0:-1]), x[-1])


def digitstring(x):
    # convert numbers into list of strings by digit place
    # 101 to ['100', '00', '1']
    x = str(x)
    if len(x) == 1:
        return [x]
    else:
        return [x[0]+'0' * (len(x) - 1)] + digitstring(x[1:])


def explainflags(x):
    flags = {
        '1': 'template having multiple segments in sequencing',
        '2': 'each segment properly aligned according to the aligner',
        '4': 'segment unmapped',
        '8': 'next segment in the template unmapped',
        '10': 'SEQ being reverse complemented',
        '20': 'SEQ of the next segment in the template being reverse complemented',
        '40': 'the first segment in the template',
        '80': 'the last segment in the template',
        '100': 'secondary alignment',
        '200': 'not passing filters, such as platform/vendor quality controls',
        '400': 'PCR or optical duplicate',
        '800': 'supplementary alignment'
    }
    flagplus = {
        '3': ['1', '2'], '5': ['1', '4'], '6': ['2', '4'], '7': ['1', '2', '4'],
        '9': ['1', '8'], 'a': ['2', '8'], 'b': ['1', '2', '8'], 'c': ['4', '8'],
        'd': ['1', '4', '8'], 'e': ['2', '4', '8'], 'f': ['1', '2', '4', '8'],
        'A': ['2', '8'], 'B': ['1', '2', '8'], 'C': ['4', '8'],
        'D': ['1', '4', '8'], 'E': ['2', '4', '8'], 'F': ['1', '2', '4', '8']
    }
    bit = hex(int(x))[2:]
    rawlist = digitstring(bit)
    bitlist = list()
    for x in rawlist:
        if x in list(flags.keys()):
            bitlist.append(x)
        else:
            bitlist.extend([i + '0' * (len(x) - 1) for i in flagplus[x[0]]])

            return {k: flags[k] for k in bitlist}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'stack the sam file'
    )
    parser.add_argument('-r', '--reference', help = 'reference fasta file')
    parser.add_argument('-s', '--sam', help = 'input sam file')
    parser.add_argument('-o', '--outprefix', help = 'output prefix')
    args = vars(parser.parse_args())

    fa = [record for record in SeqIO.parse(args['reference'], "fasta")]

    samcol = [
        'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR',
        'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL'
    ]


    sam = pd.read_table(
        args['sam'], sep='\t', header=None
    )
    sam = sam[sam.columns[0:11]]
    sam.columns = samcol

    sam = sam.loc[sam['CIGAR'] != '*',]
    sam.reset_index(drop=True, inplace=True)
    sam['CIGARLIST'] = sam['CIGAR'].map(lambda x: cigarre.findall(x))

    samstart = sam['POS'].min() - 1
    samend = sam['POS'].max() + sam['SEQ'].map(len).max()
    refseq = str(fa[0].seq).upper()
    # refseq = str(fa[0].seq[samstart:samend]).upper()
    reflen = len(refseq)
    refli = list(refseq)

    # start writing
    fseq = open(args['outprefix'] + '.seq', 'w')
    fseq.write(refseq + '\n')
    fstack = open(args['outprefix'] + '.stack', 'w')
    fstack.write(refseq + '\n')
    seqstat = dict()
    seqstat['A'] = [0] * len(refseq)
    seqstat['T'] = [0] * len(refseq)
    seqstat['G'] = [0] * len(refseq)
    seqstat['C'] = [0] * len(refseq)
    seqstat['N'] = [0] * len(refseq)
    seqstat['Deletion'] = [0] * len(refseq)
    seqstat['Other'] = [0] * len(refseq)

    for idx in sam.index:
        flag = explainflags(sam['FLAG'][idx])
        seq = sam['SEQ'][idx]
        cigar = sam['CIGARLIST'][idx]
        seqstart = sam['POS'][idx] - 1

        if ''.join(cigar).find('I') >= 0:
            continue
        nowseqsite = 0
        nowrefsite = seqstart - samstart
        outseq = ' ' * nowrefsite

        for x in cigar:
            cnum, ctype = splitcigar(x)
            if ctype == 'M':
                outseq = outseq + seq[nowseqsite: nowseqsite + cnum]
                nowrefsite += cnum
            elif ctype == 'S':
                pass
            elif ctype == 'D':
                outseq = outseq + ' ' * cnum * 2
                nowrefsite += cnum * 2
            # elif ctype == 'I':
            #     outseq = outseq + seq[nowseqsite: nowseqsite + cnum]
            #     stacks = [x[0:nowrefsite + 1] + ' ' * cnum + x[nowrefsite:] for x in stacks]
            #     nowrefsite += cnum
            nowseqsite += cnum
        fseq.write(outseq + '\n')
        inseq = False           # check whether covered by sequence
        for i in range(min(len(outseq), len(refseq))):
            if not inseq:
                if outseq[i] != ' ':
                    inseq = True
                else:
                    next
            if outseq[i] in ['A', 'T', 'G', 'C', 'N']:
                seqstat[outseq[i]][i] += 1
            elif outseq[i] == ' ':
                seqstat['Deletion'][i] += 1
            else:
                seqstat['Other'][i] += 1
        seqli = list(outseq)
        comp = list(map(lambda x: x[0] == x[1], zip(refli[0:len(seqli)], seqli)))
        seqchange = list(map(lambda x: '.' if x[0] else x[1], zip(comp, seqli)))
        result = ''.join(seqchange)
        fstack.write(result + '\n')
    fseq.close()
    fstack.close()
    fstat = open(args['outprefix'] + '.stat', 'w')
    fstat.write('seq\tA\tT\tG\tC\tN\tDeletion\tOther\n')
    fstat.write(
        '\n'.join(
            '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}'.format(
                refseq[i], seqstat['A'][i], seqstat['T'][i],
                seqstat['G'][i], seqstat['C'][i], seqstat['N'][i],
                seqstat['Deletion'][i], seqstat['Other'][i]
            ) for i in range(len(refseq))
        )
    )
    fstat.close()
