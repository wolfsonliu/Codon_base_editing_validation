#! /usr/bin/env python3

import sys
import argparse

parser = argparse.ArgumentParser(
    description = 'stat sites'
)
parser.add_argument(
    '-i', '--input',
    help='sequence per line, first line as reference'
)
parser.add_argument(
    '-o', '--output', help='tab separated file'
)
args = vars(parser.parse_args())

__stdin__ = sys.stdin
if args['input'] is not None:
    sys.stdin = open(args['input'], 'r')

__stdout__ = sys.stdout
if args['output'] is not None:
    sys.stdout = open(args['output'], 'w')

fin = sys.stdin
fout = sys.stdout
refseq = fin.readline()
count = {
    'A': [0] * len(refseq),
    'T': [0] * len(refseq),
    'G': [0] * len(refseq),
    'C': [0] * len(refseq),
    'D': [0] * len(refseq),
    'O': [0] * len(refseq)
}

for seq in fin:
    for i in range(len(seq)):
        if seq[i] in ['A', 'T', 'G', 'C']:
            count[seq[i]][i] += 1
        elif seq[i] == ' ':
            count['D'][i] += 1
        else:
            count['O'][i] += 1
fin.close()

countsum = [
    max(sum([count[k][i] for k in count]),1) for i in range(len(refseq))
]

rate = {
    'A': [count['A'][i]/countsum[i] for i in range(len(refseq))],
    'T': [count['T'][i]/countsum[i] for i in range(len(refseq))],
    'G': [count['G'][i]/countsum[i] for i in range(len(refseq))],
    'C': [count['C'][i]/countsum[i] for i in range(len(refseq))],
    'D': [count['D'][i]/countsum[i] for i in range(len(refseq))],
    'O': [count['O'][i]/countsum[i] for i in range(len(refseq))]
}

fout.write('seq\tA\tRate.A\tT\tRate.T\tG\tRate.G\tC\tRate.C\tDel\tRate.Del\tOther\tRate.Other\tSum\n')
fout.write(
    '\n'.join(
        '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}'.format(
            refseq[i],
            count['A'][i], rate['A'][i],
            count['T'][i], rate['T'][i],
            count['G'][i], rate['G'][i],
            count['C'][i], rate['C'][i],
            count['D'][i], rate['D'][i],
            count['O'][i], rate['O'][i],
            countsum[i]
        ) for i in range(len(refseq))
    )
)
fout.close()
####################
