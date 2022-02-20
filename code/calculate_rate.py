#! /usr/bin/env python3

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'calculate rate'
    )
    parser.add_argument('-i', '--input', help = 'input file')
    parser.add_argument('-a', '--aim', help = 'aim edit')
    parser.add_argument('-o', '--outprefix', help = 'output prefix')

    args = vars(parser.parse_args())

    result = dict()
    result['total_reads'] = 0
    result['cover_reads'] = 0
    result['aimedit'] = 0
    result['deletion'] = 0
    result['wildtype'] = 0
    result['otheredit'] = 0
    with open(args['input'], 'r') as f:
        refseq = f.readline().strip()
        aimsite = int((len(refseq) - 1)/2)
        for x in f:
            result['total_reads'] += 1
            x = x.rstrip()
            if len(x) <= aimsite + 1: continue
            xs = x.lstrip()
            mapseq = refseq[len(x) - len(xs):len(x)]
            msite = aimsite - (len(x) - len(xs))
            if msite < 0:
                # jump reads not cover aim site
                continue
            if xs.find(' ') > -1:
                # reads with deletion
                result['deletion'] += 1
                continue
            result['cover_reads'] += 1 # reads cover aim site
            blist = list(map(lambda x,y: x == y, list(mapseq), list(xs)))
            if sum(blist) == len(mapseq):
                # reads of wild type
                result['wildtype'] += 1
                continue
            else:
                if sum(blist) == len(mapseq) - 1 and xs[msite] == args['aim']:
                    # reads of aim edit
                    result['aimedit'] += 1
                    continue
                else:
                    # reads with other edits
                    result['otheredit'] += 1

    with open(args['outprefix'] + '.txt', 'w') as f:
        f.write('\n'.join([', '.join([str(x), str(y)]) for x,y in result.items()]))

################################################################################
