import argparse
import collections
import re

def parse_args():

    parser = argparse.ArgumentParser(description='Finds SNPs in list A in B.')

    parser.add_argument('--A',
                        dest='inputA',
                        metavar='FILE',
                        help='list of SNPs A')

    parser.add_argument('--B',
                        dest='inputB',
                        metavar='FILE',
                        help='list of SNPs B')

    parser.add_argument('--out',
                        dest='output',
                        metavar='FILE',
                        help='list of SNPs B that is in A')

    parser.add_argument('--isnot',
                        dest="isnot",
                        action="store_true",
                        help="list of SNPs B that is NOT in A")

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':

    opts  = parse_args()
    fileA = opts.inputA
    fileB = opts.inputB
    isnot = opts.isnot
    outfile = opts.output

    listA_dict = collections.defaultdict(lambda: False)
    listB_found = collections.defaultdict(lambda: False)
    with open(fileA) as fin:
        for line in fin:
            if not line.startswith("chr"):
                continue
            arr = line.strip().split("\t")
            listA_dict[arr[0]] = True

    with open(outfile, 'w') as fout:
        with open(fileB) as fin:
            if isnot:
                # find snps in B that are NOT in A
                for line in fin:
                    if not line.startswith("chr"):
                        continue
                    arr = line.strip().split("\t")
                    if not listA_dict[arr[0]] and not listB_found[arr[0]]:
                        fout.write(line)
                        listB_found[arr[0]] = True
            else:
                # find snps in B that are in A
                for line in fin:
                    if not line.startswith("chr"):
                        continue
                    arr = line.strip().split("\t")
                    if listA_dict[arr[0]] and not listB_found[arr[0]]:
                        fout.write(line)
                        listB_found[arr[0]] = True
                