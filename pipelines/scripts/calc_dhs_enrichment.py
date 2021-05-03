import os, sys
import argparse
import collections
import numpy as np
import scipy.stats as ss

SNPRES_FIELDS = ['rsid', 'chrom', 'pos', 'logp', 'maf']
class SNPRes(collections.namedtuple('_SNPRes', SNPRES_FIELDS)):
    __slots__ = ()

# need background first
# read list of snp_ids
# read file with dhs regions
# calc dhs
# - is it multi tissue? different script?
# write output with enrichment and p-value and nº SNPs in DHS and fraction

def find_annotated(res_dict, dhs_file, isannotated=False):
    dhs = open(dhs_file)
    line = dhs.readline()
    prev_chrm = 0
    nannot = 0
    nannot_type = collections.defaultdict(int)

    sorted_res_dict = dict()
    for chrm in range(1,23):
        sorted_res_dict[chrm] = sorted(res_dict[chrm])

    while line:
        arr = line.rstrip().split("\t")
        if arr[0][3:] == "X" or arr[0][3:] == "Y":
            line = dhs.readline()
            continue
        chrm = int(arr[0][3:])
        start = int(arr[1])
        end = int(arr[2])
        if isannotated:
            atype = arr[9]
        if chrm != prev_chrm:
            remaining = sorted_res_dict[chrm]
            checked = 0
        if len(remaining) == 0:
            ## No more SNPs in this chromosome, just continue reading the DHS file
            line = dhs.readline()
        else:
            for pos in remaining:
                if pos < start:
                    checked += 1
                    remaining = sorted_res_dict[chrm][checked:]
                    continue # go to next SNP
                elif pos > end:
                    line = dhs.readline()
                    break # go to next DHS line, keep checking the remaining results
                else:
                    # this is an annotated SNP
                    checked += 1
                    remaining = sorted_res_dict[chrm][checked:]
                    nannot += 1
                    if isannotated:
                        nannot_type[atype] += 1
                    continue # go to next SNP
        prev_chrm = chrm
    dhs.close()
    return nannot, nannot_type

def tejaas(filepath):
    res = list()
    with open(filepath, 'r') as mfile:
        for line in mfile:
            if not line.startswith("chr"):
                continue
            arr   = line.strip().split("\t")
            rsid  = arr[0]
            chrom = int(arr[1])
            pos   = int(arr[2])
            maf   = float(arr[3])
            q     = float(arr[4])
            mu    = float(arr[5])
            sigma = float(arr[6])
            p     = float(arr[7])
            if sigma == 0:
                continue
            logp  = np.log10(p) if p != 0 else pvalue( (q - mu) / sigma)
            res.append(SNPRes(rsid=rsid, chrom=chrom, pos=pos, logp=-logp, maf=maf))
    return res


def parse_args():

    parser = argparse.ArgumentParser(description='Filter VCF file.')

    parser.add_argument('--in',
                        dest='input',
                        metavar='FILE',
                        help='list of SNPs')

    parser.add_argument('--dhs',
                        dest='dhsfile',
                        help='DHS regions file')

    parser.add_argument('--bg',
                        dest='bg',
                        metavar='FILE',
                        help='background file')

    parser.add_argument('--out',
                        dest='output',
                        metavar='FILE',
                        help='enrichment output file')

    parser.add_argument('--annotated',
                        dest='annotated',
                        action='store_true',
                        default=False,
                        help="is the DHS file annotated (DHS index or multi_tissue?)")

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':

    opts         = parse_args()
    eqtl_file    = opts.input
    dhs_file     = opts.dhsfile
    outputfile   = opts.output
    isannotated  = opts.annotated
    master_bg_file = opts.bg

    dhs_frac_rand = dict()
    dhs_frac_type_rand = dict()
    type_master_list = ["Cancer / epithelial","Cardiac","Digestive","Lymphoid","Musculoskeletal",\
                    "Myeloid / erythroid","Neural","Organ devel. / renal","Placental / trophoblast",\
                    "Primitive / embryonic","Pulmonary devel.", \
                    "Renal / cancer","Stromal A","Stromal B","Tissue invariant","Vascular / endothelial"]    

    if os.path.exists(master_bg_file):
        dhs_frac_type_rand = collections.defaultdict(dict)
        with open(master_bg_file) as ifile:
            next(ifile)
            for line in ifile:
                arr = line.strip().split("\t")
                if arr[0] == "all":
                    dhs_frac_rand = float(arr[2])
                else:
                    if isannotated:
                        dhs_frac_type_rand[arr[0]] = float(arr[2])
    else:
        print("Background file does not exist", master_bg_file)
        raise
    
    tissues_eval = list()
    transeqtls = dict()
    dhs_annotated = dict()
    dhs_annotated_type = dict()
    enrichment = collections.defaultdict(dict)
    pval_binom = collections.defaultdict(dict)
    
    
    if not os.path.exists(eqtl_file):
        print("File does not exist", eqtl_file)
        raise
    trans_eqtls = tejaas(eqtl_file)
    
    dhs_annotated = 0
    min_trans_eqtl = 0
    if len(trans_eqtls) > min_trans_eqtl:
        res_dict = dict()
        for chrm in range(1, 23):
            res_dict[chrm] = list()
        for x in trans_eqtls:
            res_dict[x.chrom].append(x.pos)
        dhs_annotated, dhs_annotated_type = find_annotated(res_dict, dhs_file, isannotated)
        pval_binom["all"] = ss.binom_test(dhs_annotated, len(trans_eqtls), dhs_frac_rand, alternative='greater')
        enrichment["all"] = float(dhs_annotated) / len(trans_eqtls) / dhs_frac_rand
        if isannotated:
            for dhs_type in type_master_list:
                pval_binom[dhs_type] = ss.binom_test(dhs_annotated_type[dhs_type], len(trans_eqtls), dhs_frac_type_rand[dhs_type], alternative='greater')
                enrichment[dhs_type] = float(dhs_annotated_type[dhs_type]) / len(trans_eqtls) / dhs_frac_type_rand[dhs_type]
        print (f"{dhs_annotated} annotated out of {len(trans_eqtls)}. Enrichment = {enrichment['all']} - {pval_binom['all']}")

    with open(outputfile, 'w') as outstream:
        outstream.write("dhs_type\ttranseqtls\tinDHS\tEnrichment\tpval_binom\n")
        nteqtl = len(trans_eqtls)
        if nteqtl > 0:
            outstream.write(f"all\t{nteqtl}\t{dhs_annotated}\t{enrichment['all']}\t{pval_binom['all']}\n")
            if isannotated:
                for dhs_type in type_master_list:
                    outstream.write(f"{dhs_type}\t{nteqtl}\t{dhs_annotated_type[dhs_type]}\t{enrichment[dhs_type]}\t{pval_binom[dhs_type]}\n")
        else:
            outstream.write(f"all\t{nteqtl}\t0\t0\t1\n")
            