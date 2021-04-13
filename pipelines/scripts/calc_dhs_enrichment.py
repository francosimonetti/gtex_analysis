import os, sys

SNPRES_FIELDS = ['rsid', 'chrom', 'pos', 'logp', 'target', 'maf']
class SNPRes(collections.namedtuple('_SNPRes', SNPRES_FIELDS)):
    __slots__ = ()

# need background first
# read list of snp_ids
# read file with dhs regions
# calc dhs
# - is it multi tissue? different script?
# write output with enrichment and p-value and nº SNPs in DHS and fraction

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

    opts = parse_args()
    allsnps_file = opts.input
    dhs_file     = opts.dhsfile
    outputfile   = opts.output
    isannotated  = opts.annotated
    master_bg_file = opts.bg

    dhs_frac_rand = dict()
    dhs_frac_type_rand = dict()

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
    
    
    filefmt = f'{resdir}/{tissue}/{trans_eqtls_file}'
    if not os.path.exists(filefmt):
        print("File does not exist", filefmt)
        raise
    trans_eqtls = tejaas(filefmt)
    
    dhs_annotated = 0
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
            