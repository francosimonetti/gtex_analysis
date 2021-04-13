import sys
import numpy as np
import collections
import os
import argparse

type_master_list = ["Cancer / epithelial","Cardiac","Digestive","Lymphoid","Musculoskeletal",\
                    "Myeloid / erythroid","Neural","Organ devel. / renal","Placental / trophoblast",\
                    "Primitive / embryonic","Pulmonary devel.", \
                    "Renal / cancer","Stromal A","Stromal B","Tissue invariant","Vascular / endothelial"]    

def read_snplist(filename, mafcutoff=0.01):
    rsidlist = list()
    maflist  = list()
    with open(filename) as instream:
        next(instream)
        for line in instream:
            arr  = line.strip().split("\t")
            rsid = arr[0]
            maf  = float(arr[1])
            if maf >= mafcutoff and maf <= (1 - mafcutoff) :
                rsidlist.append(rsid)
                maflist.append(maf)
    return rsidlist, maflist

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

def annotated_random(gwrsids, nchoose, dhs_file, isannotated=False):
    chooseidx = np.sort(np.random.choice(len(gwrsids), nchoose, replace = False))
    res_dict = dict()
    for chrm in range(1, 23):
        res_dict[chrm] = list()
    for idx in chooseidx:
        var_id = gwrsids[idx]
        info = var_id.split('_')
        chrm = int(info[0][3:])
        bppos = int(info[1])
        res_dict[chrm].append(bppos)
    nannot, nannot_type = find_annotated(res_dict, dhs_file, isannotated)
    return nannot, nannot_type

def sample_rand_bg(snps_list, dhs_file, nchoose = 20000, niter = 20, isannotated=False):
    nannot_rand = list()
    nannot_type_list = list()
    nannot_type_array = list()
    print(f'Iteration', end="")
    for k in range(niter):
        nannot_k, nannot_type_k = annotated_random(snps_list, nchoose, dhs_file, isannotated)
        print(f' {k}', end="")
        nannot_rand.append(nannot_k)
        nannot_type_array.append([nannot_type_k[k] if k in nannot_type_k else 0 for k in type_master_list ])
    print("")
    frac_rand = np.mean(nannot_rand) / nchoose
    frac_rand_type = np.mean(np.array(nannot_type_array), axis=0) / nchoose
    return frac_rand, dict(zip(type_master_list, frac_rand_type))

def parse_args():

    parser = argparse.ArgumentParser(description='Filter VCF file.')

    parser.add_argument('--in',
                        dest='input',
                        metavar='FILE',
                        help='list of SNPs')

    parser.add_argument('--dhs',
                        dest='dhsfile',
                        help='DHS regions file')

    parser.add_argument('--out',
                        dest='output',
                        metavar='FILE',
                        help='background')

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
    master_bg_output_file = opts.output
    isannotated = opts.annotated
    

    # read all variants ids
    snps_list, mafs = read_snplist(allsnps_file)            
    dhs_frac_rand = dict()
    dhs_frac_type_rand = dict()

    # outdir = os.path.join(resdir, "dhs_enrichments")
    # if not os.path.exists(outdir):
    #     os.makedirs(outdir)

    print("Creating background file")
    with open(master_bg_output_file, 'w') as ofile:
        ofile.write("dhs_type\tn_snps\tdhs_frac_rand\n")
        dhs_frac_rand, dhs_frac_type_rand = sample_rand_bg(snps_list, dhs_file, nchoose = 20000, niter = 50, isannotated = isannotated)
        Nsnps = len(snps_list)
        print (f'Fraction of annotated SNPs: global - all - {dhs_frac_rand:7.4f}')
        ofile.write(f"all\t{Nsnps}\t{dhs_frac_rand}\n")
        if isannotated:
            for dhs_type in type_master_list:
                print(f"{dhs_type}\t{Nsnps * dhs_frac_type_rand[dhs_type]}\t{dhs_frac_type_rand[dhs_type]}")
                ofile.write(f"{dhs_type}\t{Nsnps * dhs_frac_type_rand[dhs_type]}\t{dhs_frac_type_rand[dhs_type]}\n")
