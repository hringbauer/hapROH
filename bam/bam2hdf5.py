import numpy as np
import h5py
import argparse
import pysam
import os
import sys

def bam2hdf5(path2bam, refHDF5, iid="", minMapQual = 30, minBaseQual = 20, outPath = ""):
    f = h5py.File(refHDF5, 'r')
    pos = f['variants/POS']
    rec = f['variants/MAP']
    ref = f['variants/REF'][:]
    alt = f['variants/ALT'][:, 0]
    assert(len(ref) == len(pos))
    assert(len(alt) == len(pos))
    assert(len(rec) == len(pos))
    # 1 here means there is just one sample
    # 2 : [ref read count, alt read count]
    ad = np.zeros((len(pos), 1, 2))
    pos = np.array(pos) - 1 # pysam uses 0-based coordinate

    major_adj = 0
    minor_adj = 0
    major_foc = 0
    minor_foc = 0
    base2index = {'A':0, 'C':1, 'G':2, 'T':3}
    sf = pysam.AlignmentFile(path2bam, "rb")
    print(f'total number of mapped reads: {sf.mapped}')
    for i, p in enumerate(pos):
        for pileupcolumn in sf.pileup('X', p-4, p+5, min_mapping_quality=minMapQual, min_base_quality=minBaseQual, truncate=True):
            basePos = pileupcolumn.pos
            rc = np.zeros(4)
            for pileupread in pileupcolumn.pileups:
                query_pos = pileupread.query_position
                seqlen = pileupread.alignment.query_alignment_length
                if not pileupread.is_del and not pileupread.is_refskip and query_pos >= 5 and query_pos < seqlen-5 :
                    baseCall = pileupread.alignment.query_sequence[query_pos]
                    rc[base2index[baseCall]] += 1
                    if basePos == p:
                        #print(f'RG: {pileupread.alignment.get_tag("RG")}')
                        if baseCall == ref[i]:
                            ad[i, 0, 0] += 1
                        elif baseCall == alt[i]:
                            ad[i, 0, 1] += 1
                        else:
                            print(f'at position {pos[i]} the read {baseCall} is neither ref {ref[i]} nor alt {alt[i]}.')

            if np.sum(rc) <= 1:
                continue # sites covered by only 1 read is not informative for estimating genotyping error
            if basePos != p:
                major_adj += np.max(rc)
                minor_adj += np.sum(rc) - np.max(rc)
            else:
                major_foc += np.max(rc)
                minor_foc += np.sum(rc) - np.max(rc)
            if np.sum(rc == 0) < 2:
                print(f'more than 2 genotype calls at {basePos}')

    print(f'number of major reads at flanking sites: {major_adj}')
    print(f'number of minor reads at flanking sites: {minor_adj}')
    print(f'number of major reads at focal sites: {major_foc}')
    print(f'number of minor reads at focal sites: {minor_foc}')

    print(f'err rate at flanking sites: {minor_adj/(minor_adj + major_adj)}')
    print(f'err rate at focal sites: {minor_foc/(minor_foc + major_foc)}')

    # finished reading bam file and we have made an estimate for genotyping error
    # now write a hdf5 file for read count at target sites
    l = len(pos)
    k = 1 # 1 bam = 1 sample
    rec = f['variants/MAP']
    gt = np.zeros((l, k, 2)) # dummy gt matrix, won't be used for inference
    dt = h5py.special_dtype(vlen=str)  # To have no problem with saving
    
    if len(outPath) != 0:
        hdf5Name = outPath
    else:
        bamFileName = os.path.basename(path2bam)
        hdf5Name = bamFileName[:bamFileName.find(".bam")] + ".hdf5"
    
    if len(iid) == 0:
        bamFileName = os.path.basename(path2bam)
        iid = bamFileName[:bamFileName.find(".bam")]
    if os.path.exists(hdf5Name):  # Do a Deletion of existing File there
        os.remove(hdf5Name)
    
    
    
    with h5py.File(hdf5Name, 'w') as f0:
    # Create all the Groups
        f_map = f0.create_dataset("variants/MAP", (l,), dtype='f')
        f_ad = f0.create_dataset("calldata/AD", (l, k, 2), dtype='i')
        f_ref = f0.create_dataset("variants/REF", (l,), dtype=dt)
        f_alt = f0.create_dataset("variants/ALT", (l,), dtype=dt)
        f_pos = f0.create_dataset("variants/POS", (l,), dtype='i')
        f_gt = f0.create_dataset("calldata/GT", (l, k, 2), dtype='i')
        f_samples = f0.create_dataset("samples", (k,), dtype=dt)

        #   Save the Data
        f_map[:] = rec
        f_ad[:] = ad
        f_ref[:] = ref.astype("S1")
        f_alt[:] = alt.astype("S1")
        f_pos[:] = pos + 1 # add one back
        f_gt[:] = gt
        f_samples[:] = np.array([iid]).astype("S10")
        print(f'saving sample as {iid}')

    # the second return value is the number of sites covered by at least 1 read
    return minor_adj/(minor_adj + major_adj), np.sum(np.sum(np.sum(ad, axis=1), axis=1) > 0), hdf5Name


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert bam file to hdf5 format that stores readcount info at target sites.')
    parser.add_argument('-b', action="store", dest="bam", type=str, required=True,
                        help="path to bam file")
    parser.add_argument('-r', action="store", dest="ref", type=str, required=True,
                        help="path to reference panel")
    parser.add_argument('-o', action="store", dest="out", type=str, required=False, default="",
                        help="path to the hdf5 file output")
    parser.add_argument('-i', action="store", dest="iid", type=str, required=False, default="",
                        help="IID of the target individual. If unspecified, will use the prefix of the bam file.")
    args = parser.parse_args()

    iid = args.iid
    if len(iid) == 0:
        bamName = os.path.basename(args.bam)
        iid = bamName[:bamName.find(".bam")]

    err, numSitesCovered, path2hdf5 = bam2hdf5(args.bam, args.ref, iid=iid, outPath=args.out)
    if numSitesCovered < 500:
        print(f'not enough sites covered to make inference...')
        sys.exit()

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")

    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom
    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom_BFGS
    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom_2d

    _, mle, low95, up95 = hapCon_chrom(iid, 'X', save=False, save_fp=False, n_ref=2504, diploid_ref=True, 
        exclude_pops=[], conPop=[], e_model="readcount_contam", p_model="SardHDF5", 
        readcounts=True, random_allele=False, post_model="Standard", path_targets = path2hdf5, 
        folder_out='/mnt/archgen/users/yilei/Data/iberian_BAM/hapCon/',
        h5_path1000g='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr',
        meta_path_ref='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv', 
        prefix_out=iid, c=np.arange(0, 0.2, 0.005), roh_in=1, roh_out=0, roh_jump=300, e_rate=err, e_rate_ref=0.0,
        max_gap=0, cutoff_post = 0.999, roh_min_l = 0.01, logfile=False)

    mle_bfgs, low95_bfgs, up95_bfgs = hapCon_chrom_BFGS(iid, 'X', save=False, save_fp=False, n_ref=2504, diploid_ref=True, 
        exclude_pops=[], conPop=[], e_model="readcount_contam", p_model="SardHDF5", 
        readcounts=True, random_allele=False, post_model="Standard", path_targets = path2hdf5, 
        folder_out='/mnt/archgen/users/yilei/Data/iberian_BAM/hapCon/',
        h5_path1000g='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr',
        meta_path_ref='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv', 
        prefix_out=iid, c=0.025, roh_in=1, roh_out=0, roh_jump=300, e_rate=err, e_rate_ref=0.0,
        max_gap=0, cutoff_post = 0.999, roh_min_l = 0.01, logfile=False)

    mle2, se = hapCon_chrom_2d(iid, 'X', save=False, save_fp=False, n_ref=2504, diploid_ref=True, 
        exclude_pops=[], conPop=[], e_model="readcount_contam", p_model="SardHDF5", 
        readcounts=True, random_allele=False, post_model="Standard", path_targets = path2hdf5,
        folder_out='/mnt/archgen/users/yilei/Data/iberian_BAM/hapCon/',
        h5_path1000g='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr',
        meta_path_ref='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv',
        prefix_out=iid, c=0.025, roh_in=1, roh_out=0, roh_jump=300, e_rate=0.01, 
        e_rate_ref=0.0, max_gap=0, cutoff_post = 0.999, roh_min_l = 0.01, logfile=False)
    

    with open(f'{iid}.hapcon.trim5.txt', 'w') as out:
        out.write(f'Number of target sites covered by at least one read: {numSitesCovered}\n')
        out.write(f'Method1: Fixing genotyping error rate\n')
        out.write(f'\tEstimated genotyping error via flanking region: {round(err, 6)}\n')
        out.write(f'\tMLE for contamination using grid search: {round(mle, 6)} ({round(low95, 6)} - {round(up95, 6)})\n')
        out.write(f'\tMLE for contamination using BFGS: {round(mle_bfgs[0], 6)} ({round(low95_bfgs[0], 6)} - {round(up95_bfgs[0], 6)})\n')
        out.write("\n")
        out.write(f'Method2: BFGS on contamination and genotyping error jointly\n')
        out.write(f'\tcontamination: {round(mle2[0], 6)}({round(mle2[0] - 1.96*se[0], 6)} - {round(mle2[0] + 1.96*se[0], 6)})\n')
        out.write(f'\tgenotyping error: {round(mle2[1], 6)}({round(mle2[1] - 1.96*se[1], 6)} - {round(mle2[1] + 1.96*se[1], 6)})\n')
        out.write("\n")

