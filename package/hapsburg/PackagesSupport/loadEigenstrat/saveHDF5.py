"""
Functions to save hdf5, and transfrom Eigenstrat to hdf5
@ Author: Harald Ringbauer, 2019, All rights reserved
"""
import numpy as np
import pandas as pd

import os as os
import h5py  # Python Package to do the HDF5.
from hapsburg.PackagesSupport.loadEigenstrat.loadEigenstrat import load_eigenstrat

import time
import pysam
from itertools import zip_longest
import sys

def save_hdf5(path, gt, ref, alt, pos, chs,
              rec, samples, ad=[], sex=[], clst=[],
              compression="gzip", gt_type="int8"):
        """Create a new HDF5 File with Input Data.
        Field Names are chosen to integrate with hapROH
        gt: Genotype data [l,k,2]
        ad: Allele depth [l,k,2]
        ref: Reference Allele [l]
        alt: Alternate Allele [l]
        pos: Position  [l]
        chs: Chromosome [l]
        m: Map position [l]
        samples: Sample IDs [k].
        Save genotype data as int8, readcount data as int16.
        ad: whether to save allele depth
        gt_type: What genotype data type save.
        sex: Vector of Sex [k]
        clst: Vector of Cluster [k]"""

        l, k, _ = np.shape(gt)  # Nr loci and Nr of Individuals

        if os.path.exists(path):  # Do a Deletion of existing File there
            os.remove(path)

        dt = h5py.special_dtype(vlen=str)  # To have no problem with saving

        with h5py.File(path, 'w') as f0:
            ### Create all the Groups
            f_map = f0.create_dataset("variants/MAP", (l,), dtype='f')
            f_ref = f0.create_dataset("variants/REF", (l,), dtype=dt)
            f_alt = f0.create_dataset("variants/ALT", (l,), dtype=dt)
            f_pos = f0.create_dataset("variants/POS", (l,), dtype='int32')
            f_ch = f0.create_dataset("variants/CH", (l,), dtype='int8')
            f_gt = f0.create_dataset("calldata/GT", (l, k, 2), dtype=gt_type, compression=compression)
            f_samples = f0.create_dataset("inds/samples", (k,), dtype=dt)
            
            # Create data types for smaples and cluster strings
            dtype_s = "S" + str(len(max(samples, key=len)))
            
            ### Save the Data
            f_map[:] = rec
            f_ref[:] = ref.astype("S1")
            f_alt[:] = alt.astype("S1")
            f_pos[:] = pos
            f_gt[:] = gt
            f_ch[:] = chs
            f_samples[:] = np.array(samples).astype(dtype_s)
            
            ### Optional Fields
            if np.shape(ad)[0]>0:
                f_ad = f0.create_dataset("calldata/AD", (l, k, 2), dtype='i')
                f_ad[:] = ad
                
            if len(sex)>0:
                f_sex = f0.create_dataset("inds/SEX", (k,), dtype=dt)
                f_sex[:] = sex.astype("S1")
                
            if len(clst)>0:
                dtype_cls = "S" + str(len(max(clst, key=len)))
                f_cls = f0.create_dataset("inds/CLS", (k,), dtype=dt)
                f_cls[:] = np.array(clst).astype(dtype_cls)
                   
            
def eigenstrat_to_hdf5(path_es = "", path_hdf5="",
                       packed=True, sep="\s+"):
    """Translates eigenstrat to hdf5 file. Saves new file.
    path_es: Path of the Eigenstrat [string]
    path_hdf5: Path of the hdf5 [string]
    """
    es = load_eigenstrat(base_path=path_es, packed=packed, sep=sep)
    df_snp = es.load_snp_df()
    df_ind = es.load_ind_df()
    
    print(f"Loading genotypes ({len(df_snp)} SNPs, {len(df_ind)} Inds)...")
    gt = es.get_geno_all()
    assert(np.shape(gt)[1]==len(df_ind))
    assert(np.shape(gt)[0]==len(df_snp))
    
    print(f"Transforming genotypes...")
    gt_new = np.zeros(np.shape(gt) + (2,), dtype="int8")
    gt_new[gt==2, :] = [1, 1]
    gt_new[gt==1, :] = [1, 0]
    gt_new[gt==3, :] = [3, 3]
    
    print(f"Saving HDF5...")
    save_hdf5(path=path_hdf5, gt=gt_new, ref=df_snp["ref"].values, alt=df_snp["alt"].values, 
              pos=df_snp["pos"].values, chs=df_snp["chr"].values, rec=df_snp["map"].values, 
              samples=df_ind["iid"].values, sex=df_ind["sex"].values, clst=df_ind["cls"].values,
              compression="gzip", gt_type="int8")
    print(f"Successfully saved data: {np.shape(gt_new)}")

def mpileup2hdf5(path2mpileup, refHDF5, iid="", s=-np.inf, e=np.inf, outPath="", output=True):
    f = h5py.File(refHDF5, 'r')
    pos = np.array(f['variants/POS'])
    ref = f['variants/REF']
    alt = f['variants/ALT']
    if len(alt.shape) == 2:
        alt = alt[:, 0].flatten()
    ref = np.array(ref).astype('str')
    alt = np.array(alt).astype('str')
    bases = np.array(['A', 'T', 'G', 'C'])
    subset1 = np.logical_and(pos >= s, pos <= e)
    subset2 = np.logical_and(np.in1d(ref, bases), np.in1d(alt, bases))
    print(f'exclude {np.sum(~subset1)} sites outside the specified region')
    print(f'exclude {np.sum(~subset2)} non-SNP sites')
    subset = np.logical_and(subset1, subset2)
    pos = pos[subset]
    rec = f['variants/MAP'][subset]
    ref = ref[subset]
    alt = alt[subset]
    # ref = f['variants/REF'][subset]
    
    # alt = f['variants/ALT']
    # is1D = len(alt.shape) == 1 # the 1240k hdf5's alt is a 2d array, while the new full 1000Genome hdf5 is 1d array
    # if is1D:
    #     alt = alt[subset]
    # else:
    #     alt = alt[subset, 0]

    assert(len(ref) == len(pos))
    assert(len(alt) == len(pos))
    assert(len(rec) == len(pos))
    # 1 here means there is just one sample
    # 2 : [ref read count, alt read count]
    ad = np.zeros((len(pos), 1, 2))

    major_adj = 0
    minor_adj = 0
    major_foc = 0
    minor_foc = 0
    l = len(pos)
    base2index = {'A':0, 'C':1, 'G':2, 'T':3}
    with open(path2mpileup) as f:
        for line in f:
            rc = np.zeros(4)
            contig, bp, _, coverage, readbases, baseQ = list(zip(*zip_longest(line.strip().split(), range(6))))[0]
            if int(coverage) == 0:
                continue

            insertion_index = readbases.find("+")
            insert = ""
            while insertion_index != -1:
                numInsertion = int(readbases[insertion_index+1])
                insert += readbases[insertion_index+2:insertion_index+2+numInsertion]
                insertion_index = readbases.find("+", insertion_index+2+numInsertion)
            
            rc[0] = readbases.count('A') + readbases.count('a') - insert.count('A') - insert.count('a')
            rc[1] = readbases.count('C') + readbases.count('c') - insert.count('C') - insert.count('c')
            rc[2] = readbases.count('G') + readbases.count('g') - insert.count('G') - insert.count('g')
            rc[3] = readbases.count('T') + readbases.count('t') - insert.count('T') - insert.count('t')
            coverage, bp = int(coverage), int(bp)
            i = np.searchsorted(pos, bp)
            if i < l and pos[i] == bp:
                # target sites
                # print(f'at target sites: {bp}, ref: {ref[i]}, alt: {alt[i]}')
                ad[i, 0, 0] = rc[base2index[ref[i]]]
                ad[i, 0, 1] = rc[base2index[alt[i]]]
                if coverage > 1:
                    major_foc += np.max(rc)
                    minor_foc += np.sum(rc) - np.max(rc)
            else:
                # sites flanking to the target sites
                if coverage > 1:
                    major_adj += np.max(rc)
                    minor_adj += np.sum(rc) - np.max(rc)

    if output:
        print(f'number of major reads at flanking sites: {int(major_adj)}')
        print(f'number of minor reads at flanking sites: {int(minor_adj)}')
        print(f'number of major reads at focal sites: {int(major_foc)}')
        print(f'number of minor reads at focal sites: {int(minor_foc)}')
        print(f'err rate at flanking sites: {minor_adj/(minor_adj + major_adj):.6f}')
        print(f'err rate at focal sites: {minor_foc/(minor_foc + major_foc):.6f}')

    # finished reading bam file and we have made an estimate for genotyping error
    # now write a hdf5 file for read count at target sites
    l = len(pos)
    k = 1 # 1 bam = 1 sample
    gt = np.zeros((l, k, 2)) # dummy gt matrix, won't be used for inference
    dt = h5py.special_dtype(vlen=str)  # To have no problem with saving
    
    if len(outPath) != 0 and not outPath.endswith("/"):
        outPath += "/"
    bamFileName = os.path.basename(path2mpileup)
    hdf5Name = outPath + bamFileName[:bamFileName.find(".mpileup")] + ".hdf5"

    if os.path.exists(hdf5Name):  # Do a Deletion of existing File there
        os.remove(hdf5Name)
    
    # iid is the sample ID name in the hdf5 file
    if len(iid) == 0:
        iid = bamFileName[:bamFileName.find(".mpileup")]
    
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
        f_pos[:] = pos
        f_gt[:] = gt
        f_samples[:] = np.array([iid]).astype("S50")
        print(f'saving sample as {iid} in {hdf5Name}')

    # the second return value is the number of sites covered by at least 1 read
    err = minor_adj/(minor_adj + major_adj)
    numSitesCovered = np.sum(np.sum(np.sum(ad, axis=1), axis=1) > 0)
    if output:
        print(f'estimated genotyping error by flanking sites: {err:.6f}')
        print(f'number of sites covered by at least one read: {numSitesCovered}, fraction covered: {numSitesCovered/len(pos):.3f}')
        print(f'hdf5 file saved to {hdf5Name}')
    return err, numSitesCovered, hdf5Name

def bam2hdf5(path2bam, refHDF5, ch="X", iid="", minMapQual=30, minBaseQual=20, s=-np.inf, e=np.inf, trim=0, outPath="", output=True):
    f = h5py.File(refHDF5, 'r')
    pos = np.array(f['variants/POS'])
    rec = f['variants/MAP']
    ref = f['variants/REF']
    alt = f['variants/ALT']
    is1D = len(alt.shape) == 1 # the 1240k hdf5's alt is a 2d array, while the new full 1000Genome hdf5 is 1d array
    if not is1D:
        alt = alt[:, 0].flatten()
    ref = np.array(ref).astype('str')
    alt = np.array(alt).astype('str')
    bases = np.array(['A', 'T', 'G', 'C'])
    subset1 = np.logical_and(pos >= s, pos <= e)
    subset2 = np.logical_and(np.in1d(ref, bases), np.in1d(alt, bases))
    print(f'exclude {np.sum(~subset1)} sites outside the specified region')
    print(f'exclude {np.sum(~subset2)} non-SNP sites')
    subset = np.logical_and(subset1, subset2)
    pos = pos[subset]
    rec = f['variants/MAP'][subset]
    ref = ref[subset]
    alt = alt[subset]

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
        for pileupcolumn in sf.pileup(ch, p-4, p+5, min_mapping_quality=minMapQual, min_base_quality=minBaseQual, truncate=True):
            basePos = pileupcolumn.pos
            rc = np.zeros(4)
            for pileupread in pileupcolumn.pileups:
                query_pos = pileupread.query_position
                seqlen = pileupread.alignment.query_alignment_length
                if not pileupread.is_del and not pileupread.is_refskip and query_pos >= trim and query_pos < seqlen - trim :
                    baseCall = pileupread.alignment.query_sequence[query_pos]
                    rc[base2index[baseCall]] += 1
                    if basePos == p:
                        #print(f'RG: {pileupread.alignment.get_tag("RG")}')
                        if baseCall == ref[i]:
                            ad[i, 0, 0] += 1
                        elif baseCall == alt[i]:
                            ad[i, 0, 1] += 1
                        # else:
                        #     print(f'at position {pos[i]} the read {baseCall} is neither ref {ref[i]} nor alt {alt[i]}.')

            if np.sum(rc) <= 1:
                continue # sites covered by only 1 read is not informative for estimating genotyping error
            if basePos != p:
                major_adj += np.max(rc)
                minor_adj += np.sum(rc) - np.max(rc)
            else:
                major_foc += np.max(rc)
                minor_foc += np.sum(rc) - np.max(rc)
            # if np.sum(rc == 0) < 2:
            #     print(f'more than 2 genotype calls at {basePos}')

    if output:
        print(f'number of major reads at flanking sites: {int(major_adj)}')
        print(f'number of minor reads at flanking sites: {int(minor_adj)}')
        print(f'number of major reads at focal sites: {int(major_foc)}')
        print(f'number of minor reads at focal sites: {int(minor_foc)}')
        print(f'err rate at flanking sites: {minor_adj/(minor_adj + major_adj):.6f}')
        print(f'err rate at focal sites: {minor_foc/(minor_foc + major_foc):.6f}')

    # finished reading bam file and we have made an estimate for genotyping error
    # now write a hdf5 file for read count at target sites
    l = len(pos)
    k = 1 # 1 bam = 1 sample
    gt = np.zeros((l, k, 2)) # dummy gt matrix, won't be used for inference
    dt = h5py.special_dtype(vlen=str)  # To have no problem with saving
    
    if len(outPath) != 0 and not outPath.endswith("/"):
        outPath += "/"
    bamFileName = os.path.basename(path2bam)
    hdf5Name = outPath + bamFileName[:bamFileName.find(".bam")] + ".hdf5"

    if os.path.exists(hdf5Name):  # Do a Deletion of existing File there
        os.remove(hdf5Name)
    
    # iid is the sample ID name in the hdf5 file
    if len(iid) == 0:
        iid = bamFileName[:bamFileName.find(".bam")]
    
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
        f_samples[:] = np.array([iid]).astype("S50")
        print(f'saving sample as {iid} in {hdf5Name}')

    # the second return value is the number of sites covered by at least 1 read
    err = minor_adj/(minor_adj + major_adj)
    numSitesCovered = np.sum(np.sum(np.sum(ad, axis=1), axis=1) > 0)
    if output:
        print(f'estimated genotyping error by flanking sites: {err:.6f}')
        print(f'number of sites covered by at least one read: {numSitesCovered}, fraction covered: {numSitesCovered/len(pos):.3f}')
        print(f'hdf5 file saved to {hdf5Name}')
    return err, numSitesCovered, hdf5Name


############################################### Experimental Function #############################################
######################################### Do as what Stephen's pileupcaller did ###################################

def mpileup2hdf5_damageAware(path2mpileup, refHDF5, iid="", s=-np.inf, e=np.inf, outPath="", output=True):
    f = h5py.File(refHDF5, 'r')
    pos = np.array(f['variants/POS'])
    ref = f['variants/REF']
    alt = f['variants/ALT']
    if len(alt.shape) == 2:
        alt = alt[:, 0].flatten()
    ref = np.array(ref).astype('str')
    alt = np.array(alt).astype('str')
    bases = np.array(['A', 'T', 'G', 'C'])
    subset1 = np.logical_and(pos >= s, pos <= e)
    subset2 = np.logical_and(np.in1d(ref, bases), np.in1d(alt, bases))
    print(f'exclude {np.sum(~subset1)} sites outside the specified region')
    print(f'exclude {np.sum(~subset2)} non-SNP sites')
    subset = np.logical_and(subset1, subset2)
    pos = pos[subset]
    rec = f['variants/MAP'][subset]
    ref = ref[subset]
    alt = alt[subset]

    assert(len(ref) == len(pos))
    assert(len(alt) == len(pos))
    assert(len(rec) == len(pos))
    # 1 here means there is just one sample
    # 2 : [ref read count, alt read count]
    ad = np.zeros((len(pos), 1, 2))

    major_adj = 0
    minor_adj = 0
    major_foc = 0
    minor_foc = 0
    l = len(pos)
    base2index = {'A':0, 'C':1, 'G':2, 'T':3}
    with open(path2mpileup) as f:
        for line in f:
            rc = np.zeros(4)
            contig, bp, refAllele, coverage, readbases, baseQ = list(zip(*zip_longest(line.strip().split(), range(6))))[0]
            if int(coverage) == 0:
                continue
            if refAllele == 'N':
                print(f'Cannot perform damage aware parsing, as the reference allele is not present in the mpileup file for {contig}: {bp}')
                sys.exit()
            insertion_index = readbases.find("+")
            insert = ""
            while insertion_index != -1:
                numInsertion = int(readbases[insertion_index+1])
                insert += readbases[insertion_index+2:insertion_index+2+numInsertion]
                insertion_index = readbases.find("+", insertion_index+2+numInsertion)
            
            rc[0] = readbases.count('A') + readbases.count('a') - insert.count('A') - insert.count('a')
            rc[1] = readbases.count('C') + readbases.count('c') - insert.count('C') - insert.count('c')
            rc[2] = readbases.count('G') + readbases.count('g') - insert.count('G') - insert.count('g')
            rc[3] = readbases.count('T') + readbases.count('t') - insert.count('T') - insert.count('t')
            rc[base2index[refAllele]] += readbases.count('.') + readbases.count(',') - insert.count('.') - insert.count(',')
            if refAllele == 'C':
                rc[1] -= readbases.count('C') + insert.count('C') # do not count any forward mapping C if the ref is C
                rc[3] -= readbases.count('T') + insert.count('T') # do not count any forward mapping T if the ref is C
            elif refAllele == 'G':
                rc[2] -= readbases.count('g') + insert.count('g') # do not count any reverse mapping G if the ref is G
                rc[0] -= readbases.count('a') + insert.count('a') # do not count any reverse mapping A if the ref is G
            coverage, bp = np.sum(rc), int(bp)
            i = np.searchsorted(pos, bp)
            if i < l and pos[i] == bp:
                # target sites
                # print(f'at target sites: {bp}, ref: {ref[i]}, alt: {alt[i]}')
                ad[i, 0, 0] = rc[base2index[ref[i]]]
                ad[i, 0, 1] = rc[base2index[alt[i]]]
                if coverage > 1:
                    major_foc += np.max(rc)
                    minor_foc += np.sum(rc) - np.max(rc)
            else:
                # sites flanking to the target sites
                if coverage > 1:
                    major_adj += np.max(rc)
                    minor_adj += np.sum(rc) - np.max(rc)

    if output:
        print(f'number of major reads at flanking sites: {int(major_adj)}')
        print(f'number of minor reads at flanking sites: {int(minor_adj)}')
        print(f'number of major reads at focal sites: {int(major_foc)}')
        print(f'number of minor reads at focal sites: {int(minor_foc)}')
        print(f'err rate at flanking sites: {minor_adj/(minor_adj + major_adj):.6f}')
        print(f'err rate at focal sites: {minor_foc/(minor_foc + major_foc):.6f}')

    # finished reading bam file and we have made an estimate for genotyping error
    # now write a hdf5 file for read count at target sites
    l = len(pos)
    k = 1 # 1 bam = 1 sample
    gt = np.zeros((l, k, 2)) # dummy gt matrix, won't be used for inference
    dt = h5py.special_dtype(vlen=str)  # To have no problem with saving
    
    if len(outPath) != 0 and not outPath.endswith("/"):
        outPath += "/"
    bamFileName = os.path.basename(path2mpileup)
    hdf5Name = outPath + bamFileName[:bamFileName.find(".mpileup")] + ".hdf5"

    if os.path.exists(hdf5Name):  # Do a Deletion of existing File there
        os.remove(hdf5Name)
    
    # iid is the sample ID name in the hdf5 file
    if len(iid) == 0:
        iid = bamFileName[:bamFileName.find(".mpileup")]
    
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
        f_pos[:] = pos
        f_gt[:] = gt
        f_samples[:] = np.array([iid]).astype("S50")
        print(f'saving sample as {iid} in {hdf5Name}')

    # the second return value is the number of sites covered by at least 1 read
    err = minor_adj/(minor_adj + major_adj)
    numSitesCovered = np.sum(np.sum(np.sum(ad, axis=1), axis=1) > 0)
    if output:
        print(f'estimated genotyping error by flanking sites: {err:.6f}')
        print(f'number of sites covered by at least one read: {numSitesCovered}, fraction covered: {numSitesCovered/len(pos):.3f}')
        print(f'hdf5 file saved to {hdf5Name}')
    return err, numSitesCovered, hdf5Name




def bamTable2hdf5(path2bamTable, refHDF5, iid="", s=-np.inf, e=np.inf, outPath="", output=True):
    f = h5py.File(refHDF5, 'r')
    pos = np.array(f['variants/POS'])
    ref = f['variants/REF']
    alt = f['variants/ALT']
    if len(alt.shape) == 2:
        alt = alt[:, 0].flatten()
    ref = np.array(ref).astype('str')
    alt = np.array(alt).astype('str')
    bases = np.array(['A', 'T', 'G', 'C'])
    subset1 = np.logical_and(pos >= s, pos <= e)
    subset2 = np.logical_and(np.in1d(ref, bases), np.in1d(alt, bases))
    print(f'exclude {np.sum(~subset1)} sites outside the specified region')
    print(f'exclude {np.sum(~subset2)} non-SNP sites')
    subset = np.logical_and(subset1, subset2)
    pos = pos[subset]
    rec = f['variants/MAP'][subset]
    ref = ref[subset]
    alt = alt[subset]
    # ref = f['variants/REF'][subset]
    
    # alt = f['variants/ALT']
    # is1D = len(alt.shape) == 1 # the 1240k hdf5's alt is a 2d array, while the new full 1000Genome hdf5 is 1d array
    # if is1D:
    #     alt = alt[subset]
    # else:
    #     alt = alt[subset, 0]

    assert(len(ref) == len(pos))
    assert(len(alt) == len(pos))
    assert(len(rec) == len(pos))
    # 1 here means there is just one sample
    # 2 : [ref read count, alt read count]
    ad = np.zeros((len(pos), 1, 2))

    major_adj = 0
    minor_adj = 0
    major_foc = 0
    minor_foc = 0
    l = len(pos)
    base2index = {'A':0, 'C':1, 'G':2, 'T':3}
    with open(path2bamTable) as f:
        for line in f:
            rc = np.zeros(4)
            contig, bp, A, C, G, T = line.strip().split()
            A, C, G, T = int(A), int(C), int(G), int(T)
            coverage = A+C+G+T
            if int(coverage) == 0:
                continue
            
            rc[0] = A
            rc[1] = C
            rc[2] = G
            rc[3] = T
            coverage, bp = int(coverage), int(bp)
            i = np.searchsorted(pos, bp)
            if i < l and pos[i] == bp:
                # target sites
                # print(f'at target sites: {bp}, ref: {ref[i]}, alt: {alt[i]}')
                ad[i, 0, 0] = rc[base2index[ref[i]]]
                ad[i, 0, 1] = rc[base2index[alt[i]]]
                if coverage > 1:
                    major_foc += np.max(rc)
                    minor_foc += np.sum(rc) - np.max(rc)
            else:
                # sites flanking to the target sites
                if coverage > 1:
                    major_adj += np.max(rc)
                    minor_adj += np.sum(rc) - np.max(rc)

    if output:
        print(f'number of major reads at flanking sites: {int(major_adj)}')
        print(f'number of minor reads at flanking sites: {int(minor_adj)}')
        print(f'number of major reads at focal sites: {int(major_foc)}')
        print(f'number of minor reads at focal sites: {int(minor_foc)}')
        print(f'err rate at flanking sites: {minor_adj/(minor_adj + major_adj):.6f}')
        print(f'err rate at focal sites: {minor_foc/(minor_foc + major_foc):.6f}')

    # finished reading bam file and we have made an estimate for genotyping error
    # now write a hdf5 file for read count at target sites
    l = len(pos)
    k = 1 # 1 bam = 1 sample
    gt = np.zeros((l, k, 2)) # dummy gt matrix, won't be used for inference
    dt = h5py.special_dtype(vlen=str)  # To have no problem with saving
    
    if len(outPath) != 0 and not outPath.endswith("/"):
        outPath += "/"
    bamFileName = os.path.basename(path2bamTable)
    hdf5Name = outPath + bamFileName[:bamFileName.find(".BamTable")] + ".hdf5"

    if os.path.exists(hdf5Name):  # Do a Deletion of existing File there
        os.remove(hdf5Name)
    
    # iid is the sample ID name in the hdf5 file
    if len(iid) == 0:
        iid = bamFileName[:bamFileName.find(".BamTable")]
    
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
        f_pos[:] = pos
        f_gt[:] = gt
        f_samples[:] = np.array([iid]).astype("S50")
        print(f'saving sample as {iid} in {hdf5Name}')

    # the second return value is the number of sites covered by at least 1 read
    err = minor_adj/(minor_adj + major_adj)
    numSitesCovered = np.sum(np.sum(np.sum(ad, axis=1), axis=1) > 0)
    if output:
        print(f'estimated genotyping error by flanking sites: {err:.6f}')
        print(f'number of sites covered by at least one read: {numSitesCovered}, fraction covered: {numSitesCovered/len(pos):.3f}')
        print(f'hdf5 file saved to {hdf5Name}')
    return err, numSitesCovered, hdf5Name