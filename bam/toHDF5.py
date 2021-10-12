import h5py
import os
import pysam
import time
import numpy as np

def mpileup2hdf5(path2mpileup, refHDF5, iid="", s=-np.inf, e=np.inf, outPath="", output=True):
    t1 = time.time()
    f = h5py.File(refHDF5, 'r')
    pos = np.array(f['variants/POS'])
    subset = np.logical_and(pos >= s, pos <= e)
    pos = pos[subset]
    rec = f['variants/MAP'][subset]
    ref = f['variants/REF'][subset]
    
    alt = f['variants/ALT']
    is1D = len(alt.shape) == 1 # the 1240k hdf5's alt is a 2d array, while the new full 1000Genome hdf5 is 1d array
    if is1D:
        alt = alt[subset]
    else:
        alt = alt[subset, 0]

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
            contig, bp, _, coverage, readbases, baseQ = line.strip().split()

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
                #print(f'at target sites: {bp}')
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
    print(f'finished reading mpileup file {path2mpileup}, takes {round(time.time()-t1, 3)}.')
    print(f'estimated genotyping error by flanking sites: {round(err, 6)}')
    print(f'number of sites covered by at least one read: {numSitesCovered}, fraction covered: {round(numSitesCovered/len(pos), 3)}')
    print(f'hdf5 file saved to {hdf5Name}')
    return err, numSitesCovered, hdf5Name

def bam2hdf5(path2bam, refHDF5, ch="X", iid="", minMapQual=30, minBaseQual=20, s=-np.inf, e=np.inf, outPath="", trim=0):
    f = h5py.File(refHDF5, 'r')
    pos = np.array(f['variants/POS'])
    subset = np.logical_and(pos >= s, pos <= e)
    pos = pos[subset]
    rec = f['variants/MAP'][subset]
    ref = f['variants/REF'][subset]
    
    alt = f['variants/ALT']
    is1D = len(alt.shape) == 1 # the 1240k hdf5's alt is a 2d array, while the new full 1000Genome hdf5 is 1d array
    if is1D:
        alt = alt[subset]
    else:
        alt = alt[subset, 0]

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
    return minor_adj/(minor_adj + major_adj), np.sum(np.sum(np.sum(ad, axis=1), axis=1) > 0), hdf5Name
