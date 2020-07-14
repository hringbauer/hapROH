"""
Functions to save hdf5, and transfrom Eigenstrat to hdf5
@ Author: Harald Ringbauer, 2019, All rights reserved
"""
import numpy as np
import pandas as pd

import os as os
import h5py  # Python Package to do the HDF5.
from hapsburg.PackagesSupport.loadEigenstrat.loadEigenstrat import load_eigenstrat

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