import numpy as np
import pandas as pd
import os
from hapsburg.PackagesSupport.pp_individual_roh_csvs import merge_called_blocks_custom
from hapsburg.preprocessing import PreProcessingHDF5, PreProcessingEigenstrat, PreProcessingEigenstratX
from hapsburg.PackagesSupport.loadEigenstrat.loadEigenstrat import load_eigenstrat


# a helper function used by all low-memory preprocessing classes
def extract_snps_hdf5_lowmem(h5, ids_ref, markers, ch, meta_path_ref, verbose=True, diploid=True):
        """Extract genotypes from h5 on ids and markers.
        If diploid, concatenate haplotypes along 0 axis.
        Extract indivuals first, and then subset to SNPs.
        Return 2D array [# haplotypes, # markers]"""
        # Important: Swap of Dimensions [loci<->individuals]
        nsample = len(ids_ref)
        ploidy = 2 if diploid else 1

        raw_gt = h5["calldata/GTbinary"][:, ids_ref, :ploidy]
        indices = markers // 8
        offset = markers % 8
        # get the offset bit from the uint8 integer at indices
        gts = np.packbits((raw_gt[indices, :, :] >> (7 - offset[:, None, None])) & 1, axis=0)
        gts = gts.reshape((-1, ploidy*nsample)).T
        # remove the second haplotype of males if ch='X'
        if ch == 'X' and diploid:
            meta_df = pd.read_csv(meta_path_ref, sep="\t")
            male_index = np.where(meta_df['sex'] == 'M')[0]
            male_index = np.nonzero(np.isin(ids_ref, male_index))[0] # reindex it based on ids_ref
            gts = gts[np.setdiff1d(np.arange(2*nsample), 2*male_index+1), :]

        if verbose:
            print(f"Extraction of {len(gts)} reference haplotypes at {len(markers)} sites complete")
        return gts

class PreProcessingHDF5_lowmem(PreProcessingHDF5):

    def load_data(self, iid="MA89", ch=6, start=-np.inf, end=np.inf):
        """Return Matrix of reference [k,l], Matrix of Individual Data [2,l],
        as well as linkage Map [l]"""

        if self.output == True:
            print(f"Loading Individual: {iid}")

        # Attach Part for the right Chromosome
        if self.h5_path1000g.endswith(".hdf5"):
            h5_path1000g = self.h5_path1000g
        else:
            h5_path1000g = self.h5_path1000g + str(ch) + ".hdf5"

        # Def Set the output folder:
        out_folder = self.set_output_folder(iid, ch)

        # Load and Merge the Data
        fs = self.load_h5(self.path_targets)
        f1000 = self.load_h5(h5_path1000g)
        i1, i2, flipped = self.merge_2hdf(fs, f1000, start, end)

        id_obs = self.get_index_iid(iid, fs, samples_field=self.samples_field)
        ids_ref, ids_con = self.get_ref_ids(f1000)

        markers = np.arange(0, len(i1))  # Which Markers to Slice out
        markers_obs = i1[markers]
        markers_ref = i2[markers]

        ### Load Target Dataset
        gts_ind = self.extract_snps_hdf5(
            fs, [id_obs], markers_obs, diploid=True, removeIncompleteHap=False)
        if self.readcounts:
            read_counts = self.extract_rc_hdf5(fs, id_obs, markers_obs)
        else:
            read_counts = []

        
        ### Flip target genotypes where flipped
        if self.flipstrand:
            if self.output:
                print(f"Flipping Ref/Alt Allele in target for {np.sum(flipped)} SNPs...")
            flip_idcs = flipped & (gts_ind[0,:]>=0) # Where Flip AND Genotype Data
            gts_ind[:,flip_idcs] = 1 - gts_ind[:,flip_idcs]
            if self.readcounts:
                read_counts[0,flipped], read_counts[1,flipped] = read_counts[1,flipped], read_counts[0,flipped]

        ### Load Reference Dataset
        gts = extract_snps_hdf5_lowmem(
            f1000, ids_ref, markers_ref, self.ch, self.meta_path_ref, verbose=self.output, diploid=self.diploid_ref)
        r_map = self.extract_rmap_hdf5(f1000, markers_ref)  # Extract LD Map
        pos = self.extract_rmap_hdf5(f1000, markers_ref, col="variants/POS")  # Extract Positions

        # get pCon: allele frequency from the specified contamination population
        gts_con = self.extract_snps_contaminationPop(f1000, ids_con, markers_ref)
        pCon = np.mean(gts_con, axis=0)

        # Do optional Processing Steps (based on boolean flags in class)
        # I need to pass the overhang to the optional_postprocessing function because it might be changed due to filtering to sites with data
        gts_ind, gts, r_map, pos, pCon, out_folder = self.optional_postprocessing(
            gts_ind, gts, r_map, pos, out_folder, pCon, read_counts)

        return gts_ind, gts, r_map, pos, pCon, out_folder

    def optional_postprocessing(self, gts_ind, gts, r_map, pos, out_folder, pCon, read_counts=[]):
        """Postprocessing steps of gts_ind, gts, r_map, and the folder,
        based on boolean fields of the class."""

        ### I commented out the following lines because the hdf5 file of the target individual only contains the markers with data
        if self.only_calls:
            called = self.markers_called(gts_ind, read_counts)
            gts_ind = gts_ind[:, called]
            r_map = r_map[called]
            pos = pos[called]
            pCon = pCon[called]
            called_indices = np.where(called)[0]
            indices = called_indices // 8
            offset = called_indices % 8
            gts = np.packbits((gts[:, indices] >> (7 - offset[None, :])) & 1, axis=1)
            
            if self.output:
                print(f"Subset to markers with data: {np.sum(called)} / {len(called)}")
                print(f"Fraction SNPs covered: {np.sum(called) / len(called):.4f}")
                
            if len(read_counts) > 0:
                read_counts = read_counts[:, called]

        if self.save == True:
            self.save_info(out_folder, r_map, pos,
                           gt_individual=gts_ind, read_counts=read_counts)

        if (self.readcounts == True) and len(read_counts) > 0:   # Switch to Readcount
            if self.output:
                print(f"Loading Readcounts...")
                print(f"Mean Readcount on markers with data: {np.mean(read_counts) * 2:.5f}")
            
            ### Downsample to target average coverage
            if not self.downsample:
                gts_ind = read_counts
            else:
                print('downsample readcounts data. This seems to make the emission model more robust to various sources of errors...')
                ################## downsample to ~1x data ################################################
                meanDepth = np.mean(np.sum(read_counts, axis=0))
                p = (1/meanDepth) * self.downsample
                sample_binom_ref = np.random.binomial(read_counts[0], p)
                sample_binom_alt = np.random.binomial(read_counts[1], p)
                gts_ind = np.zeros_like(read_counts)
                gts_ind[0] = sample_binom_ref
                gts_ind[1] = sample_binom_alt


        ### Shuffle Target Allele     
        if (self.random_allele == True) and (self.readcounts == False):     
            if self.output:
                print("Shuffling phase of target...")
            gts_ind = self.destroy_phase_func(gts_ind)
        
        return gts_ind, gts, r_map, pos, pCon, out_folder


class PreProcessingEigenstrat_lowmem(PreProcessingEigenstrat):
    def optional_postprocessing(self, gts_ind, gts, r_map, pos, out_folder, read_counts=[]):
        """Postprocessing steps of gts_ind, gts, r_map, and the folder,
        based on boolean fields of the class."""

        if self.only_calls:
            called = self.markers_called(gts_ind, read_counts)
            gts_ind = gts_ind[:, called]
            r_map = r_map[called]
            pos = pos[called]
            called_indices = np.where(called)[0]
            indices = called_indices // 8
            offset = called_indices % 8
            gts = np.packbits((gts[:, indices] >> (7 - offset[None, :])) & 1, axis=1)
            
            
            if self.output:
                print(f"Subset to markers with data: {np.shape(gts)[1]} / {len(called)}")
                print(f"Fraction SNPs covered: {np.shape(gts)[1] / len(called):.4f}")
                
            if len(read_counts) > 0:
                read_counts = read_counts[:, called]

        if self.save == True:
            self.save_info(out_folder, r_map, pos,
                           gt_individual=gts_ind, read_counts=read_counts)

        if (self.readcounts == True) and len(read_counts) > 0:   # Switch to Readcount
            if self.output:
                print(f"Loading Readcounts...")
                print(f"Mean Readcount on markers with data: {np.mean(read_counts) * 2:.5f}")
            gts_ind = read_counts
        
        ### Shuffle Target Allele     
        if (self.random_allele == True) and (self.readcounts == False):     
            if self.output == True:
                print("Shuffling phase of target...")
            gts_ind = self.destroy_phase_func(gts_ind)
        
        return gts_ind, gts, r_map, pos, out_folder

    def load_data(self, iid="MA89", ch=6):
        """Return Matrix of reference [k,l], Matrix of Individual Data [2,l],
        as well as linkage Map [l] and the output folder.
        Save the loaded data if self.save==True
        Various modifiers in class fields (check also PreProcessingHDF5)"""
        self.ch = ch  # To remember Chromosome
        if self.output:
            print(f"Loading Individual: {iid}")

        out_folder = self.set_output_folder(iid, ch)

        ### Load Reference HDF5 object
        h5_path1000g = self.get_1000G_path(self.h5_path1000g, ch)
        f1000 = self.load_h5(h5_path1000g)

        ### Load Eigenstrat object
        es = load_eigenstrat(base_path=self.path_targets, 
                             sep=self.sep, packed=self.packed)

        markers_obs, markers_ref, flipped = self.merge_es_hdf5(es, f1000)
        
        id_obs = self.es_get_index_iid(es, iid)
        ids_ref = self.get_ref_ids(f1000) # the second placeholder is here because my modified get_ref_ids 

        r_map = self.extract_rmap_hdf5(f1000, markers_ref)  # Extract LD Map
        pos = self.extract_rmap_hdf5(f1000, markers_ref, col="variants/POS")  # Extract Positions

        if len(self.segM) > 0:
            r_map, markers_obs, markers_ref = self.get_segment(
                r_map, markers_obs, markers_ref)
            
        ######################################################
        ### Extract genotypes: 1) For target 2) For Ref (where target covered)
        ### Extraction for target (no RC here)
        gts_ind = self.extract_snps_es(es, id_obs, markers_obs)
        
        ### Produce empty readcounts
        read_counts = []
        if self.readcounts:
            read_counts = self.to_read_counts(gts_ind)
        
        ### Subset to markers called in target
        if self.only_calls:
            called = self.markers_called(gts_ind, read_counts)
            gts_ind = gts_ind[:, called]
            r_map = r_map[called]
            pos = pos[called]
            flipped=flipped[called]   # Where to flip markers
            markers_ref = markers_ref[called] # Subset to IIDs to actually extract from Ref
             
            if len(read_counts) > 0:
                read_counts = read_counts[:, called]
            if self.output:
                print(f"Reduced to markers with data: {np.sum(called)} / {len(called)}")
                print(f"Fraction SNPs covered: {np.sum(called) / len(called):.4f}")
            self.only_calls=False # To avoid doing it in final step again
        
        #### Extract the reference genotype
        gts = extract_snps_hdf5_lowmem(f1000, ids_ref, markers_ref, self.ch, self.meta_path_ref, verbose=self.output, diploid=self.diploid_ref)
                
        ### Flip target genotypes where flipped
        if self.flipstrand:
            if self.output:
                print(f"Flipping Ref/Alt Alleles in target for {np.sum(flipped)} SNPs...")
            flip_idcs = flipped & (gts_ind[0,:]>=0) # Where Flip AND Genotype Data
            gts_ind[:,flip_idcs] = 1 - gts_ind[:,flip_idcs]
            
            if len(read_counts)>0:
                read_counts[0,flipped], read_counts[1,flipped] = read_counts[1,flipped], read_counts[0,flipped]

        ### Do optional Processing Steps (only covered, destroy phase, save)
        gts_ind, gts, r_map, pos, out_folder = self.optional_postprocessing(
            gts_ind, gts, r_map, pos, out_folder, read_counts)

        return gts_ind, gts, r_map, pos, out_folder


class PreProcessingEigenstratX_lowmem(PreProcessingEigenstrat_lowmem, PreProcessingEigenstratX):
    def set_output_folder(self, iid, ch="X"):
        return PreProcessingEigenstratX.set_output_folder(self, iid, ch)
    
    def get_1000G_path(self, h5_path1000g, ch="X"):
        return PreProcessingEigenstratX.get_1000G_path(self, h5_path1000g, ch)
    
    def es_get_index_iid(self, es, iid):
        return PreProcessingEigenstratX.es_get_index_iid(self, es, iid)
    
    def extract_snps_es(self, es, id, markers):
        return PreProcessingEigenstratX.extract_snps_es(self, es, id, markers)
    
    def load_data(self, iid="MA89", ch="X"):
        return PreProcessingEigenstrat_lowmem.load_data(self, iid, ch)
    
    # def optional_postprocessing(self, gts_ind, gts, r_map, pos, out_folder, read_counts=[]):
    #     return super(PreProcessingEigenstrat_lowmem, self).optional_postprocessing(gts_ind, gts, r_map, pos, out_folder, read_counts)


def load_preprocessing_lowmem(p_model="Eigenstrat", conPop=[], save=True, output=True):
    """Factory method to load the Transition Model.
    Return"""

    if p_model == "Eigenstrat":   ### Detects automatically what mode to use
        p_obj = PreProcessingEigenstrat_lowmem(save=save, output=output,
                                        packed=-1, sep=r"\s+") # -1: Unknown 
    elif p_model == "EigenstratX":
        p_obj = PreProcessingEigenstratX_lowmem(save=save, output=output,
                                         packed=-1, sep=r"\s+")
        
    elif (p_model == "HDF5") or (p_model == "SardHDF5"):
        p_obj = PreProcessingHDF5_lowmem(conPop, save=save, output=output)
        
    else:
        raise NotImplementedError(f"Preprocessing Model string {p_model} not found.")

    return p_obj