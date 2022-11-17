import sys
sys.path.insert(0, '/mnt/archgen/users/yilei/tools/hapROH/package')
from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind  # Need this import
indname=sys.argv[1]
hapsb_ind(iid=indname, chs=range(3,4), processes=1, 
          path_targets='/home/xiaowen_jia/LAI/genotype/NG10', 
          h5_path1000g='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/maf5_auto/maf5_chr', 
          meta_path_ref='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv', 
          folder_out='', prefix_out='', 
          e_model="haploid", p_model="Eigenstrat",
          random_allele=True, readcounts=False,
          delete=False, logfile=True, combine=True)