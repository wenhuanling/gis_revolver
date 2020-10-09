library(tidyverse)
library(revolver)
library(vcfR)
library(maftools)
library(parallel)
library(evoverse.datasets)

# the location of the loci file
loci_location = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S03_PyClone/output/output_20201001/S02_PyClone/random_1500mutations/pyclone_results/output"
# the location of the maf 
maf_location = "/mnt/projects/zhaiww1/planet/data_storage/DNA_data/hg38/funcotator_results"

# to load the masterfile for reading
master_file = read_csv("/mnt/projects/zhaiww1/liver_cancer_data/shared_group_folder/masterfiles/hg38/LiverTCR_HG38_masterfile_DNA_Logisitic_current.csv")
# the list of patient with normal tmb
normal_tmb_patient_master_list = read_tsv("/mnt/projects/zhaiww1/liver_cancer_data/shared_group_folder/masterfiles/hg38/non-problematic-DNA-pat.list", col_names = F)
# filter out all the normal tmb patients, simplifying the data frame
master_file = master_file %>% filter(Unified_ID %in% normal_tmb_patient_master_list$X1) %>% select(Unified_ID, DNA_lib)
# load the driver gene list
driver_gene_list = scan("/mnt/projects/kayaan/computational_genomics/liver_integrative_analysis/Driver_conclusion/concensus_1349/withkoreanindels/concensus_1349_koreanindels.txt", "cha")

# find all the loci files in the location
loci_list = data.frame(fn_dna=list.files(loci_location, pattern = "loci.tsv",recursive = T)) %>%
  separate(fn_dna, c("dna_lib", "tables", "loci"), sep = "/", remove = F)
# combine and get the unified id according to the respective files
pat_id = cbind(loci_list, master_file[match(loci_list$dna_lib, master_file$DNA_lib),1])
pat_id = pat_id %>% select(Unified_ID)

ccf_list = list()
# to load all the ccf information/loci file into 1 list using function
ccf_list = mclapply(pat_id$Unified_ID, revolver_input, mc.cores = 10)

ccf_list = Reduce(rbind, ccf_list)
write_tsv(ccf_list, "/mnt/projects/whling/cbp/Planet/revolver/input_for_revolver.tsv")

my_cohort = revolver_cohort(revolver_input, CCF_parser = revolver::CCF_parser, ONLY.DRIVER = T, 
                            MIN.CLUSTER.SIZE = 0, annotation = "HG38 cohort")
# to find out the duplicated driver genes and variantID
my_cohort_dup_variant = revolver:::get_duplicates(my_cohort) %>% dplyr::select(variantID) %>% unique()
my_cohort_drivers_dups = Stats_drivers(my_cohort) %>% filter(N_tot == 1) %>% pull(variantID)
# to remove the duplicated driver events
my_cohort_no_dup = revolver:::remove_drivers(my_cohort, my_cohort_drivers_dups) 
# to remove duplicated variantID
#my_cohort_no_dup = revolver:::remove_(my_cohort_no_dup, my_cohort_drivers_dups$variantID)
#### to fit the revolver cohort for plotting of various plots ####
# compute the number of mutation trees in the cohort
my_cohort_mutation_tree = compute_mutation_trees(my_cohort)
my_cohort_clone_tree = compute_clone_trees(my_cohort_no_dup) 
# fit the cohort_mutation tree
my_cohort_fit = revolver_fit(my_cohort_mutation_tree)
# cluster the fitted cohort
my_cohort_cluster = revolver_cluster(my_cohort_fit)
# compute jackknife with cluster fit
my_cohort_jackknife = revolver_jackknife(my_cohort_cluster)

# function to load all the ccf information into 1 dataframe
revolver_input = function(pat_id){
  # filtered out only the input patient id
  master_file_filtered = master_file %>% filter(Unified_ID == pat_id)
  # find all the maf files
  sv_vcf_list_tmp <- data.frame(fn_dna=list.files(maf_location, pattern = ".MUTECT2.maf",recursive = T)) %>% 
    separate(fn_dna, c("dna_lib", "mutect2", "maf"), sep = "\\.", remove = F) %>% filter(dna_lib %in% master_file_filtered$DNA_lib)
  maf_list = list()
  maf_list = lapply(setNames(sv_vcf_list_tmp$fn_dna %>% as.character(), sv_vcf_list_tmp$dna_lib), function(name){
    print(name)
    maf = read_tsv(paste0(maf_location, "/", name), comment = "#") %>% 
      unite(chr_pos, c("Chromosome", "Start_Position"), sep = ":", remove = F)
    #maf_no_col = ncol(maf)
    name = gsub(".MUTECT2.maf", "", name)
    #return(maf_col = data.frame(pat_id, name, maf_no_col))
    return(maf %>% mutate(sample_name = name) %>% select(Hugo_Symbol, chr_pos, sample_name))
  })
  #return(maf_list)
  print("rbind maflist")
  maf_list = Reduce(rbind,maf_list)
  # list out all the loci
  print("load list of loci")
  loci_list_tmp = loci_list %>% filter(dna_lib %in% master_file_filtered$DNA_lib)
  # # read the loci files
  print("read loci files")
  print(loci_list_tmp$fn_dna)
  loci = read_tsv(paste0(loci_location, "/", loci_list_tmp$fn_dna)) %>%
    separate(mutation_id, c("unknown", "chr_pos"), sep = "_", remove = F)
  # # grab all the dna libs (unique)
  loci_dna = unique(loci$sample_id)
  # # cbind the hugo_symbol to loci dataframe
  loci = cbind(loci, maf_list[match(loci$chr_pos, maf_list$chr_pos), 1]) %>% unite(ccf_pre, c("sample_id", "cellular_prevalence"), sep = ":", remove = F)
  loci = loci %>% dplyr::select(-cellular_prevalence_std, -variant_allele_frequency, -cellular_prevalence) %>% spread(sample_id,ccf_pre) %>%
    unite(CCF, loci_dna, sep = ";") %>% mutate(patientID = pat_id)
  # # to filter whether the samples are is.clonal (normally 0 are clonal clone from cluster, highest mean and size) and is.driver
  loci = loci %>% mutate(is.clonal = if_else(loci$cluster_id == 0, T, F)) %>%
    mutate(is.driver = if_else(Hugo_Symbol %in% driver_gene_list, T, F)) %>%
    unite(Misc, c("patientID", "chr_pos"),sep = ":", remove = F) %>% rename(variantID = Hugo_Symbol) %>% mutate(cluster = as.character(cluster_id)) %>%
    select(Misc, patientID, variantID, cluster, is.driver, is.clonal, CCF)
  # return the final loci
  return(loci)
}

#revolver_input = read_tsv("/media/GIS-Aquila-Cluster/mnt/projects/whling/cbp/Planet/revolver/input_for_revolver.tsv")
#revolver_input = revolver_input %>% mutate(cluster = as.character(cluster))
# to form the revolver cohort for the analysis
