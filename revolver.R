library(tidyverse)
library(revolver)
library(vcfR)
library(maftools)
library(parallel)
library(lubridate)
library(evoverse.datasets)

rm(list=ls())

#### all the global parameters ####
# the location of the loci.tsv files
loci_location = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S03_PyClone/output/output_20201001/S02_PyClone/random_1500mutations/pyclone_results/output"
# the location of the MAF files
maf_location = "/mnt/projects/zhaiww1/planet/data_storage/DNA_data/hg38/funcotator_results"
# set today's date
today_date = today()

# load the driver gene list
driver_gene_list = scan("/mnt/projects/kayaan/computational_genomics/liver_integrative_analysis/Driver_conclusion/concensus_1349/withkoreanindels/concensus_1349_koreanindels.txt", "cha")

# the list of patient with normal TMB
normal_tmb_patient_master_list = read_tsv("/mnt/projects/zhaiww1/liver_cancer_data/shared_group_folder/masterfiles/hg38/non-problematic-DNA-pat.list", col_names = F)
# non_normal_tmb = read_tsv("/mnt/projects/zhaiww1/planet/data_storage/DNA_data/hg38/work_folder/Scripts/masterfile/LiverTCR_HG38_master_file_tmb_200_samples.tsv")

# to load the master file (current version) and all the normal tmb patients, 
# simplifying the data frame to only the Unified ID and DNA libs for easy reference
master_file = read_csv("/mnt/projects/zhaiww1/liver_cancer_data/shared_group_folder/masterfiles/hg38/LiverTCR_HG38_masterfile_DNA_Logisitic_current.csv") %>%
  filter(Unified_ID %in% normal_tmb_patient_master_list$X1) %>% select(Unified_ID, DNA_lib)
# master_file = read_csv("/mnt/projects/zhaiww1/liver_cancer_data/shared_group_folder/masterfiles/hg38/LiverTCR_HG38_masterfile_DNA_Logisitic_current.csv") %>% 
#   filter(!DNA_lib %in% non_normal_tmb$DNA_lib) %>% select(Unified_ID, DNA_lib)

##### loci.tsv files ####
# find all the loci files in the location
loci_list = data.frame(fn_dna=list.files(loci_location, pattern = "loci.tsv",recursive = T)) %>%
  separate(fn_dna, c("dna_lib", "tables", "loci"), sep = "/", remove = F)
# combine and get the unified id according to the loci_list
pat_id = cbind(loci_list, master_file[match(loci_list$dna_lib, master_file$DNA_lib),1]) %>% select(Unified_ID)
# testing out which patient is givine problem
#error_pat = c("B023","C018","C019","A011","HEP270","B024","PB_95613","PB_01647","PB_67021","F008","F009" ,"HEP241" ,"HEP235")
# problematic patient
error_pat = "C018"
# filter out the problematic patient
pat_id = pat_id %>% filter(!Unified_ID %in% error_pat)

##### Main section of code for generating the revolver input ####
# create empty list 
ccf_list = list()
# to load all the ccf information/loci file into 1 list using function: revolver_input, using mclapply using 10 cores
ccf_list = mclapply(pat_id$Unified_ID, revolver_input, mc.cores = 10)
# combine the list into 1 single dataframe
ccf_list = Reduce(rbind, ccf_list)
# write into a tsv file
write_tsv(ccf_list, paste0("/mnt/projects/whling/cbp/Planet/revolver/input_for_revolver_fix_",today_date,".tsv"))

summary(ccf_list)

##### the main revolver function ####
# to create the revolver cohort
# options: Only driver genes set to TRUE
# min cluster size: set to 0
my_cohort = revolver_cohort(ccf_list, CCF_parser = revolver::CCF_parser, ONLY.DRIVER = T, 
                            MIN.CLUSTER.SIZE = 0, annotation = "HG38 cohort normal TMB")
# to find out the duplicated driver genes and variantID
my_cohort_variant_dups = revolver:::get_duplicates(my_cohort) %>% dplyr::select(patientID) %>% unique()
my_cohort_driver_dups = Stats_drivers(my_cohort) %>% filter(N_tot == 1) %>% pull(variantID)
# to remove the duplicated driver events
my_cohort_no_dup = revolver:::remove_drivers(my_cohort, my_cohort_driver_dups, check = T) 
# to remove duplicated variantID
my_cohort_no_dup = revolver:::remove_patients(my_cohort_no_dup, my_cohort_variant_dups$patientID)

#### to fit the revolver cohort for plotting of various plots ####
# compute the number of mutation/clone trees in the cohort
# currently the clone tree function seems to prompt error when running it
my_cohort_mutation_tree = compute_mutation_trees(my_cohort)
#my_cohort_clone_tree = compute_clone_trees(my_cohort_no_dup) 
# fit the cohort_mutation tree
my_cohort_fit = revolver_fit(my_cohort_mutation_tree)
# find the cluster fit using the input from my_cohort_fit
my_cohort_cluster = revolver_cluster(my_cohort_fit)
# compute jackknife fit with my_cohort_cluster
my_cohort_jackknife = revolver_jackknife(my_cohort_cluster)

#### generate and print out the figures ####
# get the list of patient ID for the plotting
pat_id_new = my_cohort_no_dup$patients
# generate and print out the figures
print_pdf(cohort = my_cohort_no_dup, phylo_tree = my_cohort_phylo_tree, pat_list = pat_id_new, 
          location = "/mnt/projects/whling/cbp/Planet/revolver/hg38_revolver_plots.pdf", fit_tree = my_cohort_fit,
          cluster = my_cohort_cluster, jackknife = my_cohort_jackknife)

#### function to load all the ccf information into 1 dataframe ####
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

#### print pdf function ####
print_pdf = function(cohort,phylo_tree, pat_list, location, fit_tree, cluster, jackknife) {
  pdf(location, width = 16.53, height = 11.69 )
  print(cohort)
  overall = plot(cohort)
  print(overall)
  driver_onccurrence = plot_drivers_occurrence(fit_tree)
  print(driver_onccurrence)
  driver_graph = plot_drivers_graph(fit_tree)
  print(driver_graph)
  driver_clonality = plot_drivers_clonality(fit_tree)
  print(driver_clonality)
  penalty = plot_penalty(fit_tree)
  print(penalty)
  det_index = plot_DET_index(fit_tree)
  print(det_index)
  cluster_cohort = plot_clusters(cluster)
  print(cluster_cohort)
  dendrogram = plot_dendrogram(jackknife)
  print(dendrogram)
  jackknife_coclustering = plot_jackknife_coclustering(jackknife)
  print(jackknife_coclustering)
  jackknife_cluster = plot_jackknife_cluster_stability(jackknife)
  print(jackknife_cluster)
  jackknife_stable = plot_jackknife_trajectories_stability(jackknife)
  print(jackknife_stable)
  for (i in pat_list) {
    print(i)
    p_tree = plot_patient_trees(phylo_tree, i)
    print(p_tree)
    p_data = plot_patient_data(phylo_tree, i)
    print(p_data)
  }
  dev.off()
}

