source(here::here("setup.R"))
source(here::here("src", "data_preprocessing", "combine_peptides_mq.R"))
source(here::here("src", "data_preprocessing", "combine_peptides_msf.R"))
source(here::here("src", "data_preprocessing", "combine_peptides_casanovo.R"))
source(here::here("src", "data_preprocessing", "additional_preprocessing_casanovo.R"))
source(here::here("src", "data_preprocessing", "preprocessing_functions.R"))


# Load metadata
mabs_seq <- readr::read_tsv(here("metadata", "mabs_cdr3_variable_full.tsv"))
sample_info <- readr::read_tsv(here::here("metadata", "detailed_sample_description.tsv"))
rawfiles_info <- readr::read_tsv(here::here("metadata", "rawfiles_description_with_blanks.tsv"), guess_max = 5e3)
imgt_genes <- readr::read_tsv(here::here("metadata", "imgt_genes.tsv"))
shared_mab_subseq <- readr::read_tsv(here::here("metadata", "mabs_shared_subseq.tsv"))



# Merge raw peptides

mq_peptides_raw <- combine_peptides_mq(here("data", "raw", "MQ"))
msf_peptides_raw <- combine_peptides_msf(here("data", "raw", "MSF"))
casanovo_peptides_raw <- combine_peptides_casanovo(here("data", "raw", "casanovo")) %>%
  additional_preprocessing_casanovo(min_search_engine_score = 0.8, min_peptide_length = 7) %>%
  rename(intensity = search_engine_score)

all_peptides_raw <- rbind(mq_peptides_raw, msf_peptides_raw, casanovo_peptides_raw)



# Preprocessing 

all_mab_peptides_preprocessed <- all_peptides_raw %>%
  keep_peptides_mapped_to_mabs(mabs_seq=mabs_seq) %>%
  annotate_contaminations(sample_info=sample_info, rawfiles_info=rawfiles_info) %>%
  annotate_cdr3_peptides(mabs_seq=mabs_seq) %>%
  annotate_imgt_peptides(imgt_genes=imgt_genes) %>%
  annotate_shared_peptides(shared_mab_subseq=shared_mab_subseq) %>%
  annotate_carryover_peptides(rawfiles_info=rawfiles_info) %>%
  select(Sample, Rawfilenumber, Sequence, Protease, run, replicate, intensity, 
         concentration, has_blood, tool, match_ig_type, is_cdr3_related, cdr3_overlap, 
         is_contamination, is_imgt_peptide, is_shared_peptide, is_carryover_peptide)



# Save preprocessed data
  
# all peptides
write_tsv(all_mab_peptides_preprocessed, 
          here::here("data", "processed", "peptides_all_TP_FP_with_blanks.tsv"))

# all peptides which passed the preprocessing
all_mab_peptides_preprocessed %>%
  filter(!(is_imgt_peptide | is_shared_peptide | is_carryover_peptide)) %>%
  filter(Sample!= "blank") %>%
  mutate(Sample = as.integer(Sample)) %>%
  write_tsv(here::here("data", "processed", "peptides_preprocessed_TP_FP.tsv"))

# TP peptides which passed the preprocessing
all_mab_peptides_preprocessed %>%
  filter(!(is_imgt_peptide | is_shared_peptide | is_carryover_peptide)) %>%
  filter(!is_contamination & (Sample!= "blank")) %>%
  mutate(Sample = as.integer(Sample)) %>%
  write_tsv(here::here("data", "processed", "peptides_preprocessed_TP.tsv"))
  
