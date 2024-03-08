library(stringr)
library(tidyverse)
library(tibble)
library(pbapply)
library(pbmcapply)
library(dplyr)
library(DECIPHER)
library(bioseq)
library(parallel)
library(doParallel)
library(foreach)
library(rlist)
library(jsonlite)
library(Biostrings)
library(data.table)
library(future)
library(future.apply)

#the purpose of this is to build the training data using the strains below.

iteration_number <- "167"

numCores <- detectCores()

cl <- makeCluster(numCores)

# Read in the DNA sequences from the fasta file

wgyale_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/wigglesworthia_glossinidia_yale.txt")
sf236_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/spiroplasma_floricola_236.txt")
amyn25_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/acidianus_manzaensis_yn25.txt")
psdsm12145_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/pseudopedobacter_saltans_dsm12145.txt")
cmic167_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/caldivirga_maquilingensis_ic167.txt")
sf2a2457t_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/shigella_flexneri_2a_2457t.txt")
lff6_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/lactobacillus_fermentum_f6.txt")
gkjs1_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/gloeobacter_kilaueensis_js1.txt")
obtav5_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/opitutaceae_bacterium_tav5.txt")
rdkctc42856_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/roseateles_depolymerans_kctc42856.txt")
afdsm7358_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/actinoplanes_friuliensis_dsm7358.txt")
godsm43160_fasta <- readDNAStringSet("D:/RStuff/codon_overlap/bacterial_genes/geodermatophilus_obscurus_dsm43160.txt")


# Create an empty list to store the genes
wgyale_gene_list <- list()
sf236_gene_list <- list()
amyn25_gene_list <- list()
psdsm12145_gene_list <- list()
cmic167_gene_list <- list()
sf2a2457t_gene_list <- list()
lff6_gene_list <- list()
gkjs1_gene_list <- list()
obtav5_gene_list <- list()
rdkctc42856_gene_list <- list()
afdsm7358_gene_list <- list()
godsm43160_gene_list <- list()


# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(wgyale_fasta)) {
  gene_seq <- wgyale_fasta[[i]]
  wgyale_gene_list[[i]] <- data.frame(wgyale_genes = as.character(gene_seq))
}

# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(sf236_fasta)) {
  gene_seq <- sf236_fasta[[i]]
  sf236_gene_list[[i]] <- data.frame(sf236_genes = as.character(gene_seq))
}

# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(amyn25_fasta)) {
  gene_seq <- amyn25_fasta[[i]]
  amyn25_gene_list[[i]] <- data.frame(amyn25_genes = as.character(gene_seq))
}

# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(psdsm12145_fasta)) {
  gene_seq <- psdsm12145_fasta[[i]]
  psdsm12145_gene_list[[i]] <- data.frame(psdsm12145_genes = as.character(gene_seq))
}

# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(cmic167_fasta)) {
  gene_seq <- cmic167_fasta[[i]]
  cmic167_gene_list[[i]] <- data.frame(cmic167_genes = as.character(gene_seq))
}

# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(sf2a2457t_fasta)) {
  gene_seq <- sf2a2457t_fasta[[i]]
  sf2a2457t_gene_list[[i]] <- data.frame(sf2a2457t_genes = as.character(gene_seq))
}

# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(lff6_fasta)) {
  gene_seq <- lff6_fasta[[i]]
  lff6_gene_list[[i]] <- data.frame(lff6_genes = as.character(gene_seq))
}

# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(gkjs1_fasta)) {
  gene_seq <- gkjs1_fasta[[i]]
  gkjs1_gene_list[[i]] <- data.frame(gkjs1_genes = as.character(gene_seq))
}

# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(obtav5_fasta)) {
  gene_seq <- obtav5_fasta[[i]]
  obtav5_gene_list[[i]] <- data.frame(obtav5_genes = as.character(gene_seq))
}

# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(rdkctc42856_fasta)) {
  gene_seq <- rdkctc42856_fasta[[i]]
  rdkctc42856_gene_list[[i]] <- data.frame(rdkctc42856_genes = as.character(gene_seq))
}

# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(afdsm7358_fasta)) {
  gene_seq <- afdsm7358_fasta[[i]]
  afdsm7358_gene_list[[i]] <- data.frame(afdsm7358_genes = as.character(gene_seq))
}

# Loop through each sequence in the fasta file and extract the gene
for (i in 1:length(godsm43160_fasta)) {
  gene_seq <- godsm43160_fasta[[i]]
  godsm43160_gene_list[[i]] <- data.frame(godsm43160_genes = as.character(gene_seq))
}



# Combine the data frames into a single data frame
wgyale_gene_df <- do.call(rbind, wgyale_gene_list)

wgyale_gene_df$wgyale_genes <- as.character(wgyale_gene_df$wgyale_genes)

sf236_gene_df <- do.call(rbind, sf236_gene_list)

sf236_gene_df$wgyale_genes <- as.character(sf236_gene_df$sf236_genes)

amyn25_gene_df <- do.call(rbind, amyn25_gene_list)

amyn25_gene_df$amyn25_genes <- as.character(amyn25_gene_df$amyn25_genes)

psdsm12145_gene_df <- do.call(rbind, psdsm12145_gene_list)

psdsm12145_gene_df$psdsm12145_genes <- as.character(psdsm12145_gene_df$psdsm12145_genes)

cmic167_gene_df <- do.call(rbind, cmic167_gene_list)

cmic167_gene_df$cmic167_genes <- as.character(cmic167_gene_df$cmic167_genes)

sf2a2457t_gene_df <- do.call(rbind, sf2a2457t_gene_list)

sf2a2457t_gene_df$sf2a2457t_genes <- as.character(sf2a2457t_gene_df$sf2a2457t_genes)

lff6_gene_df <- do.call(rbind, lff6_gene_list)

lff6_gene_df$lff6_genes <- as.character(lff6_gene_df$lff6_genes)

gkjs1_gene_df <- do.call(rbind, gkjs1_gene_list)

gkjs1_gene_df$gkjs1_genes <- as.character(gkjs1_gene_df$gkjs1_genes)

obtav5_gene_df <- do.call(rbind, obtav5_gene_list)

obtav5_gene_df$obtav5_genes <- as.character(obtav5_gene_df$obtav5_genes)

rdkctc42856_gene_df <- do.call(rbind, rdkctc42856_gene_list)

rdkctc42856_gene_df$rdkctc42856_genes <- as.character(rdkctc42856_gene_df$rdkctc42856_genes)

afdsm7358_gene_df <- do.call(rbind, afdsm7358_gene_list)

afdsm7358_gene_df$afdsm7358_genes <- as.character(afdsm7358_gene_df$afdsm7358_genes)

godsm43160_gene_df <- do.call(rbind, godsm43160_gene_list)

godsm43160_gene_df$godsm43160_genes <- as.character(godsm43160_gene_df$godsm43160_genes)


# Filter genes that are >= 300 in length
wgyale_gene_filtered <- subset(wgyale_gene_df, nchar(wgyale_genes) >= 350)

wgyale_df <- data.frame(genes = wgyale_gene_filtered$wgyale_genes)

wgyale_df <- wgyale_df %>% filter(!str_detect(genes, "[MN]"))

sf236_gene_filtered <- subset(sf236_gene_df, nchar(sf236_genes) >= 350)

sf236_df <- data.frame(genes = sf236_gene_filtered$sf236_genes)

sf236_df <- sf236_df %>% filter(!str_detect(genes, "[MN]"))

amyn25_gene_filtered <- subset(amyn25_gene_df, nchar(amyn25_genes) >= 350)

amyn25_df <- data.frame(genes = amyn25_gene_filtered$amyn25_genes)

amyn25_df <- amyn25_df %>% filter(!str_detect(genes, "[MN]"))

psdsm12145_gene_filtered <- subset(psdsm12145_gene_df, nchar(psdsm12145_genes) >= 350)

psdsm12145_df <- data.frame(genes = psdsm12145_gene_filtered$psdsm12145_genes)

psdsm12145_df <- psdsm12145_df %>% filter(!str_detect(genes, "[MN]"))

cmic167_gene_filtered <- subset(cmic167_gene_df, nchar(cmic167_genes) >= 350)

cmic167_df <- data.frame(genes = cmic167_gene_filtered$cmic167_genes)

cmic167_df <- cmic167_df %>% filter(!str_detect(genes, "[MN]"))

sf2a2457t_gene_filtered <- subset(sf2a2457t_gene_df, nchar(sf2a2457t_genes) >= 350)

sf2a2457t_df <- data.frame(genes = sf2a2457t_gene_filtered$sf2a2457t_genes)

sf2a2457t_df <- sf2a2457t_df %>% filter(!str_detect(genes, "[MN]"))

lff6_gene_filtered <- subset(lff6_gene_df, nchar(lff6_genes) >= 350)

lff6_df <- data.frame(genes = lff6_gene_filtered$lff6_genes)

lff6_df <- lff6_df %>% filter(!str_detect(genes, "[MN]"))

gkjs1_gene_filtered <- subset(gkjs1_gene_df, nchar(gkjs1_genes) >= 350)

gkjs1_df <- data.frame(genes = gkjs1_gene_filtered$gkjs1_genes)

gkjs1_df <- gkjs1_df %>% filter(!str_detect(genes, "[MN]"))

obtav5_gene_filtered <- subset(obtav5_gene_df, nchar(obtav5_genes) >= 350)

obtav5_df <- data.frame(genes = obtav5_gene_filtered$obtav5_genes)

obtav5_df <- obtav5_df %>% filter(!str_detect(genes, "[MN]"))

rdkctc42856_gene_filtered <- subset(rdkctc42856_gene_df, nchar(rdkctc42856_genes) >= 350)

rdkctc42856_df <- data.frame(genes = rdkctc42856_gene_filtered$rdkctc42856_genes)

rdkctc42856_df <- rdkctc42856_df %>% filter(!str_detect(genes, "[MN]"))

afdsm7358_gene_filtered <- subset(afdsm7358_gene_df, nchar(afdsm7358_genes) >= 350)

afdsm7358_df <- data.frame(genes = afdsm7358_gene_filtered$afdsm7358_genes)

afdsm7358_df <- afdsm7358_df %>% filter(!str_detect(genes, "[MN]"))

godsm43160_gene_filtered <- subset(godsm43160_gene_df, nchar(godsm43160_genes) >= 350)

godsm43160_df <- data.frame(genes = godsm43160_gene_filtered$godsm43160_genes)

godsm43160_df <- godsm43160_df %>% filter(!str_detect(genes, "[MN]"))


# make combined dataframe with all genes
# change rbind to be either one strain df, or combine multiple df files.

combined_genes <- rbind(godsm43160_df)

# Set the number of rows you want to select. This will define how large dataset is prior to any downstream filtering based on group size.
num_rows <- 450000

# if using equal counts of overlaps, set the desired count per overlap group.
# (note, the min group must be at least this much, or it won't work). There are 98 groups, before adding in neg controls. So calculate accordingly.
target_num <- 120

# Use sample to select a random subset of row indices
selected_indices <- sample(1:nrow(combined_genes), num_rows, replace = TRUE)

# Subset the data frame based on the selected indices
selected_df <- data.frame(selected_sequences = combined_genes[selected_indices, ])

# Slide a random section from the sequence that is sequence_size nucleotides in length

select_sequence <- function(input_sequence, sequence_size) {

  x <- c("TAA", "TGA", "TAG")
  input_seq <- as.character(input_sequence)
  max_start_pos = nchar(as.character(input_seq)) - sequence_size - 6
  div_by_3_nums <- 1:nchar(input_seq)
  div_by_3_nums <- div_by_3_nums[div_by_3_nums %% 3 == 0]
  div_by_3_nums <- div_by_3_nums[div_by_3_nums <= max_start_pos]
  start_seq <- sample(div_by_3_nums, 1)
  end_seq <- start_seq + (sequence_size - 6)
  slid_seqs <- substr(input_seq, start_seq + 1, end_seq)
  new_sequence <- tolower(paste("ATG", slid_seqs, sample(x, 1), sep = ""))
  return(new_sequence)
}


##### Apply the select_sequence function over all genes, and convert it into a character df

# Set up future plan for parallel execution
plan(multisession, workers = 11)  # Adjust the number of workers as needed

# One common approach is to use the number of seconds since the Epoch (1970-01-01 00:00:00 UTC)
random_seed <- as.integer(Sys.time())

# Set the seed for reproducibility; otherwise, use random_seed
set.seed(random_seed)

# Record the seed for future reference
cat("The seed used for this session: ", random_seed, "\n")

# Use future_mapply with the future backend and ensure RNG kind is set
system.time(df_ran_selected <- as.character(
  future_mapply(select_sequence, selected_df$selected_sequences, MoreArgs = list(sequence_size = 105), future.seed=TRUE)
))

df_ran_selected <- data.frame(pos_ran_selected = as.character(unlist(df_ran_selected)))

df_ran_selected$length <- nchar(as.character(df_ran_selected$pos_ran_selected))

#reverse an input nucleotide sequence

reversed_nt <- function(input_nt_seq){
  
  nt_split_as_string <- toString(input_nt_seq)
  
  nt_split <- strsplit(nt_split_as_string, NULL)[[1]]
  
  reverse_nt <- paste(rev(nt_split), collapse="")
  
  reverse_nt
  
}


#reverse complement function using the Biostrings package

clusterEvalQ(cl, library(Biostrings))  # Load bioseq on each worker

reverse_complement <- function(input_nt_seq){
  
  sequence_as_dna_string <- DNAString(input_nt_seq)
  
  rev_complement <- tolower(toString(reverseComplement(sequence_as_dna_string)))
  return(rev_complement)
}


#apply reverse_complement to the pos_ran_selected

system.time(df_ran_selected$pos_seq_reverse_complement <- unlist(parLapply(cl,
                                                               df_ran_selected$pos_ran_selected,
                                                               reverse_complement)))

# Note, for this function randomly selects a sequence with stop codon, but may contain more than 1 stop codon
# However, the net effect is longer sequences overall, even with more with greater than 1 stop codon.

search_for_stop <- function(nt_sequence) {
  # Check if the sequence contains at least one stop codon
  if (!grepl("tag|tga|taa", nt_sequence)) {
    return("Null")
  }
  
  repeat {
    nt_sequence_new_tag <- sub("tag.*", "", nt_sequence)
    nt_sequence_new_tga <- sub("tga.*", "", nt_sequence)
    nt_sequence_new_taa <- sub("taa.*", "", nt_sequence)
    
    # Create a list of sequences with their corresponding stop codons
    sequences <- list(
      tag = paste(nt_sequence_new_tag, "tag", sep = ""),
      tga = paste(nt_sequence_new_tga, "tga", sep = ""),
      taa = paste(nt_sequence_new_taa, "taa", sep = "")
    )
    
    # Generate a random index to select one of the sequences
    random_index <- sample(1:3, 1)
    
    # Check if the selected sequence is shorter than nt_sequence
    if (nchar(sequences[[random_index]]) < nchar(nt_sequence)) {
      return(sequences[[random_index]])
    }
  }
}

#apply search_for_stop to the pos_seq_reverse_complement

system.time(df_ran_selected$seq_to_stop <- unlist(parLapply(cl,
                                                df_ran_selected$pos_seq_reverse_complement,
                                                search_for_stop)))


### appending

# now we need to filter by those sequences that do not have an internal stop codon 
# using a new function based on the test above.
test_for_internal_codon <- function(nt_sequence) {
  
  # Extract codons from the sequence
  list_of_seq_codons <- regmatches(nt_sequence, gregexpr(".{3}", nt_sequence))[[1]]
  
  # Filter for stop codons
  stop_codons_count <- sum(list_of_seq_codons %in% c("tga", "tag", "taa"))
  
  # Check the number of stop codons and return the result
  if (stop_codons_count == 1) {
    return("1 stop codon")
  } else {
    return("More than 1 stop codon")
  }
  
}

# apply test_for_internal_codon to the pos_seq_reverse_complement

df_ran_selected$stop_codon_count_for <- unlist(parLapply(cl,
                                                         tolower(df_ran_selected$pos_ran_selected), 
                                                         test_for_internal_codon))

# filter the df based on the stop_codon_count_rev to remove those rows that have more than one stop codon in the appended_sequence column

df_ran_selected <- df_ran_selected %>% filter(stop_codon_count_for == ("1 stop codon"))

# function that will append a leader sequence to the terminal overlap sequence
append_lead_nt_to_overlap_new <- function(terminal_nt_sequence, input_source_for_leader, sequence_length){
  
  nucleotides <- c("a", "t", "g", "c")

  leader_test_0 <- terminal_nt_sequence
  leader_test_1 <- paste(sample(nucleotides, 1), leader_test_0, sep = "", collapse = "")
  leader_test_2 <- paste(sample(nucleotides, 1), leader_test_1, sep = "", collapse = "")
  leader_test_3 <- paste(sample(nucleotides, 1), leader_test_2, sep = "", collapse = "")
  
  
  if (nchar(leader_test_0) %% 3 == 0){
    terminal_nt_sequence_new <- leader_test_0
  } else if (nchar(leader_test_1) %% 3 == 0){
    terminal_nt_sequence_new <- leader_test_1
  } else if (nchar(leader_test_2) %% 3 == 0){
    terminal_nt_sequence_new <- leader_test_2
  } else if (nchar(leader_test_3) %% 3 == 0){
    terminal_nt_sequence_new <- leader_test_3
  }

  input_seq <- as.character(terminal_nt_sequence_new)
  source_for_leader <- tolower(as.character(sample(input_source_for_leader[[1]], 1)))
  max_start_pos <- nchar(as.character(source_for_leader)) - (sequence_length) - 12
  div_by_3_nums <- 1:nchar(source_for_leader)
  div_by_3_nums <- div_by_3_nums[div_by_3_nums %% 3 == 0]
  div_by_3_nums <- div_by_3_nums[div_by_3_nums <= max_start_pos]
  start_seq <- sample(div_by_3_nums, 1)
  end_seq <- start_seq - nchar(input_seq) + sequence_length
  slide_seqs <- substr(source_for_leader, start_seq + 1, end_seq)
  new_sequence <- paste("atg", slide_seqs, input_seq, sep = "")
  return(new_sequence)
  
}


## Apply append_lead_nt_to_overlap_new to the seq_to_stop
## this number should be the total desired nucleotide length, minus 3 to account for the stop codon already present in the overlap.

df_ran_selected$rev_ran_selected <- unlist(pbmapply(append_lead_nt_to_overlap_new,
                                                           df_ran_selected$seq_to_stop,
                                                           combined_genes,
                                                           102))


#apply test_for_internal_codon to the new rev sequence

df_ran_selected$stop_codon_count_rev <- unlist(parLapply(cl,
                                                         df_ran_selected$rev_ran_selected, 
                                                         test_for_internal_codon))

#filter the df based on the stop_codon_count_rev to remove those rows that have more than one stop codon in the appended_sequence column

new_df_filtered <- df_ran_selected %>% filter(stop_codon_count_for == ("1 stop codon"))
new_df_filtered <- new_df_filtered %>% filter(stop_codon_count_rev == ("1 stop codon"))


overlap_nt_df_forward <- data.frame(tolower(new_df_filtered$pos_ran_selected))
overlap_nt_df_reverse <- data.frame(tolower(new_df_filtered$rev_ran_selected))

overlap_nt_df <- data.frame(overlap_nt_df_forward, 
                            overlap_nt_df_reverse)

colnames(overlap_nt_df) <- c("nt_sequence_1", "nt_sequence_2")

reversed_nt(overlap_nt_df$nt_sequence_2[1])

#adding back in the seq_to_rev_stop overlap component

overlap_nt_df$seq_to_rev_stop <- new_df_filtered$seq_to_stop

#add a column that includes the number of overlapping nucleotides

overlap_nt_df$count_seq_to_rev_stop <- pblapply(as.character(overlap_nt_df$seq_to_rev_stop), nchar)

#filtering by count_seq_to_rev_stop value to eliminate the very low overlaps

overlap_nt_df <- overlap_nt_df[overlap_nt_df$count_seq_to_rev_stop >= 1, ]

# Group by 'count_seq_to_rev_stop' and count each
frequency_df <- overlap_nt_df %>%
  group_by(count_seq_to_rev_stop) %>%
  summarize(Counts = n())

# View the resulting dataframe
print(frequency_df)

# Modify the codon_overlap_groups function to return nchar value directly
codon_overlap_groups <- function(count_seq_to_rev_stop){
  return(as.character(count_seq_to_rev_stop))
}

# Assign groups based on the nchar value
overlap_nt_df$group <- unlist(pblapply(overlap_nt_df$count_seq_to_rev_stop[1:nrow(overlap_nt_df)],
                                                codon_overlap_groups))

##### NOTE: when you do not want to make uniforms, comment out the lines below

# Count the number of rows for each unique group
group_counts <- table(overlap_nt_df$group)

# Identify the group with the minimum count
min_group <- names(which.min(group_counts))
#target_num <- min(group_counts)
#target_num <- 120

# Normalize function
normalize_group <- function(group_num, target_num) {
  num_group <- sum(overlap_nt_df$group == group_num)
  if (num_group > target_num) {
    return(sample(which(overlap_nt_df$group == group_num), num_group - target_num))
  } else {
    return(integer(0)) # empty integer vector
  }
}

# Apply the function to each unique group
all_drop_indices <- unlist(lapply(names(group_counts), function(x) normalize_group(x, target_num)))

# Drop the rows from the dataframe
overlap_nt_df_balanced <- overlap_nt_df[-all_drop_indices,]

# Visualize groups
hist(as.numeric(unlist(overlap_nt_df_balanced$count_seq_to_rev_stop)))

overlap_nt_df_balanced_for_hist <- as.numeric(unlist(overlap_nt_df_balanced$count_seq_to_rev_stop))
hist(overlap_nt_df_balanced_for_hist,
     breaks=seq(from=min(overlap_nt_df_balanced_for_hist),
                to=max(overlap_nt_df_balanced_for_hist)))

hist(as.numeric(unlist(overlap_nt_df_balanced$group)))
nrow(overlap_nt_df)
nrow(overlap_nt_df_balanced)

# test forward and reverse sequence length, should both be equal

nchar(as.character(overlap_nt_df$forward_nt_pos))
nchar(as.character(overlap_nt_df$reverse_nt_pos))


#when we do not want to make uniforms, comment out the lines below

overlap_nt_df <- overlap_nt_df_balanced

#translate the nt sequences in overlap_nt_df columns, and add these as new columns

translate_nt_sequence <- function(nt_sequence){

  #note, this function uses the seq_translate function from the bioseq library
  aa_sequence <- toString(seq_translate(dna(toString(nt_sequence))))
  return(aa_sequence)
  
}

#adding columns of aa sequences

clusterEvalQ(cl, library(bioseq))  # Load bioseq on each worker


system.time(overlap_nt_df$aa_sequence_1 <- unlist(parLapply(cl,
                                                overlap_nt_df$nt_sequence_1, 
                                                translate_nt_sequence)))

system.time(overlap_nt_df$aa_sequence_2 <- unlist(parLapply(cl,
                                                overlap_nt_df$nt_sequence_2, 
                                                translate_nt_sequence)))

######### end to be integrated

########### generating negative controls
#random_nt_sequence_neg_control function, modified to only output the forward and reverse negative control sequences
#as opposed to the additional sequences that are generated with the positive control function

# Set the number of rows you want to select
#num_rows_neg <- sum(overlap_nt_df$count_seq_to_rev_stop == 4) * 200000 # set this value to determine the num of negs

# Set the number of rows you want to select
num_rows_neg <- 2 #use the commented out line above for unknown generation; given this is for knowns, setting this to 2 as we don't need neg controls

# Use sample to select a random subset of row indices
selected_indices_neg <- sample(1:nrow(combined_genes), num_rows_neg, replace = TRUE)

# Subset the data frame based on the selected indices
selected_df_neg <- data.frame(selected_indices_neg = combined_genes[selected_indices_neg, ])


# Slide a random section from the sequence that is 300 nucleotides in length
# apply the select_sequence function over all genes, and convert it into a character df

neg_df_ran_selected <- as.character(pbmapply(select_sequence, selected_df_neg$selected_indices_neg, 105))

neg_df_ran_selected <- data.frame(neg_ran_selected = as.character(unlist(neg_df_ran_selected)))

#convert these to aa, and add this as another column, and then convert the aa to a number and add that column

neg_df_ran_selected$aa_sequence <- unlist(parLapply(cl,
                                                    neg_df_ran_selected$neg_ran_selected, 
                                                    translate_nt_sequence))

colnames(neg_df_ran_selected) <- c("random_sequence_neg_control",
                                   "neg_control_aa")

##################### 
#building the neg control df before combining with the positive control df
#these numbers serve as variables based on the row number of the negative control data frame
number_of_rows_neg_div_2 <- nrow(neg_df_ran_selected)/2
number_of_rows_neg_div_2_plus_1 <- (nrow(neg_df_ran_selected)/2) + 1
number_of_rows_neg <- nrow(neg_df_ran_selected)

#these steps will divide the neg control data frame and build it with side-by-side nt, aa, and aa_to_number columns

#nt columns
neg_control_df <- data.frame(neg_df_ran_selected$random_sequence_neg_control[1:number_of_rows_neg_div_2])
colnames(neg_control_df) <- c("nt_sequence_1")
neg_control_df$nt_sequence_2 <- neg_df_ran_selected$random_sequence_neg_control[number_of_rows_neg_div_2_plus_1:number_of_rows_neg]
#aa columns
neg_control_df$aa_sequence_1 <- neg_df_ran_selected$neg_control_aa[1:number_of_rows_neg_div_2]
neg_control_df$aa_sequence_2 <- neg_df_ran_selected$neg_control_aa[number_of_rows_neg_div_2_plus_1:number_of_rows_neg]
#aa_to_number columns
neg_control_df$aa_sequence_1_as_num <- neg_df_ran_selected$neg_control_aa_as_num[1:number_of_rows_neg_div_2]
neg_control_df$aa_sequence_2_as_num <- neg_df_ran_selected$neg_control_aa_as_num[number_of_rows_neg_div_2_plus_1:number_of_rows_neg]
#adding null values for seq_to_rev_stop column (this is a column from the positive control df we are trying to match)
neg_control_df$seq_to_rev_stop <- c("null")
#adding "0" values for count_seq_to_rev_stop column (this is a column from the positive control df we are trying to match)
neg_control_df$count_seq_to_rev_stop <- c("0")

# add groups to neg_control_df

neg_control_df$group <- unlist(pblapply(neg_control_df$count_seq_to_rev_stop[1:nrow(neg_control_df)],
                                    codon_overlap_groups))

#combine the positive control df and the negative control df

pos_neg_df <- rbind(overlap_nt_df, neg_control_df)

#now export the pos_neg_df to csv
#needed to use this apply function on the df to permit exporting, as it yielded an error otherwise.
pos_neg_df_export <- apply(pos_neg_df, c(1, 2), as.character)

write.csv(pos_neg_df_export, paste("D:/RStuff/codon_overlap/genes_pos_neg_df_", 
                                   iteration_number, ".csv", sep = ""), row.names = TRUE)


# Define hist variable as numeric, and unlist
pos_neg_df_for_hist <- as.numeric(unlist(pos_neg_df$count_seq_to_rev_stop))

# Histogram with bins of width 1
hist(pos_neg_df_for_hist, 
     breaks=seq(from=min(pos_neg_df_for_hist), 
                to=max(pos_neg_df_for_hist), by=1))

nrow(pos_neg_df)

#### test for GC content of nt sequences

calculate_gc_frequency <- function(input_string) {
  
  converted_input_string <- toString(input_string)
  
  gc_frequency <- (stringr::str_count(converted_input_string, "g") + stringr::str_count(converted_input_string, "c"))/nchar(converted_input_string)
  
  return(as.numeric(gc_frequency))
}

#test calculate_gc_frequency

test_pos <- calculate_gc_frequency(pos_neg_df$nt_sequence_1[1:nrow(pos_neg_df)/2])
test_neg <- calculate_gc_frequency(pos_neg_df$nt_sequence_1[(1 + nrow(pos_neg_df)/2):nrow(pos_neg_df)])
mean(test_pos)
mean(test_neg)


pos_neg_df <- data.frame(apply(pos_neg_df, c(1, 2), as.character))

write.csv(pos_neg_df, paste("D:/RStuff/codon_overlap/genes_pos_neg_df_", 
                            iteration_number, ".csv", sep = ""), row.names = TRUE)
#############

pos_neg_df$stop_codon_count_for <- unlist(parLapply(cl,
                                                    pos_neg_df$nt_sequence_1,
                                                    test_for_internal_codon))

pos_neg_df$stop_codon_count_rev <- unlist(parLapply(cl,
                                                    pos_neg_df$nt_sequence_2,
                                                    test_for_internal_codon))

#filter the df based on the stop_codon_count_rev to remove those rows that have more than one stop codon in the appended_sequence column

new_df_filtered <- pos_neg_df %>% filter(stop_codon_count_for == ("1 stop codon"))
new_df_filtered <- new_df_filtered %>% filter(stop_codon_count_rev == ("1 stop codon"))
new_df_filtered <- new_df_filtered %>% filter(count_seq_to_rev_stop > ("0")) #use this line only to filter out the neg controls if setting up known dataset

######################################################################
#### NOTE ######## ~~~~~~~~
### Comment out this code if not attempting to produce only unknowns.
# pos_neg_df <- new_df_filtered
# 
# filtered_df <- pos_neg_df %>% filter(group == 0)
# 
# # Randomly select 15000 rows from the filtered DataFrame
# set.seed(as.integer(Sys.time()))
# pos_neg_df <- sample_n(filtered_df, 30000)
######################################################################

pos_neg_df <- new_df_filtered

replace_letters <- function(str) {
  new_str <- gsub("a", "B", str)
  new_str <- gsub("g", "J", new_str)
  new_str <- gsub("t", "O", new_str)
  new_str <- gsub("c", "U", new_str)
  return(new_str)
}


pos_neg_df$for_parsing <- unlist(pblapply(pos_neg_df$seq_to_rev_stop, replace_letters))

add_space_between_letters <- function(str) {
  new_str <- gsub("(.)", "\\1 ", str)
  return(new_str)
}

pos_neg_df$aa_sequence_1_spaces <- pblapply(pos_neg_df$aa_sequence_1, add_space_between_letters)
pos_neg_df$aa_sequence_2_spaces <- pblapply(pos_neg_df$aa_sequence_2, add_space_between_letters)
pos_neg_df$for_parsing_spaces <- unlist(pblapply(pos_neg_df$for_parsing, add_space_between_letters))

replace_n <- function(x) {
  x[grep("n", x, ignore.case = TRUE)] <- "No overlap"
  return(x)
}

pos_neg_df$for_parsing_spaces <- unlist(pblapply(pos_neg_df$for_parsing_spaces, replace_n))

concatenate_with_space <- function(str1, str2) {
  new_str <- paste(str1, str2, sep = "")
  return(new_str)
}

pos_neg_df$concat_pos_neg <- mapply(concatenate_with_space, pos_neg_df$aa_sequence_1_spaces, pos_neg_df$aa_sequence_2_spaces)

# Write the column to a CSV file
write.csv(pos_neg_df$concat_pos_neg, paste("D:/RStuff/codon_overlap/genes_concat_pos_neg_", 
                                           iteration_number, ".csv", sep = ""), row.names = TRUE)

write.csv(pos_neg_df$for_parsing_spaces, paste("D:/RStuff/codon_overlap/genes_for_parsing_pos_neg_", 
                                               iteration_number, ".csv", sep = ""), row.names = TRUE)

# convert the list to a JSON string
concat_pos_neg <- toJSON(pos_neg_df$concat_pos_neg)

# write the JSON string to a file
write(concat_pos_neg, file = paste("D:/RStuff/codon_overlap/genes_concat_pos_neg_", 
                                   iteration_number, ".json", sep = ""))

# convert the list to a JSON string
overlap_pos_neg <- toJSON(pos_neg_df$for_parsing_spaces)

# write the JSON string to a file
write(overlap_pos_neg, file = paste("D:/RStuff/codon_overlap/genes_overlap_", 
                                    iteration_number, ".json", sep = ""))






