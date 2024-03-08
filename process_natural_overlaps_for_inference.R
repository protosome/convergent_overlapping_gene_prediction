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

# Define the path to the file
aa_seq_file_path <- "D:/RStuff/codon_overlap/Codon Overlap Examples from Paper, supplementary/aa_seq_df_convergent.csv"

# Import the CSV file
aa_seq_file_path_df <- read.csv(aa_seq_file_path)

aa_seq_df <- aa_seq_file_path_df[,-1]

# Display the first few rows of the dataframe
head(aa_seq_file_path_df)

filtered_aa_seq_df <- aa_seq_df %>% filter(overlapping_length < 102)

process_sequences <- function(aa_seq_1, aa_seq_2) {
  # Keep only the final 34 amino acids of each sequence
  aa_seq_1_trimmed <- substr(aa_seq_1, nchar(aa_seq_1) - 33, nchar(aa_seq_1))
  aa_seq_2_trimmed <- substr(aa_seq_2, nchar(aa_seq_2) - 33, nchar(aa_seq_2))
  
  # Replace the first amino acid with an 'M'
  aa_seq_1_processed <- paste0("M", substr(aa_seq_1_trimmed, 2, 34))
  aa_seq_2_processed <- paste0("M", substr(aa_seq_2_trimmed, 2, 34))
  
  # Add an asterisk after each sequence
  aa_seq_1_processed <- paste0(aa_seq_1_processed, "*")
  aa_seq_2_processed <- paste0(aa_seq_2_processed, "*")
  
  # Concatenate the two sequences
  concatenated_seq <- paste0(aa_seq_1_processed, aa_seq_2_processed)
  
  # Add a space between every character
  final_seq <- paste(unlist(strsplit(concatenated_seq, "")), collapse = " ")
  
  # Return the processed sequence
  return(final_seq)
}

#test the processing function
test <- process_sequences(filtered_aa_seq_df$aa_for_gene_1[1], filtered_aa_seq_df$aa_for_gene_2[1])

nchar(test)

#use function to add new column to filtered dataframe

filtered_aa_seq_df$process_concatenated <- unlist(mapply(process_sequences, filtered_aa_seq_df$aa_for_gene_1, filtered_aa_seq_df$aa_for_gene_2, SIMPLIFY = FALSE))

filtered_aa_seq_df$process_concatenated[1]

# Specify the path where you want to save the CSV file
file_path_to_save <- "D:/RStuff/codon_overlap/draft_paper/processed_sequences.csv"

# Save the dataframe to the CSV file
write.csv(filtered_aa_seq_df, file_path_to_save, row.names = FALSE)















