import tensorflow as tf
import numpy as np
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import LayerNormalization, MultiHeadAttention
from tensorflow.keras.layers.experimental.preprocessing import TextVectorization
import matplotlib.pyplot as plt
import pandas as pd
import csv
import json
import nltk
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.Seq import CodonTable
from tensorflow.keras.losses import Loss
import tensorflow_addons as tfa
import itertools as it
import pandas as pd
import dask.dataframe as dd
from dask.diagnostics import ProgressBar

### FOR TOKENIZER, OPTIMALLY USE THE TEST DATASET SINCE IT IS SMALL
# to test the trained model
with open("D:/RStuff/codon_overlap/concat_pos_neg_114.json") as concat:
    concat_test_model = json.load(concat)

with open("D:/RStuff/codon_overlap/overlap_114.json") as overlap:
    overlap_test_model = json.load(overlap)
    
### Dataset to test the model
### This is the Genes dataset taken from actual DNA sequences

with open("D:/RStuff/codon_overlap/genes_concat_pos_neg_167.json") as concat:
    genes_concat_test_model = json.load(concat)

with open("D:/RStuff/codon_overlap/genes_overlap_167.json") as overlap:
    genes_overlap_test_model = json.load(overlap)

# Tokenization and text vectorization
max_length = 103 #note, this should be the max len of the overlap, not input sequences
vocab_size = 35


# We must use the same tokenizer dataset as that used to train the original model
def tokenize(sentences):
    tokenizer = tf.keras.preprocessing.text.Tokenizer(
        num_words=vocab_size, filters='')
    tokenizer.fit_on_texts(sentences)
    return tokenizer


concat_tokenizer = tokenize(concat_test_model)
overlap_tokenizer = tokenize(overlap_test_model)

def vectorize(tokenizer, sentences):
    seqs = tokenizer.texts_to_sequences(sentences)
    return tf.keras.preprocessing.sequence.pad_sequences(seqs, maxlen=max_length, padding='post')


concat_test_model_np = np.array(concat_test_model)
overlap_test_model_np = np.array(overlap_test_model)


concat_test_model_np_tokenizer = tokenize(concat_test_model_np)
overlap_test_model_np_tokenizer = tokenize(overlap_test_model_np)


concat_test_model_np_concat_vectorized = vectorize(concat_test_model_np_tokenizer, concat_test_model_np)
overlap_test_model_np_overlap_vectorized = vectorize(overlap_test_model_np_tokenizer, overlap_test_model_np)

#these np variables have been updated to convert the genes datasets.
concat_test_model_np = np.array(genes_concat_test_model)
overlap_test_model_np = np.array(genes_overlap_test_model)



####################################################################################
### Combined function to predict the overlapping dna sequence for two aa's

def predict_overlapping_sequence(trained_model, aa_sequences):
    
    def translate_to_dna_with_all_options(aa_sequence: str) -> list:
    
    
    # Codons and their frequencies for each amino acid based on the E. coli table
        back_translation_code_with_all_options = {
            'A': [('GCG', 0.27), ('GCT', 0.26), ('GCC', 0.26), ('GCA', 0.21)],
            'C': [('TGC', 0.53), ('TGT', 0.47)],
            'D': [('GAT', 0.63), ('GAC', 0.37)],
            'E': [('GAA', 0.68), ('GAG', 0.32)],
            'F': [('TTT', 0.58), ('TTC', 0.42)],
            'G': [('GGC', 0.35), ('GGT', 0.32), ('GGG', 0.25), ('GGA', 0.08)],
            'H': [('CAT', 0.56), ('CAC', 0.44)],
            'I': [('ATT', 0.48), ('ATC', 0.39), ('ATA', 0.14)],
            'K': [('AAA', 0.74), ('AAG', 0.26)],
            'L': [('CTG', 0.43), ('CTT', 0.13), ('CTC', 0.13), ('TTA', 0.14), ('CTA', 0.07), ('TTG', 0.13)],
            'M': [('ATG', 1.00)],
            'N': [('AAC', 0.60), ('AAT', 0.40)],
            'P': [('CCG', 0.52), ('CCA', 0.19), ('CCT', 0.16), ('CCC', 0.13)],
            'Q': [('CAG', 0.66), ('CAA', 0.34)],
            'R': [('CGT', 0.36), ('CGC', 0.36), ('CGG', 0.11), ('AGA', 0.08), ('AGG', 0.05), ('CGA', 0.04)],
            'S': [('AGC', 0.24), ('TCC', 0.24), ('TCT', 0.17), ('TCG', 0.15), ('TCA', 0.14), ('AGT', 0.15)],
            'T': [('ACC', 0.36), ('ACA', 0.28), ('ACG', 0.25), ('ACT', 0.11)],
            'V': [('GTG', 0.46), ('GTT', 0.28), ('GTC', 0.15), ('GTA', 0.11)],
            'W': [('TGG', 1.00)],
            'Y': [('TAT', 0.59), ('TAC', 0.41)],
            '*': [('TAA', 0.61), ('TGA', 0.30), ('TAG', 0.09)]
        }
    
        def choose_codon_based_on_frequency(codons):
            """
            Randomly choose a codon based on the frequency distribution.
            
            Args:
            - codons (list): A list of tuples where each tuple contains a codon and its frequency
            
            Returns:
            - list: A list containing the chosen codon
            """
            # Extract codon names and their frequencies
            codon_names = [codon for codon, _ in codons]
            codon_freqs = [freq for _, freq in codons]
            
            # Normalize the frequencies to ensure they sum up to 1
            total_frequency = sum(codon_freqs)
            normalized_freqs = [freq/total_frequency for freq in codon_freqs]
            
            # Randomly select a codon based on the frequency distribution
            chosen_codon = np.random.choice(codon_names, p = normalized_freqs)
            
            return [chosen_codon]
    
        # For each amino acid in the sequence, choose a codon based on its frequency
        list_of_list_of_codons = [choose_codon_based_on_frequency(back_translation_code_with_all_options[aa]) for aa in aa_sequence]
        
        # Combine the codons to form the nucleotide sequence
        list_of_combinations = [''.join(combination) for combination in it.product(*list_of_list_of_codons)]
        
        return list_of_combinations

    def translate(model, tokenizer, input_text):
        input_vectorized = vectorize(tokenizer, [input_text])
        prediction = model.predict(input_vectorized)[0]
        predicted_indices = np.argmax(prediction, axis=-1)
        #note that, we use the same tokenizer from the training data so that the output has the same tokens assigned
        predicted_text = overlap_tokenizer.sequences_to_texts([predicted_indices])[
            0].strip()
        return predicted_text


    def convert_chars_for_translated_overlap(string):
        string = string.replace("no overlap", "#", 1)  # Replace first occurrence of "no overlap" with "#"
        string = string.replace("overlap", "#", 1)  # Replace first occurrence of "overlap" with "#"
        string = string.replace("no", "#", 1)  # Replace first occurrence of "no" with "#"
        if string == "#":
            return string
        new_string = ""
        for char in string:
            if char.lower() == "b":
                new_string += "A"
            elif char.lower() == "j":
                new_string += "G"
            elif char.lower() == "o":
                new_string += "T"
            elif char.lower() == "u":
                new_string += "C"
            else:
                new_string += char
        return new_string.replace(" ", "")

    # Running the above sequence through the model, using the tokenizer run on the data from the original model training process
    translated_overlap = translate(trained_model, concat_test_model_np_tokenizer, aa_sequences)
    
    overlap_output_from_model = convert_chars_for_translated_overlap(translated_overlap)
    
    # Generate the rc to include in the next small df
    reverse_complement_overlap = Seq(overlap_output_from_model)
    reverse_complement_overlap = str(reverse_complement_overlap.reverse_complement())
    
    # Create a DataFrame with two rows: original and reversed sequences
    # Note, the rc sequence actually comes first here, followed by the model output overlap sequence
    output_forward_reverse = pd.DataFrame({"overlap_sequence": [reverse_complement_overlap, overlap_output_from_model]})
    
    # Splitting and collapsing the sequences
    sequence_list = aa_sequences.split()
    collapsed_sequence = ''.join(sequence_list)
    
    # Splitting the sequence at the asterisk and keeping the asterisk
    split_sequences_with_asterisk = [seq + '*' for seq in collapsed_sequence.split('*') if seq]
    
    # Creating a dataframe from the sequences
    df_sequences = pd.DataFrame(split_sequences_with_asterisk, columns=["amino_acid_sequence"])
    
    # Adding the back translated dna sequences to the data frame.
    df_sequences['nt_sequence'] = df_sequences['amino_acid_sequence'].apply(lambda aa_seq: translate_to_dna_with_all_options(aa_seq)[0])
    
    df_sequences["overlap_seq"] = output_forward_reverse
    
    #here, we are attaching the overlap sequence to the known coding sequence, generated above. It is
    #added at the location based on the length of the overlap_seq
    def modify_sequence(df):
        # Check if 'overlap_seq' column exists in the DataFrame
        if 'overlap_seq' not in df.columns:
            raise ValueError("DataFrame must contain 'overlap_seq' column.")
    
        # Calculate the length of the overlap sequence
        df['overlap_length'] = df['overlap_seq'].apply(len)
    
        # Modify the sequences
        df['modified_sequence'] = df.apply(lambda row: row['nt_sequence'][:-row['overlap_length']] + row['overlap_seq'], axis=1)
    
        return df
    
    modify_sequence(df_sequences)
    
    # Custom function to translate using Seq
    def translate_with_seq(sequence):
        coding_dna = Seq(sequence)
        return str(coding_dna.translate())
    
    # Apply the custom function to the 'nt_sequence' column
    df_sequences['translated_sequence'] = df_sequences['modified_sequence'].apply(translate_with_seq)
    
    #print(compare_aa)
    return(df_sequences)

################################################
# Using the predict_overlapping_sequence function
################################################

model_to_predict = tf.keras.models.load_model('codon_overlap_tf_transformer_saved_model/my_model_new_dataset_20240111_1')

formatted_aa_seq = "M A G S I V K E L H P S D H V V I R A L P S S R H V S S A R L E Q Q * M Q V I C D P L T G L N S S E A S T A L G L G I A G G E V E G V L S * "

#directly calling the function once
predict_overlapping_sequence(model_to_predict, formatted_aa_seq)

testing_df = predict_overlapping_sequence(model_to_predict, formatted_aa_seq)

testing_df

#function to call the prediction function multiple times to account for the variability in output
#if matches are acheived in the aa sequences, it provides the dataframe. Else, it states that there is
#no predicted significant overlap

def find_matching_sequence(model_to_predict, formatted_aa_seq, max_attempts, _counter=[0]):
    _counter[0] += 1
    print(f"Parsing row {_counter[0]}...")
    
    attempts = 0
    match_found = False
    matching_dataframe = None

    while attempts < max_attempts:
        attempts += 1

        try:
            result_df = predict_overlapping_sequence(model_to_predict, formatted_aa_seq)
        except CodonTable.TranslationError as e:
            print(f"Error: {e}")
             # Print the "overlap_length" at the end of each iteration
            print(f"Attempt {attempts}... Overlap Length: {matching_dataframe['overlap_length'].iloc[0] if matching_dataframe is not None else 'N/A'}")
            continue

        # Check if the "nt_sequence" and "translated_sequence" columns match in both rows
        if (result_df['amino_acid_sequence'].iloc[0] == result_df['translated_sequence'].iloc[0] and
            result_df['amino_acid_sequence'].iloc[1] == result_df['translated_sequence'].iloc[1]):
            match_found = True
            matching_dataframe = result_df
            break

    if match_found:
         # Print the "overlap_length" at the end of each iteration
        print(f"Attempt {attempts}... Overlap Length: {matching_dataframe['overlap_length'].iloc[0] if matching_dataframe is not None else 'N/A'}")
        return matching_dataframe
    else:
        print("There is no predicted significant overlap")
        return None



formatted_aa_seq = "M E D L R S D R Q P E F T Q I D C E L C F A D G E K V K I F I E K L * M E K A H E E S T I V N G V I N G K V K G G F T V E L H G I R A F L * " 

testing_df = find_matching_sequence(model_to_predict, formatted_aa_seq, 100)

testing_df

#############
###### Start building dataframe of variable sequences by pulling the forward sequence with overlap attached
#############

# Initialize an empty DataFrame to store the results
result_df_collection = pd.DataFrame()

# The number of times you want to run the function
num_iterations = 1000

# The amino acid sequence you are using for each prediction
formatted_aa_seq = "M I Y Y T I L S K L A S D A E K T Q T G L D K A T G L V R S E L G S * M I P Q L P S S L R T S P V A L S N P V W V F S A S E A S F D R I V * "

for i in range(num_iterations):
    # Run the function and get the resulting DataFrame
    testing_df = find_matching_sequence(model_to_predict, formatted_aa_seq, 100)
    
    if testing_df is not None:
        # Extract the first row of the 'modified_sequence' column
        first_row_modified_seq = testing_df['modified_sequence'].iloc[0]
        
        # Append this to result_df_collection
        result_df_collection = result_df_collection.append({'modified_sequence': first_row_modified_seq}, ignore_index=True)

# The result_df_collection DataFrame will have all the first row 'modified_sequence' from each run.
print(result_df_collection)

###### END Attempting to collect variations


######################
#attempting to evaluate performance
######################

def convert_chars(string):
    if string == "No overlap":
        return "none"
    string = string.replace(" ", "")  # Remove spaces from the string
    new_string = ""
    for char in string:
        if char.lower() == "b":
            new_string += "A"
        elif char.lower() == "j":
            new_string += "G"
        elif char.lower() == "o":
            new_string += "T"
        elif char.lower() == "u":
            new_string += "C"
        else:
            new_string += char
    return new_string

def convert_chars_for_translated_overlap(string):
    string = string.replace("no overlap", "#", 1)  # Replace first occurrence of "no overlap" with "#"
    string = string.replace("overlap", "#", 1)  # Replace first occurrence of "overlap" with "#"
    string = string.replace("no", "#", 1)  # Replace first occurrence of "no" with "#"
    if string == "#":
        return string
    new_string = ""
    for char in string:
        if char.lower() == "b":
            new_string += "A"
        elif char.lower() == "j":
            new_string += "G"
        elif char.lower() == "o":
            new_string += "T"
        elif char.lower() == "u":
            new_string += "C"
        else:
            new_string += char
    return new_string.replace(" ", "")

def translate(model, tokenizer, input_text):
    input_vectorized = vectorize(tokenizer, [input_text])
    prediction = model.predict(input_vectorized)[0]
    predicted_indices = np.argmax(prediction, axis=-1)
    #note that, we use the same tokenizer from the training data so that the output has the same tokens assigned
    predicted_text = overlap_tokenizer.sequences_to_texts([predicted_indices])[
        0].strip()
    return predicted_text

def align_sequences(seq1, seq2):
    # Perform the alignment using the Needleman-Wunsch algorithm
    alignments = pairwise2.align.globalxx(seq1, seq2)
    
    # Find the highest alignment score
    max_score = float('-inf')
    for a in alignments:
        score = a[2]
        if score > max_score:
            max_score = score
    
    return max_score

def normalize_score(row, col_name, align_score_col):
    # Get the predicted sequence from the specified column of the row
    prediction = row[col_name]
    
    # Check if the known overlap is "none", meaning no known overlap
    if row['converted_overlap'] == "none":
        # Count the number of "#" characters in the predicted sequence
        num_hashes = prediction.count("#")
        
        # Calculate the number of characters in the predicted sequence that are not "#"
        num_extra_chars = len(prediction) - num_hashes
        
        # If the predicted sequence contains only "#" characters, it's a perfect non-overlap prediction
        # and the score is set to the maximum value of 1
        if num_hashes > 0 and num_extra_chars == 0:
            return 1
        
        # If there are other characters in the predicted sequence besides "#", a penalty is calculated
        if num_extra_chars > 0:
            # The penalty is determined based on the proportion of extra characters in the predicted sequence
            penalty = 1 - (num_extra_chars / len(prediction))
            
            # The penalty is adjusted to ensure:
            # - A minimum score of 0.01 if there's at least one "#" character
            # - A score of 0 if there's no "#" character
            # - The score does not exceed 1
            return max(0.01 if num_hashes > 0 else 0, min(1, round(penalty, 2)))
        
        # If none of the above conditions are met, the score is 0
        return 0
    
    # Check if there's a known overlap (i.e., the converted_overlap value is not empty)
    elif len(row['converted_overlap']) > 0:
        # Calculate the alignment score normalized by the length of the known overlap
        # The score is bounded between 0 and 1
        return min(1, max(0, round(row[align_score_col] / len(row['converted_overlap']), 2)))
    
    # Default case: If none of the above conditions are met, the score is 0
    return 0


# Testing the predictions
test_df = pd.DataFrame(concat_test_model_np, columns=['concat'])
test_df['overlap'] = overlap_test_model_np
test_df['converted_overlap'] = test_df['overlap'].apply(convert_chars)
test_df['length_actual'] = np.where(test_df['converted_overlap'] == "none", 1, test_df['converted_overlap'].str.len())

#### Check the concat column for a letter X

test_df = test_df[~test_df['concat'].str.contains('x', case=False, na=False)]

#### In this code, we set test_df to be only those rows with length_actual equal to 1 (no known overlap)
(test_df['length_actual'] == 1).sum()
test_df = test_df[test_df['length_actual'] == 1]

# Model translations
test_df['translated_overlap'] = test_df['concat'].apply(lambda x: translate(model_to_predict, concat_test_model_np_tokenizer, x))
test_df['converted_translated_overlap'] = test_df['translated_overlap'].apply(convert_chars_for_translated_overlap)
test_df['alignment_score'] = test_df.apply(lambda row: align_sequences(row['converted_overlap'], row['converted_translated_overlap']), axis=1)
test_df['length_1'] = test_df['converted_translated_overlap'].str.len()

# Normalize alignment_score
test_df['normalized_alignment_score'] = test_df.apply(lambda row: normalize_score(row, 'converted_translated_overlap', 'alignment_score'), axis=1)

# Define a function to find the matching overlap_seq value
def find_matching_overlap_seq(model, sequence, max_attempts):
    result_df = find_matching_sequence(model, sequence, max_attempts)
    if result_df is not None:
        return result_df['overlap_seq'].iloc[1] #note, this is 1 and not 0 to get the sequence consistent with the concat translation
    else:
        return '#'

# Add a new column to your test_df to test with find_matching_overlap_seq function
_counter = 0
test_df['matching_overlap_seq_2'] = test_df['concat'].apply(lambda x: find_matching_overlap_seq(model_to_predict, x, 100))
test_df['alignment_score_2'] = test_df.apply(lambda row: align_sequences(row['converted_overlap'], row['matching_overlap_seq_2']), axis=1)
test_df['length_2'] = test_df['matching_overlap_seq_2'].str.len()

# Normalize alignment_score_2
test_df['normalized_alignment_score_2'] = test_df.apply(lambda row: normalize_score(row, 'matching_overlap_seq_2', 'alignment_score_2'), axis=1)

#save/checkpoint

test_df.to_csv('D:/Users/jkm78/Spyder/test_df_98_overlap_prediction_from_uniform_knowns.csv', index=False)

test_df = pd.read_csv('D:/Users/jkm78/Spyder/test_df_98_overlap_prediction_from_uniform_knowns.csv')













