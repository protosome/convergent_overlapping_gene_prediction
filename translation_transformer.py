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
from tensorflow.keras.losses import Loss
import tensorflow_addons as tfa


### FOR TOKENIZER, OPTIMALLY USE THE TEST DATASET SINCE IT IS SMALL
# to test the trained model
with open("D:/RStuff/codon_overlap/concat_pos_neg_114.json") as concat:
    concat_test_model = json.load(concat)

with open("D:/RStuff/codon_overlap/overlap_114.json") as overlap:
    overlap_test_model = json.load(overlap)

#genes dataset with 1,104,568 size. This is the dataset with a consistent number of overlaps generated for genes with GC content ranging
#from ~23 to 75%
    
with open("D:/RStuff/codon_overlap/concat_24.json") as concat:
    concat = json.load(concat)

with open("D:/RStuff/codon_overlap/overlap_24.json") as overlap:
    overlap = json.load(overlap)



#identify available GPU
print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))


# Tokenization and text vectorization
max_length = 103 #note, this should be the max len of the overlap, not input sequences
vocab_size = 35


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


concat_vectorized = vectorize(concat_tokenizer, concat)
overlap_vectorized = vectorize(overlap_tokenizer, overlap)

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(
    concat_vectorized, overlap_vectorized, test_size=0.0125, random_state=42)


# Transformer model
embedding_dim = 36
num_heads = 24
ffn_dim = 36
    
class TransformerBlock(tf.keras.layers.Layer):
    def __init__(self, dropout_rate=0.1, **kwargs):
        super(TransformerBlock, self).__init__(**kwargs)
        self.att = MultiHeadAttention(num_heads=num_heads, key_dim=embedding_dim)
        self.dropout1 = tf.keras.layers.Dropout(dropout_rate)  # Dropout after MultiHeadAttention
        self.norm1 = LayerNormalization(epsilon=1e-6)
        
        self.ffn = tf.keras.Sequential([
            tf.keras.layers.Dense(ffn_dim, activation='relu'),
            tf.keras.layers.Dense(embedding_dim)
        ])
        self.dropout2 = tf.keras.layers.Dropout(dropout_rate)  # Dropout after Feed Forward Network
        self.norm2 = LayerNormalization(epsilon=1e-6)

    def call(self, inputs, training):
        attn_output = self.att(inputs, inputs)
        attn_output = self.dropout1(attn_output, training=training)  # Apply dropout
        out1 = self.norm1(inputs + attn_output)
        
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)  # Apply dropout
        return self.norm2(out1 + ffn_output)


# Transformer model
class PositionalEncoding(tf.keras.layers.Layer):
    def __init__(self, position, d_model):
        super(PositionalEncoding, self).__init__()
        self.pos_encoding = self.positional_encoding(position, d_model)

    def get_angles(self, position, i, d_model):
        angle_rates = 1 / np.power(10000, (2 * (i // 2)) / np.float32(d_model))
        return position * angle_rates

    def positional_encoding(self, position, d_model):
        angle_rads = self.get_angles(
            np.arange(position)[:, np.newaxis],
            np.arange(d_model)[np.newaxis, :],
            d_model
        )

        angle_rads[:, 0::2] = np.sin(angle_rads[:, 0::2])
        angle_rads[:, 1::2] = np.cos(angle_rads[:, 1::2])

        pos_encoding = angle_rads[np.newaxis, ...]

        return tf.cast(pos_encoding, dtype=tf.float32)

    def call(self, inputs):
        return inputs + self.pos_encoding[:, :tf.shape(inputs)[1], :]


# transformer model, with additional blocks (and ability to add blocks)
def build_transformer_model(num_blocks=4):  
    inputs = tf.keras.layers.Input(shape=(max_length,))
    x = tf.keras.layers.Embedding(vocab_size, embedding_dim)(inputs)
    x = PositionalEncoding(max_length, embedding_dim)(x)

    # Adding multiple Transformer blocks
    for _ in range(num_blocks):
        x = TransformerBlock()(x, training=True)

    x = tf.keras.layers.Dense(vocab_size, activation='softmax')(x)
    
    return tf.keras.Model(inputs=inputs, outputs=x)

model = build_transformer_model()

model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy', 
              metrics=['accuracy'])


model.summary() # This model has 521,963 parameters

# This model achieved 97.2% val accuracy
history = model.fit(X_train, y_train, batch_size=32, epochs=5,
          validation_data=(X_test, y_test))


#plot test and validation accuracy history
plt.plot(history.history['accuracy'], label='accuracy')
plt.plot(history.history['val_accuracy'], label = 'val_accuracy')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.ylim([0.5, 1])
plt.legend(loc='lower right')

model.save('codon_overlap_tf_transformer_saved_model/my_model_new_dataset_20240111_1')  

my_model_combined_data_1 = tf.keras.models.load_model('codon_overlap_tf_transformer_saved_model/my_model_new_dataset_20240111_1')

