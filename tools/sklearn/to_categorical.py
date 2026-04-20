import argparse
import json
import warnings

import h5py
import numpy as np
import pandas as pd

warnings.simplefilter("ignore")


def _get_longest_sequence_length(fasta_file):
    max_len = 0
    max_id = None
    for name in fasta_file.keys():
        seq_len = len(fasta_file[name])
        if seq_len > max_len:
            max_len = seq_len
            max_id = name

    return max_len, max_id


def encode_dna_sequences(fasta_path, padding, outfile, outfile_matrix):
    from galaxy_ml.preprocessors import GenomeOneHotEncoder
    import pyfaidx

    seq_length = None
    fasta_file = pyfaidx.Fasta(fasta_path)
    if padding:
        seq_length, max_id = _get_longest_sequence_length(fasta_file)
        print("Longest sequence is %s with length %d" % (max_id, seq_length))
    print("Padding: {}".format(padding))
    X = np.arange(len(fasta_file.keys())).reshape(-1, 1)
    genome_encoder = GenomeOneHotEncoder(
        fasta_path=fasta_path, seq_length=seq_length, padding=padding
    )
    genome_encoder.fit(X)
    encoded_dna_sequences = genome_encoder.transform(X)
    flatted_enc_seqs = encoded_dna_sequences.flatten().reshape(encoded_dna_sequences.shape[0], -1)
    np.savetxt(outfile, np.asarray(flatted_enc_seqs, dtype=int), fmt="%d", delimiter="\t")
    with h5py.File(outfile_matrix, "w") as handle:
        handle.create_dataset("data", data=encoded_dna_sequences, compression="gzip")


def seq_to_kmers(sequence, k=3):
    sequence = sequence.upper()
    return [sequence[idx: idx + k] for idx in range(len(sequence) - k + 1)]


def build_kmer_vocabulary(sequences, k):
    vocabulary = {"<PAD>": 0, "<UNK>": 1}
    for sequence in sequences:
        for kmer in seq_to_kmers(sequence, k):
            if kmer not in vocabulary:
                vocabulary[kmer] = len(vocabulary)

    if len(vocabulary) == 2:
        raise ValueError(
            "No DNA k-mers were generated. Check that k is not longer than all sequences."
        )

    return vocabulary


def encode_sequence_kmers(sequence, vocabulary, k):
    return [vocabulary.get(kmer, vocabulary["<UNK>"]) for kmer in seq_to_kmers(sequence, k)]


def pad_encoded_sequences(encoded_sequences, pad_value=0):
    max_len = max(len(sequence) for sequence in encoded_sequences)
    return [sequence + [pad_value] * (max_len - len(sequence)) for sequence in encoded_sequences]


def encode_dna_kmers(fasta_path, k, outfile, outfile_vocab):
    import pyfaidx

    if k < 1:
        raise ValueError("k-mer size must be at least 1.")

    fasta_file = pyfaidx.Fasta(fasta_path)
    sequences = [str(fasta_file[name]).upper() for name in fasta_file.keys()]
    vocabulary = build_kmer_vocabulary(sequences, k)
    encoded_sequences = [
        encode_sequence_kmers(sequence, vocabulary, k) for sequence in sequences
    ]
    padded_sequences = np.asarray(
        pad_encoded_sequences(encoded_sequences, pad_value=vocabulary["<PAD>"]),
        dtype=int,
    )
    np.savetxt(outfile, padded_sequences, fmt="%d", delimiter="\t")
    with open(outfile_vocab, "w") as handle:
        json.dump(vocabulary, handle, indent=4, sort_keys=False)
        handle.write("\n")


def encode_labels(infile, input_header, outfile, num_classes=None):
    from keras.utils import to_categorical

    header = "infer" if input_header else None
    input_vector = pd.read_csv(infile, sep="\t", header=header)
    output_matrix = to_categorical(input_vector, num_classes=num_classes)
    np.savetxt(outfile, np.asarray(output_matrix, dtype=int), fmt="%d", delimiter="\t")


def main(args):
    task_type = args.encoder_task_type
    num_classes = args.num_classes
    header = "infer" if args.labels_header == "booltrue" else None
    padding = True if args.padding == "booltrue" else False
    sequence_encoding = args.sequence_encoding
    kmer_size = args.kmer_size

    if task_type == "label_encoder":
        encode_labels(args.labels_path, header, args.outfile, num_classes=num_classes)
    elif task_type == "dna_encoder":
        if sequence_encoding == "one_hot":
            encode_dna_sequences(args.fasta_path, padding, args.outfile, args.outfile_matrix)
        elif sequence_encoding == "kmer":
            encode_dna_kmers(args.fasta_path, kmer_size, args.outfile, args.outfile_vocab)
        else:
            raise ValueError("Unsupported DNA sequence encoding: %s" % sequence_encoding)
    else:
        raise ValueError("Unsupported encoder type: %s" % task_type)


if __name__ == "__main__":
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-l", "--labels_path", dest="labels_path")
    aparser.add_argument("-d", "--labels_header", dest="labels_header", default=False)
    aparser.add_argument("-t", "--encoder_task_type", dest="encoder_task_type", required=True)
    aparser.add_argument("-y", "--num_classes", dest="num_classes", type=int, default=None)
    aparser.add_argument("-p", "--padding", dest="padding", default="boolfalse")
    aparser.add_argument("-s", "--sequence_encoding", dest="sequence_encoding", default="one_hot")
    aparser.add_argument("-k", "--kmer_size", dest="kmer_size", type=int, default=3)
    aparser.add_argument("-f", "--fasta_path", dest="fasta_path")
    aparser.add_argument("-o", "--outfile", dest="outfile", required=True)
    aparser.add_argument("-m", "--outfile_matrix", dest="outfile_matrix")
    aparser.add_argument("-v", "--outfile_vocab", dest="outfile_vocab")
    args = aparser.parse_args()

    print(args)

    main(
        args
    )
