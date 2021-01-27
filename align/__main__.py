from align import algs

# Load fasta file names for alignments
def load_pairs(pairs_txt_file):
    # Read in pairs
    with open(pairs_txt_file) as pf:
        # Read in each line, remove whitespace, and make list of tuples representing each pair
        pairs_files = [tuple(line.strip().split()) for line in pf]
    return pairs_files

# Load sequences given fasta file names
def load_sequences(fasta_file_names):
    pa = algs.PairwiseAligner("BLOSUM50", -3, -1)
    sequence_pairs = [tuple([pa.load_fasta(f1), pa.load_fasta(f2)]) for (f1, f2) in fasta_file_names]
    return sequence_pairs

# Question 1
def do_question_1(allpairs):
    sw = algs.SmithWaterman("BLOSUM50", -11, -3)
    sw.load_scoring_matrix()
    scores = []
    for (s1, s2) in allpairs:
        s = sw.score(s1, s2)
        scores.append(s)
    return scores

# Evaluate alignment algorithms using a set of "true" and "false" alignments
def main():
    # Read in names of fasta files for true & false alignments
    pospairs_files = load_pairs("scoring_matrices/Pospairs.txt")
    negpairs_files = load_pairs("scoring_matrices/Negpairs.txt")
    # Read in sequences from those fasta files
    pospairs = load_sequences(pospairs_files)
    negpairs = load_sequences(negpairs_files)
    allpairs = pospairs + negpairs
    # Question 1
    q1_alignment_scores = do_question_1(allpairs)
    with open("part2_results_data/q1.txt", "w") as q1:
        q1.writelines("%s\n" % score for score in q1_alignment_scores)


main()