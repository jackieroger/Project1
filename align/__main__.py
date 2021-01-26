from align import algs

# Load "true" and "false" alignments
def load_pairs(pairs_file):
    # Read in pairs
    with open(pairs_file) as pf:
        # Read in each line, remove whitespace, and make list of tuples representing each pair
        pairs_files = [tuple(line.strip().split()) for line in pf]
    return pairs_files



# Evaluate alignment algorithms using a set of "true" and "false" alignments
def main():
    # Read in true & false alignments
    pospairs = load_pairs("scoring_matrices/Pospairs.txt")
    negpairs = load_pairs("scoring_matrices/Negpairs.txt")

main()