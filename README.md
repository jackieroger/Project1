# Project 1 - Sequence Alignment

![BuildStatus](https://github.com/jackieroger/Project1/workflows/HW1/badge.svg?event=push)

### Main
To run the code in align/\_\_main\_\_.py which contains the code for part 2, run the following command from the root directory of this project:
```
python -m align
```

### Testing
To run the unit tests for part 1, run the following command from the root directory of this project:
```
python -m pytest test/*
```

## Class methods & attributes

#### load_fasta()
- Takes in a string representing a path to a fasta file
- Returns a string representing the sequence contained in the file

#### load_scoring_matrix()
- Takes in no parameters
- Loads the appropriate scoring matrix based on the name of the scoring matrix that was passed in during object initialization
- Creates attributes: **score_matrix** (2d numpy array representing the user-indicated score matrix) and **residue_indices** (dictionary mapping between indices and their column/row in the score matrix)

#### align()
- Takes in 2 sequences (strings) to be aligned
- Returns an alignment (2-element list of strings representing the aligned version of each of the 2 sequences)
- Creates attributes: **alignment_score** (int representing alignment score for the 2 inputted sequences) and **alignment** (list containing the aligned versions of seq1 and seq2)

#### score()
- This method is similar to align(), but just scores an alignment without doing a traceback to find the actual alignment
- Takes in 2 sequences (strings) to be aligned
- Returns an alignment score (int)
- Creates attribute: **alignment_score** (int representing alignment score for the 2 inputted sequences)

## Example usage

Example using smith-waterman
```
from align import algs

# Create sw object
sw = algs.SmithWaterman("BLOSUM50", -3, -1)

# Load scoring matrix
sw.load_scoring_matrix()

# Load sequence data
seq1 = sw.load_fasta("test_data/prot-0004.fa")
seq2 = sq.load_fasta("test_data/prot-0915.fa")

# Do alignment
sw_alignment = sw.align(seq1, seq2) # run alignment
sw_alignment_score = sw.alignment_score

# Alternatively, just score alignment (without doing traceback to get alignment)
sw_alignment_score = sw.score(seq1, seq2)
```
