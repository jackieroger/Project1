# Project 1 - Sequence Alignment

![BuildStatus](https://github.com/jackieroger/Project1/workflows/HW1/badge.svg?event=push)

### Testing for part 1
To run the unit tests for part 1, run the following command from the root directory of this project:
```
python -m pytest test/*
```

### Responses (both parts) and code for part 2
These are all contained in **Jackie_Roger_BMI203_HW1**. Additional copies of the plots are in the **plots** directory and additional copies of some of the results data are in the **part2_results_data** directory.

## Class methods & attributes

#### \_\_init\_\_()
- This is called automatically when a new class object is created
- Takes in match/mismatch score matrix type (string), gap opening penalty (int), and gap extension penalty (int)
- Returns class object (either SmithWaterman or NeedlemanWunsch)

#### load_scoring_matrix()
- Takes in no parameters
- Loads the appropriate scoring matrix based on the name of the scoring matrix that was passed in during object initialization
- Creates attributes: **score_matrix** (2d numpy array representing the user-indicated score matrix) and **residue_indices** (dictionary mapping between indices and their column/row in the score matrix)

#### load_fasta()
- Takes in a string representing a path to a fasta file
- Returns a string representing the sequence contained in the file

#### align()
- Takes in 2 sequences (strings) to be aligned
- Returns an alignment (2-element list of strings representing the aligned version of each of the 2 sequences)
- Creates attributes: **alignment_score** (int representing alignment score for the 2 inputted sequences) and **alignment** (list containing the aligned versions of seq1 and seq2)
- Optional parameter: for needleman-wunsch alignment, add overlap=True to do an overlap alignment

#### score()
- This method is similar to align(), but just scores an alignment without doing a traceback to find the actual alignment
- Takes in 2 sequences (strings) to be aligned
- Returns an alignment score (int)
- Creates attribute: **alignment_score** (int representing alignment score for the 2 inputted sequences)
- Optional parameter: for needleman-wunsch alignment, add overlap=True to do an overlap alignment

## Example usage

Example using smith-waterman (run from project root directory)
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
sw_alignment = sw.align(seq1, seq2)
sw_alignment_score = sw.alignment_score

# Alternatively, just score alignment (without doing traceback to get alignment)
sw_alignment_score = sw.score(seq1, seq2)
```

### For the curious user...

If you'd like to explore the alignment score matrices generated when align() is called, the class attributes (2d numpy arrays) representing them are called:
- **m**: match/mismatch matrix
- **ix**: gap matrix for seq1
- **iy**: gap matrix for seq2