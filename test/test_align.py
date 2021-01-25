import pytest
import numpy as np
from align import algs

@pytest.fixture
def some_relevant_data():
	return np.ones(10)

# For three different fasta files, check that the correct fasta file is inputted
# and the correct fasta sequence is outputted
def test_fasta_io():
	pa = algs.PairwiseAligner("PAM100", -10, -3)
	seq1 = pa.load_fasta("prot-0004.fa")
	correct_seq1 = "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"
	seq2 = pa.load_fasta("prot-0915.fa")
	correct_seq2 = "MKKATCLTDDQRWQSVLARDPNADGEFVFAVRTTGIFCRPSCRARHALRENVSFYANASEALAAGFRPCKRCQPDKANPRQHRLDKITHACR"
	# This one has lowercase letters so it checks the lower to upper conversion too
	seq3 = pa.load_fasta("prot-0819.fa")
	correct_seq3 = "ALLSFERKYRVRGGTLIGGDLFDFWVGPYFVGFFGVSAIFFIFLGVSLIGYAASQGPTWDPFAISINPPDLKYGLGAAPLLEGGFWQAITVCALGAFISWMLREVEISRKLGIGWHVPLAFCVPIFMFCVLQVFRPLLLGSWGHAFPYGILSHLDWVNNFGYQYLNWHYNPGHMSSVSFLFVNAMALGLHGGLILSVANPGDGDKVKTAEHENQYFRDVVGYSIGALSIHRLGLFLASNIFLTGAFGTIASGPFXXXXXXXXXXXXLDIPFWS"
	# Check that they match
	assert seq1 == correct_seq1
	assert seq2 == correct_seq2
	assert seq3 == correct_seq3

# Helper function to test scoring matrix I/O that hardcodes the correct one
def make_correct_scoring_matrix():
	# This is the entire BLOSUM50 matrix
	row0 = [5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0, -2, -1, -1, -5]
	row1 = [-2, 7, -1, -2, -4, 1, 0, -3, 0, -4, -3, 3, -2, -3, -3, -1, -1, -3, -1, -3, -1, 0, -1, -5]
	row2 = [-1, -1, 7, 2, -2, 0, 0, 0, 1, -3, -4, 0, -2, -4, -2, 1, 0, -4, -2, -3, 4, 0, -1, -5]
	row3 = [-2, -2, 2, 8, -4, 0, 2, -1, -1, -4, -4, -1, -4, -5, -1, 0, -1, -5, -3, -4, 5, 1, -1, -5]
	row4 = [-1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -2, -5]
	row5 = [-1, 1, 0, 0, -3, 7, 2, -2, 1, -3, -2, 2, 0, -4, -1, 0, -1, -1, -1, -3, 0, 4, -1, -5]
	row6 = [-1, 0, 0, 2, -3, 2, 6, -3, 0, -4, -3, 1, -2, -3, -1, -1, -1, -3, -2, -3, 1, 5, -1, -5]
	row7 = [0, -3, 0, -1, -3, -2, -3, 8, -2, -4, -4, -2, -3, -4, -2, 0, -2, -3, -3, -4, -1, -2, -2, -5]
	row8 = [-2, 0, 1, -1, -3, 1, 0, -2, 10, -4, -3, 0, -1, -1, -2, -1, -2, -3, 2, -4, 0, 0, -1, -5]
	row9 = [-1, -4, -3, -4, -2, -3, -4, -4, -4, 5, 2, -3, 2, 0, -3, -3, -1, -3, -1, 4, -4, -3, -1, -5]
	row10 = [-2, -3, -4, -4, -2, -2, -3, -4, -3, 2, 5, -3, 3, 1, -4, -3, -1, -2, -1, 1, -4, -3, -1, -5]
	row11 = [-1, 3, 0, -1, -3, 2, 1, -2, 0, -3, -3, 6, -2, -4, -1, 0, -1, -3, -2, -3, 0, 1, -1, -5]
	row12 = [-1, -2, -2, -4, -2, 0, -2, -3, -1, 2, 3, -2, 7, 0, -3, -2, -1, -1, 0, 1, -3, -1, -1, -5]
	row13 = [-3, -3, -4, -5, -2, -4, -3, -4, -1, 0, 1, -4, 0, 8, -4, -3, -2, 1, 4, -1, -4, -4, -2, -5]
	row14 = [-1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -2, -5]
	row15 = [1, -1, 1, 0, -1, 0, -1, 0, -1, -3, -3, 0, -2, -3, -1, 5, 2, -4, -2, -2, 0, 0, -1, -5]
	row16 = [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 2, 5, -3, -2, 0, 0, -1, 0, -5]
	row17 = [-3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1, 1, -4, -4, -3, 15, 2, -3, -5, -2, -3, -5]
	row18 = [-2, -1, -2, -3, -3, -1, -2, -3, 2, -1, -1, -2, 0, 4, -3, -2, -2, 2, 8, -1, -3, -2, -1, -5]
	row19 = [0, -3, -3, -4, -1, -3, -3, -4, -4, 4, 1, -3, 1, -1, -3, -2, 0, -3, -1, 5, -4, -3, -1, -5]
	row20 = [-2, -1, 4, 5, -3, 0, 1, -1, 0, -4, -4, 0, -3, -4, -2, 0, 0, -5, -3, -4, 5, 2, -1, -5]
	row21 = [-1, 0, 0, 1, -3, 4, 5, -2, 0, -3, -3, 1, -1, -4, -1, 0, -1, -2, -2, -3, 2, 5, -1, -5]
	row22 = [-1, -1, -1, -1, -2, -1, -1, -2, -1, -1, -1, -1, -1, -2, -2, -1, 0, -3, -1, -1, -1, -1, -1, -5]
	row23 = [-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1]
	sm = np.array([row0, row1, row2, row3, row4, row5, row6, row7, row8, row9, row10, row11, row12, row13, row14, row15, row16, row17, row18, row19, row20, row21, row22, row23])
	return(sm)

# For a given scoring matrix type, check that the correct matrix is loaded and saved
def test_scoring_matrix_io():
	pa = algs.PairwiseAligner("BLOSUM50", -5, -1)
	pa.load_scoring_matrix()
	# Call helper to hardcode correct matrix
	correct_scoring_matrix = make_correct_scoring_matrix()
	# Also as a bonus, make sure PAM250 loads without error because this one was problematic
	# (no need for assert statement for this one because if it doesn't run successfully,
	# this test automatically fails with error anyway)
	pa_bonus = algs.PairwiseAligner("PAM250", -5, -1)
	pa_bonus.load_scoring_matrix()
	# Check that they match (BLOSUM50)
	assert np.array_equal(pa.score_matrix, correct_scoring_matrix)

# For 4 different sequences, check that the alignment of each with itself is correct
def test_identical():
	# Sequences
	seq1 = "ZEAL"
	seq2 = "WEIRD"
	seq3 = "VLLDGNGEVVQNGGTYYLLPQVWAQGGGVQLAKTGEETCPLTVVQSPNELSDGKPIRIES"
	seq4 = "RLRSAFIPDDDKVRIGFAYAPKCAPSPWWTVVEGLSVKLSEDESTQFDYPFKFEQVSDQLHSYKLLYCEGKHEKCASIGINRDQKGYRRLVVTEDYPLTVVLKKDE"
	# Correct alignments to check against
	correct_align1 = [seq1, seq1]
	correct_align2 = [seq2, seq2]
	correct_align3 = [seq3, seq3]
	correct_align4 = [seq4, seq4]
	# Test smith-waterman
	sw = algs.SmithWaterman("BLOSUM62", -13, -4)
	sw.load_scoring_matrix()
	# Alignments
	sw_align1 = sw.align(seq1, seq1)
	sw_align2 = sw.align(seq2, seq2)
	sw_align3 = sw.align(seq3, seq3)
	sw_align4 = sw.align(seq4, seq4)
	# Test needleman-wunsch
	nw = algs.NeedlemanWunsch("PAM250", -9, -2)
	nw.load_scoring_matrix()
	# Alignments
	nw_align1 = nw.align(seq1, seq1)
	nw_align2 = nw.align(seq2, seq2)
	nw_align3 = nw.align(seq3, seq3)
	nw_align4 = nw.align(seq4, seq4)
	# Check that they match
	assert sw_align1 == correct_align1
	assert sw_align2 == correct_align2
	assert sw_align3 == correct_align3
	assert sw_align4 == correct_align4
	assert nw_align1 == correct_align1
	assert nw_align2 == correct_align2
	assert nw_align3 == correct_align3
	assert nw_align4 == correct_align4

def test_alignment_score():
	assert True
