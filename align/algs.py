import numpy as np

# Parent alignment class with methods shared between both smith-waterman & needleman-wunsch
class PairwiseAligner:

	# Initializes PairwiseAligner class with the following parameters:
	# score matrix, gap opening penalty, & gap extension penalty
	def __init__(self, score_matrix_type, gap_opening_penalty, gap_extension_penalty):
		self.score_matrix_type = score_matrix_type
		self.gap_opening_penalty = gap_opening_penalty
		self.gap_extension_penalty = gap_extension_penalty

	# Takes in the name of a fasta file, finds that file, reads in the contents,
	# and returns the sequence contained in the file
	def load_fasta(self, fasta_name):
		seq = ""
		fasta_file = "sequences/" + fasta_name
		with open(fasta_file) as ff:
			for line in ff:
				# Ignore header
				if not line.startswith(">"):
					# Concatenate each line of sequences
					seq += line.strip().upper()
		return seq

	# Helper function for load_scoring_matrix() that makes a dictionary of residue
	# column indices for the scoring matrix
	def make_residue_indices(self, residue_string):
		residue_list = residue_string.split()
		self.residue_indices = {res: int(ind) for ind, res in enumerate(residue_list)}

	# Load the score matrix that was designated in init
	def load_scoring_matrix(self):
		score_matrix_file = "scoring_matrices/" + self.score_matrix_type + ".mat"
		with open(score_matrix_file) as smf:
			# Read in each line and remove whitespace
			lines = [line.strip() for line in smf]
			# Only keep lines that start with alphanumeric or "-" for negative numbers
			lines = [l for l in lines if l[0].isalpha() or l[0].isnumeric() or l.startswith("-")]
			# Index amino acid columns
			self.make_residue_indices(lines[0])
			# Add each score matrix row to a numpy array representing the score matrix
			score_matrix = np.loadtxt(l for l in lines[1:])
			# Save as a matrix of ints
			self.score_matrix = score_matrix.astype(np.int)

	# Make each alignment matrix
	# Description of each matrix:
	# - 3 matrices for cell values: match/mismatch (m), gap in seq1 (ix), gap in seq2 (iy))
	# - 6 matrices for pointers: pl_m, pl_ix, pl_iy, pt_m, pt_ix, pt_iy
	# - the pl matrices are location pointers & the pt matrices are type pointers
	# - the possible values for the pt matrices are m, ix and iy
	# - the possible values for the pl matrices are diag, up, and left
	# (for example, if pl_m[2,2]=="diag" & pl_m[2,2]=="ix", that means that the value for m
	# at cell [2,2] comes from the value for ix at [1,1], so there's a pointer from m[2,2]
	# to ix[1,1])
	# Cells of m, ix, and iy are initialized as -inf
	# Cells of pl_m, pl_ix, pl_iy, pt_m, pt_ix, pt_iy are initialized as j (for jackie :))
	# Visually, seq1 is along the top of the matrix and seq2 is along the left of the matrix
	def make_alignment_matrices(self):
		# Make matrices
		self.m = np.full((self.num_rows, self.num_cols), float("-inf"), dtype=float)
		self.ix = np.full((self.num_rows, self.num_cols), float("-inf"), dtype=float)
		self.iy = np.full((self.num_rows, self.num_cols), float("-inf"), dtype=float)
		self.pl_m = np.full((self.num_rows, self.num_cols), "j", dtype="object")
		self.pl_ix = np.full((self.num_rows, self.num_cols), "j", dtype="object")
		self.pl_iy = np.full((self.num_rows, self.num_cols), "j", dtype="object")
		self.pt_m = np.full((self.num_rows, self.num_cols), "j", dtype="object")
		self.pt_ix = np.full((self.num_rows, self.num_cols), "j", dtype="object")
		self.pt_iy = np.full((self.num_rows, self.num_cols), "j", dtype="object")
		# Set up pointers for left column of ix and top row of iy
		self.pl_ix[1:, 0] = "up"
		self.pt_ix[1:, 0] = "ix"
		self.pl_iy[0, 1:] = "left"
		self.pt_iy[0, 1:] = "iy"

	# Initialize alignment matrices
	# Passed to child classes
	def init_alignment_matrices(self):
		pass

	# Helper function for populate_alignment_matrices() that evaluates the neighbors for a given cell
	# and updates that cell with the correct score and pointers
	# The sw_m flag indicates if 0 should be added to the m neighbors list (for smith-waterman)
	def update_cell(self, i, j, neighbors, matrix, pt_matrix, pl_matrix, location, sw_m=False):
		# Sort neighbors (decreasing), get the max score (0th tuple), update matrix, and save pointers accordingly
		if sw_m == True:
			neighbors.append(("j", 0)) # back to j for jackie :) since no pointer needed
		# If there is a tie, the order is: diag > left > up
		neighbors.sort(key=lambda n: (n[1], n[0]), reverse=True)
		matrix[i, j] = neighbors[0][1]
		pt_matrix[i, j] = neighbors[0][0]
		pl_matrix[i, j] = location

	# Fill in each alignment matrix
	def populate_alignment_matrices(self, sw=False):
		# Iterate through each row (0th row is already set up)
		for i in range(1, self.num_rows):
			# Iterate through each column (0th column is already set up)
			for j in range(1, self.num_cols):
				# Find the value for m at that cell
				# s is the match/mismatch score
				s = self.score_matrix[self.residue_indices[self.seq1[j - 1]], self.residue_indices[self.seq2[i - 1]]]
				# m_neighbors is a list of tuples containing the neighbor name
				# and the score associated with choosing that neighbor
				m_neighbors = [("m", self.m[i - 1, j - 1] + s), ("ix", self.ix[i - 1, j - 1] + s),
							 ("iy", self.iy[i - 1, j - 1] + s)]
				# For smith-waterman, add 0 as an option for the value of the cell in m
				if sw == True:
					self.update_cell(i, j, m_neighbors, self.m, self.pt_m, self.pl_m, "diag", sw_m=True)
				else:
					self.update_cell(i, j, m_neighbors, self.m, self.pt_m, self.pl_m, "diag")
				# Find the value for ix at that cell
				ix_neighbors = []
				ix_neighbors.append(("m", self.m[i - 1, j] + self.gap_opening_penalty + self.gap_extension_penalty))
				ix_neighbors.append(("ix", self.ix[i - 1, j] + self.gap_extension_penalty))
				self.update_cell(i, j, ix_neighbors, self.ix, self.pt_ix, self.pl_ix, "up")
				# Find the value for ix at that cell
				iy_neighbors = []
				iy_neighbors.append(("m", self.m[i, j - 1] + self.gap_opening_penalty + self.gap_extension_penalty))
				iy_neighbors.append(("iy", self.iy[i, j - 1] + self.gap_extension_penalty))
				self.update_cell(i, j, iy_neighbors, self.iy, self.pt_iy, self.pl_iy, "left")

	# Calls method to populate alignment matrices with flag for NW or SW
	# Passed to child classes
	def call_populate_alignment_matrices(self):
		pass

	# Get alignment score between two sequences using alignment matrices
	# Passed to child classes
	def get_alignment_score(self):
		pass

	# Get traceback start
	# Passed to child classes
	def get_traceback_start(self):
		pass

	# Check traceback end condition
	# Passed to child classes
	def check_traceback_end_condition(self, i, j):
		pass

	# Get alignment between two sequences by tracing back through alignment matrices
	def get_alignment(self):
		# Make a dictionary to map from pointers to matrices (for both type & location)
		type_map = {}
		type_map["m"] = [self.m, self.pt_m, self.pl_m]
		type_map["ix"] = [self.ix, self.pt_ix, self.pl_ix]
		type_map["iy"] = [self.iy, self.pt_iy, self.pl_iy]
		type_map["j"] = []
		loc_map = {"diag": (-1,-1), "up": (-1,0), "left": (0,-1)}
		# Get starting point
		[start_i, start_j, start_type] = self.get_traceback_start()
		i = start_i
		j = start_j
		curr_type = start_type
		# Save current positions in sequences
		seq1_curr_pos = j - 1
		seq2_curr_pos = i - 1
		# Do traceback
		self.alignment = ["", ""]
		while self.check_traceback_end_condition(i, j) == False:
			# Update alignment
			if curr_type == "m":
				self.alignment[0] += self.seq1[seq1_curr_pos]
				self.alignment[1] += self.seq2[seq2_curr_pos]
				seq1_curr_pos -= 1
				seq2_curr_pos -= 1
			elif curr_type == "ix":
				self.alignment[0] += "-"
				self.alignment[1] += self.seq2[seq2_curr_pos]
				seq2_curr_pos -= 1
			elif curr_type == "iy":
				self.alignment[0] += self.seq1[seq1_curr_pos]
				self.alignment[1] += "-"
				seq1_curr_pos -= 1
			# Find moves to next cell location
			next_cell_location = loc_map[type_map[curr_type][2][i, j]]
			# Find & update type of next cell location
			curr_type = type_map[curr_type][1][i, j]
			# Update location
			i += next_cell_location[0]
			j += next_cell_location[1]
		# Reverse the alignment strings since they were generated backwards
		self.alignment[0] = self.alignment[0][::-1]
		self.alignment[1] = self.alignment[1][::-1]

	# Takes in 2 sequences and returns an alignment
	# The returned alignment is a 2-element list of strings corresponding to the aligned seq1 & seq2
	def align(self, seq1, seq2):
		# Save sequences
		self.seq1 = seq1
		self.seq2 = seq2
		# Make & initialize alignment matrices (visually, seq1 is on top and seq 2 is on left)
		self.num_rows = len(self.seq2) + 1
		self.num_cols = len(self.seq1) + 1
		self.make_alignment_matrices()
		self.init_alignment_matrices()
		# Fill in alignment matrices
		self.call_populate_alignment_matrices()
		# Get alignment score
		self.get_alignment_score()
		# Get alignment
		self.get_alignment()
		# Return the alignment, which is a 2 element list [alignment of seq1, alignment of seq2]
		# where dashes (-) represent gaps and underscores (_) represent offsets for local (SW) alignment
		return self.alignment

	# Takes in 2 sequences and calculates an alignment score without doing a traceback or generating alignment output
	# (useful for part 2 of the assignment)
	def score(self, seq1, seq2):
		# Save sequences
		self.seq1 = seq1
		self.seq2 = seq2
		# Make alignment & initialize alignment matrices
		self.num_rows = len(self.seq2) + 1
		self.num_cols = len(self.seq1) + 1
		self.make_alignment_matrices()
		self.init_alignment_matrices()
		# Fill in alignment matrices
		self.call_populate_alignment_matrices()
		# Get alignment score
		self.get_alignment_score()


# Child alignment class for smith-waterman alignments
class SmithWaterman(PairwiseAligner):

	# Initialize alignment matrices
	# Inherited from parent class
	def init_alignment_matrices(self):
		# Initialize values for m
		# (other values are already initialized)
		self.m[0, :] = 0
		self.m[:, 0] = 0

	# Calls method to populate alignment matrices with flag for NW or SW
	# Inherited from parent class
	def call_populate_alignment_matrices(self):
		self.populate_alignment_matrices(sw=True)

	# Get alignment score between two sequences using alignment matrices
	# Inherited from parent class
	def get_alignment_score(self):
		# Find maximum m value
		self.alignment_score = np.amax(self.m)

	# Get traceback start
	# Inherited from parent class
	def get_traceback_start(self):
		# Find locations of maximum m value(s)
		max_m_indices = np.where(self.m == np.amax(self.m))
		max_m_locations = list(zip(max_m_indices[0], max_m_indices[1]))
		# Take the 0th (arbitrary, like life)
		start_i = max_m_locations[0][0]
		start_j = max_m_locations[0][1]
		return [start_i, start_j, "m"]

	# Check traceback end condition
	# Inherited from parent class
	def check_traceback_end_condition(self, i, j):
		return (self.m[i,j] == 0)

# Child alignment class for needleman-wunsch alignments
class NeedlemanWunsch(PairwiseAligner):

	# Initialize alignment matrices
	# Inherited from parent class
	def init_alignment_matrices(self):
		# Initialize values for top left corner of m, left column of ix, and top row of iy
		# (other values are already initialized)
		self.m[0, 0] = 0
		self.ix[:, 0] = [self.gap_opening_penalty + self.gap_extension_penalty * i for i in range(self.num_rows)]
		self.iy[0, :] = [self.gap_opening_penalty + self.gap_extension_penalty * j for j in range(self.num_cols)]

	# Calls method to populate alignment matrices with flag for NW or SW
	# Inherited from parent class
	def call_populate_alignment_matrices(self):
		self.populate_alignment_matrices()

	# Get alignment score between two sequences using alignment matrices
	# Inherited from parent class
	def get_alignment_score(self):
		# Get values in the bottom right cell
		m_val = self.m[self.num_rows-1, self.num_cols-1]
		ix_val = self.ix[self.num_rows-1, self.num_cols-1]
		iy_val = self.iy[self.num_rows-1, self.num_cols-1]
		# Take max
		self.alignment_score = max(m_val, ix_val, iy_val)

	# Get traceback start
	# Inherited from parent class
	def get_traceback_start(self):
		# start i & j are the bottom right cell for needleman-wunsch
		start_i = self.num_rows - 1
		start_j = self.num_cols - 1
		# start_options is a list of tuples with the type & value for the bottom right cell
		start_options = []
		start_options.append(("m", self.m[self.num_rows-1, self.num_cols-1]))
		start_options.append(("ix", self.iy[self.num_rows - 1, self.num_cols - 1]))
		start_options.append(("iy", self.ix[self.num_rows - 1, self.num_cols - 1]))
		# Sort (decreasing), get the max score (0th tuple), and save the start type for the traceback
		start_options.sort(key=lambda n: n[1], reverse=True)
		start_type = start_options[0][0]
		return [start_i, start_j, start_type]

	# Check traceback end condition
	# Inherited from parent class
	def check_traceback_end_condition(self, i, j):
		return (i == 0 and j == 0)
