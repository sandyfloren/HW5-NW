# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._scores_mat = None

        # Init matrices for backtrace procedure
        self._back = None
  
        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Handle case of empty strings
        if seqA == "":
            if seqB == "":
                return 0., "", ""
            else:
                score = self.gap_open + self.gap_extend * (len(seqB)-1)
                return score, "-" * len(seqB), seqB
        elif seqB == "":
            score = self.gap_open + self.gap_extend * (len(seqA)-1)
            return score, seqA, "-" * len(seqA)

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # Initialize matrix private attributes for use in alignment
        # Create matrices for alignment scores, gaps, and backtracing
        self._scores_mat = np.empty((len(seqA)+1, len(seqB)+1))

        self._back = np.empty((len(seqA)+1, len(seqB)+1), dtype=object)

        
        # Global alignment 
        # Assign first row and column based on gap opening and extension penalty
        self._scores_mat[:, 0] = [0., *[self.gap_open + x * self.gap_extend for x in range(1, len(seqA)+1) ]]

        self._scores_mat[0, :] = [0., *[self.gap_open + x * self.gap_extend for x in range(1, len(seqB)+1) ]]

        self._back[0] = [(0, 0), *[(0, x) for x in range(len(seqB))]]
        self._back[:, 0] = [(0, 0), *[(x, 0) for x in range(len(seqA))]]
        

        # Helper function to compute (sim(seqA[i], '-'), sim('-', seqB[j]), sim(seqA[i], seqB[j]))
        def similarities(i, j):
            # c b
            # a 

            a_last = self._back[i, j-1]
            b_last = self._back[i-1, j]

            # Check if a gap is being extended or opened
            if a_last == (i, j-2): # Gap extension
                a_gap = self.gap_extend
            else:
                a_gap = self.gap_open + self.gap_extend

            if b_last == (i-2, j):  # Gap extension
                b_gap = self.gap_extend
            else:
                b_gap = self.gap_open + self.gap_extend

            # Compute similarity between seqA[i] and seqB[j]
            sim_i_j = self.sub_dict[(self._seqA[i-1], self._seqB[j-1])]

            return a_gap, b_gap, sim_i_j



        # Compute alignment score for each position in matrix other than first row and column 
        for i in range(1, len(seqA) + 1):     
            for j in range(1, len(seqB) + 1):
                sims = similarities(i, j)

                scores = [
                    self._scores_mat[i, j-1] + sims[0],
                    self._scores_mat[i-1, j] + sims[1],
                    self._scores_mat[i-1, j-1] + sims[2]
                    ]

                self._scores_mat[i, j] = max(scores)
                best = np.argmax(scores)

                if best == 0:
                    self._back[i, j] = (i, j-1)
                elif best == 1:
                    self._back[i, j] = (i-1, j)
                else:
                    self._back[i, j] = (i-1, j-1)
        		    
        self.alignment_score = self._scores_mat[i, j]
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Start from bottom right   
        i = len(self._seqA)
        j = len(self._seqB)
        
        while i > 0 or j > 0:
            i_prev, j_prev = self._back[i, j]

            if i_prev == i - 1: 
                self.seqA_align = self._seqA[i-1] + self.seqA_align
            else:
                self.seqA_align = '-' + self.seqA_align

            if j_prev == j - 1:
                self.seqB_align = self._seqB[j-1] + self.seqB_align
            else:
                self.seqB_align = '-' + self.seqB_align

            i, j = i_prev, j_prev
        return (self.alignment_score, self.seqA_align, self.seqB_align)
    

def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
