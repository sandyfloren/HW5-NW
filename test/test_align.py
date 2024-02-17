# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    # Note for TAs: it's not clear what is meant here by "asserting that you # have correctly filled out the your 3 alignment matrices." 
    # My NW implementation uses fewer matrices and a different strategy than 
    # the one outlined in the video, so here I just test the output of the # alignment for correctness.
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    aligner = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
    
    score, seqA, seqB = aligner.align(seq1, seq2)

    assert seqA == 'MYQR', 'seq1 alignment is incorrect'
    assert seqB == 'M-QR', 'seq2s alignment is incorrect'


def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    aligner = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
    
    score, seqA, seqB = aligner.align(seq3, seq4)

    assert score == 17, 'Alignment score for seq3 and seq4 should be 17'




