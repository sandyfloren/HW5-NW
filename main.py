# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    aligner = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)


    hs_gg_score, hs_gg_a, hs_gg_b = aligner.align(hs_seq, gg_seq)
    hs_mm_score, hs_mm_a, hs_mm_b = aligner.align(hs_seq, mm_seq)
    hs_br_score, hs_br_a, hs_br_b = aligner.align(hs_seq, br_seq)
    hs_tt_score, hs_tt_a, hs_tt_b = aligner.align(hs_seq, tt_seq)


    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print(hs_gg_a)
    print(hs_gg_b)
    print(hs_gg_score)
    print()
    print(hs_mm_a)
    print(hs_mm_b)
    print(hs_mm_score)
    print()
    print(hs_br_a)
    print(hs_br_b)
    print(hs_br_score)
    print()
    print(hs_tt_a)
    print(hs_tt_b)
    print(hs_tt_score)
    print()

if __name__ == "__main__":
    main()
