from align import *
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq

# sequnces
seq1 = str(SeqIO.read("data/COX1 Human CDS.txt", 'fasta').seq)
seq2 = str(SeqIO.read("data/COX1 Mouse CDS.txt", 'fasta').seq)
mode = "global"
# my function
scheme = {
    "match": 1,
    "mismatch": -1,
    "gap": -2
}
my_align = pairwise_alignment(
    seq1, seq2, mode=mode, align=False, scheme=scheme)

# biopython
seq1 = Seq(seq1)
seq2 = Seq(seq2)

aligner = Align.PairwiseAligner()
aligner.mode = mode

aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -2
aligner.extend_gap_score = -2

alignments = aligner.align(seq1, seq2)

print(alignments[0].score)


# comparing scores

if my_align["score"] == alignments[0].score:
    print("Result matched.")
