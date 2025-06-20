from align import *
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq


huamn = str(SeqIO.read("data/COX1 Human CDS.txt", 'fasta').seq)
yeast = str(SeqIO.read("data/COX1 Mouse CDS.txt", 'fasta').seq)


my_align = pairwise_alignment(huamn, yeast, mode="global")


# seq1 = Seq(huamn)
# seq2 = Seq(yeast)

# aligner = Align.PairwiseAligner()
# aligner.mode = "local"

# aligner.match_score = 1
# aligner.mismatch_score = -1
# aligner.open_gap_score = -2
# aligner.extend_gap_score = -2

# alignments = aligner.align(seq1, seq2)

# print(alignments[0])
# print(alignments[0].score)
