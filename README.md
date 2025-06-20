# DNA sequence alignment

This projects aim to align two sequences either globally or locally, based on Needlmen-Wunsch and Smith-Waterman algorithm implemented in python using libraries like numpy and pandas.

**Motivation**:
"I started this project because i wanted to start working on bioinfromtics and related prjects and to deepen my understanding of fundamental tools used in it like sequence alignment."

## Testing

The project is tested in a unique way. The results are compared against the results of Biopython module. in `test.py`
Several files have been provided for testing purposes in the `data` folder.

## Features

1. The code alows for costum scoring. But in case scoreing not provided it assumes m = 1, mm = -1 and gap = -2
2. The alignment is also printed in the consol for accurate visualization
3. For all sequence alignmennnt a single function operties `paiwise_alignment(seq1, seq2, mode)` and returns a python `dict` dictionary. Including keys like:

```python
    alignment = {
        "seq1": str
        "seq2": str
        "align_seq1": str
        "align_seq2": str
        "score": int
        "stats": dict
        "matrix": numpy.ndarray
    }

    stat = {
        "match": int,
        "mismatch": int,
        "seq1_gap": int,
        "seq2_gap": int,
        "identity": float
    }

```

4. in the `paiwise_alignment(seq1, seq2, mode, align=True)` can specify an align argument which decides whether the alignment will be printed or not.
5. After each time the function is run irrespective of align argument it prints a table which gives stastics pertaining to the alignment.
6. For shorter sequences i.e. below length of 20 a function called `heat_map(s_matrix, seq1: str, seq2: str)` return a heatmap for the scoring matrix.

## Usage

Using the code is very simple, `from align.py import pairwise_alignment` and specify your sequences in the funciton.

Within the repo **data** folder is provided with some sample sequences to run in the **test.py** .

## Algorithm Details

This project implements two classic sequence alignment algorithms:

- **Needleman-Wunsch (Global Alignment):**  
  A dynamic programming algorithm that finds the optimal alignment between two sequences over their entire length. It initializes the scoring matrix with cumulative gap penalties and traces back from the bottom-right cell to build the full alignment.

- **Smith-Waterman (Local Alignment):**  
  Also based on dynamic programming, but focuses on finding the highest-scoring local region between two sequences. It initializes the first row and column with zeros to allow alignments to start anywhere, sets negative scores to zero during matrix filling, and traces back from the highest scoring cell until a zero score is reached.

Both algorithms use customizable scoring schemes for matches, mismatches, and gaps, and the traceback reconstructs the optimal alignment according to the scoring matrix.

## sample input and output

### Input

```python
seq1 = "AGCTACGATCGA"
seq2 = "AGCTGCGATA"
my_align = pairwise_alignment(seq1, seq2, mode="global")
```

### Output

```text
Seq A    1  A G C T A C G A T C G A
            | | | |   | | | |     |
Seq B    1  A G C T G C G A T - - A


           stats   value
Length of align:       12
          Score:        4
          Match:        9
     Mismatches:        1
    Gap in Seq1:        0
    Gap in Seq2:        2
     Total Gaps:        2
   Identity (%):       75
```

### For larger sequences you might expact an output

(Here sequence is not shown since it will hinder the visuals)

```test
           stats   value
Length of align:     1548
          Score:      859
          Match:     1208
     Mismatches:      331
    Gap in Seq1:        6
    Gap in Seq2:        3
     Total Gaps:        9
   Identity (%):       78
```

## Potential improvement and future prospects

1. The implementation works of linear gap panelites and does not consider afine gap penalties.
2. Only align two sequenes pairwise (no multiple sequence alignment)
