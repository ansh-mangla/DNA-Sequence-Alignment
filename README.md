````markdown
# DNA Sequence Alignment

This project aligns two sequences either **globally** or **locally**, based on the **Needleman-Wunsch** and **Smith-Waterman** algorithms implemented in Python using libraries like `numpy` and `pandas`.

## Motivation
> I started this project to begin working on bioinformatics projects and to deepen my understanding of fundamental tools, such as sequence alignment.

---

## Testing

The project is tested by comparing its results against the **Biopython** module in `test.py`. Several files are provided for testing purposes in the `data` folder.

---

## Features

1. **Custom scoring:** You can provide your own scoring scheme. If not provided, the defaults are: match = 1, mismatch = -1, gap = -2.  
2. **Console visualization:** The alignment is printed in the console for clear visualization.  
3. **Single function interface:** All sequence alignments are performed using:

```python
pairwise_alignment(seq1, seq2, mode, align=True)
````

It returns a Python dictionary with the following structure:

```python
alignment = {
    "seq1": str,
    "seq2": str,
    "align_seq1": str,
    "align_seq2": str,
    "score": int,
    "stats": dict,
    "matrix": numpy.ndarray
}

stats = {
    "match": int,
    "mismatch": int,
    "seq1_gap": int,
    "seq2_gap": int,
    "identity": float
}
```

4. The `align` argument in `pairwise_alignment(seq1, seq2, mode, align=True)` determines whether the alignment is printed in the console.
5. After each run, a table is printed showing statistics related to the alignment.
6. For shorter sequences (length < 20), a function called `heat_map(s_matrix, seq1: str, seq2: str)` generates a heatmap for the scoring matrix.

---

## Usage

Import the function and specify your sequences:

```python
from align import pairwise_alignment

seq1 = "AGCTACGATCGA"
seq2 = "AGCTGCGATA"
my_align = pairwise_alignment(seq1, seq2, mode="global")
```

The repository’s **data** folder contains sample sequences to use in `test.py`.

---

## Algorithm Details

### Needleman-Wunsch (Global Alignment)

A dynamic programming algorithm that finds the optimal alignment over the entire length of two sequences.

* Initializes the scoring matrix with cumulative gap penalties.
* Traces back from the bottom-right cell to build the full alignment.

### Smith-Waterman (Local Alignment)

Also based on dynamic programming but focuses on finding the highest-scoring local region between two sequences.

* First row and column initialized to zeros to allow local alignments.
* Negative scores are set to zero during matrix filling.
* Traces back from the highest-scoring cell until zero is reached.

Both algorithms allow **customizable scoring schemes** for matches, mismatches, and gaps.

---

## Sample Input and Output

### Input

```python
seq1 = "AGCTACGATCGA"
seq2 = "AGCTGCGATA"
my_align = pairwise_alignment(seq1, seq2, mode="global")
```

### Output (short sequences)

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

### Output (larger sequences)

*(Sequence not shown for readability)*

```text
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

---

## Potential Improvements and Future Prospects

1. Currently uses **linear gap penalties**; affine gap penalties could improve alignment scoring.
2. Supports **only pairwise sequence alignment**; extending to multiple sequence alignment could be a future improvement.

```
