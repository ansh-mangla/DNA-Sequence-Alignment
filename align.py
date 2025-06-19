import numpy as np
import pandas as pd


def pairwise_alinger(seq1, seq2, scheme={"match": 1, "mismatch": -1, "gap": -2}):
    s_matrix = np.zeros((len(seq2)+1, len(seq1)+1))

    n = len(seq1)
    l = len(seq2)

    s_matrix[0, 1:] = np.arange(1, n + 1) * scheme["gap"]
    s_matrix[1:, 0] = np.arange(1, l + 1) * scheme["gap"]

    # print(s_matrix)

    for j in range(1, len(seq1)+1):
        for i in range(1, len(seq2)+1):

            diago = s_matrix[i-1, j-1] + \
                (scheme["match"] if seq1[j - 1] ==
                 seq2[i-1] else scheme["mismatch"])
            top = s_matrix[i-1, j] + scheme["gap"]
            left = s_matrix[i, j-1] + scheme["gap"]

            options = [diago, top, left]
            s_matrix[i, j] = max(options)
    # print(s_matrix)

    # the scoring matrix to visualization
    # s_matrix_df = pd.DataFrame(
    #     s_matrix,
    #     index=[seq2[n-1] if n != 0 else 0 for n in range(l+1)],
    #     columns=[seq1[n-1] if n != 0 else 0 for n in range(n+1)]
    # ).astype(int)
    # print(s_matrix_df)

    return s_matrix


def traceback(s_matrix, seq1, seq2, scheme={"match": 1, "mismatch": -1, "gap": -2}):
    i = len(seq2)
    j = len(seq1)
    align_seq1 = []
    align_seq2 = []
    while i != 0 or j > 0:

        current = s_matrix[i][j]

        diago = (
            s_matrix[i - 1, j - 1]
            + (scheme["match"] if i and j and seq1[j - 1]
               == seq2[i - 1] else scheme["mismatch"])
        )
        top = s_matrix[i - 1, j] + scheme["gap"]
        left = s_matrix[i, j - 1] + scheme["gap"]

        if current == diago and i and j:  # match or a mismatch
            align_seq1.append(seq1[j - 1])
            align_seq2.append(seq2[i - 1])
            i -= 1
            j -= 1
        elif current == top and i:  # gap in seq1
            align_seq1.append("-")
            align_seq2.append(seq2[i - 1])
            i -= 1
        else:  # left # gap in seq2
            align_seq1.append(seq1[j - 1])
            align_seq2.append("-")
            j -= 1

    align_seq1.reverse()
    align_seq2.reverse()
    return "".join(align_seq1), "".join(align_seq2), int(s_matrix[-1, -1])


def print_alignmetn(align_seq1, align_seq2, score):
    bar = []
    for n1, n2 in zip(align_seq1, align_seq2):
        if n1 == n2:
            bar.append("|")
        else:
            bar.append(" ")
    print(f"{" ".join(align_seq1)}\n{" ".join(bar)}\n{" ".join(align_seq2)}")
    print(f"score: {score}")


def global_aling(seq1, seq2, scheme={"match": 1, "mismatch": -1, "gap": -2}):
    matrix = pairwise_alinger(seq1, seq2, scheme)
    align1, align2, score = traceback(matrix, seq1, seq2)
    print_alignmetn(align1, align2, score)


scheme = {
    "match": 1,
    "mismatch": -1,
    "gap": -2
}

global_aling("AGCTACGATCGA", "AGCTGCGATA", scheme)
