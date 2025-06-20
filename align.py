""" align two sequences either globally or locally, based on Needlmen-Wunsch and Smith-Waterman algorithm implemented in python using libraries like numpy and pandas"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# generte scoring matrix


def get_scoring_matrix(seq1, seq2, mode="global", scheme={"match": 1, "mismatch": -1, "gap": -2}):
    """Initilise and populate the scoring matrix based on the mode and sequnces.
    Returns: scoring matrix"""
    s_matrix = np.zeros((len(seq2)+1, len(seq1)+1))

    n = len(seq1)
    l = len(seq2)

    if mode == "global":
        s_matrix[0, 1:] = np.arange(1, n + 1) * scheme["gap"]
        s_matrix[1:, 0] = np.arange(1, l + 1) * scheme["gap"]

    for j in range(1, len(seq1)+1):
        for i in range(1, len(seq2)+1):

            diago = s_matrix[i-1, j-1] + \
                (scheme["match"] if seq1[j - 1] ==
                 seq2[i-1] else scheme["mismatch"])
            top = s_matrix[i-1, j] + scheme["gap"]
            left = s_matrix[i, j-1] + scheme["gap"]

            options = [diago, top, left]

            if mode == "global":
                s_matrix[i, j] = max(options)
            else:
                s_matrix[i, j] = max(options) if max(options) > 0 else 0

    return s_matrix
    n = len(seq1)
    l = len(seq2)
    s_matrix = np.zeros((l+1, n+1))
    # print(s_matrix)

    for j in range(1, len(seq1)+1):
        for i in range(1, len(seq2)+1):

            diago = s_matrix[i-1, j-1] + \
                (scheme["match"] if seq1[j - 1] ==
                    seq2[i-1] else scheme["mismatch"])
            top = s_matrix[i-1, j] + scheme["gap"]
            left = s_matrix[i, j-1] + scheme["gap"]

            options = [diago, top, left]
            s_matrix[i, j] = max(options) if max(options) > 0 else 0
    # print(s_matrix)

    # the scoring matrix to visualization
    s_matrix_df = pd.DataFrame(
        s_matrix,
        index=[seq2[n-1] if n != 0 else 0 for n in range(l+1)],
        columns=[seq1[n-1] if n != 0 else 0 for n in range(n+1)]
    ).astype(int)
    print(s_matrix_df)

    return s_matrix


def print_s_matrix(s_matrix: np.ndarray, seq1: str, seq2: str):
    """print the scoring matrix"""

    n = len(seq1)
    l = len(seq2)
    s_matrix_df = pd.DataFrame(
        s_matrix,
        index=[seq2[n-1] if n != 0 else 0 for n in range(l+1)],
        columns=[seq1[n-1] if n != 0 else 0 for n in range(n+1)]
    ).astype(int)
    print(s_matrix_df)


def traceback(s_matrix: np.ndarray, seq1: str, seq2: str, mode="global", scheme={"match": 1, "mismatch": -1, "gap": -2}):
    """
    Takes a scoring matrix and sequnces that are being aligned and traces back the to get the alignment itself. 
    Returns: Alignment 
    """
    i = len(seq2)
    j = len(seq1)
    align_seq1 = []
    align_seq2 = []

    if mode == "global":
        i = len(seq2)
        j = len(seq1)
    else:
        flat_index = s_matrix.argmax()
        i, j = np.unravel_index(flat_index, s_matrix.shape)

    score = int(s_matrix[i][j])

    stat = {
        "match": int(0),
        "mismatch": int(0),
        "seq1_gap": int(0),
        "seq2_gap": int(0),
        "identity": float(0)
    }

    while i > 0 or j > 0:

        current = s_matrix[i][j]

        diago = (
            s_matrix[i - 1, j - 1]
            + (scheme["match"] if i and j and seq1[j - 1]
               == seq2[i - 1] else scheme["mismatch"])
        )
        top = s_matrix[i - 1, j] + scheme["gap"]
        left = s_matrix[i, j - 1] + scheme["gap"]

        # stope when the local sequence alinged / matched
        if mode == "local" and current == 0:
            break

        # actual traceback
        if current == diago:  # match or a mismatch
            if s_matrix[i - 1, j - 1] + scheme["match"] == current:
                stat['match'] += 1
            else:
                stat['mismatch'] += 1
            align_seq1.append(seq1[j - 1])
            align_seq2.append(seq2[i - 1])
            i -= 1
            j -= 1
        elif current == top:  # gap in seq1
            stat['seq1_gap'] += 1
            align_seq1.append("-")
            align_seq2.append(seq2[i - 1])
            i -= 1
        else:  # left # gap in seq2
            stat['seq2_gap'] += 1
            align_seq1.append(seq1[j - 1])
            align_seq2.append("-")
            j -= 1

    align_seq1.reverse()
    align_seq2.reverse()

    stat["identity"] = round(
        (stat["match"] / len("".join(align_seq1))) * 100, 3)
    alignment = {
        "seq1": seq1,
        "seq2": seq2,
        "align_seq1": "".join(align_seq1),
        "align_seq2": "".join(align_seq2),
        "score": score,
        "stats": stat,
        "matrix": s_matrix
    }
    return alignment


# print the alignment itself
def print_alignment(alignment: dict, mode="global"):
    """
    Prints the actual alignment to the terminal
    """
    bar = []

    align_seq1 = alignment["align_seq1"]
    align_seq2 = alignment["align_seq2"]

    for n1, n2 in zip(align_seq1, align_seq2):
        if n1 == n2:
            bar.append("|")
        else:
            bar.append(" ")

    space_1 = ""
    for s in range(len(str(len(align_seq1)))):
        space_1 += " "
    space_1 += " "

    s_p_1 = alignment['seq1'].find(
        align_seq1.replace(" ", "").replace("-", ""))
    s_p_2 = alignment['seq2'].find(
        align_seq2.replace(" ", "").replace("-", ""))

    k = 50  # number of bases i want in each line
    print("\n\n")
    for n in range(0, len(align_seq1), k):

        space_seq_1 = str(n+1) + space_1[len(str(n)):]
        space_seq_2 = str(n+1) + space_1[len(str(n)):]

        if mode == "local":
            space_seq_1 = str(n + s_p_1) + space_1[len(str(n + s_p_1)):]
            space_seq_2 = str(n + s_p_2) + space_1[len(str(n + s_p_2)):]

        print("Seq A    " + space_seq_1 + " ".join(align_seq1[n:n+k]))
        print("         " + space_1 + " ".join(bar[n:n+k]))
        print("Seq B    " + space_seq_2 + " ".join((align_seq2[n:n+k])))
    print("\n")

# print the stats of the alignment


def get_stats(alignment: dict):
    """
    Prints the stastics related to the alignment 

    Arguments 
    alignment = {
        "seq1": str
        "seq2": str
        "align_seq1": str
        "align_seq2": str
        "score": int
        "stats": dict
        "matrix": numpy.ndarray
    }

    stats = {
        "match": int,
        "mismatch": int,
        "seq1_gap": int,
        "seq2_gap": int,
        "identity": float
    }
    Return: None 

    """
    align_seq1 = alignment["align_seq1"]
    align_seq2 = alignment["align_seq2"]
    score = alignment["score"]
    stats = alignment["stats"]
    identity = round((stats["match"] / len(align_seq1)) * 100, 3)

    string_of_stat = pd.DataFrame({
        "stats ":
            ["Length of align: ", "Score: ", "Match: ", "Mismatches: ", "Gap in Seq1: ",
                "Gap in Seq2: ", "Total Gaps: ", "Identity (%): "],
        "value ":
        map(int, [len(align_seq1), score, stats['match'], stats['mismatch'], stats["seq1_gap"],
                  stats["seq2_gap"], stats["seq1_gap"] + stats["seq2_gap"], identity]),
    }).to_string(index=False)
    print(string_of_stat)
    print("\n")

# main function


def pairwise_alignment(seq1: str, seq2: str, mode: str, scheme={"match": 1, "mismatch": -1, "gap": -2}, align=True):
    """
    Align 2 DNA sequences either globlly or locally 

    Arguments:
    seq1: str
    seq2: str
    mode: str either "global" or "local"
    scheme ={"match": 1, "mismatch": -1, "gap": -2}
    align = True 

    Returns:
    {
        "seq1": str
        "seq2": str
        "align_seq1": str
        "align_seq2": str
        "score": int
        "stats": dict
        "matrix": numpy.ndarray
    }

    stats = {
        "match": int,
        "mismatch": int,
        "seq1_gap": int,
        "seq2_gap": int,
        "identity": float
    }

    """
    matrix = get_scoring_matrix(seq1, seq2, mode, scheme)
    alignment = traceback(matrix, seq1, seq2, mode)
    if align == True:
        print_alignment(alignment, mode)
    get_stats(alignment)
    return alignment


# only for seq below 20 bases
def heat_map(s_matrix, seq1: str, seq2: str):
    """
    Generates a heatmap of scoring matrices for sequences below length 20. 

    Arugment:
    s_matrix: numpy.array
    seq1: str
    seq2: str 

    Return: None 

    """

    if len(seq1) > 20 or len(seq2) > 20:
        print("Alignment too large to be displayed")
        return None
    else:
        pass

    plt.figure(figsize=(7, 5))

    rows = ['-'] + [n for n in seq2]
    cols = ['-'] + [n for n in seq1]

    matrix = pd.DataFrame(s_matrix, index=rows, columns=cols)
    ax = sns.heatmap(matrix, cmap="viridis", annot=True if len(seq1) < 25 else False,
                     fmt=".0f", linewidths=0.5, square=True,
                     cbar_kws={"label": "alignment score"})
    ax.set_xlabel("Sequence 1")
    ax.set_ylabel("Sequence 2")
    ax.set_title("Global‑alignment scoring matrix", pad=14)
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)
    plt.tight_layout()
    plt.show()
