import numpy


SEQ_LENGTH = 70
SEQ_COUNT = 100

alphabet = ["A", "C", "G", "U"]

scores = [
    # A  C  G  U
    [0, 0, 0, 1],  # A
    [0, 0, 1, 0],  # C
    [0, 1, 0, 1],  # G
    [1, 0, 1, 0],  # U
]


#
# scores = [
#     # A  C  G  U
#     [0, 0, 0, 1],  # A
#     [0, 0, 1, 0],  # C
#     [0, 1, 0, 0],  # G
#     [1, 0, 0, 0],  # U
# ]


def generate_random_RNA(SEQ_LENGTH):
    seq = ""
    # A C G U
    probabilities = [0.25, 0.25, 0.25, 0.25]
    for i in range(SEQ_LENGTH):
        seq += numpy.random.choice(alphabet, size=1, p=probabilities)[0]
    return seq


def nussinov(seq):
    score_matrix = [[0.0 for _ in xrange(len(seq))] for _ in xrange(len(seq))]
    # iterate over score_matrix diagonal by diagonal
    for k in reversed(xrange(1, len(seq))):
        i = 0
        j = len(seq) - k
        for _ in xrange(k):
            score_matrix[i][j] = score_matrix[i + 1][j - 1] + scores[alphabet.index(seq[i])][alphabet.index(seq[j])]
            score_matrix[i][j] = max(score_matrix[i][j], score_matrix[i + 1][j], score_matrix[i][j - 1])
            for t in xrange(i, j):
                score_matrix[i][j] = max(score_matrix[i][j], score_matrix[i][t] + score_matrix[t + 1][j])
            i += 1
            j += 1

    best_score = score_matrix[0][len(seq) - 1]
    # print_matrix(score_matrix)
    return best_score


def print_matrix(matrix):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print '\n'.join(table)


def main():
    sequences = []
    sequence_scores = []
    for i in range(SEQ_COUNT):
        new_seq = generate_random_RNA(SEQ_LENGTH)
        # new_seq = "GGGAAAUCC"
        # new_seq = "GGGACCUUCC"
        # new_seq = "GGUUCCUUCCCAA"
        nussinov_score = nussinov(new_seq)
        # print(nussinov_score)
        sequence_scores.append(nussinov_score)
        sequences.append(new_seq)

    print "Average Score:", sum(sequence_scores) / len(sequence_scores)


if __name__ == "__main__":
    main()
