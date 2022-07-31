#!/usr/bin/env python
import sys
import string
import random

#### INSTRUCTIONS FOR USE:
# call program as follows: ./gibbs.py <Motif Length> <Data File>
# make sure the gibbs.py is marked as executable: chmod +x gibbs.py
from math import log

import numpy

alphabet = ['A', 'G', 'C', 'T']

EPSILON = 0.000000000001
MIN_ITER_COUNT = 1500
MAX_ITER_COUNT = 2000


#### GibbsSampler:
#### 	INPUTS:	S - list of sequences
####		L - length of motif
####	OUTPUT:	PWM - 4xL list with frequencies of each base at each position
####                  Order of bases should be consistent with alphabet variable
def GibbsSampler(S, L):
    PWM = []
    for i in range(len(alphabet)):
        PWM.append([0.0] * L)

    ######### ADD YOUR CODE HERE ######
    # Initialize the motif position in each sequence
    A = []
    for i in xrange(len(S)):
        A.append(random.randint(0, len(S[i]) - L))

    # Until convergence:
    converge_count = 0
    iter_count = 0
    while True:
        if iter_count > MAX_ITER_COUNT:
            break
        if converge_count == 3 and iter_count > MIN_ITER_COUNT:
            break
        iter_count += 1
        prev_PWM = []
        prev_A = A[:]
        # copy PWM to prev_PWM
        for i in xrange(len(PWM)):
            prev_PWM.append([])
            for j in xrange(len(PWM[i])):
                prev_PWM[i].append(PWM[i][j])
        # re-estimate the position weight matrix (PWM) from all the motifs except one,
        i = random.randint(0, len(S) - 1) # pick a random sequence
        for k in xrange(L):
            # frequencies of each nucleotide at Kth position of motif
            frequencies = [1.0 for _ in range(len(alphabet))]
            for j in xrange(len(S)):
                if j == i:
                    continue
                motif_start_position = A[j]
                curr_position = motif_start_position + k
                curr_nucleotide = S[j][curr_position]
                frequencies[alphabet.index(curr_nucleotide)] += 1
            normalizer = sum(frequencies)
            for freq in xrange(len(frequencies)):
                frequencies[freq] /= normalizer
                PWM[freq][k] = frequencies[freq]

        # sample a new motif position for the i-th sequence
        Z = [0.0] * (len(S[i]) - L)
        for j in xrange(len(Z)):
            Z[j] = 1.0
            for k in xrange(L):
                curr_pos = j + k
                curr_nucleotide = S[i][curr_pos]
                Z[j] *= PWM[alphabet.index(curr_nucleotide)][k]
        normalizer = sum(Z)
        for j in xrange(len(Z)):
            Z[j] /= normalizer
        # sample A[i] from Z probability distribution
        A[i] = numpy.random.choice(range(len(Z)), size=1, p=Z)[0]

        # compare PWM to prev_PWM to determine convergence
        curr_convergence = True
        for i in xrange(len(PWM)):
            for j in xrange(len(PWM[i])):
                if abs(PWM[i][j] - prev_PWM[i][j]) > EPSILON:
                    curr_convergence = False
                    break
        if curr_convergence:
            converge_count += 1
        else:
            converge_count = 0
    # print "Converged after %d iterations" % iter_count

    ######### END OF YOUR CODE HERE #####
    return PWM


###### YOUR OWN FUNCTIONS HERE
# optional -- feel free to add your own functions if you want to


###### END OF YOUR FUNCTIONS

def main():
    random.seed(42)
    L = int(sys.argv[1])
    datafile = sys.argv[2]
    S = readdata(datafile)
    P = GibbsSampler(S, L)
    # for i in range(L):
    #     print "%-5d " % (i + 1),
    # print ""
    #
    # for j in range(len(alphabet)):
    #     print " %s " % alphabet[j],
    #     for i in range(L):
    #         print " %5.3f" % P[j][i],
    #     print ""
    alphabet_indices = [0, 2, 1, 3]
    for i in range(len(alphabet)):
        print " %s " % alphabet[i],
    print ""

    for i in range(L):
        for j in alphabet_indices:
            print " %5.3f" % P[j][i],
        print ""


def readdata(file):
    data = [];
    for line in open(file, 'r'):
        data.append(line[0:-1])
    return data


main()

