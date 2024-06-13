# License MIT 2024 Ahmad Retha
# Linear Time Hamming Search

import sys
from argparse import ArgumentParser


def report(position, k_mismatches, t, m):
    """ Report (print out) the position, number of mismatches and the pattern found in the text
    position: The start position
    k_mismatches: The number of mismatches
    t: The text
    m: The length of the pattern
    """
    print(f'{position},{k_mismatches}\t{t[position:position+m]}')


def ALPHA(A, c):
    """ Returns the index (position) of letter c in the Alphabet A or -1 if not found
    A: The alphabet
    c: The letter
    """
    return A.find(c)


def initialize_P(A, P, p):
    """ Fills in the Pattern table P in the preprocessing stage
    A: The alphabet
    P: The initialized pattern table
    p: The pattern
    """
    for i, c in enumerate(p):
        j = ALPHA(A, c)
        P[j] = P[j] | (1 << i)


def initialize_T(A, T, t, m):
    """ Fills in the text table T before the search stage proper
    A: The alphabet
    T: The initialized text table
    t: The text
    m: The length of the pattern 
    """
    for i, c in enumerate(t):
        if i >= m-2:
            break
        j = ALPHA(A, c)
        if j == -1:
            continue
        T[j] = T[j] | (1 << i)


def popcount(n):
    """ Returns count of set bits in n
    n: A number
    """
    return n.bit_count()


def search(A, sigma, P, T, t, n, m, k):
    """ Performs the search stage of the LTHS algorithm
    A: The alphabet
    sigma: The size of the alphabet
    P: The pattern table
    T: The text table
    t: The text
    n: The length of the text
    m: The length of the pattern
    k: K-mismatches threshold
    """
    for i in range(m-1, n):
        x = ALPHA(A, t[i])

        if x != -1:
            T[x] = T[x] | (1 << m-1)

        s = 0
        for j in range(sigma):
            diff = (P[j] ^ T[j]) & P[j]
            s = s + popcount(diff)
            T[j] = T[j] >> 1

        if s <= k:
            report(i-m+1, s, t, m)


def LTHS(A, p, k, t):
    """ Performs the Linear Time Hamming Search algorithm
    A: The alphabet, e.g 'ACGT'
    p: The pattern (needle)
    k: K-mismatches threshold (k < |p|)
    t: The text to search (haystack)
    """
    sigma = len(A)
    m = len(p)
    n = len(t)

    P = [0] * sigma
    initialize_P(A, P, p)

    T = [0] * sigma
    initialize_T(A, T, t, m)

    search(A, sigma, P, T, t, n, m, k)


if __name__ == '__main__':
    parser = ArgumentParser(
        prog='LTHS: Linear Time Hamming Search',
        description='Search for a pattern in a text with k-mismatches under the Hamming distance model',
        epilog='License MIT 2024 Ahmad Retha')
    parser.add_argument('-A', '--alphabet', help='Alphabet string, e.g. ACGT', required=True)
    parser.add_argument('-p', '--pattern', help='The pattern (needle)', required=True)
    parser.add_argument('-k', '--k-mismatches', help='K-mismatches threshold', type=int, required=True)
    parser.add_argument('-t', '--text', help='The text to search (haystack)', required=True)
    args = parser.parse_args()

    if len(args.pattern) > len(args.text):
        sys.exit('Pattern is longer than text')

    if args.k_mismatches < 0 or args.k_mismatches >= len(args.pattern):
        sys.exit('Invalid k-mismatches value')

    LTHS(args.alphabet, args.pattern, args.k_mismatches, args.text)
