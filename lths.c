// License MIT 2024 Ahmad Retha
// LTHS: Linear Time Hamming Search
//
// In this implementation you do not need to provide an alphabet argument to the
// program. It accepts an ASCII pattern and text and uses a lookup array to get
// the index of characters in constant time.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lths.h"

void printUsage() {
    printf("LTHS: Linear Time Hamming Search\n\n");
    printf("Usage: ./lths <pattern> <k-mismatches> <text>\n");
    printf("Search for a pattern in a text with k-mismatches under the Hamming distance model\n\n");
    printf("arguments:\n");
    printf("\t<PATTERN>\tThe pattern (needle)\n");
    printf("\t<K_MISMATCHES>\tK-mismatches threshold\n");
    printf("\t<TEXT>\t\tThe text to search (haystack)\n\n");
    printf("License MIT 2024 Ahmad Retha\n");
}

void report(char* t, unsigned int i, unsigned int m, unsigned int s) {
    char x[m+1];
    unsigned int pos = i - m + 1;
    strncpy(x, &t[pos], m);
    x[m] = '\0';
    printf("%d,%d\t%s\n", pos, s, x);
}

void search(char* ALPHA, unsigned int sigma, WORD* P, WORD* T, char* t, unsigned int n, unsigned int m, unsigned int k) {
    //do the rest of the positions i=m-1..n-1
    unsigned int c, h, i, j, s;
    for (i = m-1; i < n; i++) {
        c = (unsigned int)t[i];
        j = (unsigned int)ALPHA[c]; //index of letter in T

        if (j != 0) {
            T[j-1] = T[j-1] | (1 << m-1);
        }

        s = 0;
        for (h = 0; h < sigma; h++) {
            WORD diff = (P[h] ^ T[h]) & P[h];
            s = s + PopCount(diff);
            T[h] = T[h] >> 1;
        }

        if (s <= k) {
            report(t, i, m, s);
        }
    }
}

int main(int argc, char** argv) {
    if (argc == 1) {
        printUsage();
        return 0;
    }

    if (argc != 4) {
        printf("Invalid arguments\n\n");
        printUsage();
        return 1;
    }

    unsigned int m = strlen(argv[1]);
    char* p = argv[1];
    unsigned int k = (unsigned int) atoi(argv[2]);
    unsigned int n = strlen(argv[3]);
    char* t = argv[3];

    if (m > n) {
        printf("Error: Pattern is longer than text\n");
        return 1;
    }

    if (m > BITSINWORD) {
        printf("Error: Pattern is too long\n");
        return 1;
    }

    if (k >= m) {
        printf("Error: Invalid k-mismatches value\n");
        return 1;
    }

    //generate ascii alphabet look up table in time O(m+|ascii|)
    //sigma is used for counting the number of unique letters assigned in p
    //ALPHA will store the index of letters not found in p as 0
    char* ALPHA = (char*) calloc(128, sizeof(char));
    unsigned int c, i, j;
    unsigned int sigma = 0;
    for (i = 0; i < m; i++) {
        c = (unsigned int)p[i]; //ascii char number
        if (ALPHA[c] == 0) {
            ALPHA[c] = ++sigma;
        }
    }

    //create pattern table P, requiring space O(|A|) and time O(m+|A|)
    WORD* P = (WORD*) calloc(sigma, WORDSIZE);
    for (i = 0; i < m; i++) {
        c = (unsigned int)p[i]; //ascii char number
        j = (unsigned int)ALPHA[c]; //index of letter in P
        if (j != 0) {
            P[j-1] = P[j-1] | (1 << i);
        }
    }

    //create text table T, requiring space O(|A|) and time O(m+|A|)
    WORD* T = (WORD*) calloc(sigma, WORDSIZE);
    for (i = 0; i < m-1; i++) {
        c = (unsigned int)t[i];
        j = (unsigned int)ALPHA[c]; //index of letter in T
        if (j != 0) {
            T[j-1] = T[j-1] | (1 << i);
        }
    }

    search(ALPHA, sigma, P, T, t, n, m, k);

    //clean up
    free(ALPHA);
    free(P);
    free(T);

    return 0;
}
