
#define WORD unsigned long long int
#define WORDSIZE sizeof(WORD)
#define BITSINWORD (WORDSIZE*8)
#define PopCount(x) __builtin_popcountll((x))

void printUsage();
void report(char* t, unsigned int i, unsigned int m, unsigned int s);
void search(char* ALPHA, unsigned int sigma, WORD* P, WORD* T, char* t, unsigned int n, unsigned int m, unsigned int k);
