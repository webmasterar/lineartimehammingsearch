# LTHS: Linear Time Hamming Search

LTHS: Linear Time Hamming Search, by Ahmad Retha. This is an algorithm for
*pattern matching* (string searching) under the *Hamming* distance model with
*k-mismatches* in *linear time*, assuming a pattern length of less than the size
of one computer word (typically up to 64 characters) and a constant alphabet size.
I am not aware of an equivalent algorithm, so I am sharing this for the benefit
of researchers in the fields of Computer Science, Bioinformatics and related fields.

Example code is provided in C and Python under the MIT License, 2024.


## Definitions and Concepts

We define an alphabet *A* as a finite non-empty set of letters of size *|A|*.
Here, we will only be describing an algorithm using a fixed alphabet of constant
size, O(1)=|A|. We define a _string_, *x*, as a linear array of letters taken
from *A*. A string has length *m=|x|*. By *x[i]*, we denote the position of a
single letter in a non-empty string by its index *i* in *x*, for *i=0,1,...,m−1*.
We follow the convention that *i=0* is the first index of a letter in the string.
A _substring_ of length *c*, where *1<=c<=m*, is a non-empty linear array of
*x*, starting at any position *i=0..m-c* and ending at *i+c-1*.

Letters in an alphabet *A* can be represented by numbers that correspond to indices.
For example, we can map the DNA alphabet to the following numbers: `{A:0, C:1, G:2, T:3}`.
Thus, letters can be used as indices referencing positions in an array. Let *X*
be an array of size *|A|*, we can reference the second element in the array by
*X[C]*. Looking up the index corresponding to a letter can be done in constant
time *O(1)* using a data structure such as an array, a hash table or a function
call. We abstract the details of that here for the sake of brevity and maintaining
simplicity in the description.

The Hamming distance, *h*, between two strings of equal length, for example *x=ABCAAB*
and *y=ABACAC*, is the number of positions where the characters do not match
between them. In this example,*h=3*, since there are _mismatches_ at positions 2,
3 and 5.

We can use the Hamming distance model for searching a text *y* of length *n>=m*,
to determining if pattern *x* can be found in it. To perform the search, we move
*x* along *y* one position at a time and compare all characters of *x* and *y*
at each new position. This naive _moving-window_ approach has quadratic time
complexity, *O(nm)*, because we compare every character in *x* against the
characters of *y* for each position from *0..n-m*.

To search for a pattern having *k-mismatches* - that is up to *k* differences
between *x* and *y* at any position where the substring of *y* is of length *m* -
we must specify *k*, a distance threshold, where *1<k<m*. If we set *k=0*, this
would become an exact pattern matching problem and there exist more suitable
linear time algorithms (such as KMP) to do that. The value of *k* depends on a
person's requirements for what they consider to be a valid match.

The LTHS algorithm has better time complexity than the naive algorithm and the
distance threshold *k* does not factor into its complexity. The algorithm uses
bitwise operations to perform its computation so computer word size factors in
to the asymptotic time and space complexity. To keep the explanation simple,
we will assume a pattern *p* of length *m<=w*, where *w* is the computer word
size. On modern CPUs, *w=64*.

We also assume operations such as assignment of values to variables and simple
mathematical operations require constant time. Bitwise operations require time
and space *[m/w]* (ceil) to perform their calculations. We assume all operations
on a single computer word require constant time to perform. In the LTHS algorithm,
when *m<=w*, one computer word is used so operations require constant time, *O(1)=[m/w]*.

In addition to pattern length and computer word size, another factor that affects
the time and space complexity is the alphabet (*A*) size. We assume a constant
alphabet *O(1)=|A|*. In the LTHS algorithm, during the preprocessing of the pattern
*p*, we store a computer word for every unique letter in *p*, so it requires space
*O(|A|[m/w])*. In our example, *m<=w*, so one computer word is used, and the
alphabet size is constant, thus space complexity will be *O(1)=O(|A|)*. In the
search stage of the LTHS algorithm, it loops over the alphabet *A* with every
character, so although it is a constant factor that does not play into the
asymptotic time complexity, in practice a large alphabet may affect performance.
LTHS works best with small alphabets such as the DNA alphabet mentioned above.
We will be using the DNA alphabet in the example below.

Algorithm LTHS is divided into two stages: a pattern preprocessing stage and a search
stage on the text. The algorithm can be adapted so the text can be provided _on-line_
one character at a time, as it becomes available.


## Preprocessing

We first create a table *P* for the pattern of size *|A|* having a space complexity
of *O(|A|[m/w])*. All values in *P* will be set to 0 to begin with. Each position
in the table represents a letter of *A* and stores a computer word marking the
positions that a letter occurs in *p* using 1s. Let *p=AATAGC*, then reading
each letter we would set the position to 1 in each word. In the end, *P* contains:


| P[Letter] |  Word  |
|-----------|--------|
| P[A]      | 001011 |
| P[C]      | 100000 |
| P[G]      | 010000 |
| P[T]      | 000100 |


Initializing the table requires time *O(|A|)=O(|A|[m/w])*, and reading the pattern
and marking the positions requires *O(m)=O(m.[m/w])* time. We use the *bitwise-OR*
operation to add 1s to the computer word of each letter. Here is the pseudo-code:

```
P <- a list of words of size |A| initialized to 0
for i = 0..m-1:
	letter x = p[i]
	P[x] = P[x] | (1 << i)
```

Table *P* is now complete and does not need to be modified.

We then create a similarly sized table *T* for the text. This will be initialized
to 0 for all letters. The contents of this table will be updated during the search
phase which is described below.


## Searching

Let *i* denote the current position in the text *t*, and currently *i=0*. 
Before we can do the search proper, we need to read the first few characters of
the text *t*. Let *t=CCAACAGTG* and *k=2*, meaning we will only report a
match with two or fewer mismatches.

In a similar manner to how we preprocessed the pattern, we preprocess the text
for the first *m-1* characters (up to position *i=m-2* inclusive). This is what
table *T* will look like after that:


| T[Letter] |  Word  |
|-----------|--------|
| T[A]      | 001100 |
| T[C]      | 010011 |
| T[G]      | 000000 |
| T[T]      | 000000 |


We are now in the search phase and we can begin to search through the rest of
the positions *i=m-1..n-1*, moving the window by one position right for each
letter in *t*.

The first thing we do is set the last position of letter *t[i]* in *T* to 1, so
now *T[A]=101100*.

```
for i = m-1..n-1:
	letter x = t[i]
	T[x] = T[x] | (1 << m-1)
	# ... code continued below  ...
```

The next step is to loop through each letter of the alphabet *A*, and using bitwise
operations determine how many differences there are between our pattern and the
text in the current window. We use the bitwise _Exclusive-OR_ operation (denoted
by *^*) to find which positions are different between *p* and the *t* for the
current letter in the window. We then use _bitwise-AND_ (denoted by &) to keep
only the 1s at positions where we expect the letter to be in *p*. This keeps the
positions that do not match for the current letter as 1.

We then count the 1s (the mismatching positions) in constant time using the bitwise
_PopCount_ operation. This value will be added to a mismatch counter *s*. We then
bitwise _Right-Shift_ (denoted by *>>*) the word in table *T* for that letter,
shifting all the 1s one position to the right, discarding the bit on the end.
Here is the pseudo-code:

```
	s = 0
	foreach letter a in A:
		diff = (P[a] ^ T[a]) & P[a]
		s = s + PopCount(diff)
		T[a] = T[a] >> 1
	# ... code continued below  ...	
```

Finally, we check the mismatch counter is within the threshold and report the start
position for the match and the distance *s*. Here is the pseudo-code:

```
	if s <= k:
		report(i-m+1, s)
```

Searching the text requires time *O(n[m/w])*, which means it requires *O(n)* time
when *m<=w* and the alphabet is constant.


## Working Example of Search Stage

We will use the same pattern *p=AATAGC* and text *t=CCAACAGTG* as mentioned above
and will assume table *P* has already been initialized in the preprocessing stage
and table *T* has been initialized as before the search stage proper. We will be
using the DNA alphabet, *A={A,C,G,T}*.

We are now ready to continue our search example. Note below that the example is
commented with the code preceded by the # symbol to help describe the process.

Let *i=m-1*, so first we begin by reading the next character of *t* and updating
*T[A]*. The table now looks like this:


| T[Letter] |  Word  |
|-----------|--------|
| T[A]      | 101100 |
| T[C]      | 010011 |
| T[G]      | 000000 |
| T[T]      | 000000 |


Next, we reset the mismatches counter *s=0* and begin looping through alphabet
*A* and assigning each letter to *a*. 

```
letter a = 'A'
# diff = (P[a] ^ T[a]) & P[a]
000011 = (001011 ^ 101100) & 001011
# s = s + PopCount(diff)
2 = 0 + PopCount(000011)
# T[a] = T[a] >> 1
010110 = 101100 >> 1
```

```
letter a = 'C'
# diff = (P[a] ^ T[a]) & P[a]
100000 = (100000 ^ 010011) & 100000
# s = s + PopCount(diff)
3 = 2 + PopCount(100000)
# T[a] = T[a] >> 1
001001 = 010011 >> 1
```

```
letter a = 'G'
# diff = (P[a] ^ T[a]) & P[a]
010000 = (010000 ^ 000000) & 010000
# s = s + PopCount(diff)
4 = 3 + PopCount(010000)
# T[a] = T[a] >> 1
000000 = 000000 >> 1
```

```
letter a = 'T'
# diff = (P[a] ^ T[a]) & P[a]
000100 = (000100 ^ 000000) & 000100
# s = s + PopCount(diff)
5 = 4 + PopCount(000100)
# T[a] = T[a] >> 1
000000 = 000000 >> 1
```

Now we check the mismatches counter and report a match only if *s<=k*. In this
case we do not report a match.

We move to the next position in *t*, *i=m*. We update T[G], so table *T* now
looks like this:


| T[Letter] |  Word  |
|-----------|--------|
| T[A]      | 010110 |
| T[C]      | 001001 |
| T[G]      | 100000 |
| T[T]      | 000000 |


We reset the mismatches counter *s=0* and loop through alphabet *A*.

```
letter a = 'A'
# diff = (P[a] ^ T[a]) & P[a]
001001 = (001011 ^ 010110) & 001011
# s = s + PopCount(diff)
2 = 0 + PopCount(001001)
# T[a] = T[a] >> 1
001011 = 010110 >> 1
```

```
letter a = 'C'
# diff = (P[a] ^ T[a]) & P[a]
100000 = (100000 ^ 001001) & 100000
# s = s + PopCount(diff)
3 = 2 + PopCount(100000)
# T[a] = T[a] >> 1
000100 = 001001 >> 1
```

```
letter a = 'G'
# diff = (P[a] ^ T[a]) & P[a]
010000 = (010000 ^ 100000) & 010000
# s = s + PopCount(diff)
4 = 3 + PopCount(010000)
# T[a] = T[a] >> 1
010000 = 100000 >> 1
```

```
letter a = 'T'
# diff = (P[a] ^ T[a]) & P[a]
000100 = (000100 ^ 000000) & 000100
# s = s + PopCount(diff)
5 = 4 + PopCount(000100)
# T[a] = T[a] >> 1
000000 = 000000 >> 1
```

Now we check the mismatches counter and find *s>k* so we do not report a match.
We move to the next position in *t*, *i=m+1*. We update T[T], so table *T* now
looks like this:


| T[Letter] |  Word  |
|-----------|--------|
| T[A]      | 001011 |
| T[C]      | 000100 |
| T[G]      | 010000 |
| T[T]      | 100000 |


We reset *s=0* and loop through alphabet *A* once again.

```
letter a = 'A'
# diff = (P[a] ^ T[a]) & P[a]
000000 = (001011 ^ 001011) & 001011
# s = s + PopCount(diff)
0 = 0 + PopCount(000000)
# T[a] = T[a] >> 1
000101 = 001011 >> 1
```

```
letter a = 'C'
# diff = (P[a] ^ T[a]) & P[a]
100000 = (100000 ^ 000100) & 100000
# s = s + PopCount(diff)
1 = 0 + PopCount(100000)
# T[a] = T[a] >> 1
000010 = 000100 >> 1
```

```
letter a = 'G'
# diff = (P[a] ^ T[a]) & P[a]
000000 = (010000 ^ 010000) & 010000
# s = s + PopCount(diff)
1 = 1 + PopCount(000000)
# T[a] = T[a] >> 1
001000 = 010000 >> 1
```

```
letter a = 'T'
# diff = (P[a] ^ T[a]) & P[a]
000100 = (000100 ^ 100000) & 000100
# s = s + PopCount(diff)
2 = 1 + PopCount(000100)
# T[a] = T[a] >> 1
010000 = 100000 >> 1
```

Now we check the mismatches counter and find *s<=k* so we report a match with
the starting position in *t*, *i-m+1*, and the number of mismatches, *s=2*.
Comparing the subsequence of *t* at this position, *AACAGT*, and pattern *p=AATAGC*,
we confirm they do in fact have two mismatches at positions 2 and 5.

The next step in searching *t* - the last position - is left as an exercise for
the reader, though we provide table *T* after reading the next character, *i=n-1*,
as a starting point:


| T[Letter] |  Word  |
|-----------|--------|
| T[A]      | 000101 |
| T[C]      | 000010 |
| T[G]      | 101000 |
| T[T]      | 010000 |


## Conclusion

The LTHS algorithm improves on the naive approach at the cost of a little space
for storing tables *P* and *T* and these are dependent on the alphabet size. The
alphabet size is a limiting factor in the algorithm from the performance perspective,
but for a small alphabet and a fairly long pattern, this algorithm is superior to
the naive one. It uses fast bitwise operations, simple addition and logical
comparison operations to update the positions of letters and check for matches.
Assuming the pattern is of length *m<=w* and a constant sized alphabet, preprocessing
the pattern requires *O(m)* time, searching the text is performed in *O(n)*, linear
time. The space complexity is *O(|A|)*.

With some clever programming, using specific operations available on some CPU
architectures, the search stage can be optimized to perform the operations of
each alphabet letter simultaneously and sum up the mismatches in fewer operations,
thus making it even faster than the generalized implementation presented here.
