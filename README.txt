BpMatch is a tool for the fast computation of a segment mutation
distance between biological sequences, based on original algorithms
derived from Claudio Felicioli master thesis.

Copyright (C) 2003-2012 Claudio Felicioli
mail: pangon at gmail dot com

BpMatch is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This particular implementation of BpMatch uses a custom modification
of the generic suffix tree library of Christian Kreibich <christian at whoop dot org>.
All the files containing Christian Kreibich original code, even if
modified by me, maintain his original header and need to be considered
as external to BpMatch.


####### compilation

To compile all the tools you can use the included makefile, executing from the bpmatch directory:
> make all

You can compile single tools with the commands:
> make dna2st
> make testSt
> make complRev
> make bpmatch
> make countrepetitions


####### file format

Genetic sequences need to be described in plain text files with the following rules
-the characters AGCTacgt are accepted as bases
-the characters XNxn are accepted as an unknown base
-any space, tabular or newline are skipped
-any number is skipped
-any line with a starting character > is skipped


####### example of use

The standard use case is to use segments of a source sequence S to cover a target sequence T
The sequences need to be found in two files: S.dna and T.dna
The firs step is to compute the complemented reversal of S:
> ./complRev S.dna S_r.dna

Then you need to compute the two modified suffix trees:
> ./dna2st S.dna S.st
> ./dna2st S_r.dna S_r.st

Finally you can execute BpMatch, in this example the minimum length of segments is 10 and the minimum number of repetitions is 3:
> ./bpmatch S.st S_r.st T.dna 10 3
The simple output is the covered ratio.

For an extended detailed output of the coverage you need to pass one more parameter:
> ./bpmatch S.st S_r.st T.dna 10 3 coverage.txt
A file named coverage.txt is created with the detailed output

