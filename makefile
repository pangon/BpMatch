VPATH = ./src
CC = g++
#CFLAGS = -Wall -O3
CFLAGS = -w -O3

dna2stObjects = dna2st.cpp lst_debug.c lst_stree.c lst_string.c lst_algorithms.c


all : dna2st testSt complRev bpmatch countrepetitions

clear :
	rm dna2st
	rm testSt
	rm complRev
	rm bpmatch
	rm countrepetitions

dna2st : $(dna2stObjects)
	$(CC) $(CFLAGS) -o $@ $^

testSt : testSt.cpp
	$(CC) $(CFLAGS) -o $@ $^

complRev : complRev.cpp
	$(CC) $(CFLAGS) -o $@ $^

bpmatch : bpmatch.cpp
	$(CC) $(CFLAGS) -o $@ $^

countrepetitions : countrepetitions.cpp
	$(CC) $(CFLAGS) -o $@ $^

