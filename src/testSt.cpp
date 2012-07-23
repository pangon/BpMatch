/*
 * Copyright (C) 2007-2011 Claudio Felicioli
 * mail: c.felicioli@1d20.net - pangon@gmail.com
 *
 * testSt is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

struct st_node;
struct st_edge {
	int count;
	int start;
	int end;
	st_node* to;
	};

struct st_node {
	int length; //usefull only for leaves
	int leafId; //==0 if not leaf
	st_edge* edges[5];
	int firstLeaf;
	int lastLeaf;
	};

enum bpmatch_utils_base {A, G, T, C, X};
bpmatch_utils_base complement[5];
char base2char[5];

void st_unserialize_node(FILE* file, st_node* node, st_node** leaves);
void st_unserialize_edge(FILE* file, st_edge* edge, st_node** leaves);
void st_free_node(st_node* node);
void recTest(st_node* node, int depth);
void st_print(st_node* node);

//for recTest
bpmatch_utils_base recString[20];
bool* startPoints[20];
int testCount;
int testHidden;
int testHiddenCount;
st_edge* lastEdge;

char st_fileName[50];
FILE* st_file;
int stringLength;
bpmatch_utils_base* s;
st_node* st_root;
st_node** st_leaves;

int main(int argc, char *argv[]) {
	if(argc!=2) {
		printf("usage: %s suffixTree\n", argv[0]);
		exit(EXIT_FAILURE);
		}
	strncpy(st_fileName, argv[1], 48);
	st_fileName[49]='\0';
	complement[C]=G;
	complement[G]=C;
	complement[A]=T;
	complement[T]=A;
	complement[X]=X;
	base2char[A]='A';
	base2char[G]='G';
	base2char[T]='T';
	base2char[C]='C';
	base2char[X]='X';
	printf("Loading %s.\n", st_fileName);
	if((st_file=fopen(st_fileName, "r"))==NULL) {
		fprintf(stderr, "error: %s opening fail\n", st_fileName);
		exit(EXIT_FAILURE);
		}
	fread(&(stringLength), sizeof(int), 1, st_file);
	s=new bpmatch_utils_base[stringLength];
	fread(s, sizeof(bpmatch_utils_base), stringLength, st_file);
	st_root=new st_node;
	st_leaves=new st_node*[stringLength+2];
	st_unserialize_node(st_file, st_root, st_leaves);
	int st_file_fd;
	if((st_file_fd=fileno(st_file))==-1) {
		fprintf(stderr, "failed getting %s file descriptor", st_fileName);
		exit(EXIT_FAILURE);
		}
	if(fsync(st_file_fd)!=0) {
		fprintf(stderr, "failed fsync of %s", st_fileName);
		exit(EXIT_FAILURE);
		}
	if(fclose(st_file)!=0) {
		fprintf(stderr, "failed fclose of %s", st_fileName);
		exit(EXIT_FAILURE);
		}
	for(int i=0;i<stringLength;i++) printf("%c", base2char[s[i]]);
	printf("\n");
//	st_print(st_root);
	printf("suffix tree unserialized.\n");
	testHidden=0;
	testHiddenCount=0;
	recTest(st_root, 0);
	delete[] st_leaves;
	st_free_node(st_root);
	delete[] s;
	printf("termination reached without errors\n");
	exit(EXIT_SUCCESS);
	}

void st_unserialize_node(FILE* file, st_node* node, st_node** leaves) {
	fread(&(node->length), sizeof(int), 1, file);
	fread(&(node->leafId), sizeof(int), 1, file);
	if(node->leafId!=0) leaves[node->leafId]=node;
	fread(&(node->firstLeaf), sizeof(int), 1, file);
	fread(&(node->lastLeaf), sizeof(int), 1, file);
	bool presence;
	for(int i=0;i<5;i++) {
		fread(&presence, sizeof(bool), 1, file);
		if(presence) {
			node->edges[i]=new st_edge;
			st_unserialize_edge(file, node->edges[i], leaves);
			}
		else node->edges[i]=NULL;
		}
	}

void st_unserialize_edge(FILE* file, st_edge* edge, st_node** leaves) {
	fread(&(edge->count), sizeof(int), 1, file);
	fread(&(edge->start), sizeof(int), 1, file);
	fread(&(edge->end), sizeof(int), 1, file);
	edge->to=new st_node;
	st_unserialize_node(file, edge->to, leaves);
	}

void st_free_node(st_node* node) {
	for(int i=0;i<5;i++) {
		if(node->edges[i]!=NULL) {
			st_free_node(node->edges[i]->to);
			delete node->edges[i];
			}
		}
	delete node;
	}

void recTest(st_node* node, int depth) {
	if(depth>=20) {
		fprintf(stderr, "too depth: %d\n", depth);
		exit(EXIT_FAILURE);
		}
	if(depth==0) {
		for(int i=0;i<20;i++) startPoints[i]=new bool[stringLength];
		}
	if(testHiddenCount<testHidden) {
		testHiddenCount++;
		recString[depth]=bpmatch_utils_base(s[lastEdge->start+testHiddenCount]);
		printf("tocontrol ");
		for(int j=0;j<=depth;j++) printf("%c", base2char[recString[j]]);
		printf(" (%d) ...", lastEdge->count);
		testCount=0;
		for(int j=0;j<stringLength;j++) {
			if(startPoints[depth-1][j]) {
				if(j+depth<stringLength) startPoints[depth][j]=(s[j+depth]==s[lastEdge->start+testHiddenCount]);
				else startPoints[depth][j]=false;
				}
			else startPoints[depth][j]=false;
			if(startPoints[depth][j]) testCount++;
			}
		if(testCount==lastEdge->count) printf("OK!\n");
		else {
			printf("WRONG! (%d) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", testCount);
			}
		recTest(node, depth+1);
		}
	else {
		for(int i=0;i<5;i++) {
			if(node->edges[i]!=NULL) {
				recString[depth]=bpmatch_utils_base(i);
				printf("tocontrol ");
				for(int j=0;j<=depth;j++) printf("%c", base2char[recString[j]]);
				printf(" (%d) ...", node->edges[i]->count);
				//metto in startPoints[depth] quelli che rimangono ammissibili da startPoints[depth-1]
				testCount=0;
				for(int j=0;j<stringLength;j++) {
					if(depth==0) startPoints[depth][j]=(s[j]==i);
					else {
						if(startPoints[depth-1][j]) {
							if(j+depth<stringLength) startPoints[depth][j]=(s[j+depth]==i);
							else startPoints[depth][j]=false;
							}
						else startPoints[depth][j]=false;
						}
					if(startPoints[depth][j]) testCount++;
					}
				if(testCount==node->edges[i]->count) printf("OK!\n");
				else {
					printf("WRONG! (%d) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", testCount);
					}
				lastEdge=node->edges[i];
				testHidden=node->edges[i]->end-node->edges[i]->start;
				testHiddenCount=0;
				if(node->edges[i]->to->leafId>0) testHidden=0;
				recTest(node->edges[i]->to, depth+1);
				}
			}
		}
	if(depth==0) {
		for(int i=0;i<20;i++) delete[] startPoints[i];
		}
	}

void st_print(st_node* node) {
	printf("NODO l=%d, leaf=%d [%d-%d]\n", node->length, node->leafId, node->firstLeaf, node->lastLeaf);
	if(node->edges[A]==NULL && node->edges[C]==NULL && node->edges[T]==NULL && node->edges[G]==NULL && node->edges[X]==NULL) printf("-leaf\n");
	for(int i=0;i<5;i++) {
		if(node->edges[i]!=NULL) printf("-*%c*[%d-%d](count=%d)\n", base2char[i], node->edges[i]->start, node->edges[i]->end, node->edges[i]->count);
		}
	for(int i=0;i<5;i++) {
		if(node->edges[i]!=NULL) st_print(node->edges[i]->to);
		}
	}
