/*
 * Copyright (C) 2007-2011 Claudio Felicioli
 * mail: c.felicioli@1d20.net - pangon@gmail.com
 *
 * countrepetitions is free software; you can redistribute it and/or modify
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

#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
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

void st_print(st_node* node);
void st_free(void);
void st_free_node(st_node* node);
void st_unserialize(void);
void st_unserialize_node(st_node* node);
void st_unserialize_edge(st_edge* edge);
void recPrint(st_node* node);
void recCount(st_node* node);
enum bpmatch_utils_base {A, G, T, C, X};
char base2char[5];
FILE* st_file;
st_node* root;
st_node** leaves;
char st_fileName[50];
double factor;
int stringLength;
int seqSum;
int seqNumber;
bpmatch_utils_base* s;

// textual output
char* s2;
int pos;
double expected;
FILE* output;
char outputFileName[50];

int main(int argc, char *argv[]) {
	if(argc!=6) {
		fprintf(stderr, "usage: %s suffixTree factor dirSeqSum dirSeqNumber outputFile\n", argv[0]);
		exit(EXIT_FAILURE);
		}
	strncpy(st_fileName, argv[1], 49);
	st_fileName[49]='\0';
	factor=atof(argv[2]);
	if(factor<0.0) {
		fprintf(stderr, "factor must be >0.0\n");
		exit(EXIT_FAILURE);
		}
	seqSum=atoi(argv[3]);
	seqNumber=atoi(argv[4]);
	strncpy(outputFileName, argv[5], 49);
	outputFileName[49]='\0';
	if((st_file=fopen(st_fileName, "r"))==NULL) {
		fprintf(stderr, "error: %s opening fail\n", st_fileName);
		exit(EXIT_FAILURE);
		}
	printf("Loading %s.\n", st_fileName);
	st_unserialize();
	printf("suffix tree unserialized.\n");
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
//	st_print(root);
//	lunghezza,sequenza,#,posizioni
	base2char[A]='A';
	base2char[G]='G';
	base2char[T]='T';
	base2char[C]='C';
	base2char[X]='X';
	s2=new char[stringLength];
	pos=0;
	output=fopen(outputFileName, "w");
//	recCount(root);
	recPrint(root);
	int outputfd;
	if((outputfd=fileno(output))==-1) {
		fprintf(stderr, "failed getting %s file descriptor", outputFileName);
		exit(EXIT_FAILURE);
		}
	if(fsync(outputfd)!=0) {
		fprintf(stderr, "failed fsync of %s", outputFileName);
		exit(EXIT_FAILURE);
		}
	if(fclose(output)!=0) {
		fprintf(stderr, "failed fclose of %s", outputFileName);
		exit(EXIT_FAILURE);
		}
	printf("done.\n");
	delete[] s2;
	st_free();
	}

void recPrint(st_node* node) {
	int position;
	for(int i=0;i<5;i++) {
		if(node->edges[i]!=NULL && node->edges[i]->count>1) {
			for(int j=0;j<=(node->edges[i]->end-node->edges[i]->start);j++) {
				s2[pos++]=base2char[s[node->edges[i]->start+j]];
				}
			s2[pos]='\0';
			expected=pow(0.25, pos)*2*(seqSum-(pos-1)*seqNumber);
			if(expected<=0.0) expected=0.0000000000000001;
			if((double(node->edges[i]->count)/expected)>=factor) {
				fprintf(output, "%d,%s,%d/%f", pos, s2, node->edges[i]->count, expected);
				for(int u=0;u<node->edges[i]->count;u++) {
					position=1+stringLength-leaves[node->edges[i]->to->firstLeaf+u]->length;
					if(position>stringLength/2) {
						position=stringLength-position-pos+1;
						fprintf(output, ",!%d", position);
						}
					else fprintf(output, ",%d", position);
					}
				fprintf(output, "\n");
				}
			recPrint(node->edges[i]->to);
			pos-=(node->edges[i]->end-node->edges[i]->start+1);
			}
		}
	}

void recCount(st_node* node) {
	for(int i=0;i<5;i++) {
		if(node->edges[i]!=NULL && node->edges[i]->count>=2) {
			for(int j=0;j<=(node->edges[i]->end-node->edges[i]->start);j++) {
				s2[pos++]=base2char[s[node->edges[i]->start+j]];
				s2[pos]='\0';
				fprintf(output, "%d,%s,%d", pos, s2, node->edges[i]->count);
				for(int u=0;u<node->edges[i]->count;u++) {
					fprintf(output, ",%d", 1+stringLength-leaves[node->edges[i]->to->firstLeaf+u]->length);
					}
				fprintf(output, "\n");
				}
			recCount(node->edges[i]->to);
			pos-=(node->edges[i]->end-node->edges[i]->start+1);
			}
		}
	}

void st_unserialize(void) {
	fread(&(stringLength), sizeof(int), 1, st_file);
	s=new bpmatch_utils_base[stringLength];
	fread(s, sizeof(bpmatch_utils_base), stringLength, st_file);
	root=new st_node;
	leaves=new st_node*[stringLength+2];
	st_unserialize_node(root);
	}

void st_unserialize_node(st_node* node) {
	fread(&(node->length), sizeof(int), 1, st_file);
	fread(&(node->leafId), sizeof(int), 1, st_file);
	if(node->leafId!=0) leaves[node->leafId]=node;
	fread(&(node->firstLeaf), sizeof(int), 1, st_file);
	fread(&(node->lastLeaf), sizeof(int), 1, st_file);
	bool presence;
	for(int i=0;i<5;i++) {
		fread(&presence, sizeof(bool), 1, st_file);
		if(presence) {
			node->edges[i]=new st_edge;
			st_unserialize_edge(node->edges[i]);
			}
		else node->edges[i]=NULL;
		}
	}

void st_unserialize_edge(st_edge* edge) {
	fread(&(edge->count), sizeof(int), 1, st_file);
	fread(&(edge->start), sizeof(int), 1, st_file);
	fread(&(edge->end), sizeof(int), 1, st_file);
	edge->to=new st_node;
	st_unserialize_node(edge->to);
	}

void st_print(st_node* node) {
	printf("NODO l=%d, leaf=%d [%d-%d]\n", node->length, node->leafId, node->firstLeaf, node->lastLeaf);
	if(node->edges[A]==NULL && node->edges[C]==NULL && node->edges[T]==NULL && node->edges[G]==NULL && node->edges[X]==NULL) printf("-leaf\n");
	for(int i=0;i<5;i++) {
		if(node->edges[i]!=NULL) printf("-edge[%d-%d](count=%d)\n", node->edges[i]->start, node->edges[i]->end, node->edges[i]->count);
		}
	for(int i=0;i<5;i++) {
		if(node->edges[i]!=NULL) st_print(node->edges[i]->to);
		}
	}

void st_free(void) {
	delete[] leaves;
	st_free_node(root);
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
