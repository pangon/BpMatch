/*
 * Copyright (C) 2007-2011 Claudio Felicioli
 * mail: c.felicioli@1d20.net - pangon@gmail.com
 *
 * dna2st is free software; you can redistribute it and/or modify
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
#include <sys/types.h>
#include <sys/stat.h>

#include "lst_debug.h"
#include "lst_structs.h"
#include "lst_stree.h"
#include "lst_string.h"
#include "lst_algorithms.h"


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
void translate2stree(LST_Node* node, st_node* node2, st_node** leaves);
int st_compute_count(st_node* node);
void st_free(st_node* root, st_node** leaves);
void st_free_node(st_node* node);
void st_serialize(st_node* root);
void st_serialize_node(st_node* node, FILE* file);
void st_serialize_edge(st_edge* edge, FILE* file);
enum bpmatch_utils_base {A, G, T, C, X};
char base2char[5];
int scanBp(FILE* file, bpmatch_utils_base* bp);
int base_cmp(void* a, void* b);
void base_copy(void* a, void* b);
char* base_print(LST_StringIndex *index);
bpmatch_utils_base* s;
FILE* gene;
int geneSize;
int trueGeneSize;
int leafCount;
char geneName[50];
char stName[50];
bool scanBp_fasta;

int main(int argc, char *argv[]) {
	if(argc!=3) {
		printf("usage: %s geneSourceFile suffixTreeTargetFile\n", argv[0]);
		exit(EXIT_FAILURE);
		}
	strncpy(geneName, argv[1], 49);
	geneName[49]='\0';
	strncpy(stName, argv[2], 49);
	stName[49]='\0';
	if((gene=fopen(geneName, "r"))==NULL) {
		fprintf(stderr, "error: %s opening fail\n", geneName);
		exit(EXIT_FAILURE);
		}
	struct stat geneInfo;
	stat(geneName, &geneInfo);
	geneSize=geneInfo.st_size;
	base2char[A]='A';
	base2char[G]='G';
	base2char[T]='T';
	base2char[C]='C';
	base2char[X]='X';
	LST_StringSet* set=lst_stringset_new();
	s=new bpmatch_utils_base[geneSize+1];
	trueGeneSize=0;
	bpmatch_utils_base bp;
	scanBp_fasta=false;
	while(scanBp(gene, &bp)!=EOF) {
		s[trueGeneSize++]=bp;
//		printf("%c", base2char[s[trueGeneSize-1]]);
		}
//	printf("\n");
	int gene_fd;
	if((gene_fd=fileno(gene))==-1) {
		fprintf(stderr, "failed getting gene file descriptor");
		exit(EXIT_FAILURE);
		}
	if(fsync(gene_fd)!=0) {
		fprintf(stderr, "failed fsync of gene");
		exit(EXIT_FAILURE);
		}
	if(fclose(gene)!=0) {
		fprintf(stderr, "failed fclose of gene");
		exit(EXIT_FAILURE);
		}
	s[trueGeneSize++]=X;
	lst_stringclass_set_defaults(base_cmp, base_copy, base_print);
	LST_String* string=lst_string_new(s, sizeof(bpmatch_utils_base), trueGeneSize-1);
	lst_stringset_add(set, string);
	LST_STree* tree=lst_stree_new(set);
//	lst_debug_print_tree(tree);
	st_node* root=new st_node;
	st_node** leaves=new st_node*[trueGeneSize+2];
	root->length=0;
	leafCount=1;
	translate2stree(tree->root_node, root, leaves);
	lst_stree_free(tree);
	st_compute_count(root);
//	st_print(root);
	st_serialize(root);
	delete[] s;
	st_free(root, leaves);
	exit(EXIT_SUCCESS);
	}

int base_cmp(void* a, void* b) {
	if(*((bpmatch_utils_base*)a)==*((bpmatch_utils_base*)b)) return(0);
	if(*((bpmatch_utils_base*)a)<*((bpmatch_utils_base*)b)) return(-1);
	return(1);
	}

void base_copy(void* a, void* b) {
	*((bpmatch_utils_base*)b)=*((bpmatch_utils_base*)a);
	}

char* base_print(LST_StringIndex *index) {
	if(index->start_index==index->string->num_items-1) return("<eos>");
	char* out=new char[100];
	for(int i=index->start_index;i<=*(index->end_index);i++) out[i-index->start_index]=base2char[s[i]];
	out[*(index->end_index)-index->start_index+1]='\0';
	return(out);
	}

void st_serialize_node(st_node* node, FILE* file) {
	fwrite(&(node->length), sizeof(int), 1, file);
	fwrite(&(node->leafId), sizeof(int), 1, file);
	fwrite(&(node->firstLeaf), sizeof(int), 1, file);
	fwrite(&(node->lastLeaf), sizeof(int), 1, file);
	bool presence;
	for(int i=0;i<5;i++) {
		if(node->edges[i]!=NULL) {
			presence=true;
			fwrite(&presence, sizeof(bool), 1, file);
			st_serialize_edge(node->edges[i], file);
			}
		else {
			presence=false;
			fwrite(&presence, sizeof(bool), 1, file);
			}
		}
	}

void st_serialize_edge(st_edge* edge, FILE* file) {
	fwrite(&(edge->count), sizeof(int), 1, file);
	fwrite(&(edge->start), sizeof(int), 1, file);
	fwrite(&(edge->end), sizeof(int), 1, file);
	st_serialize_node(edge->to, file);
	}

void st_serialize(st_node* root) {
	printf("serializing suffix tree...\n");
	FILE* st_file=fopen(stName, "w");
	fwrite(&(trueGeneSize), sizeof(int), 1, st_file);
	fwrite(s, sizeof(bpmatch_utils_base), trueGeneSize, st_file);
	st_serialize_node(root, st_file);
	int st_file_fd;
	if((st_file_fd=fileno(st_file))==-1) {
		fprintf(stderr, "failed getting st_file file descriptor");
		exit(EXIT_FAILURE);
		}
	if(fsync(st_file_fd)!=0) {
		fprintf(stderr, "failed fsync of st_file");
		exit(EXIT_FAILURE);
		}
	if(fclose(st_file)!=0) {
		fprintf(stderr, "failed fclose of st_file");
		exit(EXIT_FAILURE);
		}
	printf("done\n");
	}

void st_free(st_node* root, st_node** leaves) {
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

void translate2stree(LST_Node* node, st_node* node2, st_node** leaves) {
	for(int i=0;i<5;i++) node2->edges[i]=NULL;
	st_node* tmp;
	for(LST_Edge* edge=node->kids.lh_first;edge;edge=edge->siblings.le_next) {
		switch(s[edge->range.start_index]) {
			case A:
				node2->edges[A]=new st_edge;
				node2->edges[A]->start=edge->range.start_index;
				node2->edges[A]->end=*(edge->range.end_index);
				tmp=new st_node;
				node2->edges[A]->to=tmp;
				tmp->length=node2->length+(node2->edges[A]->end-node2->edges[A]->start)+1;
				translate2stree(edge->dst_node, tmp, leaves);
				break;
			case C:
				node2->edges[C]=new st_edge;
				node2->edges[C]->start=edge->range.start_index;
				node2->edges[C]->end=*(edge->range.end_index);
				tmp=new st_node;
				node2->edges[C]->to=tmp;
				tmp->length=node2->length+(node2->edges[C]->end-node2->edges[C]->start)+1;
				translate2stree(edge->dst_node, tmp, leaves);
				break;
			case T:
				node2->edges[T]=new st_edge;
				node2->edges[T]->start=edge->range.start_index;
				node2->edges[T]->end=*(edge->range.end_index);
				tmp=new st_node;
				node2->edges[T]->to=tmp;
				tmp->length=node2->length+(node2->edges[T]->end-node2->edges[T]->start)+1;
				translate2stree(edge->dst_node, tmp, leaves);
				break;
			case G:
				node2->edges[G]=new st_edge;
				node2->edges[G]->start=edge->range.start_index;
				node2->edges[G]->end=*(edge->range.end_index);
				tmp=new st_node;
				node2->edges[G]->to=tmp;
				tmp->length=node2->length+(node2->edges[G]->end-node2->edges[G]->start)+1;
				translate2stree(edge->dst_node, tmp, leaves);
				break;
			case X:
				node2->edges[X]=new st_edge;
				node2->edges[X]->start=edge->range.start_index;
				node2->edges[X]->end=*(edge->range.end_index);
				tmp=new st_node;
				node2->edges[X]->to=tmp;
				tmp->length=node2->length+(node2->edges[X]->end-node2->edges[X]->start)+1;
				translate2stree(edge->dst_node, tmp, leaves);
				break;
			}
		}
	if(!(node->kids.lh_first)) {
		node2->leafId=leafCount++;
		leaves[node2->leafId]=node2;
		node2->firstLeaf=node2->leafId;
		node2->lastLeaf=node2->leafId;
		}
	else {
		node2->leafId=0;
		node2->firstLeaf=0;
		node2->lastLeaf=0;
		for(int i=0;i<5;i++) {
			if(node2->edges[i]!=NULL) {
				if(node2->firstLeaf==0 || node2->firstLeaf>node2->edges[i]->to->firstLeaf) node2->firstLeaf=node2->edges[i]->to->firstLeaf;
				if(node2->lastLeaf==0 || node2->lastLeaf<node2->edges[i]->to->lastLeaf) node2->lastLeaf=node2->edges[i]->to->lastLeaf;
				}
			}
		}
	}

int st_compute_count(st_node* node) {
	int sum=0;
	for(int i=0;i<5;i++) {
		if(node->edges[i]!=NULL) {
			node->edges[i]->count=st_compute_count(node->edges[i]->to);
			sum+=node->edges[i]->count;
			}
		}
	if(sum==0) return(1);
	return(sum);
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

int scanBp(FILE* file, bpmatch_utils_base* bp) {
	char c;
	if(fscanf(file, "%c", &c)==EOF) return EOF;

	if(scanBp_fasta || c=='>') {
			scanBp_fasta=(c!='\n');
			return(scanBp(file, bp));
			}

	if(c=='\r' || c=='\n' || c==' ' || c=='\t' || c=='0' || c=='1' || c=='2' || c=='3' || c=='4' || c=='5' || c=='6' || c=='7' || c=='8' || c=='9') return(scanBp(file, bp));
	switch(c) {
		case 'A':
		case 'a':
			*bp=A;
			break;
		case 'G':
		case 'g':
			*bp=G;
			break;
		case 'T':
		case 't':
			*bp=T;
			break;
		case 'C':
		case 'c':
			*bp=C;
			break;
			break;
		case 'N':
		case 'n':
		case 'X':
		case 'x':
			*bp=X;
			break;
		default:
			fprintf(stderr, "failed traducing char %c to base.\n", c);
			exit(EXIT_FAILURE);
		}
	return(1);
	}
