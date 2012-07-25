/*
 * bpmatch calculate, from sequences S and T, the maximum coverage of T
 * using only subsequences and complemented reversed subsequences of S,
 * with minimum length "l", with at least "minRep" occurrences in S,
 * possibly overlapped, and, in such a maximum coverage, minimize the
 * number of subsequences used.
 * version 1.6
 *
 * Copyright (C) 2003-2012 Claudio Felicioli
 * mail: c.felicioli@1d20.net - pangon@gmail.com
 *
 * bpmatch is free software; you can redistribute it and/or modify
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
#include <sys/types.h>
#include <sys/stat.h>
#define DEBUG (false)


enum bpmatch_utils_base {A, G, T, C, X};

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

struct st_search {
	st_node* current_node;
	st_edge* current_edge;
	int current_depth;
	int current_edge_depth;
	int current_count;
	st_node* root;
	bpmatch_utils_base* string;
	};

bpmatch_utils_base complement[5];
char base2char[5];

void intestation();
void st_unserialize_node(FILE* file, st_node* node, st_node** leaves);
void st_unserialize_edge(FILE* file, st_edge* edge, st_node** leaves);
void st_print(st_node* node);
void st_free_node(st_node* node);
st_node* st_match(st_node* node, bpmatch_utils_base base, bpmatch_utils_base* string);
int scanBp(FILE* file, bpmatch_utils_base* bp);
void epilogo();
bool composite_match(st_search* suffix_tree_a, st_search* suffix_tree_b, bpmatch_utils_base base);
bool composite_match_target_check(st_search* suffix_tree_a, st_search* suffix_tree_b, st_search* suffix_tree_c, st_search* suffix_tree_d, bpmatch_utils_base base);
int single_match(st_search* suffix_tree, bpmatch_utils_base base);
void initialize_new_match(st_search* match_struct);

int last_valid_direct_count;
int last_valid_reverse_count;
int last_valid_direct_count_target;
int last_valid_reverse_count_target;

char st_fileName[50];
char rst_fileName[50];
char st2_fileName[50];
char rst2_fileName[50];
char target_fileName[50];
char output_fileName[100];
int l;
int minRep;
int minRep2;
FILE* st_file;
FILE* rst_file;
FILE* st2_file;
FILE* rst2_file;
FILE* target_file;
FILE* output_file;
int stringLength;
bpmatch_utils_base* s;
bpmatch_utils_base* rs;
bpmatch_utils_base* t;
bpmatch_utils_base* rt;
st_node* st_root;
st_node** st_leaves;
st_node* rst_root;
st_node** rst_leaves;
st_node* st2_root;
st_node** st2_leaves;
st_node* rst2_root;
st_node** rst2_leaves;
int at;
bool fileoutput;
bool scanBp_fasta;
bool check_target_repetitions;

int main(int argc, char *argv[]) {
	if(argc!=6 && argc!=7 && argc!=8 && argc!=9) {
		printf("usage: %s sourceGeneSuffixTree reverseSourceGeneSuffixTree TargetGene minimumLength minimumRepetition [outputFile]\n", argv[0]);
		printf("       %s sourceGeneSuffixTree reverseSourceGeneSuffixTree TargetGeneSuffixTree reverseTargetGeneSuffixTree minimumLength minimumSourceRepetition minimumTargetRepetition [outputFile]\n", argv[0]);
		exit(EXIT_FAILURE);
		}

	check_target_repetitions=(argc==8 || argc==9);

	strncpy(st_fileName, argv[1], 49);
	st_fileName[49]='\0';
	strncpy(rst_fileName, argv[2], 49);
	rst_fileName[49]='\0';

	fileoutput=(argc==7 || argc==9);

	if(check_target_repetitions) {
		strncpy(st2_fileName, argv[3], 49);
		st2_fileName[49]='\0';
		strncpy(rst2_fileName, argv[4], 49);
		rst2_fileName[49]='\0';

		l=atoi(argv[5]);
		minRep=atoi(argv[6]);
		minRep2=atoi(argv[7]);

		if(fileoutput) strncpy(output_fileName, argv[8], 99);
		}
	else {
		strncpy(target_fileName, argv[3], 49);
		target_fileName[49]='\0';

		l=atoi(argv[4]);
		minRep=atoi(argv[5]);

		if(fileoutput) strncpy(output_fileName, argv[6], 99);
		}

	if(fileoutput) {
		output_fileName[99]='\0';
		if((output_file=fopen(output_fileName, "w"))==NULL) {
			fprintf(stderr, "error: %s opening fail\n", output_fileName);
			exit(EXIT_FAILURE);
			}
		}

//	intestation();
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
//	printf("Loading %s.\n", st_fileName);
	if((st_file=fopen(st_fileName, "r"))==NULL) {
		fprintf(stderr, "error: %s opening fail\n", st_fileName);
		exit(EXIT_FAILURE);
		}
	if(fread(&(stringLength), sizeof(int), 1, st_file)!=1) {
		fprintf(stderr, "failed fread of stringLength\n");
		exit(EXIT_FAILURE);
		}
	s=new bpmatch_utils_base[stringLength];
	if((int)fread(s, sizeof(bpmatch_utils_base), stringLength, st_file)!=stringLength) {
		fprintf(stderr, "failed fread of source string\n");
		exit(EXIT_FAILURE);
		}
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
//	for(int i=0;i<stringLength;i++) printf("%c", base2char[s[i]]);
//	printf("\n");
//	st_print(st_root);
//	printf("suffix tree unserialized.\n");
//	printf("Loading %s.\n", rst_fileName);
	if((rst_file=fopen(rst_fileName, "r"))==NULL) {
		fprintf(stderr, "error: %s opening fail\n", rst_fileName);
		exit(EXIT_FAILURE);
		}
	if(fread(&(stringLength), sizeof(int), 1, rst_file)!=1) {
		fprintf(stderr, "failed fread of (reversed) stringLength\n");
		exit(EXIT_FAILURE);
		}
	rs=new bpmatch_utils_base[stringLength];
	if((int)fread(rs, sizeof(bpmatch_utils_base), stringLength, rst_file)!=stringLength) {
		fprintf(stderr, "failed fread of (reversed) source string\n");
		exit(EXIT_FAILURE);
		}
	rst_root=new st_node;
	rst_leaves=new st_node*[stringLength+2];
	st_unserialize_node(rst_file, rst_root, rst_leaves);
	int rst_file_fd;
	if((rst_file_fd=fileno(rst_file))==-1) {
		fprintf(stderr, "failed getting %s file descriptor", rst_fileName);
		exit(EXIT_FAILURE);
		}
	if(fsync(rst_file_fd)!=0) {
		fprintf(stderr, "failed fsync of %s", rst_fileName);
		exit(EXIT_FAILURE);
		}
	if(fclose(rst_file)!=0) {
		fprintf(stderr, "failed fclose of %s", rst_fileName);
		exit(EXIT_FAILURE);
		}
//	for(int i=0;i<stringLength;i++) printf("%c", base2char[rs[i]]);
//	printf("\n");
//	st_print(rst_root);
//	printf("suffix tree unserialized.\n");

	int trueStringLength;
	if(check_target_repetitions) {
//		printf("Loading %s.\n", st2_fileName);
		if((st2_file=fopen(st2_fileName, "r"))==NULL) {
			fprintf(stderr, "error: %s opening fail\n", st2_fileName);
			exit(EXIT_FAILURE);
			}
		if(fread(&(stringLength), sizeof(int), 1, st2_file)!=1) {
			fprintf(stderr, "failed fread of stringLength\n");
			exit(EXIT_FAILURE);
			}
		t=new bpmatch_utils_base[stringLength];
		if((int)fread(t, sizeof(bpmatch_utils_base), stringLength, st2_file)!=stringLength) {
			fprintf(stderr, "failed fread of source string\n");
			exit(EXIT_FAILURE);
			}
		trueStringLength=stringLength;
		st2_root=new st_node;
		st2_leaves=new st_node*[stringLength+2];
		st_unserialize_node(st2_file, st2_root, st2_leaves);
		int st2_file_fd;
		if((st2_file_fd=fileno(st2_file))==-1) {
			fprintf(stderr, "failed getting %s file descriptor", st2_fileName);
			exit(EXIT_FAILURE);
			}
		if(fsync(st2_file_fd)!=0) {
			fprintf(stderr, "failed fsync of %s", st2_fileName);
			exit(EXIT_FAILURE);
			}
		if(fclose(st2_file)!=0) {
			fprintf(stderr, "failed fclose of %s", st2_fileName);
			exit(EXIT_FAILURE);
			}
//		printf("Loading %s.\n", rst2_fileName);
		if((rst2_file=fopen(rst2_fileName, "r"))==NULL) {
			fprintf(stderr, "error: %s opening fail\n", rst2_fileName);
			exit(EXIT_FAILURE);
			}
		if(fread(&(stringLength), sizeof(int), 1, rst2_file)!=1) {
			fprintf(stderr, "failed fread of (reversed) stringLength\n");
			exit(EXIT_FAILURE);
			}
		rt=new bpmatch_utils_base[stringLength];
		if((int)fread(rt, sizeof(bpmatch_utils_base), stringLength, rst2_file)!=stringLength) {
			fprintf(stderr, "failed fread of (reversed) source string\n");
			exit(EXIT_FAILURE);
			}
		rst2_root=new st_node;
		rst2_leaves=new st_node*[stringLength+2];
		st_unserialize_node(rst2_file, rst2_root, rst2_leaves);
		int rst2_file_fd;
		if((rst2_file_fd=fileno(rst2_file))==-1) {
			fprintf(stderr, "failed getting %s file descriptor", rst2_fileName);
			exit(EXIT_FAILURE);
			}
		if(fsync(rst2_file_fd)!=0) {
			fprintf(stderr, "failed fsync of %s", rst2_fileName);
			exit(EXIT_FAILURE);
			}
		if(fclose(rst2_file)!=0) {
			fprintf(stderr, "failed fclose of %s", rst2_fileName);
			exit(EXIT_FAILURE);
			}
		}
	else {
		if((target_file=fopen(target_fileName, "r"))==NULL) {
			fprintf(stderr, "error: %s opening fail\n", target_fileName);
			exit(EXIT_FAILURE);
			}
		struct stat target_fileInfo;
		stat(target_fileName, &target_fileInfo);
		stringLength=target_fileInfo.st_size;
		t=new bpmatch_utils_base[stringLength+1];
		trueStringLength=0;
		bpmatch_utils_base bp;
		scanBp_fasta=false;
		while(scanBp(target_file, &bp)!=EOF) t[trueStringLength++]=bp;
		t[trueStringLength++]=X;
		int target_file_fd;
		if((target_file_fd=fileno(target_file))==-1) {
			fprintf(stderr, "failed getting gene file descriptor");
			exit(EXIT_FAILURE);
			}
		if(fsync(target_file_fd)!=0) {
			fprintf(stderr, "failed fsync of gene");
			exit(EXIT_FAILURE);
			}
		if(fclose(target_file)!=0) {
			fprintf(stderr, "failed fclose of gene");
			exit(EXIT_FAILURE);
			}
		}

//	printf("\n");
//	printf("*************************************************\n");
//	printf("matching - l=%d - minRep=%d\n", l, minRep);
//	printf("*************************************************\n");
	bool end=false;
	bool found;
	int shiftTo=0;
	int maxPrefix, maxSuffix;
	int state=0;
	at=0;
	int i;
	int coverage=0;
	st_search* search_source_direct=new st_search;
	st_search* search_source_reverse=new st_search;
	st_search* search_target_direct=new st_search;
	st_search* search_target_reverse=new st_search;
	search_source_direct->root=st_root;
	search_source_direct->string=s;
	search_source_reverse->root=rst_root;
	search_source_reverse->string=rs;
	if(check_target_repetitions) {
		search_target_direct->root=st2_root;
		search_target_direct->string=t;
		search_target_reverse->root=rst2_root;
		search_target_reverse->string=rt;
		}
	while(!end) {
		found=false;
		if(state==0) {
			//case 0
			if(DEBUG) printf("case 0\n");
			//searching for direct seq
			i=l-1;
			//using reverse from actual+l-1 to actual (backward parsing)
			initialize_new_match(search_source_direct);
			initialize_new_match(search_source_reverse);
			if(check_target_repetitions) {
				initialize_new_match(search_target_direct);
				initialize_new_match(search_target_reverse);
				while(i>=0 && at+i<trueStringLength && composite_match_target_check(search_source_direct, search_source_reverse, search_target_direct, search_target_reverse, complement[t[at+i]])) i--;
				}
			else {
				while(i>=0 && at+i<trueStringLength && composite_match(search_source_direct, search_source_reverse, complement[t[at+i]])) i--;
				}
			if(i>=0) {
				//not found
				shiftTo=at+i+1;
				if(DEBUG) printf("not found [%d -> %d]\n", at, shiftTo);
				}
			else {
				//direct sequence found
				i=0;
				//using direct from actual (forward parsing)
				initialize_new_match(search_source_direct);
				initialize_new_match(search_source_reverse);
				if(check_target_repetitions) {
					initialize_new_match(search_target_direct);
					initialize_new_match(search_target_reverse);
					while(at+i<trueStringLength && composite_match_target_check(search_source_direct, search_source_reverse, search_target_direct, search_target_reverse, t[at+i])) i++;
					}
				else {
					while(at+i<trueStringLength && composite_match(search_source_direct, search_source_reverse, t[at+i])) i++;
					}
				found=true;
				shiftTo=at+i;
				if(DEBUG) printf("***************************************************found [%d -> %d]\n", at, shiftTo);
				}
			}
		else {
			//case 1
			if(DEBUG) printf("case 1\n");
			//searching for direct (eventually ovelapping) seq
			maxSuffix=0;
			maxPrefix=0;
			//using direct from actual (forward parsing)
			initialize_new_match(search_source_direct);
			initialize_new_match(search_source_reverse);
			if(check_target_repetitions) {
				initialize_new_match(search_target_direct);
				initialize_new_match(search_target_reverse);
				while(at+maxSuffix<trueStringLength && composite_match_target_check(search_source_direct, search_source_reverse, search_target_direct, search_target_reverse, t[at+maxSuffix])) maxSuffix++;
				}
			else {
				while(at+maxSuffix<trueStringLength && composite_match(search_source_direct, search_source_reverse, t[at+maxSuffix])) maxSuffix++;
				}
			if(maxSuffix>=l) {
				//direct sequence found
				found=true;
				shiftTo=at+maxSuffix;
				if(DEBUG) printf("1**************************************************found [%d -> %d]\n", at, shiftTo);
				}
			else {
				if(DEBUG) printf("compute maxPrefix\n");
				//using reverse from actual (backward parsing)
				initialize_new_match(search_source_direct);
				initialize_new_match(search_source_reverse);
				if(check_target_repetitions) {
					initialize_new_match(search_target_direct);
					initialize_new_match(search_target_reverse);
					while(maxPrefix<l && at-maxPrefix>=0 && composite_match_target_check(search_source_direct, search_source_reverse, search_target_direct, search_target_reverse, complement[t[at-maxPrefix]])) maxPrefix++;
					}
				else {
					while(maxPrefix<l && at-maxPrefix>=0 && composite_match(search_source_direct, search_source_reverse, complement[t[at-maxPrefix]])) maxPrefix++;
					}
				maxPrefix--;
				while(!found && maxPrefix+maxSuffix>=l) {
					if(DEBUG) printf("*ite %d-%d\n", maxPrefix, maxSuffix);
					i=maxSuffix-l;
					//using direct from actual (forward parsing)
					initialize_new_match(search_source_direct);
					initialize_new_match(search_source_reverse);
					if(check_target_repetitions) {
						initialize_new_match(search_target_direct);
						initialize_new_match(search_target_reverse);
						while(at+i<trueStringLength && composite_match_target_check(search_source_direct, search_source_reverse, search_target_direct, search_target_reverse, t[at+i])) i++;
						}
					else {
						while(at+i<trueStringLength && composite_match(search_source_direct, search_source_reverse, t[at+i])) i++;
						}
					if(i==maxSuffix) {
						//direct sequence found
						found=true;
						shiftTo=at+maxSuffix;
						}
					else {
						if(DEBUG) printf("error at %d\n", i);
						maxSuffix=i;
						}
					}
				if(!found) {
					//not found
					shiftTo=at+1;
					}
				}
			}

		if(found) {
			if(fileoutput) {
				fprintf(output_file, "%i,", at);
				if(state==1) fprintf(output_file, "*");
				for(int j=at;j<shiftTo;j++) fprintf(output_file, "%c", base2char[t[j]]);
				if(check_target_repetitions) fprintf(output_file, ",direct-reverse:%d-%d,targetrepetitionsdirect-reverse:%d-%d\n", last_valid_direct_count, last_valid_reverse_count, last_valid_direct_count_target, last_valid_reverse_count_target);
				else fprintf(output_file, ",direct-reverse:%d-%d\n", last_valid_direct_count, last_valid_reverse_count);
				}
			coverage+=shiftTo-at;
			state=1;
			}
		else state=0;

		at=shiftTo;
		if(at>=trueStringLength) end=true;
		}

//	printf("\n\ncoverage: %d/%d\n", coverage, trueStringLength-1);
	printf("%f\n", (double(coverage))/((double)(trueStringLength-1)));

	if(fileoutput) {
		fprintf(output_file, "\ncoverage: %d/%d\n", coverage, trueStringLength-1);
		fprintf(output_file, "%f\n", (double(coverage))/((double)(trueStringLength-1)));
		if(fclose(output_file)!=0) {
			fprintf(stderr, "failed fclose of output");
			exit(EXIT_FAILURE);
			}
		}

//	printf("\ncoverage: %d/%d\n", coverage, trueStringLength, subSeqId);
//	printf("%f (%d sequences used)\n\n", double(coverage)/double(trueStringLength-1), subSeqId);
///	printf("%f", double(coverage)/double(trueStringLength-1));
	delete[] st_leaves;
	st_free_node(st_root);
	delete[] rst_leaves;
	st_free_node(rst_root);
	delete[] s;
	delete[] rs;
	delete[] t;
//	printf("termination reached without errors\n");
	exit(EXIT_SUCCESS);
	}

void initialize_new_match(st_search* match_struct) {
	match_struct->current_node=match_struct->root;
	match_struct->current_depth=0;
	match_struct->current_edge_depth=0;
	match_struct->current_count=1;
	}

void intestation() {
	if(system("clear")==-1) {
		fprintf(stderr, "failed system call\n");
		exit(EXIT_FAILURE);
		}
	printf("/*\n");
	printf(" * bpmatch calculate, from sequences S and T, the maximum coverage of T\n");
	printf(" * using only subsequences and complemented reversed subsequences of S,\n");
	printf(" * with minimum length \"l\", with at least \"minRep\" occurrences in S,\n");
	printf(" * possibly overlapped, and, in such a maximum coverage, minimize the\n");
	printf(" * number of subsequences used.\n");
	printf(" * version 1.6\n");
	printf(" *\n");
	printf(" * Copyright (C) 2003-2012 Claudio Felicioli\n");
	printf(" * mail: c.felicioli@1d20.net - pangon@gmail.com\n");
	printf(" *\n");
	printf(" * bpmatch is free software; you can redistribute it and/or modify\n");
	printf(" * it under the terms of the GNU General Public License as published by\n");
	printf(" * the Free Software Foundation; either version 3 of the License, or\n");
	printf(" * (at your option) any later version.\n");
	printf(" *\n");
	printf(" * This program is distributed in the hope that it will be useful,\n");
	printf(" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
	printf(" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n");
	printf(" * GNU General Public License for more details.\n");
	printf(" *\n");
	printf(" * You should have received a copy of the GNU General Public License\n");
	printf(" * along with this program. If not, see <http://www.gnu.org/licenses/>.\n");
	printf("*/\n\n");
	}

void st_unserialize_node(FILE* file, st_node* node, st_node** leaves) {
	if(fread(&(node->length), sizeof(int), 1, file)!=1) {
		fprintf(stderr, "failed st_unserialize_node of length\n");
		exit(EXIT_FAILURE);
		}
	if(fread(&(node->leafId), sizeof(int), 1, file)!=1) {
		fprintf(stderr, "failed st_unserialize_node of leafId\n");
		exit(EXIT_FAILURE);
		}
	if(node->leafId!=0) leaves[node->leafId]=node;
	if(fread(&(node->firstLeaf), sizeof(int), 1, file)!=1) {
		fprintf(stderr, "failed st_unserialize_node of firstLeaf\n");
		exit(EXIT_FAILURE);
		}
	if(fread(&(node->lastLeaf), sizeof(int), 1, file)!=1) {
		fprintf(stderr, "failed st_unserialize_node of lastLeaf\n");
		exit(EXIT_FAILURE);
		}
	bool presence;
	for(int i=0;i<5;i++) {
		if(fread(&presence, sizeof(bool), 1, file)!=1) {
			fprintf(stderr, "failed st_unserialize_node of presence\n");
			exit(EXIT_FAILURE);
			}
		if(presence) {
			node->edges[i]=new st_edge;
			st_unserialize_edge(file, node->edges[i], leaves);
			}
		else node->edges[i]=NULL;
		}
	}

void st_unserialize_edge(FILE* file, st_edge* edge, st_node** leaves) {
	if(fread(&(edge->count), sizeof(int), 1, file)!=1) {
		fprintf(stderr, "failed st_unserialize_edge of count\n");
		exit(EXIT_FAILURE);
		}
	if(fread(&(edge->start), sizeof(int), 1, file)!=1) {
		fprintf(stderr, "failed st_unserialize_edge of start\n");
		exit(EXIT_FAILURE);
		}
	if(fread(&(edge->end), sizeof(int), 1, file)!=1) {
		fprintf(stderr, "failed st_unserialize_edge of end\n");
		exit(EXIT_FAILURE);
		}
	edge->to=new st_node;
	st_unserialize_node(file, edge->to, leaves);
	}

bool composite_match(st_search* suffix_tree_a, st_search* suffix_tree_b, bpmatch_utils_base base) {
	suffix_tree_a->current_count=single_match(suffix_tree_a, base);
	suffix_tree_b->current_count=single_match(suffix_tree_b, base);
	if(suffix_tree_a->current_count+suffix_tree_b->current_count>=minRep) {
		last_valid_direct_count=suffix_tree_a->current_count;
		last_valid_reverse_count=suffix_tree_b->current_count;
		return(true);
		}
	return(false);
	}

bool composite_match_target_check(st_search* suffix_tree_a, st_search* suffix_tree_b, st_search* suffix_tree_c, st_search* suffix_tree_d, bpmatch_utils_base base) {
	suffix_tree_a->current_count=single_match(suffix_tree_a, base);
	suffix_tree_b->current_count=single_match(suffix_tree_b, base);
	suffix_tree_c->current_count=single_match(suffix_tree_c, base);
	suffix_tree_d->current_count=single_match(suffix_tree_d, base);
	if(suffix_tree_a->current_count+suffix_tree_b->current_count>=minRep && suffix_tree_c->current_count+suffix_tree_d->current_count>=minRep2) {
		last_valid_direct_count=suffix_tree_a->current_count;
		last_valid_reverse_count=suffix_tree_b->current_count;
		last_valid_direct_count_target=suffix_tree_c->current_count;
		last_valid_reverse_count_target=suffix_tree_d->current_count;
		return(true);
		}
	return(false);
	}

int single_match(st_search* suffix_tree, bpmatch_utils_base base) {
	if(suffix_tree->current_count==0) return(0);
	if(DEBUG) printf("match %c [node leafId=%d]\n", base2char[base], suffix_tree->current_node->leafId);
	if(base==X) return(0);
	if(suffix_tree->current_node->leafId==0) {
		//internal node
		if(suffix_tree->current_depth<suffix_tree->current_edge_depth) {
			suffix_tree->current_depth++;
			if(suffix_tree->string[suffix_tree->current_edge->start+suffix_tree->current_depth]!=base) return(0);
			return(suffix_tree->current_count);
			}
		else {
			if(suffix_tree->current_node->edges[base]==NULL) return(0);
			if(DEBUG) printf("OK-%d (not leaf)\n", suffix_tree->current_node->edges[base]->count);
			suffix_tree->current_edge=suffix_tree->current_node->edges[base];
			suffix_tree->current_edge_depth=suffix_tree->current_edge->end-suffix_tree->current_edge->start;
			suffix_tree->current_depth=0;
			suffix_tree->current_node=suffix_tree->current_edge->to;
			return(suffix_tree->current_edge->count);
			}
		}
	else {
		//leaf
		suffix_tree->current_depth++;
		if(DEBUG) printf("LEAF: %c\n", base2char[suffix_tree->string[suffix_tree->current_edge->start+suffix_tree->current_depth]]);
		if(suffix_tree->string[suffix_tree->current_edge->start+suffix_tree->current_depth]!=base) return(0);
		if(DEBUG) printf("OK-1 (leaf)\n");
		return(1);
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

void st_free_node(st_node* node) {
	for(int i=0;i<5;i++) {
		if(node->edges[i]!=NULL) {
			st_free_node(node->edges[i]->to);
			delete node->edges[i];
			}
		}
	delete node;
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
