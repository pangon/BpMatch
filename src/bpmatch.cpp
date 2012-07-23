/*
 * bpmatch calculate, from sequences S and T, the maximum coverage of T
 * using only subsequences and complemented reversed subsequences of S,
 * with minimum length "l", with at least "minRep" occurrences in S,
 * possibly overlapped, and, in such a maximum coverage, minimize the
 * number of subsequences used.
 * version 1.5
 *
 * Copyright (C) 2003-2011 Claudio Felicioli
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

void intestation();
void st_unserialize_node(FILE* file, st_node* node, st_node** leaves);
void st_unserialize_edge(FILE* file, st_edge* edge, st_node** leaves);
void st_print(st_node* node);
void st_free_node(st_node* node);
st_node* st_match(st_node* node, bpmatch_utils_base base, bpmatch_utils_base* string);
int scanBp(FILE* file, bpmatch_utils_base* bp);
void epilogo();

//need for match
st_edge* st_lastEdge;
int st_matchTmp;
int tmp_depth;

char st_fileName[50];
char rst_fileName[50];
char target_fileName[50];
char output_fileName[100];
int l;
int minRep;
FILE* st_file;
FILE* rst_file;
FILE* target_file;
FILE* output_file;
int stringLength;
bpmatch_utils_base* s;
bpmatch_utils_base* rs;
bpmatch_utils_base* t;
st_node* st_root;
st_node** st_leaves;
st_node* rst_root;
st_node** rst_leaves;
int at;
bool fileoutput;
bool scanBp_fasta;

int main(int argc, char *argv[]) {
	if(argc!=6 && argc!=7) {
		printf("usage: %s sourceGeneSuffixTree reverseSourceGeneSuffixTree TargetGene minimumLength minimumRepetition [outputFile]\n", argv[0]);
		exit(EXIT_FAILURE);
		}
	strncpy(st_fileName, argv[1], 48);
	st_fileName[49]='\0';
	strncpy(rst_fileName, argv[2], 48);
	rst_fileName[49]='\0';
	strncpy(target_fileName, argv[3], 48);
	target_fileName[49]='\0';

	fileoutput=(argc==7);
	if(fileoutput) {
		strncpy(output_fileName, argv[6], 98);
		output_fileName[99]='\0';
		if((output_file=fopen(output_fileName, "w"))==NULL) {
			fprintf(stderr, "error: %s opening fail\n", output_fileName);
			exit(EXIT_FAILURE);
			}
		}

	l=atoi(argv[4]);
	minRep=atoi(argv[5]);
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
//	for(int i=0;i<stringLength;i++) printf("%c", base2char[s[i]]);
//	printf("\n");
//	st_print(st_root);
//	printf("suffix tree unserialized.\n");
//	printf("Loading %s.\n", rst_fileName);
	if((rst_file=fopen(rst_fileName, "r"))==NULL) {
		fprintf(stderr, "error: %s opening fail\n", rst_fileName);
		exit(EXIT_FAILURE);
		}
	fread(&(stringLength), sizeof(int), 1, rst_file);
	rs=new bpmatch_utils_base[stringLength];
	fread(rs, sizeof(bpmatch_utils_base), stringLength, rst_file);
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
	if((target_file=fopen(target_fileName, "r"))==NULL) {
		fprintf(stderr, "error: %s opening fail\n", target_fileName);
		exit(EXIT_FAILURE);
		}
	struct stat target_fileInfo;
	stat(target_fileName, &target_fileInfo);
	stringLength=target_fileInfo.st_size;
	t=new bpmatch_utils_base[stringLength+1];
	int trueStringLength=0;
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
//	printf("\n");
//	printf("*************************************************\n");
//	printf("matching - l=%d - minRep=%d\n", l, minRep);
//	printf("*************************************************\n");
	bool end=false;
	bool found, foundR;
	int shiftTo, shiftToR;
	bool dontdo=false;
	bool dontdoR=false;
	int maxPrefix, maxSuffix;
	int state=0;
	at=0;
	int i;
	st_node* tmp_node;
//	int subSeqId=0;
	int coverage=0;
	while(!end) {
		found=false;
		foundR=false;
		if(state==0) {
			//case 0
			if(DEBUG) printf("case 0\n");
			if(dontdo) dontdo=false;
			else {
				//searching for direct seq
				if(DEBUG) printf("searching for direct seq...\n");
				i=l-1;
				//using reverse from actual+l-1 to actual (backward parsing)
				tmp_node=rst_root;
				tmp_depth=0;
				st_matchTmp=0;
				while(i>=0 && at+i<trueStringLength && (tmp_node=st_match(tmp_node, complement[t[at+i]], rs))!=NULL) i--;
				if(i>=0) {
					//not found
					shiftTo=at+i+1;
					if(DEBUG) printf("not found [%d -> %d]\n", at, shiftTo);
					}
				else {
					//direct sequence found
					i=0;
					//using direct from actual (forward parsing)
					tmp_node=st_root;
					tmp_depth=0;
					st_matchTmp=0;
					while(at+i<trueStringLength && (tmp_node=st_match(tmp_node, t[at+i], s))!=NULL) i++;
					found=true;
					shiftTo=at+i;
					if(DEBUG) printf("***************************************************found [%d -> %d]\n", at, shiftTo);
					}
				}
			if(dontdoR) dontdoR=false;
			else {
				//searching for reverse seq
				if(DEBUG) printf("searching for reverse seq...\n");
				i=l-1;
				//using direct from actual+l-1 to actual (backward parsing)
				tmp_node=st_root;
				tmp_depth=0;
				st_matchTmp=0;
				while(i>=0 && at+i<trueStringLength && (tmp_node=st_match(tmp_node, complement[t[at+i]], s))!=NULL) i--;
				if(i>=0) {
					//not found
					shiftToR=at+i+1;
					if(DEBUG) printf("not found [%d -> %d]\n", at, shiftToR);
					}
				else {
					//reverse sequence found
					i=0;
					//using reverse from actual (forward parsing)
					tmp_node=rst_root;
					tmp_depth=0;
					st_matchTmp=0;
					while(at+i<trueStringLength && (tmp_node=st_match(tmp_node, t[at+i], rs))!=NULL) i++;
					foundR=true;
					shiftToR=at+i;
					if(DEBUG) printf("##################################################found [%d -> %d]\n", at, shiftToR);
					}
				}
			}
		else {
			//case 1
			if(DEBUG) printf("case 1\n");
			//searching for direct (eventually ovelapping) seq
			if(DEBUG) printf("searching for direct (eventually ovelapping) seq\n");
			maxSuffix=0;
			maxPrefix=0;
			//using direct from actual (forward parsing)
			tmp_node=st_root;
			tmp_depth=0;
			st_matchTmp=0;
			while(at+maxSuffix<trueStringLength && (tmp_node=st_match(tmp_node, t[at+maxSuffix], s))!=NULL) maxSuffix++;
			if(maxSuffix>=l) {
				//direct sequence found
				found=true;
				shiftTo=at+maxSuffix;
				if(DEBUG) printf("1**************************************************found [%d -> %d]\n", at, shiftTo);
				}
			else {
				if(DEBUG) printf("compute maxPrefix\n");
				//using reverse from actual (backward parsing)
				tmp_node=rst_root;
				tmp_depth=0;
				st_matchTmp=0;
				while(maxPrefix<l && at-maxPrefix>=0 && (tmp_node=st_match(tmp_node, complement[t[at-maxPrefix]], rs))!=NULL) maxPrefix++;
				maxPrefix--;
				int ite=0;
				while(!found && maxPrefix+maxSuffix>=l) {
					if(DEBUG) printf("*ite %d-%d\n", maxPrefix, maxSuffix);
					if(ite++>8) exit(1);
					i=maxSuffix-l;
					//using direct from actual (forward parsing)
					tmp_node=st_root;
					tmp_depth=0;
					st_matchTmp=0;
					while(at+i<trueStringLength && (tmp_node=st_match(tmp_node, t[at+i], s))!=NULL) i++;
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
			//searching for reverse (eventually ovelapping) seq
			if(DEBUG) printf("searching for reverse (eventually ovelapping) seq\n");
			maxSuffix=0;
			maxPrefix=0;
			//using reverse from actual (forward parsing)
			tmp_node=rst_root;
			tmp_depth=0;
			st_matchTmp=0;
			while(at+maxSuffix<trueStringLength && (tmp_node=st_match(tmp_node, t[at+maxSuffix], rs))!=NULL) maxSuffix++;
			if(maxSuffix>=l) {
				//direct sequence found
				foundR=true;
				shiftToR=at+maxSuffix;
				if(DEBUG) printf("1############################################found [%d -> %d]\n", at, shiftTo);
				}
			else {
				if(DEBUG) printf("compute maxPrefix\n");
				//using direct from actual (backward parsing)
				tmp_node=st_root;
				tmp_depth=0;
				st_matchTmp=0;
				while(maxPrefix<l && at-maxPrefix>=0 && (tmp_node=st_match(tmp_node, complement[t[at-maxPrefix]], s))!=NULL) maxPrefix++;
				maxPrefix--;
				while(!foundR && maxPrefix+maxSuffix>=l) {
					if(DEBUG) printf("#ite\n");
					i=maxSuffix-l;
					//using reverse from actual (forward parsing)
					tmp_node=rst_root;
					tmp_depth=0;
					st_matchTmp=0;
					while(at+i<trueStringLength && (tmp_node=st_match(tmp_node, t[at+i], rs))!=NULL) i++;
					if(i==maxSuffix) {
						//direct sequence found
						foundR=true;
						shiftToR=at+maxSuffix;
						}
					else maxSuffix=i;
					}
				if(!foundR) {
					//not found
					shiftToR=at+1;
					}
				}
			}
		if(!found && !foundR) {
			if(shiftToR==shiftTo) at=shiftTo;
			else {
				if(shiftToR<shiftTo) {
					at=shiftToR;
					dontdo=true;
					}
				else {
					at=shiftTo;
					dontdoR=true;
					}
				}
			state=0;
			}
		if(found && foundR) {
			if(shiftToR<shiftTo) {
				if(fileoutput) {
					fprintf(output_file, "%i,", at);
					if(state==1) fprintf(output_file, "*");
					for(int j=at;j<shiftTo;j++) fprintf(output_file, "%c", base2char[t[j]]);
					fprintf(output_file, ",direct\n");
					}
/*				printf("%d,direct,", subSeqId++);
				for(int j=at;j<shiftTo;j++) printf("%c", base2char[t[j]]);
				if(state==1) printf(",<=%d,%d\n", at, shiftTo-1);
				else printf(",%d,%d\n", at, shiftTo-1);*/
				coverage+=shiftTo-at;
				at=shiftTo;
				}
			else {
				if(fileoutput) {
					fprintf(output_file, "%i,", at);
					if(state==1) fprintf(output_file, "*");
					for(int j=at;j<shiftToR;j++) fprintf(output_file, "%c", base2char[t[j]]);
					fprintf(output_file, ",reverse\n");
					}
/*				printf("%d,reverse,", subSeqId++);
				for(int j=at;j<shiftToR;j++) printf("%c", base2char[t[j]]);
				if(state==1) printf(",<=%d,%d\n", at, shiftToR-1);
				else printf(",%d,%d\n", at, shiftToR-1);*/
				coverage+=shiftToR-at;
				at=shiftToR;
				}
			state=1;
			}
		if(found && !foundR) {
			if(fileoutput) {
				fprintf(output_file, "%i,", at);
				if(state==1) fprintf(output_file, "*");
				for(int j=at;j<shiftTo;j++) fprintf(output_file, "%c", base2char[t[j]]);
				fprintf(output_file, ",direct\n");
				}
/*			printf("%d,direct,", subSeqId++);
			for(int j=at;j<shiftTo;j++) printf("%c", base2char[t[j]]);
			if(state==1) printf(",<=%d,%d\n", at, shiftTo-1);
			else printf(",%d,%d\n", at, shiftTo-1);*/
			coverage+=shiftTo-at;
			at=shiftTo;
			state=1;
			}
		if(!found && foundR) {
			if(fileoutput) {
				fprintf(output_file, "%i,", at);
				if(state==1) fprintf(output_file, "*");
				for(int j=at;j<shiftToR;j++) fprintf(output_file, "%c", base2char[t[j]]);
				fprintf(output_file, ",reverse\n");
				}
/*			printf("%d,reverse,", subSeqId++);
			for(int j=at;j<shiftToR;j++) printf("%c", base2char[t[j]]);
			if(state==1) printf(",<=%d,%d\n", at, shiftToR-1);
			else printf(",%d,%d\n", at, shiftToR-1);*/
			coverage+=shiftToR-at;
			at=shiftToR;
			state=1;
			}
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

void intestation() {
	system("clear");
	printf("/*\n");
	printf(" * bpmatch calculate, from sequences S and T, the maximum coverage of T\n");
	printf(" * using only subsequences and complemented reversed subsequences of S,\n");
	printf(" * with minimum length \"l\", with at least \"minRep\" occurrences in S,\n");
	printf(" * possibly overlapped, and, in such a maximum coverage, minimize the\n");
	printf(" * number of subsequences used.\n");
	printf(" * version 1.5\n");
	printf(" *\n");
	printf(" * Copyright (C) 2003-2011 Claudio Felicioli\n");
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

//st_edge* st_lastEdge;
//int st_matchTmp;
//int tmp_depth;

st_node* st_match(st_node* node, bpmatch_utils_base base, bpmatch_utils_base* string) {
	if(DEBUG) printf("match %c [node leafId=%d]\n", base2char[base], node->leafId);
//	if(DEBUG) if(node->lastLeaf-node->firstLeaf<20) st_print(node);
	if(base==X) return(NULL);
	if(node->leafId==0) {
		//internal node
		if(st_matchTmp<tmp_depth) {
			st_matchTmp++;
			if(string[st_lastEdge->start+st_matchTmp]!=base) return(NULL);
			return(node);
			}
		else {
			if(node->edges[base]==NULL) return(NULL);
			if(node->edges[base]->count<minRep) return(NULL);
			if(DEBUG) printf("OK-%d (not leaf)\n", node->edges[base]->count);
			st_lastEdge=node->edges[base];
			tmp_depth=st_lastEdge->end-st_lastEdge->start;
			st_matchTmp=0;
			return(node->edges[base]->to);
			}
		}
	else {
		//leaf
		st_matchTmp++;
		if(1<minRep) return(NULL);
		if(DEBUG) printf("LEAF: %c\n", base2char[string[st_lastEdge->start+st_matchTmp]]);
		if(string[st_lastEdge->start+st_matchTmp]!=base) return(NULL);
		if(DEBUG) printf("OK-1 (leaf)\n");
		return(node);
		}
//	if(DEBUG) printf("%c found (rep %d)\n", base2char[base], node->edges[base]->count);
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
