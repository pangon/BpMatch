/*
 * Copyright (C) 2007-2011 Claudio Felicioli
 * mail: c.felicioli@1d20.net - pangon@gmail.com
 *
 * complRev is free software; you can redistribute it and/or modify
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


enum bpmatch_utils_base {A, G, T, C, X};

int scanBp(FILE* file, bpmatch_utils_base* bp);

bpmatch_utils_base complement[5];
char base2char[5];
char inputName[50];
char outputName[50];
FILE* input;
FILE* output;
bool scanBp_fasta;

int main(int argc, char *argv[]) {
	if(argc!=3) {
		printf("usage: %s input output\n", argv[0]);
		exit(EXIT_FAILURE);
		}
	strncpy(inputName, argv[1], 48);
	inputName[49]='\0';
	strncpy(outputName, argv[2], 48);
	outputName[49]='\0';
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
	if((input=fopen(inputName, "r"))==NULL) {
		fprintf(stderr, "error: %s opening fail\n", inputName);
		exit(EXIT_FAILURE);
		}
	struct stat inputInfo;
	stat(inputName, &inputInfo);
	int l=inputInfo.st_size;
	bpmatch_utils_base s[l];
	bpmatch_utils_base s2[l];
	bpmatch_utils_base bp;
	int trueL=0;
	scanBp_fasta=false;
	while(scanBp(input, &bp)!=EOF) s[trueL++]=bp;
	int target_file_fd;
	if((target_file_fd=fileno(input))==-1) {
		fprintf(stderr, "failed getting gene file descriptor");
		exit(EXIT_FAILURE);
		}
	if(fsync(target_file_fd)!=0) {
		fprintf(stderr, "failed fsync of gene");
		exit(EXIT_FAILURE);
		}
	if(fclose(input)!=0) {
		fprintf(stderr, "failed fclose of gene");
		exit(EXIT_FAILURE);
		}
	for(int i=0;i<trueL;i++) s2[i]=complement[s[trueL-1-i]];
	if((output=fopen(outputName, "w"))==NULL) {
		fprintf(stderr, "error: %s opening fail\n", outputName);
		exit(EXIT_FAILURE);
		}
	for(int i=0;i<trueL;i++) fprintf(output, "%c", base2char[s2[i]]);
	fprintf(output, "\n");
	if((target_file_fd=fileno(output))==-1) {
		fprintf(stderr, "failed getting gene file descriptor");
		exit(EXIT_FAILURE);
		}
	if(fsync(target_file_fd)!=0) {
		fprintf(stderr, "failed fsync of gene");
		exit(EXIT_FAILURE);
		}
	if(fclose(output)!=0) {
		fprintf(stderr, "failed fclose of gene");
		exit(EXIT_FAILURE);
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
