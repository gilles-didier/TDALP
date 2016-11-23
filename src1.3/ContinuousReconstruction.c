/*
    'tdalp' computes the time-dependent-asymmetric-linear-parsimony parametric reconstruction of a phylogenetic tree.

    Copyright (C) 2016  Gilles DIDIER

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Utils.h"
#include "Tree.h"
#include "StateTree.h"
#include "ContinuousParsimony.h"
#include "DrawTreeCairo.h"
#include "DrawTreePSTricks.h"
#include "DrawTreeGeneric.h"
#include "DrawTreeResultGeneric.h"

#define OUTPUT "output"

#define MIN_TIME 0.001

#define SIZE_BUFFER_CHAR 300
#define HELPMESSAGE "NAME\n\ttdalp - compute time-dependent-asymmetric-linear-parsimony parametric reconstruction\n\t\nSYNOPSIS\n\ttdalp [OPTIONS] <inputFile> <outputFile>\n\nDESCRIPTION\n\tCompute the time-dependent-asymmetric-linear-parsimony parametric reconstruction of the tree contained in <inputFile> (it must be in Newick format). It outputs the parametric reconstruction in text format in the file <outputFile> and in graphic format in the file <outputFile> with the corresponding extension (i.e .pdf, .png etc.)\n\n\tOptions are\n\t-s <file>\n\t\tload the file containing the character states in '.csv' format : a series of lines '<name of the node> <state value>'\n\t-t <type>\n\t\tset the type of function applied on the branch length of the tree: \n\t\t\t-t i -> inverse\n\t\t\t-t u -> identity\n\t-f <format>\n\t\tset the graphic format of the output\n\t\t\t-f 1 -> pdf\n\t\t\t-f 2 -> postscript\n\t\t\t-f 3 -> png\n\t\t\t-f 4 -> svg\n\t\t\t-f 5 -> pstricks\n\t-x <file>\n\t\toutput the detail of the parametric reconstruction of the root in the file <file>\n\t-c <r1> <g1> <b1> <r2> <g2> <b2>\n\t\tset the color scale going for the color (r1,g1,b1) to the color (r2,g2,b2) (rgb code of colors have components in the range [0,1])\n\t-h\n\t\tdisplay help\n"

double fUn(double t) {
	return 1;
}

double fInv(double t) {
	return 1/t;
}

int main(int argc, char **argv) {		
	char inputFileNameTree[SIZE_BUFFER_CHAR], inputFileNameState[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], outputFileNameRoot[SIZE_BUFFER_CHAR], 
	option[256], type = 'i', format = '0';
	int i;
	FILE *fi, *fo, *fs;
	TypeInfoDrawTreeGeneric info;
	
	info.param.start = (TypeRGB) {.red = 1., .green = 1., .blue = 0.};
	info.param.end = (TypeRGB) {.red = 1., .green = 0., .blue = 0.};
	inputFileNameState[0] = '\0';
	outputFileNameRoot[0] = '\0';
	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%s", inputFileNameState) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -s");
		}
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%s", outputFileNameRoot) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -s");
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &type) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a character is required after option -t");
		}
		if(option['c']) {
			option['c'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.start.red)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.start.green)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.start.blue)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.end.red)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.end.green)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.end.blue)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a character is required after option -f");
		}
		if(option['h']) {
			printf("%s\n", HELPMESSAGE);
			exitProg(ExitOk, NULL);
		}
	}
	if (i>=argc || sscanf(argv[i++], "%s", inputFileNameTree) != 1) exitProg(ErrorArgument, "input file tree");
	if (i>=argc || sscanf(argv[i++], "%s", outputFileName) != 1) exitProg(ErrorArgument, "output file");
	if((fi = fopen(inputFileNameTree, "r"))) {
		TypeTree *tree;
		TypeResult *res;
		int n;
		tree = readTree(fi);
		fixTreeTime(tree);
		fillState(tree);
		if(strlen(inputFileNameState)>0 && (fs = fopen(inputFileNameState, "r")))
			setState(fs, tree);
		switch(type) {
			case 'i':
				res = getResult(tree, fInv);
				break;
			case 'u':
			default:
				res = getResult(tree, fUn);
		}
		if((fo = fopen(outputFileName, "w"))) {
			fprintStateTree(fo, tree, display_time_both);
			fprintf(fo, "\n\n");
			fprintResultList(fo, res, tree);
			fclose(fo);
		} else
			exitProg(ErrorWriting, outputFileName);
		if(outputFileNameRoot[0] != '\0' && format != '0') {
			char *tmp, outputFileNameG[SIZE_BUFFER_CHAR];
			if((tmp = strrchr(outputFileNameRoot, '.')) != NULL)
				tmp[0] = '\0';
				info.param.vmax = res[0].val[0];
				info.param.vmin = res[0].val[0];
				for(i=0; i<tree->size; i++) {
					int j;
					for(j=0; j<=res[i].sizePar; j++) {
						if(res[i].val[j]>info.param.vmax)
							info.param.vmax = res[i].val[j];
						if(res[i].val[j]<info.param.vmin)
							info.param.vmin = res[i].val[j];
					}
				}
				if(info.param.vmax==info.param.vmin)
					info.param.vmax++;
				info.param.vscale = 1/(info.param.vmax-info.param.vmin);
				switch(format) {
				case '1':
					sprintf(outputFileNameG, "%s.pdf", outputFileNameRoot);
					setFunctPDF(&(info.funct));
					drawSpecialResultFileGeneric(outputFileNameG, &(res[tree->root]), &info);
					break;
				case '2':
					sprintf(outputFileNameG, "%s.ps", outputFileNameRoot);
					setFunctPS(&(info.funct));
					drawSpecialResultFileGeneric(outputFileNameG, &(res[tree->root]), &info);
					break;
				case '3':
					sprintf(outputFileNameG, "%s.png", outputFileNameRoot);
					setFunctPNG(&(info.funct));
					drawSpecialResultFileGeneric(outputFileNameG, &(res[tree->root]), &info);
					break;
				case '4':
					sprintf(outputFileNameG, "%s.svg", outputFileNameRoot);
					setFunctSVG(&(info.funct));
					drawSpecialResultFileGeneric(outputFileNameG, &(res[tree->root]), &info);
					break;
				case '5':
					sprintf(outputFileNameG, "%s.tex", outputFileNameRoot);
					setFunctPSTricks(&(info.funct));
					drawSpecialResultFileGeneric(outputFileNameG, &(res[tree->root]), &info);
					break;
				default:
				;
			}
		}
		if(format != '0') {
			char *tmp, outputFileNameG[SIZE_BUFFER_CHAR];
			if((tmp = strrchr(outputFileName, '.')) != NULL)
				tmp[0] = '\0';
			bltoabsTime(tree);
			if(tree->minTime == NO_TIME || tree->minTime == 0.)
				tree->minTime = tree->time[tree->root]*0.9;
			if(tree->maxTime == NO_TIME) {
				tree->maxTime = 0.;
				for(n=0; n<tree->size; n++)
					if(tree->time[n]>tree->maxTime)
						tree->maxTime = tree->time[n];
			}
			switch(format) {
				case '1':
					sprintf(outputFileNameG, "%s.pdf", outputFileName);
					setFunctPDF(&(info.funct));
					drawTreeResultFileGeneric(outputFileNameG, tree, res, &info);
					break;
				case '2':
					sprintf(outputFileNameG, "%s.ps", outputFileName);
					setFunctPS(&(info.funct));
					drawTreeResultFileGeneric(outputFileNameG, tree, res, &info);
					break;
				case '3':
					sprintf(outputFileNameG, "%s.png", outputFileName);
					setFunctPNG(&(info.funct));
					drawTreeResultFileGeneric(outputFileNameG, tree, res, &info);
					break;
				case '4':
					sprintf(outputFileNameG, "%s.svg", outputFileName);
					setFunctSVG(&(info.funct));
					drawTreeResultFileGeneric(outputFileNameG, tree, res, &info);
					break;
				case '5':
					sprintf(outputFileNameG, "%s.tex", outputFileName);
					setFunctPSTricks(&(info.funct));
					drawTreeResultFileGeneric(outputFileNameG, tree, res, &info);
					break;
				default:
				;
			}
		}
		for(n=0; n<tree->size; n++)
			freeResult(res[n]);
		free((void*)res);
		if(tree->info != NULL) {
			free((void*)tree->info);
			tree->info = NULL;
		}
		freeTree(tree);
	} else
		exitProg(ErrorReading, inputFileNameTree);
	return 0;
}

