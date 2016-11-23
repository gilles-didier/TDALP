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




#ifndef ContinuousParsimonyF
#define ContinuousParsimonyF
#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"
#include "Fraction.h"
#include "StateTree.h"

typedef struct COEFF {
	double a, b, c, d;
} TypeCoeff;

typedef struct FRACTION_INTER {
	TypeFraction inf, sup;
} TypeFractionInter;

typedef struct COLUMN {
	int sizeVal;
	TypeCoeff *coeff;
	double *boundVal;
} TypeColumn;

typedef double TypeFunctTime(double);

/*having sizePar == -1 means that the state is known*/
typedef struct F {
	int sizePar;
	TypeFraction *boundPar;
	TypeColumn *col;
} TypeF;


typedef struct RESULT {
	int sizePar;
	TypeFraction *boundPar;
	double *val;
} TypeResult;


void fillParametricContinuousParsimony(TypeTree *tree, TypeFunctTime phi);
TypeResult *getResult(TypeTree *tree, TypeFunctTime phi);
TypeResult getResultFromF(TypeF *f, int unknown);
TypeResult getResultFromH(TypeF *h, TypeResult *a, int unknown);
TypeResult *fillResult(TypeF *funcFH, TypeTree *tree, TypeFunctTime phi);
void freeF(TypeF *f);
TypeListDouble *getErrors(TypeResult *res, TypeIndexStateList *sta, TypeTree *tree);
TypeFractionInter getBestInter(TypeResult *res, TypeIndexStateList *sta, TypeTree *tree);
void freeResult(TypeResult r);


/*print result as flat list*/	
void fprintResultList(FILE *f, TypeResult *res, TypeTree *tree);
void fprintFunc(FILE *f, TypeF *func);
void fprintFuncGnuplot(FILE *f, TypeF *func);
void fprintFuncDebug(FILE *f, TypeF *func);
void fprintResult(FILE *f, TypeResult *r);
void fprintResultDebug(FILE *f, TypeResult *r);
void fprintStateTreeResult(FILE *f, TypeResult *res, TypeTree *tree);
void sprintResult(char *s, TypeResult *r);
void sprintResultDebug(char *s,  TypeResult *r);
void fprintResultDebug(FILE *f, TypeResult *r);
void fprintColDebug(FILE *f, TypeColumn *c);

#endif
