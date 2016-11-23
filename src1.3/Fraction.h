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




#ifndef FractionF
#define FractionF
#include <stdlib.h>
#include <stdio.h>

typedef struct FRACTION {
	double num, den;
} TypeFraction;

TypeFraction fract(double num, double den);
int cmpFract(TypeFraction a, TypeFraction b);
double fractToPercent(TypeFraction a);
void cannonizeFract(TypeFraction *fr);
void fprintFracLatex(FILE *f, TypeFraction frac);
void fprintFrac(FILE *f, TypeFraction frac);
void sprintFrac(char *s, TypeFraction frac);
void sprintFracLatex(char *s, TypeFraction frac);
void fprintfracDebug(FILE *f, double n, double d);

#endif
