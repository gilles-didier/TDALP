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




#ifndef DrawTreeResultGenericF
#define DrawTreeResultGenericF
#include <stdio.h>
#include "Tree.h"
#include "DrawTreeGeneric.h"
#include "ContinuousParsimony.h"

#ifdef __cplusplus
extern "C" {
#endif

void drawTreeResultFileGeneric(char *filename, TypeTree *tree, TypeResult *r, TypeInfoDrawTreeGeneric *info);
void drawSpecialResultFileGeneric(char *filename, TypeResult *r, TypeInfoDrawTreeGeneric *info);
void drawScaleVPst(double x, double y, double size, TypeInfoDrawTreeGeneric *info);


#ifdef __cplusplus
}
#endif

#endif
