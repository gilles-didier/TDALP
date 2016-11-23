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




#include <limits.h>
#include <string.h>
#include <math.h>
#include "ContinuousParsimony.h"

#define RADIUS 0.5
#define XOFFSET 0.0
#define YOFFSET 0.0
#define XOFF 0.2
#define YOFF -0.4
#define SIZEBUF 10
#define FZERO {0,1}
#define FINFTY {1,0}
#define FNINFTY {-1,0}
#define INFTY 1E20
#define SIZE_BUFFER_CHAR 2000
#define WIDTH "stringWidth"

static void getColumn(TypeColumn *colF, TypeColumn *colH, int jp, int jm, double timeFact);
static void FtoH(TypeF *f, TypeF *h, double timeFact, int known);
static void fillResultRec(int n, TypeResult *res, TypeF *funcFH, TypeTree *tree, TypeFunctTime phi);
static void freeCol(TypeColumn c);
static void freeFuncF(TypeF f);

void compactResult(TypeResult *r) {
	int i, ind = 0;
	for(i=0; i<r->sizePar; i++) {
		if(r->val[i+1] != r->val[i]) {
			r->val[ind] = r->val[i];
			r->boundPar[ind++] = r->boundPar[i];
			cannonizeFract(&(r->boundPar[ind-1]));
		}
	}
	r->val[ind] = r->val[r->sizePar];
	r->sizePar = ind;
	r->boundPar = (TypeFraction*) realloc((void*)r->boundPar, r->sizePar*sizeof(TypeFraction));
	r->val = (double*) realloc((void*)r->val, (r->sizePar+1)*sizeof(double));
}

void freeResult(TypeResult r) {
	if(r.boundPar != NULL)
		free((void*)r.boundPar);
	if(r.val != NULL)
		free((void*)r.val);
}

void freeCol(TypeColumn c) {
	if(c.boundVal != NULL)
		free((void*)c.boundVal);
	if(c.coeff != NULL)
		free((void*)c.coeff);
}

void freeFuncF(TypeF f) {
	if(f.boundPar != NULL)
		free((void*)f.boundPar);
	if(f.col != NULL) {
		int i;
		for(i=0; i<=f.sizePar; i++)
			freeCol(f.col[i]);
		free((void*)f.col);
	}
}
/*	
TypeResult getResultFromF(TypeF *f, int unknown) {
	TypeResult r;
	if(unknown) {
		int i, sizeBuf = SIZEBUF;
		r.boundPar = (TypeFraction*) malloc(sizeBuf*sizeof(TypeFraction));
		r.val = (double*) malloc((sizeBuf+1)*sizeof(double));
		r.sizePar = 0;
		for(i=0; i<=f->sizePar; i++) {
			TypeFraction finf, fsup;
			int j;
			if(i>0)
				finf = f->boundPar[i-1];
			else 
				finf = fract(0,1);			
			if(i<f->sizePar)
				fsup = f->boundPar[i];
			else
				fsup = fract(1,0);
			for(j=1; j<=f->col[i].sizeVal && cmpFract(fract(f->col[i].coeff[j].b, f->col[i].coeff[j].a), finf) <= 0; j++);
			for(;j<=f->col[i].sizeVal && cmpFract(fract(f->col[i].coeff[j].b, f->col[i].coeff[j].a), fsup) < 0; j++) {
				if(r.sizePar >= sizeBuf) {
					sizeBuf += SIZEBUF;
					r.boundPar = (TypeFraction*) realloc((void*) r.boundPar, sizeBuf*sizeof(TypeFraction));
					r.val = (double*) realloc((void*) r.val, (sizeBuf+1)*sizeof(double));
				}
				r.val[r.sizePar] = f->col[i].boundVal[j-1];
				r.boundPar[r.sizePar++] = fract(f->col[i].coeff[j].b, f->col[i].coeff[j].a);
			}
			r.val[r.sizePar] = f->col[i].boundVal[j-1];
			if(i<f->sizePar) {
				if(r.sizePar >= sizeBuf) {
					sizeBuf += SIZEBUF;
					r.boundPar = (TypeFraction*) realloc((void*) r.boundPar, sizeBuf*sizeof(TypeFraction));
					r.val = (double*) realloc((void*) r.val, (sizeBuf+1)*sizeof(double));
				}
				r.boundPar[r.sizePar++] = f->boundPar[i];
			}
		}
	} else {
		r.boundPar = NULL;
		r.val = (double*) malloc(sizeof(double));
		r.sizePar = 0;
		r.val[0] = f->col[0].boundVal[0];
	}
	return r;
}
*/

TypeResult getResultFromF(TypeF *f, int unknown) {
	TypeResult r;
	if(unknown) {
		int i, sizeBuf = SIZEBUF;
		r.boundPar = (TypeFraction*) malloc(sizeBuf*sizeof(TypeFraction));
		r.val = (double*) malloc((sizeBuf+1)*sizeof(double));
		r.sizePar = 0;
		for(i=0; i<=f->sizePar; i++) {
			TypeFraction finf, fsup;
			int j, jmin, jmax;
			if(i>0)
				finf = f->boundPar[i-1];
			else 
				finf = fract(0,1);			
			if(i<f->sizePar)
				fsup = f->boundPar[i];
			else
				fsup = fract(1,0);
			for(jmin=0; jmin<=f->col[i].sizeVal && (-finf.num*f->col[i].coeff[jmin].a+finf.den*f->col[i].coeff[jmin].b)<=0; jmin++);
			for(jmax=f->col[i].sizeVal+1; jmax>=0 && (-fsup.num*f->col[i].coeff[jmax-1].a+fsup.den*f->col[i].coeff[jmax-1].b)>=0; jmax--);
			for(j=jmin; j<=jmax; j++) {
				if(r.sizePar >= sizeBuf) {
					sizeBuf += SIZEBUF;
					r.boundPar = (TypeFraction*) realloc((void*) r.boundPar, sizeBuf*sizeof(TypeFraction));
					r.val = (double*) realloc((void*) r.val, (sizeBuf+1)*sizeof(double));
				}				
				r.val[r.sizePar] = f->col[i].boundVal[j-1];
				if(j<jmax)
					r.boundPar[r.sizePar++] = fract(f->col[i].coeff[j].b, f->col[i].coeff[j].a);
			}
			if(i<f->sizePar) 
				r.boundPar[r.sizePar++] = f->boundPar[i];
		}
		r.boundPar = (TypeFraction*) realloc((void*) r.boundPar, r.sizePar*sizeof(TypeFraction));
		r.val = (double*) realloc((void*) r.val, (r.sizePar+1)*sizeof(double));
		compactResult(&r);
	} else {
		r.boundPar = NULL;
		r.val = (double*) malloc(sizeof(double));
		r.sizePar = 0;
		r.val[0] = f->col[0].boundVal[0];
	}
	return r;
}

TypeResult getResultFromH(TypeF *h, TypeResult *a, int unknown) {
	TypeResult r;
	if(unknown) {
		TypeFraction minPar;
		int ih, ia, sizeBuf = SIZEBUF;
		r.boundPar = (TypeFraction*) malloc(sizeBuf*sizeof(TypeFraction));
		r.val = (double*) malloc((sizeBuf+1)*sizeof(double));
		r.sizePar = 0;
		ih = 0; ia = 0;
		do {
			if(r.sizePar >= sizeBuf) {
				sizeBuf += SIZEBUF;
				r.boundPar = (TypeFraction*) realloc((void*) r.boundPar, sizeBuf*sizeof(TypeFraction));
				r.val = (double*) realloc((void*) r.val, (sizeBuf+1)*sizeof(double));
			}
			if(h->col[ih].sizeVal == 0 || (a->val[ia] >= h->col[ih].boundVal[0] && a->val[ia] <= h->col[ih].boundVal[h->col[ih].sizeVal-1])) {
				r.val[r.sizePar] = a->val[ia];
			} else {
				if(a->val[ia] < h->col[ih].boundVal[0]) {
					r.val[r.sizePar] = h->col[ih].boundVal[0];
				} else {
					r.val[r.sizePar] = h->col[ih].boundVal[h->col[ih].sizeVal-1];
				}
			}
			minPar = fract(1,0);
			if(ih<h->sizePar && cmpFract(h->boundPar[ih],minPar)<0)
				minPar = h->boundPar[ih];
			if(ia<a->sizePar && cmpFract(a->boundPar[ia],minPar)<0)
				minPar = a->boundPar[ia];
			if(ih<h->sizePar && cmpFract(h->boundPar[ih],minPar)<=0)
				ih++;
			if(ia<a->sizePar && cmpFract(a->boundPar[ia],minPar)<=0)
				ia++;
			if(cmpFract(fract(1,0), minPar) != 0)
				r.boundPar[r.sizePar++] = minPar;
		} while(cmpFract(fract(1,0),minPar) != 0);		
	} else {
		r.sizePar = 0;
		r.boundPar = NULL;
		r.val = (double*) malloc(sizeof(double));
		r.val[0] = h->col[0].boundVal[0];
	}
	compactResult(&r);
	return r;	
}


void getColumn(TypeColumn *colF, TypeColumn *colH, int jp, int jm, double timeFact) {
	int j;
	colH->sizeVal = jm-jp-1;
	if(jp > 0)
		colH->sizeVal++;
	if(jm <= colF->sizeVal)
		colH->sizeVal++;
	colH->boundVal = (double*) malloc(colH->sizeVal*sizeof(double));
	colH->coeff = (TypeCoeff*) malloc((colH->sizeVal+1)*sizeof(TypeCoeff));
	colH->sizeVal = 0;
	if(jp > 0) {
		colH->coeff[colH->sizeVal].a = timeFact;
		colH->coeff[colH->sizeVal].b = 0.;
		colH->coeff[colH->sizeVal].c = colF->boundVal[jp-1]*(timeFact-colF->coeff[jp].a)+colF->coeff[jp].c;
		colH->coeff[colH->sizeVal].d = colF->boundVal[jp-1]*colF->coeff[jp].b+colF->coeff[jp].d;
		colH->boundVal[colH->sizeVal++] = colF->boundVal[jp-1];
	}
	if(jm <= colF->sizeVal) {
		for(j=jp; j<jm; j++) {
			colH->coeff[colH->sizeVal] = colF->coeff[j];
			colH->boundVal[colH->sizeVal++] = colF->boundVal[j];
		}
		colH->coeff[colH->sizeVal].a = 0.;
		colH->coeff[colH->sizeVal].b = timeFact;
		colH->coeff[colH->sizeVal].c = -colF->boundVal[jm-1]*colF->coeff[jm].a+colF->coeff[jm].c;
		colH->coeff[colH->sizeVal].d = colF->boundVal[jm-1]*(colF->coeff[jm].b-timeFact)+colF->coeff[jm].d;
	} else {
		for(j=jp; j<colF->sizeVal; j++) {
			colH->coeff[colH->sizeVal] = colF->coeff[j];
			colH->boundVal[colH->sizeVal++] = colF->boundVal[j];
		}
		colH->coeff[colH->sizeVal] = colF->coeff[colF->sizeVal];
	}
}


void FtoH(TypeF *f, TypeF *h, double timeFact, int unknown) {
	int sizeParBuf, i;
	TypeFraction boundp, boundm;
	if(unknown) {
		sizeParBuf = 2*f->sizePar+2;
		h->sizePar = 0;
		h->boundPar = (TypeFraction*) malloc(sizeParBuf*sizeof(TypeFraction));
		h->col = (TypeColumn*) malloc((sizeParBuf+1)*sizeof(TypeColumn));
		for(i=0; i<=f->sizePar; i++) {
			TypeFraction finf, fsup;
			int jp, jm, sjp, sjm, ejp, ejm;
			if(i>0)
				finf = f->boundPar[i-1];
			else 
				finf = fract(0,1);			
			if(i<f->sizePar)
				fsup = f->boundPar[i];
			else
				fsup = fract(1,0);
			if(f->col[i].coeff[0].a<timeFact) {
				sjp=0; ejp=0;
			} else {
				for(sjp=1; sjp<=f->col[i].sizeVal && cmpFract(fract(f->col[i].coeff[sjp].b, f->col[i].coeff[sjp].a-timeFact), finf)<=0; sjp++);
				for(ejp=sjp; ejp<=f->col[i].sizeVal && cmpFract(fract(f->col[i].coeff[ejp].b, f->col[i].coeff[ejp].a-timeFact), fsup)<0 && cmpFract(fract(f->col[i].coeff[ejp].b, f->col[i].coeff[ejp].a-timeFact), finf)>0; ejp++);
			}
			if(f->col[i].coeff[f->col[i].sizeVal].b<timeFact) {
				sjm=f->col[i].sizeVal+1; ejm=f->col[i].sizeVal+1;
			} else {
				for(ejm=f->col[i].sizeVal; ejm>0 && cmpFract(fract(f->col[i].coeff[ejm-1].b-timeFact, f->col[i].coeff[ejm-1].a), fsup)>=0; ejm--);
				for(sjm=ejm; sjm>0 && cmpFract(fract(f->col[i].coeff[sjm-1].b-timeFact, f->col[i].coeff[sjm-1].a), finf)>0 && cmpFract(fract(f->col[i].coeff[sjm-1].b-timeFact, f->col[i].coeff[sjm-1].a), fsup)<0; sjm--);
			}
			jp = sjp;
			jm = sjm;
			getColumn(&(f->col[i]), &(h->col[h->sizePar]), jp, jm, timeFact);
			while(jp < ejp || jm < ejm) {
				if(jp<ejp)
					boundp = fract(f->col[i].coeff[jp].b, f->col[i].coeff[jp].a-timeFact);
				else
					boundp = fsup;
				if(jm<ejm)
					boundm = fract(f->col[i].coeff[jm].b-timeFact, f->col[i].coeff[jm].a);
				else
					boundm = fsup;
				if(h->sizePar >= sizeParBuf) {
					sizeParBuf += f->sizePar+2;
					h->boundPar = (TypeFraction*) realloc((void*) h->boundPar, sizeParBuf*sizeof(TypeFraction));
					h->col = (TypeColumn*) realloc((void*) h->col, (sizeParBuf+1)*sizeof(TypeColumn));
				}
				if(cmpFract(boundp, boundm)<=0)
					h->boundPar[h->sizePar] = boundp;
				else
					h->boundPar[h->sizePar] = boundm;
				if(cmpFract(boundp, h->boundPar[h->sizePar])<=0)
					jp++;
				if(cmpFract(boundm, h->boundPar[h->sizePar])<=0)
					jm++;
				h->sizePar++;
				getColumn(&(f->col[i]), &(h->col[h->sizePar]), jp, jm, timeFact);
			}
			if(i<f->sizePar) {
				if(h->sizePar >= sizeParBuf) {
					sizeParBuf += f->sizePar+2;
					h->boundPar = (TypeFraction*) realloc((void*) h->boundPar, sizeParBuf*sizeof(TypeFraction));
					h->col = (TypeColumn*) realloc((void*) h->col, (sizeParBuf+1)*sizeof(TypeColumn));
				}
				h->boundPar[h->sizePar] = f->boundPar[i];
				h->sizePar++; 
			}		
		}
		h->boundPar = (TypeFraction*) realloc((void*) h->boundPar, h->sizePar*sizeof(TypeFraction));
		h->col = (TypeColumn*) realloc((void*) h->col, (h->sizePar+1)*sizeof(TypeColumn));
	} else {
		h->sizePar = f->sizePar;
		h->boundPar = (TypeFraction*) malloc(h->sizePar*sizeof(TypeFraction));
		h->col = (TypeColumn*) malloc((h->sizePar+1)*sizeof(TypeColumn));
		for(i=0; i<=h->sizePar; i++) {
			h->col[i].sizeVal = 1;
			h->col[i].boundVal = (double*) malloc(sizeof(double));
			h->col[i].boundVal[0] = f->col[i].boundVal[0];
			h->col[i].coeff = (TypeCoeff*) malloc(2*sizeof(TypeCoeff));
			h->col[i].coeff[0].a = timeFact;
			h->col[i].coeff[0].b = 0;
			h->col[i].coeff[0].c = f->col[i].boundVal[0]*(timeFact-f->col[i].coeff[0].a)+f->col[i].coeff[0].c;
			h->col[i].coeff[0].d = f->col[i].boundVal[0]*f->col[i].coeff[0].b+f->col[i].coeff[0].d;
			h->col[i].coeff[1].a = 0;
			h->col[i].coeff[1].b = timeFact;
			h->col[i].coeff[1].c = -f->col[i].boundVal[0]*f->col[i].coeff[0].a+f->col[i].coeff[0].c;
			h->col[i].coeff[1].d = f->col[i].boundVal[0]*(-timeFact+f->col[i].coeff[0].b)+f->col[i].coeff[0].d;
			if(i<h->sizePar)
				h->boundPar[i] = f->boundPar[i];
		}
	}
}

void fillResultRec(int n, TypeResult *res, TypeF *funcFH, TypeTree *tree, TypeFunctTime phi) {
	int c;
	for(c=tree->node[n].child; c >= 0; c=tree->node[c].sibling) {
		res[c] = getResultFromH(&(funcFH[c]), &(res[n]), isUnknown(c, tree));
		fillResultRec(c, res, funcFH, tree, phi);
	}
}

TypeResult *fillResult(TypeF *funcFH, TypeTree *tree, TypeFunctTime phi) {
	TypeResult *res;
	res = (TypeResult*) malloc(tree->size*sizeof(TypeResult));
	res[tree->root] = getResultFromF(&(funcFH[tree->root]), isUnknown(tree->root, tree));
	fillResultRec(tree->root, res, funcFH, tree, phi);
	return res;
}

void fillInfoParametric(int n, TypeF *funcFH, TypeTree *tree, TypeFunctTime phi) {
	int c, nchild;
	nchild = 0;
	for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling) {
		fillInfoParametric(c, funcFH, tree, phi);
		nchild++;
	}
	if(nchild>0) {
		int k, sizePar;
		TypeF *h;
		h = (TypeF*) malloc(nchild*sizeof(TypeF));
		k = 0;
		sizePar = 0;
		for(c=tree->node[n].child; c!= NOSUCH; c=tree->node[c].sibling) {
			FtoH(&(funcFH[c]), &(h[k]), phi(tree->time[c]), isUnknown(c, tree));
			freeFuncF(funcFH[c]);
			funcFH[c] = h[k];
			sizePar += h[k].sizePar;
			k++;
		}
		if(isUnknown(n, tree)) {
			TypeFraction minPar;
			int *curPar, *curVal;
			curPar = (int*) malloc(nchild*sizeof(int));
			curVal = (int*) malloc(nchild*sizeof(int));
			funcFH[n].sizePar = 0;
			funcFH[n].boundPar = (TypeFraction*) malloc(sizePar*sizeof(TypeFraction));
			funcFH[n].col = (TypeColumn*) malloc((sizePar+1)*sizeof(TypeColumn));
			for(k=0; k<nchild; k++)
				curPar[k] = 0; 
			do {
				double minVal;
				int sizeVal = 0;
				for(k=0; k<nchild; k++)
					sizeVal += h[k].col[curPar[k]].sizeVal;
				funcFH[n].col[funcFH[n].sizePar].coeff = (TypeCoeff*) malloc((sizeVal+1)*sizeof(TypeCoeff));
				funcFH[n].col[funcFH[n].sizePar].boundVal = (double*) malloc(sizeVal*sizeof(double));
				funcFH[n].col[funcFH[n].sizePar].sizeVal = 0;
				for(k=0; k<nchild; k++)
					curVal[k] = 0;
				do {
					funcFH[n].col[funcFH[n].sizePar].coeff[funcFH[n].col[funcFH[n].sizePar].sizeVal].a = 0;
					funcFH[n].col[funcFH[n].sizePar].coeff[funcFH[n].col[funcFH[n].sizePar].sizeVal].b = 0;
					funcFH[n].col[funcFH[n].sizePar].coeff[funcFH[n].col[funcFH[n].sizePar].sizeVal].c = 0;
					funcFH[n].col[funcFH[n].sizePar].coeff[funcFH[n].col[funcFH[n].sizePar].sizeVal].d = 0;
					for(k=0; k<nchild; k++) {
						funcFH[n].col[funcFH[n].sizePar].coeff[funcFH[n].col[funcFH[n].sizePar].sizeVal].a += h[k].col[curPar[k]].coeff[curVal[k]].a;
						funcFH[n].col[funcFH[n].sizePar].coeff[funcFH[n].col[funcFH[n].sizePar].sizeVal].b += h[k].col[curPar[k]].coeff[curVal[k]].b;
						funcFH[n].col[funcFH[n].sizePar].coeff[funcFH[n].col[funcFH[n].sizePar].sizeVal].c += h[k].col[curPar[k]].coeff[curVal[k]].c;
						funcFH[n].col[funcFH[n].sizePar].coeff[funcFH[n].col[funcFH[n].sizePar].sizeVal].d += h[k].col[curPar[k]].coeff[curVal[k]].d;
					}
					minVal = INFTY;
					for(k=0; k<nchild; k++)
						if(curVal[k]<h[k].col[curPar[k]].sizeVal && h[k].col[curPar[k]].boundVal[curVal[k]]<minVal)
							minVal = h[k].col[curPar[k]].boundVal[curVal[k]];
					for(k=0; k<nchild; k++)
						if(curVal[k]<h[k].col[curPar[k]].sizeVal && h[k].col[curPar[k]].boundVal[curVal[k]]<=minVal)
							curVal[k]++;
					if(minVal<INFTY) {
						funcFH[n].col[funcFH[n].sizePar].boundVal[funcFH[n].col[funcFH[n].sizePar].sizeVal] = minVal;
						funcFH[n].col[funcFH[n].sizePar].sizeVal++;
					}
				} while(minVal<INFTY);
				funcFH[n].col[funcFH[n].sizePar].coeff = (TypeCoeff*) realloc((void*)funcFH[n].col[funcFH[n].sizePar].coeff, (funcFH[n].col[funcFH[n].sizePar].sizeVal+1)*sizeof(TypeCoeff));
				if(funcFH[n].col[funcFH[n].sizePar].sizeVal>0)
					funcFH[n].col[funcFH[n].sizePar].boundVal = (double*) realloc((void*)funcFH[n].col[funcFH[n].sizePar].boundVal, funcFH[n].col[funcFH[n].sizePar].sizeVal*sizeof(double));
				else {
					free((void*)funcFH[n].col[funcFH[n].sizePar].boundVal);
					funcFH[n].col[funcFH[n].sizePar].boundVal = NULL;
				}
				minPar = fract(1,0);
				for(k=0; k<nchild; k++)
					if(curPar[k]<h[k].sizePar && cmpFract(h[k].boundPar[curPar[k]],minPar)<0)
						minPar = h[k].boundPar[curPar[k]];
				if(cmpFract(minPar, fract(1,0)) < 0) {
					funcFH[n].boundPar[funcFH[n].sizePar] = minPar;
					funcFH[n].sizePar++;
					for(k=0; k<nchild; k++)
						if(curPar[k]<h[k].sizePar && cmpFract(h[k].boundPar[curPar[k]],minPar) <= 0)
							curPar[k]++;
				}
			} while(cmpFract(minPar, fract(1,0)) < 0);
			if(funcFH[n].sizePar>0)
				funcFH[n].boundPar = (TypeFraction*) realloc((void*)funcFH[n].boundPar, funcFH[n].sizePar*sizeof(TypeFraction));
			else {
				free((void*) funcFH[n].boundPar);
				funcFH[n].boundPar = NULL;
			}
			funcFH[n].col = (TypeColumn*) realloc((void*)funcFH[n].col, (funcFH[n].sizePar+1)*sizeof(TypeColumn));
			free((void*) curPar);
			free((void*) curVal);			
		} else {
			int *curPar;
			TypeFraction minPar;
			funcFH[n].sizePar = 0;
			funcFH[n].boundPar = (TypeFraction*) malloc(sizePar*sizeof(TypeFraction));
			funcFH[n].col = (TypeColumn*) malloc((sizePar+1)*sizeof(TypeColumn));
			curPar = (int*) malloc(nchild*sizeof(int));
			for(k=0; k<nchild; k++)
				curPar[k] = 0;
			do {
				funcFH[n].col[funcFH[n].sizePar].coeff = (TypeCoeff*) malloc(sizeof(TypeCoeff));
				funcFH[n].col[funcFH[n].sizePar].sizeVal = 1;
				funcFH[n].col[funcFH[n].sizePar].boundVal = (double*) malloc(sizeof(double));
				funcFH[n].col[funcFH[n].sizePar].boundVal[0] = ((double*)tree->info)[n];
				funcFH[n].col[funcFH[n].sizePar].coeff[0].a = 0;
				funcFH[n].col[funcFH[n].sizePar].coeff[0].b = 0;
				funcFH[n].col[funcFH[n].sizePar].coeff[0].c = 0.;
				funcFH[n].col[funcFH[n].sizePar].coeff[0].d = 0.;			
				for(k=0; k<nchild; k++) {
					int j;
					for(j=0; j<h[k].col[curPar[k]].sizeVal && ((double*)tree->info)[n]>h[k].col[curPar[k]].boundVal[j]; j++);
					funcFH[n].col[funcFH[n].sizePar].coeff[0].c += h[k].col[curPar[k]].coeff[j].c-h[k].col[curPar[k]].coeff[j].a*((double*)tree->info)[n];
					funcFH[n].col[funcFH[n].sizePar].coeff[0].d += h[k].col[curPar[k]].coeff[j].d+h[k].col[curPar[k]].coeff[j].b*((double*)tree->info)[n];
				}
				minPar = fract(1,0);
				for(k=0; k<nchild; k++)
					if(curPar[k]<h[k].sizePar && cmpFract(h[k].boundPar[curPar[k]],minPar)<0)
						minPar = h[k].boundPar[curPar[k]];
				if(cmpFract(minPar, fract(1,0)) < 0) {
					funcFH[n].boundPar[funcFH[n].sizePar] = minPar;
					funcFH[n].sizePar++;
					for(k=0; k<nchild; k++)
						if(curPar[k]<h[k].sizePar && cmpFract(h[k].boundPar[curPar[k]],minPar) <= 0)
							curPar[k]++;
				}
			} while(cmpFract(minPar, fract(1,0)) < 0);
			free((void*)curPar);
			funcFH[n].col = (TypeColumn*) realloc((void*)funcFH[n].col, (funcFH[n].sizePar+1)*sizeof(TypeColumn));
		}
		free((void*) h);
	} else {
		if(isUnknown(n, tree)) {
			funcFH[n].sizePar = 0;
			funcFH[n].boundPar = NULL;
			funcFH[n].col = (TypeColumn*) malloc(sizeof(TypeColumn));
			funcFH[n].col[0].sizeVal = 0;
			funcFH[n].col[0].boundVal = NULL;
			funcFH[n].col[0].coeff = (TypeCoeff*) malloc(sizeof(TypeCoeff));
			funcFH[n].col[0].coeff[0].a = 0;
			funcFH[n].col[0].coeff[0].b = 0;
			funcFH[n].col[0].coeff[0].c = 0.;
			funcFH[n].col[0].coeff[0].d = 0.;			
		} else {
			funcFH[n].sizePar = 0;
			funcFH[n].boundPar = NULL;
			funcFH[n].col = (TypeColumn*) malloc(sizeof(TypeColumn));
			funcFH[n].col[0].sizeVal = 0;
			funcFH[n].col[0].boundVal = (double*) malloc(sizeof(double));
			funcFH[n].col[0].boundVal[0] = ((double*)tree->info)[n];
			funcFH[n].col[0].coeff = (TypeCoeff*) malloc(sizeof(TypeCoeff));
			funcFH[n].col[0].coeff[0].a = 0;
			funcFH[n].col[0].coeff[0].b = 0;
			funcFH[n].col[0].coeff[0].c = 0.;
			funcFH[n].col[0].coeff[0].d = 0.;			
		}
	}
}

void fillParametricContinuousParsimony(TypeTree *tree, TypeFunctTime phi) {
	TypeF *funcFH;
	int n;
	funcFH = (TypeF*) malloc(tree->size*sizeof(TypeF));
	fillInfoParametric(tree->root, funcFH, tree, phi);
	for(n=0; n<tree->size; n++)
		freeFuncF(funcFH[n]);
	free((void*)funcFH);
}

TypeResult *getResult(TypeTree *tree, TypeFunctTime phi) {
	TypeF *funcFH;
	TypeResult *res;
	int n;
	funcFH = (TypeF*) malloc(tree->size*sizeof(TypeF));
	fillInfoParametric(tree->root, funcFH, tree, phi);
	res = fillResult(funcFH, tree, phi);
	for(n=0; n<tree->size; n++)
		freeFuncF(funcFH[n]);
	free((void*)funcFH);
	return res;
}

void fprintFuncGnuplot(FILE *f, TypeF *func) {
	int p, v;
	fprintf(f, "splot ");
	if(func->sizePar>0) {
		fprintf(f, "(x<(%lf/%lf))?(", func->boundPar[0].num, func->boundPar[0].den);
		if(func->col[0].sizeVal>0) {
			fprintf(f, "(y<(%lf))?-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].boundVal[0], func->col[0].coeff[0].a, func->col[0].coeff[0].b, func->col[0].coeff[0].c, func->col[0].coeff[0].d);
			for(v=1; v<func->col[0].sizeVal; v++)
				fprintf(f, ":((y<(%lf))?-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].boundVal[v], func->col[0].coeff[v].a, func->col[0].coeff[v].b, func->col[0].coeff[v].c, func->col[0].coeff[v].d);
			fprintf(f, ":-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].coeff[func->col[0].sizeVal].a, func->col[0].coeff[func->col[0].sizeVal].b, func->col[0].coeff[func->col[0].sizeVal].c, func->col[0].coeff[func->col[0].sizeVal].d);
			for(v=1; v<func->col[0].sizeVal; v++)
				fprintf(f, ")");
		} else {
			fprintf(f, "-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].coeff[0].a, func->col[0].coeff[0].b, func->col[0].coeff[0].c, func->col[0].coeff[0].d);
		}
		fprintf(f, ")");
		for(p=1; p<func->sizePar; p++) {
			fprintf(f, ":((x<(%lf/%lf))?", func->boundPar[p].num, func->boundPar[p].den);
			if(func->col[0].sizeVal>0) {
				fprintf(f, "(y<(%lf))?-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].boundVal[0], func->col[p].coeff[0].a, func->col[p].coeff[0].b, func->col[p].coeff[0].c, func->col[p].coeff[0].d);
				for(v=1; v<func->col[p].sizeVal; v++)
					fprintf(f, ":((y<(%lf))?-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].boundVal[v], func->col[p].coeff[v].a, func->col[p].coeff[v].b, func->col[p].coeff[v].c, func->col[p].coeff[v].d);
				fprintf(f, ":-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].coeff[func->col[p].sizeVal].a, func->col[p].coeff[func->col[p].sizeVal].b, func->col[p].coeff[func->col[p].sizeVal].c, func->col[p].coeff[func->col[p].sizeVal].d);
				for(v=1; v<func->col[p].sizeVal; v++)
					fprintf(f, ")");
			} else {
				fprintf(f, "-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].coeff[0].a, func->col[p].coeff[0].b, func->col[p].coeff[0].c, func->col[p].coeff[0].d);
			}
		}
		fprintf(f, ":(");
		if(func->col[0].sizeVal>0) {
			fprintf(f, "(y<(%lf))?-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].boundVal[0], func->col[p].coeff[0].a, func->col[p].coeff[0].b, func->col[p].coeff[0].c, func->col[p].coeff[0].d);
			for(v=1; v<func->col[p].sizeVal; v++)
				fprintf(f, ":((y<(%lf))?-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].boundVal[v], func->col[p].coeff[v].a, func->col[p].coeff[v].b, func->col[p].coeff[v].c, func->col[p].coeff[v].d);
			fprintf(f, ":-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].coeff[func->col[p].sizeVal].a, func->col[p].coeff[func->col[p].sizeVal].b, func->col[p].coeff[func->col[p].sizeVal].c, func->col[p].coeff[func->col[p].sizeVal].d);
			for(v=1; v<func->col[p].sizeVal; v++)
				fprintf(f, ")");
		} else {
			fprintf(f, "-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].coeff[0].a, func->col[p].coeff[0].b, func->col[p].coeff[0].c, func->col[p].coeff[0].d);
		}
		fprintf(f, ")");
		for(p=1; p<func->sizePar; p++)
			fprintf(f, ")");
	} else {
		if(func->col[0].sizeVal>0) {
			fprintf(f, "(y<(%lf))?-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].boundVal[0], func->col[0].coeff[0].a, func->col[0].coeff[0].b, func->col[0].coeff[0].c, func->col[0].coeff[0].d);
			for(v=1; v<func->col[0].sizeVal; v++)
				fprintf(f, ":((y<(%lf))?-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].boundVal[v], func->col[0].coeff[v].a, func->col[0].coeff[v].b, func->col[0].coeff[v].c, func->col[0].coeff[v].d);
			fprintf(f, ":-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].coeff[func->col[0].sizeVal].a, func->col[0].coeff[func->col[0].sizeVal].b, func->col[0].coeff[func->col[0].sizeVal].c, func->col[0].coeff[func->col[0].sizeVal].d);
			for(v=1; v<func->col[0].sizeVal; v++)
				fprintf(f, ")");
		} else {
			fprintf(f, "-%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].coeff[0].a, func->col[0].coeff[0].b, func->col[0].coeff[0].c, func->col[0].coeff[0].d);
		}
	}
}

void fprintfl(FILE *f, double a, double b) {
	if(a>0.0)
		fprintf(f, "%.1lf\\gamma", a);
	if(a>0.0 && b>0.0)
		fprintf(f, "+");
	if(b>0.0)
		fprintf(f, "%.1lf", b);
}

void fprintfrac(FILE *f, int n, int d) {
	if(d==0) {
		fprintf(f, " \\infty ");
	} else {
		if(d==1 || n ==0) {
			fprintf(f, " %d ", n);
		} else {
			fprintf(f, " \\frac{%d}{%d} ", n, d);
		}
	}
}

void fprintfs(FILE *f, TypeCoeff c) {
	if(c.a>0 || c.b>0) {
		if(c.a>0 && c.b>0)
			fprintf(f, "(");
		if(c.a>0) {
			fprintf(f, "-");
			if(c.a!=1)
				fprintf(f, "%.2lf", c.a);
			fprintf(f, "\\gamma");
		}
		if(c.a>0 && c.b>0)
			fprintf(f, "+");
		if(c.b>0)
			fprintf(f, "%.2lf", c.b);
		if(c.a>0 && c.b>0)
			fprintf(f, ")");
		fprintf(f, " x");
	}
	if(c.c>0. || c.d>0.) {
			fprintf(f, "+");
		if(c.c>0.0) {
			if(c.c!=1.)
				fprintf(f, "%.2lf", c.c);
			fprintf(f, "\\gamma");
		}
		if(c.c>0.0 && c.d>0.0)
			fprintf(f, "+");
		if(c.d!=0.0)
			fprintf(f, "%.2lf", c.d);
	}
}

void fprintFunc(FILE *f, TypeF *func) {
	int p, v;
	fprintf(f, "$$\\begin{array}{");
	for(p=0; p<func->sizePar; p++)
		fprintf(f, "cc");
	fprintf(f, "c");
	fprintf(f, "}\n");
	for(p=0; p<func->sizePar; p++) {
		fprintf(f, "& ");
		fprintfrac(f, func->boundPar[p].num, func->boundPar[p].den);
		fprintf(f, "& ");
	}
	fprintf(f, "\\\\ \n");
	fprintf(f, "\\begin{array}{rc}\n");
	fprintf(f, "& ");
	fprintfs(f, func->col[0].coeff[0]);
	fprintf(f, "\\\\\n ");
	for(v=0; v<func->col[0].sizeVal; v++) {
		fprintf(f, "%.2lf & \\raisedrule[0.15em]{0.5pt}\\\\\n", func->col[0].boundVal[v]);
	fprintf(f, "& ");
	fprintfs(f, func->col[0].coeff[v+1]);
	fprintf(f, "\\\\\n ");
	}
	fprintf(f, "\\end{array} \n");
	for(p=1; p<=func->sizePar; p++) {
		fprintf(f, "& | &");
		fprintf(f, "\\begin{array}{rc}\n");
		fprintf(f, "& ");
		fprintfs(f, func->col[p].coeff[0]);
		fprintf(f, "\\\\\n ");
		for(v=0; v<func->col[p].sizeVal; v++) {
			fprintf(f, "%.2lf & \\raisedrule[0.15em]{0.5pt}\\\\\n", func->col[p].boundVal[v]);
		fprintf(f, "& ");
		fprintfs(f, func->col[p].coeff[v+1]);
		fprintf(f, "\\\\\n ");
		}
		fprintf(f, "\\end{array} \n");
	}
	fprintf(f, "\\end{array}$$");
}

void fprintptDebug(FILE *f, TypeCoeff c, double x) {
	fprintf(f, "(%.2lf, %.2lf)", -c.a*x+c.c, c.b*x+c.d);
}
void fprintfsDebug(FILE *f, TypeCoeff c) {
	if(c.a>0 || c.b>0) {
		if(c.a>0 && c.b>0)
			fprintf(f, "(");
		if(c.a>0) {
			fprintf(f, "-");
			if(c.a!=1)
				fprintf(f, "%.2lf", c.a);
			fprintf(f, "g");
		}
		if(c.a>0 && c.b>0)
			fprintf(f, "+");
		if(c.b>0 && (c.b != 1 || c.a > 0))
			fprintf(f, "%.2lf", c.b);
		if(c.a>0 && c.b>0)
			fprintf(f, ")");
		fprintf(f, "x");
	}
	if(c.c!=0. || c.d!=0.) {
		if(c.c>0.0) {
			fprintf(f, "+");
			if(c.c!=1.)
				fprintf(f, "%.2lf", c.c);
			fprintf(f, "g");
		} else {
			if(c.c<0.0) {
				if(c.c!=-1.) {
					fprintf(f, "%.2lfg", c.c);
				} else {
					fprintf(f, "-g");
				}
			}
		}
		if(c.d>0.0) {
			fprintf(f, "+");
			fprintf(f, "%.2lf", c.d);
		} else {
			if(c.d<0.0) {
				fprintf(f, "%.2lf", c.d);
			}
		}
	}
}

/*void fprintfracDebug(FILE *f, double n, double d) {
	if(d==0) {
		fprintf(f, "inf");
	} else {
		if(d==1 || n ==0) {
			fprintf(f, "%.2lf", n);
		} else {
			fprintf(f, "%.2lf/%.2lf", n, d);
		}
	}
}*/

void fprintColDebug(FILE *f, TypeColumn *c) {
	int v;
	fprintfsDebug(f, c->coeff[0]);
	for(v=1; v<=c->sizeVal; v++) {
		fprintf(f, "\n%.2lf   ", c->boundVal[v-1]);
		fprintptDebug(f, c->coeff[v-1], c->boundVal[v-1]);
		fprintf(f,"\n");
		fprintfsDebug(f, c->coeff[v]);
	}
	fprintf(f,"\n");
}

void fprintFuncDebug(FILE *f, TypeF *func) {
	int p, v;
	fprintf(f, "%d ", 0);
	fprintfsDebug(f, func->col[0].coeff[0]);
	for(v=1; v<=func->col[0].sizeVal; v++) {
		fprintf(f, "%d ", v);
		fprintf(f, "\n%.2lf   ", func->col[0].boundVal[v-1]);
		fprintptDebug(f, func->col[0].coeff[v-1], func->col[0].boundVal[v-1]);
		fprintf(f,"\n");
		fprintfsDebug(f, func->col[0].coeff[v]);
	}
	for(p=1; p<=func->sizePar; p++) {
		fprintf(f, "\n*********");
		fprintfracDebug(f, func->boundPar[p-1].num, func->boundPar[p-1].den);
		fprintf(f, "*********\n");
		fprintfsDebug(f, func->col[p].coeff[0]);
		for(v=1; v<=func->col[p].sizeVal; v++) {
			fprintf(f, "\n%.2lf   ", func->col[p].boundVal[v-1]);
			fprintptDebug(f, func->col[p].coeff[v-1], func->col[p].boundVal[v-1]);
			fprintf(f,"\n");
			fprintfsDebug(f, func->col[p].coeff[v]);
		}
	}
}

void fprintStateTreeResult(FILE *f, TypeResult *res, TypeTree *tree) {
	int n;
	char **myName, **myComment, **tmpComment;
	
	myName = (char**) malloc(tree->size*sizeof(char*));
	myComment = (char**) malloc(tree->size*sizeof(char*));
	tmpComment = tree->comment;
	for(n=0; n<tree->size; n++) {
		char buffer[500000];
		if(isUnknown(n, tree)) {
			myName[n] = (char*) malloc((strlen(UNKNOWN_STRING)+1)*sizeof(char));
			strcpy(myName[n], UNKNOWN_STRING);
		} else {
			sprintf(buffer, "%.3lf", ((double*)tree->info)[n]);
			myName[n] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
			strcpy(myName[n], buffer);
		}
		sprintResult(buffer, &(res[n]));
		myComment[n] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
		strcpy(myComment[n], buffer);
	}
	tmpComment = tree->comment;
	tree->comment = myComment;
	fprintTree(f, tree, display_name);
	tree->comment = tmpComment;
	for(n=0; n<tree->size; n++) {
		free((void*)myName[n]);
		free((void*)myComment[n]);
	}
	free((void*)myName);
	free((void*)myComment);
}


void sprintResult(char *s, TypeResult *r) {
	int i;
	char *tmp;
	tmp = s;
	tmp += sprintf(tmp, "\\hbox{$\\scriptstyle\\{%.3lf", r->val[0]);
	for(i=1; i<=r->sizePar; i++)
		tmp += sprintf(tmp, "\\stackrel{\\frac{%lf}{%lf}}{,}%.3lf", r->boundPar[i-1].num, r->boundPar[i-1].den, r->val[i]);
	tmp += sprintf(tmp, "\\}$}");
} 

void sprintResultDebug(char *s,  TypeResult *r) {
	int i;
	char *tmp;
	tmp = s;
	tmp += sprintf(tmp, "%.3lf", r->val[0]);
	for(i=1; i<r->sizePar; i++)
			tmp += sprintf(tmp, " >%lf/%lf< %.3lf", r->boundPar[i-1].num, r->boundPar[i-1].den, r->val[i]);
}

void fprintResultDebug(FILE *f, TypeResult *r) {
	int i;
	fprintf(f, "%.3lf", r->val[0]);
	for(i=1; i<=r->sizePar; i++)
		fprintf(f, " >%lf/%lf< %.3lf", r->boundPar[i-1].num, r->boundPar[i-1].den, r->val[i]);
}

TypeFraction *mergeFractList(TypeFraction *l, int *size, TypeResult *res) {
	int i, j, ind;
	TypeFraction *new;
	new = (TypeFraction*) malloc((*size+res->sizePar)*sizeof(TypeFraction));
	i = 0; j = 0; ind = 0;
	while(i<*size || j<res->sizePar) {
		TypeFraction min = fract(1,0);
		if(i<*size && cmpFract(l[i], min)<=0)
			min = l[i];
		if(j<res->sizePar && cmpFract(res->boundPar[j], min)<=0)
			min = res->boundPar[j];
		if(i<*size && cmpFract(min, l[i])<=0)
			i++;
		if(j<res->sizePar && cmpFract(min, res->boundPar[j])<=0)
			j++;
		new[ind++] = min;
	}
	*size = ind;
	new = (TypeFraction*) realloc((void*)new, ind*sizeof(TypeFraction));
	return new;
}



	
/*print result as flat list*/	
void fprintResultList(FILE *f, TypeResult *res, TypeTree *tree) {
	int n;
	if(tree->size<=0)
		return;
	for(n=0; n<tree->size; n++) {
		fprintIdentTimeComment(f, n, tree, display_name);
		fprintf(f, "\t");
		fprintResultDebug(f, &(res[n]));
		fprintf(f, "\n");
	}
}

TypeListDouble *getErrors(TypeResult *res, TypeIndexStateList *sta, TypeTree *tree) {
	int i, k, size;
	TypeFraction *tot;
	TypeListDouble *err;
	double min;
	err = (TypeListDouble*) malloc(sizeof(TypeListDouble));
	tot = (TypeFraction*) malloc(res[sta->index[0]].sizePar*sizeof(TypeFraction));
	for(size=0; size<res[sta->index[0]].sizePar; size++)
		tot[size] = res[sta->index[0]].boundPar[size];
	for(k=1; k<sta->size; k++) {
		TypeFraction *tmp;
		tmp = mergeFractList(tot, &size, &(res[sta->index[k]]));
		free((void*) tot);
		tot = tmp;
	}
	err->size = size+1;
	err->val = (double*) malloc(err->size*sizeof(double));

	for(i=0; i<size; i++)
		err->val[i] = 0.;
	for(k=0; k<sta->size; k++) {
		int i, j;
		j = 0;
		for(i=0; i<=size; i++) {
			err->val[i] += fabs(sta->state[k]-res[sta->index[k]].val[j]);
			if(cmpFract(res[sta->index[k]].boundPar[j], tot[i]) <= 0)
				j++;
		}
	}
	min = err->val[0];
	for(i=1; i<=size; i++)
		if(err->val[i]<min) {
			min = err->val[i];
		}
	return err;
}
		
TypeFractionInter getBestInter(TypeResult *res, TypeIndexStateList *sta, TypeTree *tree) {
	int i, k, imin, size;
	TypeFraction *tot;
	TypeListDouble *err;
	TypeFractionInter inter;
	double min;
	err = (TypeListDouble*) malloc(sizeof(TypeListDouble));
	tot = (TypeFraction*) malloc(res[sta->index[0]].sizePar*sizeof(TypeFraction));
	for(size=0; size<res[sta->index[0]].sizePar; size++)
		tot[size] = res[sta->index[0]].boundPar[size];
	for(k=1; k<sta->size; k++) {
		TypeFraction *tmp;
		tmp = mergeFractList(tot, &size, &(res[sta->index[k]]));
		free((void*) tot);
		tot = tmp;
	}
	err->size = size+1;
	err->val = (double*) malloc(err->size*sizeof(double));

	for(i=0; i<size; i++)
		err->val[i] = 0.;
	for(k=0; k<sta->size; k++) {
		int i, j;
		j = 0;
		for(i=0; i<=size; i++) {
			err->val[i] += fabs(sta->state[k]-res[sta->index[k]].val[j]);
			if(cmpFract(res[sta->index[k]].boundPar[j], tot[i]) <= 0)
				j++;
		}
	}
	min = err->val[0]; imin = 0;
	for(i=1; i<=size; i++)
		if(err->val[i]<min) {
			min = err->val[i];
			imin = i;
		}
	if(imin > 0)
		inter.inf = tot[imin-1];
	else
		inter.inf = fract(0,1);
	if(imin < size)
		inter.sup = tot[imin];
	else
		inter.sup = fract(1,0);
	free((void*) tot);
	return inter;
}


