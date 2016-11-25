#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <locale.h>
#include "Utils.h"
#include "DrawTreeResultGeneric.h"

#define NB_MIN 4
#define OFFSET 10
#define TICK_LENGTH 10
#define FONT_NAME "Helvetica"
#define FONT_SIZE 9.
#define LABEL_SEP 10.
#define MAX_STRING_SIZE 500
#define STANDARD_WIDTH 1000
#define CHAR_WIDTH 10
#define HALF_PI asin(1)

static double drawNodeGeneric(int n, int parent, TypeTree *tree, TypeResult *r,  TypeInfoDrawTreeGeneric *param);
static void drawTreeGeneric(TypeTree *tree, TypeResult *r, TypeInfoDrawTreeGeneric *param);
static double convex(double s, double e, double v);

double convex(double s, double e, double v) {
	return (1.-v)*s+v*e;
}

TypeRGB getRGB(TypeRGB start, TypeRGB end, double v) {
	TypeRGB res = {.red = convex(start.red, end.red, v), .green = convex(start.green, end.green, v), .blue = convex(start.blue, end.blue, v)};
	return res;
}

void drawResultPie(double x, double y, TypeResult *r, TypeInfoDrawTreeGeneric *info) {
	x++;
	if(r->sizePar == 0) {
		info->funct.fillWedge(getRGB(info->param.start, info->param.end, (r->val[0]-info->param.vmin)*info->param.vscale), x, y, 0, HALF_PI, &(info->param));
	} else {
		int i;
		info->funct.fillWedge(getRGB(info->param.start, info->param.end, (r->val[0]-info->param.vmin)*info->param.vscale), x, y, 0, HALF_PI-(asin(((double)r->boundPar[0].den)/(sqrt(pow((double)r->boundPar[0].num, 2.)+pow((double)r->boundPar[0].den, 2.))))), &(info->param));
		for(i=0; i<r->sizePar-1; i++)
			info->funct.fillWedge(getRGB(info->param.start, info->param.end, (r->val[i+1]-info->param.vmin)*info->param.vscale), x, y, HALF_PI-(asin(((double)r->boundPar[i].den)/(sqrt(pow((double)r->boundPar[i].num, 2.)+pow((double)r->boundPar[i].den, 2.))))), HALF_PI-(asin(((double)r->boundPar[i+1].den)/(sqrt(pow((double)r->boundPar[i+1].num, 2.)+pow((double)r->boundPar[i+1].den, 2.))))), &(info->param));
		info->funct.fillWedge(getRGB(info->param.start, info->param.end, (r->val[r->sizePar]-info->param.vmin)*info->param.vscale), x, y, HALF_PI-(asin(((double)r->boundPar[r->sizePar-1].den)/(sqrt(pow((double)r->boundPar[r->sizePar-1].num, 2.)+pow((double)r->boundPar[r->sizePar-1].den, 2.))))), HALF_PI, &(info->param));
	}
	info->funct.drawWedge(x, y, 0, HALF_PI, &(info->param));
}

void drawResultPieSpecial(double x, double y, TypeResult *r, TypeInfoDrawTreeGeneric *info) {
	int i;
	char buffer[MAX_STRING_SIZE];
	double a;
	x++;
	if(r->sizePar == 0) {
		info->funct.fillWedge(getRGB(info->param.start, info->param.end, (r->val[0]-info->param.vmin)*info->param.vscale), x, y, 0, HALF_PI, &(info->param));
	} else {
		int i;
		info->funct.fillWedge(getRGB(info->param.start, info->param.end, (r->val[0]-info->param.vmin)*info->param.vscale), x, y, 0, HALF_PI-(asin(((double)r->boundPar[0].den)/(sqrt(pow((double)r->boundPar[0].num, 2.)+pow((double)r->boundPar[0].den, 2.))))), &(info->param));
		for(i=0; i<r->sizePar-1; i++)
			info->funct.fillWedge(getRGB(info->param.start, info->param.end, (r->val[i+1]-info->param.vmin)*info->param.vscale), x, y, HALF_PI-(asin(((double)r->boundPar[i].den)/(sqrt(pow((double)r->boundPar[i].num, 2.)+pow((double)r->boundPar[i].den, 2.))))), HALF_PI-(asin(((double)r->boundPar[i+1].den)/(sqrt(pow((double)r->boundPar[i+1].num, 2.)+pow((double)r->boundPar[i+1].den, 2.))))), &(info->param));
		info->funct.fillWedge(getRGB(info->param.start, info->param.end, (r->val[r->sizePar]-info->param.vmin)*info->param.vscale), x, y, HALF_PI-(asin(((double)r->boundPar[r->sizePar-1].den)/(sqrt(pow((double)r->boundPar[r->sizePar-1].num, 2.)+pow((double)r->boundPar[r->sizePar-1].den, 2.))))), HALF_PI, &(info->param));
	}
	if(r->sizePar>0) {
		a = -(HALF_PI-(asin(((double)r->boundPar[0].den)/(sqrt(pow((double)r->boundPar[0].num, 2.)+pow((double)r->boundPar[0].den, 2.))))))/2.;
		sprintf(buffer, "%.2lf", r->val[0]);
		info->funct.drawTextAngle(x+(info->param.radius+info->param.tickLength)*cos(a),y+(info->param.radius+info->param.tickLength)*sin(a), a, buffer, "l", &(info->param));
		a = -HALF_PI+(asin(((double)r->boundPar[0].den)/(sqrt(pow((double)r->boundPar[0].num, 2.)+pow((double)r->boundPar[0].den, 2.)))));
		info->funct.drawDottedLine(x, y, x+info->param.radius*cos(a),y+info->param.radius*sin(a), &(info->param));
		for(i=1; i<r->sizePar; i++) {
			a = -(HALF_PI-(asin(((double)r->boundPar[i-1].den)/(sqrt(pow((double)r->boundPar[i-1].num, 2.)+pow((double)r->boundPar[i-1].den, 2.)))))+HALF_PI-(asin(((double)r->boundPar[i].den)/(sqrt(pow((double)r->boundPar[i].num, 2.)+pow((double)r->boundPar[i].den, 2.))))))/2.;
			sprintf(buffer, "%.2lf", r->val[i]);
			info->funct.drawTextAngle(x+(info->param.radius+info->param.tickLength)*cos(a),y+(info->param.radius+info->param.tickLength)*sin(a), a, buffer, "l", &(info->param));
			a = -HALF_PI+(asin(((double)r->boundPar[i].den)/(sqrt(pow((double)r->boundPar[i].num, 2.)+pow((double)r->boundPar[i].den, 2.)))));
			info->funct.drawDottedLine(x, y, x+info->param.radius*cos(a),y+info->param.radius*sin(a), &(info->param));
		}
		a = -(2*HALF_PI-(asin(((double)r->boundPar[r->sizePar-1].den)/(sqrt(pow((double)r->boundPar[r->sizePar-1].num, 2.)+pow((double)r->boundPar[r->sizePar-1].den, 2.))))))/2.;
		sprintf(buffer, "%.2lf", r->val[r->sizePar]);
		info->funct.drawTextAngle(x+(info->param.radius+info->param.tickLength)*cos(a),y+(info->param.radius+info->param.tickLength)*sin(a), a, buffer, "l", &(info->param));
	} else {
		a = -HALF_PI/2.;
		sprintf(buffer, "%.2lf", r->val[0]);
		info->funct.drawTextAngle(x+(info->param.radius+info->param.tickLength)*cos(a),y+(info->param.radius+info->param.tickLength)*sin(a), a, buffer, "l", &(info->param));
	}
	info->funct.drawWedge(x, y, 0, HALF_PI, &(info->param));
}

#define GRAD_SEP 8
void drawScaleVPst(double x, double y, double size, TypeInfoDrawTreeGeneric *info) {
    int flag = 0;
    double start, end, step, cur;
    info->funct.drawLine(x, y+GRAD_SEP, x+size, y+GRAD_SEP, &(info->param));
    if((info->param.vmax-info->param.vmin) <= 0.)
		return;
    step = pow(10., floor(log10(info->param.vmax-info->param.vmin)));
    if((info->param.vmax-info->param.vmin)/step<NB_MIN)
		step /= 2.;
    start = step*ceil(info->param.vmin/step);
    end = step*floor(info->param.vmax/step);
	flag = step<1.;
	for(cur = start; cur<=end; cur += step) {
		char tmp[500];
        info->funct.drawLine(x+(cur-info->param.vmin)*info->param.vscale*size, y+GRAD_SEP, x+(cur-info->param.vmin)*info->param.vscale*size, y+GRAD_SEP+info->param.tickLength, &(info->param));
		if(flag)
			sprintf(tmp, "%.1lf", cur);
		else
			sprintf(tmp, "%.0lf", cur);
        info->funct.drawText(x+(cur-info->param.vmin)*info->param.vscale*size, y+GRAD_SEP+info->param.tickLength, tmp, "t", &(info->param));
	}
	info->funct.fillGradient(info->param.start, info->param.end, x, y-1, x+size, y+GRAD_SEP-2, &(info->param));
}

void drawTreeResultFileGeneric(char *filename, TypeTree *tree, TypeResult *r, TypeInfoDrawTreeGeneric *info) {
	int i;
	double *timeSave;
	
	timeSave = tree->time;
	tree->time = (double*) malloc(tree->size*sizeof(double));
	if(tree->time != NULL) {
		int hasUnknown = 0;
		for(i=0; i<tree->size; i++) {
			tree->time[i] = timeSave[i];
			if(timeSave[i] == NO_TIME)
				hasUnknown = 1;
		}
		if(hasUnknown)
			fillUnknownTimes(info->param.tmin, info->param.tmax,  tree);
	} else {
		for(i=0; i<tree->size; i++)
			tree->time[i] = 1.;
		bltoabsTime(tree);
	}
	info->param.tmin = tree->minTime;
    info->param.tmax = tree->maxTime;
    info->param.vmax = r[0].val[0];
	info->param.vmin = r[0].val[0];
    for(i=0; i<tree->size; i++) {
		int j;
		for(j=0; j<=r[i].sizePar; j++) {
			if(r[i].val[j]>info->param.vmax)
				info->param.vmax = r[i].val[j];
			if(r[i].val[j]<info->param.vmin)
				info->param.vmin = r[i].val[j];
		}
	}
	if(info->param.vmax==info->param.vmin)
		info->param.vmax++;
	info->param.vscale = 1/(info->param.vmax-info->param.vmin);
	info->funct.start(filename, tree, &(info->param));
	info->param.labelSep += info->param.radius;
	info->param.scale = (info->param.width-info->param.labelWidth-info->param.xoffset-info->param.labelSep)/((info->param.tmax-info->param.tmin));
	drawTreeGeneric(tree, r, info);
	drawScaleGeneric(info->param.xoffset, info->param.leafSep*(countLeaves(tree)+1), info);
	drawScaleVPst(info->param.width/4, info->param.leafSep*(countLeaves(tree)+2.1), info->param.width/2, info);
//	drawScaleVPst(info->param.xoffset+info->param.width-info->param.labelWidth-info->param.labelSep, info->param.leafSep*(countLeaves(tree)+1), info->param.labelWidth, info);
	info->funct.end(&(info->param));
	free((void*) tree->time);
	tree->time = timeSave;
}

void drawSpecialResultFileGeneric(char *filename, TypeResult *r, TypeInfoDrawTreeGeneric *info) {
	info->funct.startStd(filename, 300, 300, &(info->param));
	info->param.radius = 250;
	info->param.tickLength = 5;
	drawResultPieSpecial(10, 290, r, info);
	info->funct.end(&(info->param));
}

void drawTreeGeneric(TypeTree *tree, TypeResult *r, TypeInfoDrawTreeGeneric *info) {
	int tmp;
	double min, max, y;
	if(tree->size<=0)
		return;
	if((tmp = tree->node[tree->root].child) >= 0) {
		min = drawNodeGeneric(tmp, tree->root, tree, r, info);
		max = min;
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling) {
			max = drawNodeGeneric(tmp, tree->root, tree, r, info);
		}
	} else {
		max = info->param.leafCur+info->param.leafSep/2.;
		min = max;
	}
	info->funct.drawLine((tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset, min, (tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset, max, &(info->param));
	y = (min+max)/2;
	info->funct.drawLine((tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, info->param.xoffset, info->param.yoffset+y, &(info->param));
	if(tree->name != NULL && tree->name[tree->root] != NULL)
		info->funct.drawText((tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, y+info->param.yoffset, tree->name[tree->root], "l", &(info->param));
	drawResultPie((tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset+0.5, y+info->param.yoffset+info->param.roffset, &(r[tree->root]), info);
}

double drawNodeGeneric(int n, int parent, TypeTree *tree, TypeResult *r, TypeInfoDrawTreeGeneric *info) {
	double min, max, y;
	if(tree->node[n].child >= 0) {
		int tmp = tree->node[n].child;
		min = drawNodeGeneric(tmp, n, tree, r, info);
		max = min;
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling)
			max = drawNodeGeneric(tmp, n, tree, r, info);
		y = (min+max)/2;
		info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+min, (tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+max, &(info->param));
		if(tree->name && tree->name[n])
			info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, y+info->param.yoffset, tree->name[n], "l", &(info->param));
			drawResultPie((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+0.5, y+info->param.yoffset+info->param.roffset, &(r[n]), info);
	} else {
		info->param.leafCur += info->param.leafSep;
		y = info->param.leafCur;
		if(tree->name && tree->name[n])
			info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, info->param.yoffset+y, tree->name[n], "l", &(info->param));
		drawResultPie((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+0.5, y+info->param.yoffset+info->param.roffset, &(r[n]), info);
	}
	info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, (tree->time[parent]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, &(info->param));
	return y;
}
















