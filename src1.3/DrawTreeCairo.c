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
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <locale.h>
#include <cairo.h>
#include <cairo-ps.h>
#include <cairo-pdf.h>
#include <cairo-svg.h>
#include "Utils.h"
#include "DrawTreeCairo.h"

#define CAIRO_FONT_SIZE 50
#define CAIRO_LINE_WIDTH 1.
#define CAIRO_RADIUS 10
#define CAIRO_OFFSET 10
#define TICKLENGTH 3
#define LABELSEP 5
#define LEAFSEP 15
#define HALF_PI asin(1)

typedef struct CAIRO_INFO {
	cairo_surface_t *surface;
	cairo_t *cr;
	char *filename;
} TypeCairoInfo;

static double getMaxLeafLabelWidthCairo(char **name, int size, cairo_t *cr);
static void setFunctCairo(TypeFunctDrawTreeGeneric *funct);
static void setParamCairo(TypeParamDrawTreeGeneric *param);

double getMaxLeafLabelWidthCairo(char **name, int size, cairo_t *cr) {
	int n;
	double max = 0.;
	if(name == NULL)
		return 0.;
	for(n=0; n<size; n++)
		if(name[n] != NULL) {
			cairo_text_extents_t extentsT;
			cairo_text_extents(cr, name[n], &extentsT);
			if(extentsT.width>max)
				max = extentsT.width;
		}
	return max;
}


void drawTextCairo(double x, double y, char *text, char *al, TypeParamDrawTreeGeneric *param) {
	cairo_t *cr = ((TypeCairoInfo*)param->info)->cr;
	cairo_save(cr);
	cairo_set_source_rgba (cr, 0., 0., 0., 1);
	cairo_font_extents_t extentsF;
	cairo_text_extents_t extentsT;
	double xs, ys;
	if(text == NULL)
		return;
	cairo_font_extents (cr, &extentsF);
	cairo_text_extents (cr, text, &extentsT);
	if(al != NULL && strlen(al)>0) {
		if(strlen(al) == 1) {
			switch(al[0]) {
				case 'B':
					xs = x-extentsT.width/2.;
					ys = y+extentsF.descent;
					break;
				case 'b':
					xs = x-extentsT.width/2.;
					ys = y;
					break;
				case 't':
					xs = x-extentsT.width/2.;
					ys = y+extentsF.descent+extentsF.ascent;
					break;
				case 'l':
					xs = x;
					ys = y+extentsF.ascent/3.;
					break;
				case 'r':
					xs = x-extentsT.width;
					ys = y+extentsF.ascent/3.;
					break;
				case 'c':
				default:
					xs = x-extentsT.width/2.;
					ys = y+extentsF.ascent/3.;
			}
		} else {
			switch(al[0]) {
				case 'B':
					ys = y+extentsF.descent;
					break;
				case 'b':
					ys = y;
					break;
				case 't':
					ys = y+extentsF.descent+extentsF.ascent;
					break;
				case 'c':
				default:
					ys = y+extentsF.ascent/3.;
			}
			switch(al[1]) {
				case 'l':
					xs = x;
					break;
				case 'r':
					xs = x-extentsT.width;
					break;
				case 'c':
				default:
					xs = x-extentsT.width/2.;
			}
		}
	} else {
		xs = x-extentsT.width/2.;
		ys = y+extentsF.ascent/3.;
	}
	cairo_move_to (cr, xs, ys);
	cairo_show_text (cr, text);
	cairo_restore(cr);
}

void drawTextAngleCairo(double x, double y, double a, char *text, char *al, TypeParamDrawTreeGeneric *param) {
	cairo_t *cr = ((TypeCairoInfo*)param->info)->cr;
	cairo_save(cr);
	cairo_set_source_rgba (cr, 0., 0., 0., 1);
	cairo_font_extents_t extentsF;
	cairo_text_extents_t extentsT;
	double xs, ys;
	if(text == NULL)
		return;
	cairo_font_extents (cr, &extentsF);
	cairo_text_extents (cr, text, &extentsT);
	cairo_set_font_size (cr, 6);
	if(al != NULL && strlen(al)>0) {
		if(strlen(al) == 1) {
			switch(al[0]) {
				case 'B':
					xs = x-extentsT.width/2.;
					ys = y+extentsF.descent;
					break;
				case 'b':
					xs = x-extentsT.width/2.;
					ys = y;
					break;
				case 't':
					xs = x-extentsT.width/2.;
					ys = y+extentsF.descent+extentsF.ascent;
					break;
				case 'l':
					xs = x;
					ys = y+extentsF.ascent/3.;
					break;
				case 'r':
					xs = x-extentsT.width;
					ys = y+extentsF.ascent/3.;
					break;
				case 'c':
				default:
					xs = x-extentsT.width/2.;
					ys = y+extentsF.ascent/3.;
			}
		} else {
			switch(al[0]) {
				case 'B':
					ys = y+extentsF.descent;
					break;
				case 'b':
					ys = y;
					break;
				case 't':
					ys = y+extentsF.descent+extentsF.ascent;
					break;
				case 'c':
				default:
					ys = y+extentsF.ascent/3.;
			}
			switch(al[1]) {
				case 'l':
					xs = x;
					break;
				case 'r':
					xs = x-extentsT.width;
					break;
				case 'c':
				default:
					xs = x-extentsT.width/2.;
			}
		}
	} else {
		xs = x-extentsT.width/2.;
		ys = y+extentsF.ascent/3.;
	}
	cairo_move_to(cr, xs, ys);
	cairo_rotate(cr, a);
	cairo_show_text(cr, text);
	cairo_restore(cr);
}

void drawLineCairo(double x1, double y1, double x2, double y2, TypeParamDrawTreeGeneric *param) {
	cairo_t *cr = ((TypeCairoInfo*)param->info)->cr;
	cairo_save(cr);
	cairo_set_source_rgba (cr, 0., 0., 0., 1);
	cairo_set_line_width (cr, CAIRO_LINE_WIDTH);
	cairo_move_to (cr, x1, y1);
	cairo_line_to (cr, x2, y2);
	cairo_stroke (cr);
	cairo_restore(cr);
}

void drawDottedLineCairo(double x1, double y1, double x2, double y2, TypeParamDrawTreeGeneric *param) {
	double dashes[] = {1.0,  /* ink */2.0,  /* skip */};
	int ndash  = sizeof (dashes)/sizeof(dashes[0]);
	cairo_t *cr = ((TypeCairoInfo*)param->info)->cr;
	cairo_save(cr);
	cairo_set_source_rgba (cr, 0., 0., 0., 1);
	cairo_set_dash (cr, dashes, ndash, -dashes[0]);
	cairo_set_line_width (cr, CAIRO_LINE_WIDTH/2);
	cairo_move_to (cr, x1, y1);
	cairo_line_to (cr, x2, y2);
	cairo_stroke (cr);
	cairo_restore(cr);
}


void fillWedgeCairo(TypeRGB rgb, double x, double y, double a, double b, TypeParamDrawTreeGeneric *param) {
	cairo_t *cr = ((TypeCairoInfo*)param->info)->cr;
	cairo_save(cr);
	cairo_set_source_rgb (cr, rgb.red, rgb.green, rgb.blue);
	cairo_set_line_width (cr, CAIRO_LINE_WIDTH);
	cairo_move_to (cr, x, y);
	cairo_line_to (cr, x+param->radius*cos(-b), y+param->radius*sin(-b));
	cairo_arc(cr, x, y, param->radius, -b, -a);
	cairo_line_to (cr, x, y);
	cairo_fill (cr);
	cairo_restore(cr);
}

void drawWedgeCairo(double x, double y, double a, double b, TypeParamDrawTreeGeneric *param) {
	cairo_t *cr = ((TypeCairoInfo*)param->info)->cr;
	cairo_save(cr);
	cairo_set_source_rgb (cr, 0, 0, 0);
	cairo_set_line_width (cr, CAIRO_LINE_WIDTH/2.);
	cairo_new_sub_path(cr);
	cairo_move_to (cr, x, y);
	cairo_line_to (cr, x+param->radius*cos(-b), y+param->radius*sin(-b));
	cairo_arc(cr, x, y, param->radius, -b, -a);
	cairo_line_to (cr, x, y);
	cairo_stroke (cr);
	cairo_restore(cr);
}

void fillGradientCairo(TypeRGB rgb0, TypeRGB rgb1, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param) {
	cairo_t *cr = ((TypeCairoInfo*)param->info)->cr;
	cairo_save(cr);
	cairo_pattern_t *linpat;
	linpat = cairo_pattern_create_linear(x0, 0, x1, 0);
	cairo_pattern_add_color_stop_rgb(linpat, 0.,  rgb0.red,  rgb0.green,  rgb0.blue);
	cairo_pattern_add_color_stop_rgb(linpat, 1.,  rgb1.red,  rgb1.green,  rgb1.blue);
	cairo_new_sub_path(cr);
	cairo_rectangle(cr, x0, y0, x1-x0, y1-y0);
	cairo_set_source(cr, linpat);
	cairo_fill(cr);
	cairo_pattern_destroy(linpat);
	cairo_restore(cr);
}

void setParamCairo(TypeParamDrawTreeGeneric *param) {
	param->xoffset = CAIRO_OFFSET;
	param->yoffset = 0;
	param->leafSep = LEAFSEP;
	param->labelSep = LABELSEP;
	param->ycenter = 3.;
	param->tickLength = TICKLENGTH;
	param->ydec = 2.;
	param->radius = CAIRO_RADIUS;
	param->roffset = CAIRO_RADIUS/2.;
	param->leafCur = 0.;
}

void setFunctCairo(TypeFunctDrawTreeGeneric *funct) {
	funct->drawText = drawTextCairo;
	funct->drawTextAngle = drawTextAngleCairo;
	funct->drawLine = drawLineCairo;
	funct->drawDottedLine = drawDottedLineCairo;
	funct->fillWedge = fillWedgeCairo;
	funct->drawWedge = drawWedgeCairo;
	funct->fillGradient = fillGradientCairo;
}
void setFunctSVG(TypeFunctDrawTreeGeneric *funct) {
	funct->startStd = startSVGStd;
	funct->start = startSVG;
	funct->end = endCairo;
	setFunctCairo(funct);
}

void startSVGStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param) {
	cairo_surface_t *surface;
	cairo_t *cr;
	cairo_font_extents_t extentsF;
	surface = cairo_svg_surface_create (filename, width, height);
	cr = cairo_create(surface);
	cairo_select_font_face (cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size (cr, 12);
	cairo_font_extents (cr, &extentsF);
	setParamCairo(param);
	param->height = height;
	param->width = width;
	param->info = (void*) malloc(sizeof(TypeCairoInfo));
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->surface = surface;
	((TypeCairoInfo*)param->info)->filename = filename;
	cairo_set_source_rgb (((TypeCairoInfo*)param->info)->cr, 0, 0, 0);
}	

void startSVG(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param) {
	cairo_surface_t *surface;
	cairo_t *cr;
	cairo_font_extents_t extentsF;
	surface = cairo_svg_surface_create (filename, 500, 500);
	cr = cairo_create(surface);
	cairo_select_font_face (cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size (cr, 12);
	cairo_font_extents (cr, &extentsF);
	param->height = LEAFSEP*(countLeaves(tree)+1)+extentsF.descent+extentsF.ascent+2*LABELSEP+TICKLENGTH;
	if(countLeaves(tree)>40)
		if(countLeaves(tree)>40)
		param->width = param->height/sqrt(2);
	else
		param->width = 2*param->height/sqrt(2);
	else
		param->width = 2*param->height/sqrt(2);
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
	surface = cairo_svg_surface_create (filename, param->width, param->height);
	cr = cairo_create (surface);
	setParamCairo(param);
	param->tmin = tree->minTime;
    param->tmax = tree->maxTime;
	param->info = (void*) malloc(sizeof(TypeCairoInfo));
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->surface = surface;
	((TypeCairoInfo*)param->info)->filename = filename;
	param->labelWidth = getMaxLeafLabelWidthCairo(tree->name, tree->size, ((TypeCairoInfo*)param->info)->cr);	
	param->scale = (param->width-param->labelWidth-param->xoffset-param->labelSep)/((param->tmax-param->tmin));
	cairo_set_source_rgb (((TypeCairoInfo*)param->info)->cr, 0, 0, 0);
}	
	
void endCairo(TypeParamDrawTreeGeneric *param) {
	cairo_surface_show_page(((TypeCairoInfo*)param->info)->surface);
	cairo_destroy(((TypeCairoInfo*)param->info)->cr);
	cairo_surface_finish(((TypeCairoInfo*)param->info)->surface);
	cairo_surface_destroy(((TypeCairoInfo*)param->info)->surface);
	free((void*)param->info);
}


void setFunctPDF(TypeFunctDrawTreeGeneric *funct) {
	funct->startStd = startPDFStd;
	funct->start = startPDF;
	funct->end = endCairo;
	setFunctCairo(funct);
}

void startPDFStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param) {
	cairo_surface_t *surface;
	cairo_t *cr;
	cairo_font_extents_t extentsF;
	surface = cairo_pdf_surface_create (filename, width, height);
	cr = cairo_create(surface);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size(cr, 12);
	cairo_font_extents(cr, &extentsF);
	param->height = height;
	param->width = width;
	param->info = malloc(sizeof(TypeCairoInfo));
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->surface = surface;
	((TypeCairoInfo*)param->info)->filename = filename;
	setParamCairo(param);
	cairo_set_source_rgb (cr, 0, 0, 0);
}

void startPDF(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param) {
	cairo_surface_t *surface;
	cairo_t *cr;
	cairo_font_extents_t extentsF;
	surface = cairo_pdf_surface_create (filename, 500, 500);
	cr = cairo_create(surface);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size(cr, 12);
	cairo_font_extents(cr, &extentsF);
	param->height = LEAFSEP*(countLeaves(tree)+1)+extentsF.descent+extentsF.ascent+2*LABELSEP+TICKLENGTH;
	if(countLeaves(tree)>40)
		param->width = param->height/sqrt(2);
	else
		param->width = 2*param->height/sqrt(2);
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
	surface = cairo_pdf_surface_create (filename, param->width, param->height);
	cr = cairo_create (surface);
	param->info = malloc(sizeof(TypeCairoInfo));
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->surface = surface;
	((TypeCairoInfo*)param->info)->filename = filename;
	setParamCairo(param);
	param->tmin = tree->minTime;
    param->tmax = tree->maxTime;
	param->labelWidth = getMaxLeafLabelWidthCairo(tree->name, tree->size, cr);	
	param->scale = (param->width-param->labelWidth-param->xoffset-param->labelSep)/((param->tmax-param->tmin));
	cairo_set_source_rgb (cr, 0, 0, 0);
}

void setFunctPS(TypeFunctDrawTreeGeneric *funct) {
	funct->startStd = startPSStd;
	funct->start = startPS;
	funct->end = endCairo;
	setFunctCairo(funct);
}

void startPSStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param) {
	cairo_surface_t *surface;
	cairo_t *cr;
	cairo_font_extents_t extentsF;
	surface = cairo_ps_surface_create (filename, width, height);
	cr = cairo_create(surface);
	cairo_select_font_face (cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size (cr, 12);
	cairo_font_extents (cr, &extentsF);
	setParamCairo(param);
	param->height = height;
	param->width = width;
	param->info = (void*) malloc(sizeof(TypeCairoInfo));
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->surface = surface;
	((TypeCairoInfo*)param->info)->filename = filename;
	cairo_set_source_rgb (((TypeCairoInfo*)param->info)->cr, 0, 0, 0);
}

void startPS(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param) {
	cairo_surface_t *surface;
	cairo_t *cr;
	cairo_font_extents_t extentsF;
	surface = cairo_ps_surface_create (filename, 500, 500);
	cr = cairo_create(surface);
	cairo_select_font_face (cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size (cr, 12);
	cairo_font_extents (cr, &extentsF);
	param->height = LEAFSEP*(countLeaves(tree)+1)+extentsF.descent+extentsF.ascent+2*LABELSEP+TICKLENGTH;
	if(countLeaves(tree)>40)
		param->width = param->height/sqrt(2);
	else
		param->width = 2*param->height/sqrt(2);
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
	surface = cairo_ps_surface_create (filename, param->width, param->height);
	cr = cairo_create (surface);
	setParamCairo(param);
	param->tmin = tree->minTime;
    param->tmax = tree->maxTime;
	param->info = (void*) malloc(sizeof(TypeCairoInfo));
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->surface = surface;
	((TypeCairoInfo*)param->info)->filename = filename;
	param->labelWidth = getMaxLeafLabelWidthCairo(tree->name, tree->size, ((TypeCairoInfo*)param->info)->cr);	
	param->scale = (param->width-param->labelWidth-param->xoffset-param->labelSep)/((param->tmax-param->tmin));
	cairo_set_source_rgb (((TypeCairoInfo*)param->info)->cr, 0, 0, 0);
}

void setFunctPNG(TypeFunctDrawTreeGeneric *funct) {
	funct->startStd = startPNGStd;
	funct->start = startPNG;
	funct->end = endPNG;
	setFunctCairo(funct);
}

void startPNGStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param) {
	cairo_surface_t *surface;
	cairo_t *cr;
	cairo_font_extents_t extentsF;
	surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, width, height);
	cr = cairo_create(surface);
	cairo_select_font_face (cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size (cr, 12);
	cairo_font_extents (cr, &extentsF);
	setParamCairo(param);
	param->height = height;
	param->width = width;
	param->info = (void*) malloc(sizeof(TypeCairoInfo));
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->surface = surface;
	((TypeCairoInfo*)param->info)->filename = filename;
	cairo_set_source_rgb (((TypeCairoInfo*)param->info)->cr, 0, 0, 0);
}

void startPNG(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param) {
	cairo_surface_t *surface;
	cairo_t *cr;
	cairo_font_extents_t extentsF;
	surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 500, 500);
	cr = cairo_create(surface);
	cairo_select_font_face (cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size (cr, 12);
	cairo_font_extents (cr, &extentsF);
	param->height = LEAFSEP*(countLeaves(tree)+1)+extentsF.descent+extentsF.ascent+2*LABELSEP+TICKLENGTH;
	if(countLeaves(tree)>40)
		param->width = param->height/sqrt(2);
	else
		param->width = 2*param->height/sqrt(2);
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
	surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, param->width, param->height);
	cr = cairo_create (surface);
	setParamCairo(param);
	param->tmin = tree->minTime;
    param->tmax = tree->maxTime;
	param->info = (void*) malloc(sizeof(TypeCairoInfo));
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->cr = cr;
	((TypeCairoInfo*)param->info)->surface = surface;
	((TypeCairoInfo*)param->info)->filename = filename;
	param->labelWidth = getMaxLeafLabelWidthCairo(tree->name, tree->size, ((TypeCairoInfo*)param->info)->cr);	
	param->scale = (param->width-param->labelWidth-param->xoffset-param->labelSep)/((param->tmax-param->tmin));
	cairo_set_source_rgb (((TypeCairoInfo*)param->info)->cr, 0, 0, 0);
}

void endPNG(TypeParamDrawTreeGeneric *param) {
	printf("status %s\n", cairo_status_to_string (cairo_surface_write_to_png (((TypeCairoInfo*)param->info)->surface, ((TypeCairoInfo*)param->info)->filename)));
	cairo_destroy (((TypeCairoInfo*)param->info)->cr);
	cairo_surface_destroy (((TypeCairoInfo*)param->info)->surface);
	free((void*)param->info);
}
