#ifndef DRAWCURVE
#define DRAWCURVE

#define MAXDCX 1000	// max vertices

#include "crease.h"

typedef enum _curve_type {
	CTYPE_UNDEF, CTYPE_TRIM, CTYPE_FOLD, CTYPE_FOLD_LEFT, CTYPE_FOLD_RIGHT
} curve_type;

typedef enum _cvertex_type {
	CVTYPE_UNDEF,
	CVTYPE_RUL, CVTYPE_RUL_LEFT, CVTYPE_RUL_RIGHT,
	CVTYPE_CREASE, CVTYPE_CREASE_START, CVTYPE_CREASE_END,
	CVTYPE_PAPER, CVTYPE_PAPER_TOP, CVTYPE_PAPER_BOTTOM, CVTYPE_PAPER_LEFT, CVTYPE_PAPER_RIGHT
} cvertex_type;

class curvedraw
{
public:
	curvedraw::curvedraw();
	curvedraw::~curvedraw();

	void init();

	curve_type ctype;

	int cvcnt;
	double cvx[MAXDCX], cvy[MAXDCX];

	// loop curve if distance between start and end points < th_d[mm]
	void closecv( double th_d );

	static int IS_drawcurve( curvedraw *dc, crease *c, int psx, int psy, int pex, int pey,
		double *ISx, double *ISy, int ISsize, int &IScnt, 
		cvertex_type *Xtype, int *Xidx, double *Xlen, int *Cidx );

	// return 1 if clockwise, -1 if ccw
	static int IS_check_cw( curvedraw *dc, crease *c, int psx, int psy, int pex, int pey,
		double *ISx, double *ISy, int ISsize, int &IScnt, 
		cvertex_type *Xtype, int *Xidx, double *Xlen, int *Cidx );
};

#endif
