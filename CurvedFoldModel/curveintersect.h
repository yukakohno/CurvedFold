#ifndef FOLDCURVE
#define FOLDCURVE

#define MAXFCX 100	// max vertices

#include "crease.h"
#include "curvedraw.h"

class curveintersect
{
public:
	curveintersect::curveintersect();
	curveintersect::~curveintersect();

	void init();

	curve_type ctype;
	int cridx;	// index of crease
	int dcidx;	// index of curvedraw

	int cvcnt;
	double cvx[MAXFCX], cvy[MAXFCX];

	// type of intersecting vertex
	cvertex_type cvtype[MAXFCX];

	// index of ruling if cvtype== CVTYPE_RUL_*
	// index of X if cvtype== CVTYPE_FCURVE_*
	// undefined if cvtype== CVTYPE_PAPER_*
	int cvidx[MAXFCX];

	// length of ruling if cvtype== CVTYPE_RUL_*
	// distance from X[cvidx] if cvtype== CVTYPE_FCURVE_*
	// distance from corner if cvtype== CVTYPE_PAPER_*
	double cvlen[MAXFCX];

	int create_drawcurve( curvedraw *dc, crease *c, int psx, int psy, int pex, int pey );
	void inverse_c();

	// loop curve if distance between start and end points < th_d[mm]
	void closecv( double th_d );
};

#endif
