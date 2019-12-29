
/*************************************************************
	GraphWindowCP.h
*************************************************************/

#include <FL/Fl.h>
#include <FL/Fl_Double_Window.h>
#include <FL/Fl_JPEG_Image.h>
#include <FL/Fl_Box.h>
#include "common.h"
#include "../CurvedFoldModel/papermodel.h"

#ifndef GRAPHWINDOWCP
#define GRAPHWINDOWCP

class GraphWindowCP : public Fl_Double_Window
{
protected:
	void draw();
	int handle(int);
	void init();

	int push_x, push_y;	// position clicked (used in dragging)

public:
	GraphWindowCP(int X, int Y, int W, int H, const char *L);
	GraphWindowCP(int X, int Y, int W, int H);

	Fl_JPEG_Image *jpg;
	int disp_axis;
	int disp_X;
	int disp_TNB;
	int disp_R;
	int disp_PLY;
	int disp_PTN;
	int disp_CP, ppos, pprm;

	int flg_addcurve; // 0:-, 1:trim, 2:fold

	void *cwin;
	papermodel *ppm;

	double ox, oy, wscl, hscl;

	double ofs_x, ofs_y;
	double ofs_a, ofs_sina, ofs_cosa;
	int flg_psize, ofs_psx, ofs_psy, ofs_pex, ofs_pey;

	int cridx, cpidx;
	double org_cpx, org_cpy;

	int ccnt;
	int cx[MAXDCX], cy[MAXDCX];

	double mt[9];

	void fixCurve();
	void fixCurve0();
};

#endif // GRAPHWINDOWCP


