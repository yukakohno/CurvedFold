#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <windows.h>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

void ControlPanel::cb_btn_optfold(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);

	int ret = ppm->optFold();
	if (This->cb_optmat->value() && ppm->tgcnt > 3) {
		ppm->optMat(CMODE_B); // B: Ü‚èü‚ÆÜ‚èŠp“x‚©‚çA3D‹Èü‚ð‹‚ß‚é
	}
	crease* c = &(This->ppm.crs[0]);

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}

void ControlPanel::cb_btn_opttr(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);

	int ret = ppm->optTorsion();

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}
