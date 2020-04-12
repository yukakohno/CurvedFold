#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <windows.h>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

void ControlPanel::cb_btn_optfold(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);

	int ret = ppm->optFold();

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
