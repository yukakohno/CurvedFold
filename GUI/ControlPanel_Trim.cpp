#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Windows.h>
#include "ControlPanel.h"

void ControlPanel::cb_cb_addcurve( Fl_Widget *wgt, void *idx )
{
	ControlPanel *This = (ControlPanel *)idx;
	int flg=0;
	if( This->cb_trimcurve->value() ){
		flg=1;
	}
	if( This->cb_foldcurve->value() ){
		flg=2;
	}
	This->gwin->flg_addcurve = flg;
	This->gwin_cp->flg_addcurve = flg;
	This->gwin_cp->ccnt=0;
}

void ControlPanel::cb_btn_proccurve( Fl_Widget *wgt, void *idx )
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.set_postproc_type( PPTYPE_ADDCREASE );
	This->ppm.postproc();
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}


