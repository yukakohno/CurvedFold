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

void ControlPanel::cb_btn_resetcurve( Fl_Widget *wgt, void *idx )
{
	ControlPanel *This = (ControlPanel *)idx;
	int cidx = This->vs_cidx->value();
	//This->ppm.resetpostproc( cidx );
	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}

void ControlPanel::cb_vs_cidx(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->gwin_gr->cidx = ((Fl_Value_Slider*)wgt)->value();
	This->gwin_gr->redraw();
}

void ControlPanel::cb_vs_otype(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	//This->ppm.hCrs_evaltype = ((Fl_Value_Slider*)wgt)->value();
	int grpopt = This->value_grpopt();
	switch( grpopt ){
		case 0:	This->ppm.hCrs_evaltype = EVTYPE_TORSION;	break;
		case 1:	This->ppm.hCrs_evaltype = EVTYPE_BETALR;	break;
		case 2:	This->ppm.hCrs_evaltype = EVTYPE_RULCROSS;	break;
	}
}

void ControlPanel::cb_btn_rectifyC(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.flg_interpolate = 0;
	This->ppm.flg_rectifyCrs = This->vs_cidx->value();
	//This->refresh(0);
	This->ppm.set_postproc_type( PPTYPE_OPTIMIZE );
	This->ppm.postproc();
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
	This->gwin->ppm = &(This->ppm);
	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
	This->setorg();
	This->ppm.flg_rectifyCrs = 0;
	This->ppm.flg_interpolate = 1;
}

void ControlPanel::cb_btn_resetC(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.flg_rectifyCrs = 0;
	int cidx = This->vs_cidx->value();
	if( 0 < cidx && cidx < This->ppm.rcrcnt ){
		This->ppm.rcrs[cidx]->flg_org=0;
	} else if( cidx < 0 && -cidx < This->ppm.lcrcnt ){
		This->ppm.lcrs[-cidx]->flg_org=0;
	}
	This->ppm.initHCrs();
	This->refresh(0);
	This->setorg();
}

void ControlPanel::cb_btn_fixC(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.flg_rectifyCrs = 0;
	//This->ppm.CP2cv0();
	//This->ppm.cv0to1();
	//This->ppm.initHCrs();
}
