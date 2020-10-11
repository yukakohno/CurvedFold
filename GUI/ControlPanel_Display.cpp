#include "ControlPanel.h"

void ControlPanel::cb_cb_dispaxis(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	GraphWindowCP *gwin_cp = This->gwin_cp;
	gwin->disp_axis = gwin_cp->disp_axis = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
	gwin_cp->redraw();
}

void ControlPanel::cb_cb_dispX(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	GraphWindowCP *gwin_cp = This->gwin_cp;
	gwin->disp_X = gwin_cp->disp_X = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
	gwin_cp->redraw();
}

void ControlPanel::cb_cb_dispTNB(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	GraphWindowCP *gwin_cp = This->gwin_cp;
	gwin->disp_TNB = gwin_cp->disp_TNB = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
	gwin_cp->redraw();
}

void ControlPanel::cb_cb_dispR(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	GraphWindowCP *gwin_cp = This->gwin_cp;
	gwin->disp_R = gwin_cp->disp_R = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
	gwin_cp->redraw();
}

void ControlPanel::cb_cb_dispRlen(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	GraphWindowCP *gwin_cp = This->gwin_cp;
	gwin->disp_maxrlen = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
	//gwin_cp->redraw();
}

void ControlPanel::cb_cb_dispAx2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	GraphWindowCP *gwin_cp = This->gwin_cp;
	gwin->disp_axis2 = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
	//gwin_cp->redraw();
}

void ControlPanel::cb_cb_dispPLY(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	GraphWindowCP *gwin_cp = This->gwin_cp;
	gwin->disp_PLY = gwin_cp->disp_PLY = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
	gwin_cp->redraw();
}

void ControlPanel::cb_cb_dispPTN(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	GraphWindowCP *gwin_cp = This->gwin_cp;
	gwin->disp_PTN = gwin_cp->disp_PTN = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
	gwin_cp->redraw();
}

void ControlPanel::cb_cb_dispCP(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	GraphWindowCP *gwin_cp = This->gwin_cp;
	gwin->disp_CP = gwin_cp->disp_CP = ((Fl_Check_Button*)wgt)->value();
	if( gwin->disp_CP ){
		gwin->pprm = gwin_cp->pprm = This->value_grpparam();
		gwin->ppos = gwin_cp->ppos = This->vs_ppos->value();
	}
	gwin->redraw();
	gwin_cp->redraw();
}

void ControlPanel::cb_cb_CurveEnd(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int val = ((Fl_Check_Button*)wgt)->value();
	if( val==0 ){
		This->ppm.flg_postproc_type = PPTYPE_X;
		This->refresh(1);
	} else {
		This->ppm.flg_postproc_type = PPTYPE_PRICURVE;
		This->refresh(1);
		This->ppm.flg_postproc_type = PPTYPE_UNDEF;
	}
}

void ControlPanel::cb_cb_dispLSMT(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	gwin->disp_LIN_SMOOTH = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
}

void ControlPanel::cb_cb_dispOFST(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	gwin->disp_POLY_OFFSET = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
}

void ControlPanel::cb_cb_dispTGT(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	gwin->disp_TGT = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
	GraphWindowCP *gwin_cp = This->gwin_cp;
	gwin_cp->disp_TGT = ((Fl_Check_Button*)wgt)->value();
	gwin_cp->redraw();
}
