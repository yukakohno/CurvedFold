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
	This->push_hist(c->Pcnt, CMODE_B, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}

void ControlPanel::cb_btn_optfold2(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
	int mode = This->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C, CMODE_R
	int prm = This->value_grpparam();	// P_CV2D, P_CV3D, P_TRSN, P_FLDA, P_RULL, P_RULR, ...

	if ( !(mode == CMODE_R && (prm == P_RULL || prm == P_RULR || prm == P_RULL1 || prm == P_RULR1)) )
	{
		return;
	}

	int minval=-1;
	double mintgap = 10000;
	for (int val = This->vs_xmang1->minimum(); val < This->vs_xmang1->maximum(); val++)
	{
		//This->vs_xmang0->value(val);
		//This->btn_R2TA0->do_callback();

		double mang = val / 180.0 * M_PI;
		ppm->re_sidx = 0;
		ppm->re_eidx = c->Xcnt;
		int ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, mang);

		if (ret == 0) {
			ppm->set_postproc_type(PPTYPE_PRICURVE);
			ppm->postproc();
			ppm->set_postproc_type(PPTYPE_UNDEF);
			if (This->cb_optmat->value() && ppm->tgcnt > 3) {
				ppm->optMat(CMODE_R);
			}
			if (mintgap > ppm->avetgap) {
				minval = val;
				mintgap = ppm->avetgap;
			}
		}
	}
	if (minval > -1) {
		This->vs_xmang0->value(minval);
		This->vs_xmang1->value(minval);
		This->btn_R2TA0->do_callback();
	}
	This->push_hist(c->Pcnt, CMODE_R, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);
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
