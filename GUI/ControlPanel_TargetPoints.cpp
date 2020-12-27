#include "ControlPanel.h"

void ControlPanel::cb_btn_addTgt(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	GraphWindow3DCF* gwin = This->gwin;
	GraphWindowCP* gwin_cp = This->gwin_cp;

	if ( This->mode_modtgt ) {
		This->mode_modtgt = false;
		gwin->disp_modTGT = gwin_cp->disp_modTGT = -2;
		((Fl_Button*)wgt)->label("add / modify / delete");
	}
	else {
		This->mode_modtgt = true;
		This->vs_tgtidx->do_callback();
		((Fl_Button*)wgt)->label("return");
	}

	gwin->redraw();
	gwin_cp->redraw();
}

void ControlPanel::cb_vs_tgtidx(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	GraphWindow3DCF* gwin = This->gwin;
	GraphWindowCP* gwin_cp = This->gwin_cp;
	papermodel* ppm = &(This->ppm);
	if (This->mode_modtgt == false) {
		return;
	}

	if (This->ppm.tgcnt == 0) {
		gwin->disp_modTGT = gwin_cp->disp_modTGT = -1;
	} else {
		int index = ((Fl_Value_Slider*)wgt)->value();
		if (index < 0 || This->ppm.tgcnt <= index) {
			return;
		}
		gwin->disp_modTGT = gwin_cp->disp_modTGT = index;
	}
	//This->vs_tgtx->value(ppm->tgx[index] - ppm->ogx[index]);
	//This->vs_tgty->value(ppm->tgy[index] - ppm->ogy[index]);
	//This->vs_tgtz->value(ppm->tgz[index] - ppm->ogz[index]);
	//This->vs_tgtn->value();

	gwin->redraw();
	gwin_cp->redraw();
}
#if 0
void ControlPanel::cb_vs_tgtxyz(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	GraphWindow3DCF* gwin = This->gwin;
	GraphWindowCP* gwin_cp = This->gwin_cp;
	papermodel* ppm = &(This->ppm);
	if (This->mode_modtgt == false) {
		return;
	}
	int index = This->vs_tgtidx->value();
	if (index < 0 || This->ppm.tgcnt <= index) {
		return;
	}
	ppm->tgx[index] = ppm->ogx[index] + This->vs_tgtx->value();
	ppm->tgy[index] = ppm->ogy[index] + This->vs_tgty->value();
	ppm->tgz[index] = ppm->ogz[index] + This->vs_tgtz->value();

	gwin->redraw();
	gwin_cp->redraw();
}

void ControlPanel::cb_vs_tgtn(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	GraphWindow3DCF* gwin = This->gwin;
	GraphWindowCP* gwin_cp = This->gwin_cp;
	papermodel* ppm = &(This->ppm);
	if (This->mode_modtgt == false) {
		return;
	}

	// TOOD

	gwin->redraw();
	gwin_cp->redraw();
}
#endif