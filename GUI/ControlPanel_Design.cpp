#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Windows.h>
#include "ControlPanel.h"

void ControlPanel::cb_cb_rectifyA(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.flg_rectifyA = ((Fl_Check_Button*)wgt)->value()!=0;
	This->ppm.set_postproc_type( PPTYPE_PRICURVE );
	This->refresh(1);
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
}

void ControlPanel::cb_cb_rectifyT(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.flg_rectifyT = ((Fl_Check_Button*)wgt)->value()!=0;
	This->ppm.set_postproc_type( PPTYPE_PRICURVE );
	This->refresh(1);
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
}

void ControlPanel::cb_cb_rectifyR(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.flg_rectifyR = ((Fl_Check_Button*)wgt)->value()!=0;
	This->ppm.set_postproc_type( PPTYPE_PPEDGE );
	//This->refresh(1);
	This->ppm.postproc();
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
}

void ControlPanel::cb_vs_rulres(Fl_Widget *wgt, void* idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	for( int i=0; i<MAX_CRS_CNT; i++){
		This->ppm.crs[i].XcntOrg = ((Fl_Value_Slider*)wgt)->value();
		This->ppm.crs[i].XspcOrg = (double)(XCNT_DEF*XSPC_DEF)/((double)This->ppm.crs[0].XcntOrg);
	}
	This->ppm.set_postproc_type( PPTYPE_PRICURVE );
	This->refresh(1);
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
}

void ControlPanel::cb_rb_gwin(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	char wgtlbl[256];	strcpy(wgtlbl, wgt->label());

	if(!strcmp(wgtlbl,"Y-DN") || !strcmp(wgtlbl,"Z-UP") || !strcmp(wgtlbl,"X-UP"))
	{
		// reset rotation: up = 0:y-up, 1:z-up, 2:x-up
		if(!strcmp(wgtlbl,"Y-DN")){
			wgt->label("Z-UP");
			gwin->resetRotation(1);
		} else if(!strcmp(wgtlbl,"Z-UP")){
			wgt->label("X-UP");
			gwin->resetRotation(2);
		} else if(!strcmp(wgtlbl,"X-UP")){
			//gwin->resetRotation(0);
			wgt->label("Y-DN");
			gwin->resetRotation(-1);
		}
	} else if(!strcmp(wgtlbl,"SCL")){	// reset scale
		gwin->resetScale();
	} else if(!strcmp(wgtlbl,"TRNS")){	// reset translation
		gwin->resetTranslation();
	}
	gwin->redraw();
}

void ControlPanel::cb_vs_ppos(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int pos = This->vs_ppos->value();
	int prm = This->value_grpparam();	//P_CV2D, P_CV3D, P_TRSN, P_FLDA
	if( prm == -1 ){
		return;
	}

	switch( prm ){
		case P_CV2D: // curv2d
			This->vs_pval->bounds( This->kmin, This->kmax );
			This->vs_pval->step( This->kstep );
			This->vs_pval->value( This->ppm.crs[0].Px2d[pos] );
			break;
		case P_CV3D: // curvature
			This->vs_pval->bounds( This->kmin, This->kmax );
			This->vs_pval->step( This->kstep );
			This->vs_pval->value( This->ppm.crs[0].Px[pos] );
			break;
		case P_TRSN: // torsion
			This->vs_pval->bounds( This->tmin, This->tmax );
			This->vs_pval->step( This->tstep );
			This->vs_pval->value( This->ppm.crs[0].Py[pos] );
			break;
		case P_FLDA: // angle
			This->vs_pval->bounds( This->amin, This->amax );
			This->vs_pval->step( This->astep );
			This->vs_pval->value( (int)(This->ppm.crs[0].Pa[pos]*180.0/M_PI) );
			break;
	}

	This->gwin->pprm = This->gwin_cp->pprm = prm;
	This->gwin->ppos = This->gwin_cp->ppos = This->vs_ppos->value();
	if( This->gwin->disp_CP ){ This->gwin->redraw(); }
	if( This->gwin_cp->disp_CP ){ This->gwin_cp->redraw(); }
}

void ControlPanel::cb_vs_pval(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int pos = This->vs_ppos->value();
	int prm = This->value_grpparam();	//P_CV2D, P_CV3D, P_TRSN, P_FLDA
	This->ppm.set_postproc_type( PPTYPE_PRICURVE );
	switch( prm ){
		case P_CV2D: // curv2d
			c->Px2d[pos] = This->vs_pval->value();
			This->refresh(1);
			break;
		case P_CV3D: // curvature
			c->Px[pos] = This->vs_pval->value();
			This->refresh(1);
			break;
		case P_TRSN: // torsion
			c->Py[pos] = This->vs_pval->value();
			This->refresh(1);
			break;
		case P_FLDA: // angle
			c->Pa[pos] = This->vs_pval->value()*M_PI/180.0;
			This->refresh(1);
			break;
	}
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
}

void ControlPanel::cb_btn_apply(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int mode = This->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C

}
