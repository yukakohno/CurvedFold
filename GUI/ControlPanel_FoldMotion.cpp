#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Windows.h>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

void ControlPanel::cb_vs_fmot( Fl_Widget *wgt, void *idx)
{
	int i, f, fcnt;
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);

#if 1
	f = This->vs_fmot->value();
#else
	if( wgt == (Fl_Widget *)This->vs_fmot ){
		f = This->vs_fmot->value();
		This->vs_fmot2->value( f );
	} else {
		f = This->vs_fmot2->value();
		This->vs_fmot->value( f );
	}
#endif
	fcnt = This->an_fcnt;
	This->rb_fix[CMODE_B]->setonly();
	This->rb_param[P_FLDA]->setonly();

	if( 0<=f && f<fcnt ){
		for( i=0; i<c->Xcnt; i++ ){
			c->k2d[i] = This->k2d_org[i] - This->k2d_dif0[i]*f;
			c->tr[i] = This->tr_org[i] - This->tr_dif0[i]*f;
			c->alpha[i] = This->alpha_org[i] - This->alpha_dif0[i]*f;
		}
	} else if( -fcnt<f ){
		for( i=0; i<c->Xcnt; i++ ){
			c->k2d[i] = This->k2d_org[i] + This->k2d_dif1[i]*f;
			c->tr[i] = This->tr_org[i] + This->tr_dif1[i]*(-f);
			c->alpha[i] = This->alpha_org[i] + This->alpha_dif1[i]*(-f);
		}
	}

	This->ppm.set_postproc_type( PPTYPE_FMOT );
	This->ppm.flg_interpolate = 0;
	This->refresh( 0 );
	This->ppm.flg_interpolate = 1;
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
}

void ControlPanel::cb_vs_ppos2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int prm = This->value_grpparam();
	int frm = This->vs_fmot2->value() + c->FM_fidx_org;
	int pos = This->vs_ppos2->value();

	switch( prm ){
		case P_FLDA:
			This->vs_pval2->bounds( This->amin, This->amax );
			This->vs_pval2->step( This->astep );
			This->vs_pval2->value( (int)(c->Pa[pos]*180.0/M_PI) );
			break;
		case P_FLDA1:
			This->vs_pval2->bounds( This->amin, This->amax );
			This->vs_pval2->step( This->astep );
			This->vs_pval2->value( (int)(c->FM_Pa[frm][pos]*180.0/M_PI) );
			break;
		case P_TRSN:
			This->vs_pval2->bounds( This->tmin, This->tmax );
			This->vs_pval2->step( This->tstep );
			This->vs_pval2->value( c->Py[pos] );
			break;
		case P_TRSN1:
			This->vs_pval2->bounds( This->tmin, This->tmax );
			This->vs_pval2->step( This->tstep );
			This->vs_pval2->value( c->FM_Pt[frm][pos] );
			break;
	}
	This->gwin->pprm = This->gwin_cp->pprm = prm;
	This->gwin->ppos = This->gwin_cp->ppos = pos;
	if( This->gwin->disp_CP ){ This->gwin->redraw(); }
	if( This->gwin_cp->disp_CP ){ This->gwin_cp->redraw(); }
}

void ControlPanel::cb_vs_pval2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);
	crease *c = &(This->ppm.crs[0]);
	int prm = This->value_grpparam();
	int pos = This->vs_ppos2->value();
	double val = This->vs_pval2->value();
	//printf("val=%f\n",val);

	switch( prm ){
		case P_FLDA:
		case P_FLDA1:
			c->Pa[pos] = val*M_PI/180.0;
			This->rb_param[P_FLDA]->setonly();
			This->rb_fix[CMODE_B]->setonly();
			This->ppm.set_postproc_type( PPTYPE_PRICURVE );
			This->refresh( 0 );
			This->ppm.set_postproc_type( PPTYPE_UNDEF );
			break;
		case P_TRSN:
		case P_TRSN1:
			c->Py[pos] = val;
			This->rb_param[P_TRSN]->setonly();
			This->rb_fix[CMODE_B]->setonly();
			This->ppm.set_postproc_type( PPTYPE_PRICURVE );
			This->refresh( 0 );
			This->ppm.set_postproc_type( PPTYPE_UNDEF );
			break;
	}
}

void ControlPanel::cb_btn_paramopt(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int divnum = This->gwin->divnum = This->vs_divnum->value();
	int frm = This->vs_fmot2->value() + c->FM_fidx_org;

	// TODO: optimize parameters

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}

void ControlPanel::cb_btn_paramset(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int prm = This->value_grpparam();
	int frm = This->vs_fmot2->value() + c->FM_fidx_org;

	memcpy( &(c->FM_Pa[frm][0]), c->Pa, sizeof(double)*c->Pcnt );
	c->flg_FM_Pa[frm] = 1;
	memcpy( &(c->FM_Pt[frm][0]), c->Py, sizeof(double)*c->Pcnt );
	c->flg_FM_Pt[frm] = 1;
#if FM_MODE==2
	memcpy( &(c->FM_tr[frm][0]), c->tr, sizeof(double)*c->Xcnt );
	memcpy( &(c->FM_alpha[frm][0]), c->alpha, sizeof(double)*c->Xcnt );
#endif
	c->updateFM( 0 );
	cb_vs_fmot2( wgt, idx);
}

void ControlPanel::cb_btn_paramreset(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int prm = This->value_grpparam();
	int frm = This->vs_fmot2->value() + c->FM_fidx_org;

	if( frm == c->FM_fidx_min || frm == c->FM_fidx_org || frm == c->FM_fidx_max ){
		return;
	}
#if 1
	c->flg_FM_Pa[frm] = 0;
	c->flg_FM_Pt[frm] = 0;
#else
	switch( prm ){
		case P_FLDA:
		case P_FLDA1:
			c->flg_FM_Pa[frm] = 0;
			break;
		case P_TRSN:
		case P_TRSN1:
			c->flg_FM_Pt[frm] = 0;
			break;
	}
#endif
	c->updateFM( 0 );
	cb_vs_fmot2( wgt, idx);
}

void ControlPanel::cb_btn_matopt(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int divnum = This->gwin->divnum = This->vs_divnum->value();
	int frm = This->vs_fmot2->value() + c->FM_fidx_org;

	This->gwin->resetObjRot();
	This->gwin->resetObjTrans();

	double mobj[16];
	int ret = This->ppm.getMat_stitch1( divnum, -1, mobj ); // ç≈èâ10ì_ÇÃÇ›égÇ§
	for( int i=0; i<16; i++ ){ This->gwin->mObject[i] = mobj[i]; }

	This->gwin->redraw();
}

void ControlPanel::cb_btn_matset(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int frm = This->vs_fmot2->value() + c->FM_fidx_org;

	double mobj[16];
	for( int i=0; i<16; i++ ){ mobj[i] = This->gwin->mObject[i]; }
	This->gwin->resetObjRot();
	This->gwin->resetObjTrans();

	mult_m44_n44( c->FM_m3[frm], mobj, 1);
	c->flg_FM_m3[frm] = 1;
	c->updateFM( 0 );

	This->vs_fmot2->do_callback();
}

void ControlPanel::cb_btn_matreset(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int frm = This->vs_fmot2->value() + c->FM_fidx_org;

	c->flg_FM_m3[frm] = 0;
	c->updateFM( 0 );

	This->vs_fmot2->do_callback();
}

void ControlPanel::update_bx_FM( crease *c, int frm )
{
	//crease *c = &(ppm.crs[0]);
	int min = (int)vs_fmot2->minimum();
	int max = (int)vs_fmot2->maximum();
	for( int j=min; j<=max; j++ ){
		int i = j + c->FM_fidx_org;
		if( c->flg_FM_Pt[i] ){
			bx_FM_Pt[i]->color(FL_BLUE);
		} else if( i==frm ) {
			bx_FM_Pt[i]->color(FL_BLACK);
		} else {
			bx_FM_Pt[i]->color(FL_WHITE);
		}
		if( c->flg_FM_Pa[i] ){
			bx_FM_Pa[i]->color(FL_RED);
		} else if( i==frm ) {
			bx_FM_Pa[i]->color(FL_BLACK);
		}else{
			bx_FM_Pa[i]->color(FL_WHITE);
		}
		if( c->flg_FM_m3[i] ){
			bx_FM_m3[i]->color(FL_GREEN);
		} else if( i==frm ) {
			bx_FM_m3[i]->color(FL_BLACK);
		} else {
			bx_FM_m3[i]->color(FL_WHITE);
		}
	}
	redraw();
}

void ControlPanel::cb_vs_fmot2( Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int flg_usecp = This->cb_usecp->value();

	int frm = This->vs_fmot2->value() + c->FM_fidx_org;
	int prm = This->value_grpparam();

	This->rb_fix[CMODE_B]->setonly();
	This->rb_param[P_FLDA]->setonly();

	This->gwin->resetObjRot();
	This->gwin->resetObjTrans();

	memcpy( c->Pa, &(c->FM_Pa[frm][0]), sizeof(double)*c->Pcnt );
	memcpy( c->Py, &(c->FM_Pt[frm][0]), sizeof(double)*c->Pcnt );

	if( flg_usecp==0 ){
		memcpy( c->k2d, This->k2d_org, sizeof(double)*c->Xcnt );
		memcpy( c->tr, &(c->FM_tr[frm][0]), sizeof(double)*c->Xcnt );
		memcpy( c->alpha, &(c->FM_alpha[frm][0]), sizeof(double)*c->Xcnt );
		This->ppm.flg_interpolate = 0;
	}
	memcpy( c->m3, &(c->FM_m3[frm][0]), sizeof(double)*16 );
	//for( int i=0; i<16; i++ ){ This->gwin->mObject[i] = c->FM_m3[frm][i]; }
#if 0
	//float *mo = This->gwin->mObject;
	printf("frm: %d, flg_flg_FM_m3 : %d\n", frm, c->flg_FM_m3[frm]);
	printf("%f\t%f\t%f\t%f\n", c->m3[0], c->m3[1], c->m3[2], c->m3[3]);
	printf("%f\t%f\t%f\t%f\n", c->m3[4], c->m3[5], c->m3[6], c->m3[7]);
	printf("%f\t%f\t%f\t%f\n", c->m3[8], c->m3[9], c->m3[10], c->m3[11]);
	printf("%f\t%f\t%f\t%f\n\n", c->m3[12], c->m3[13], c->m3[14], c->m3[15]);
#endif

	This->update_bx_FM( c, frm );

	This->ppm.set_postproc_type( PPTYPE_FMOT );
	//This->ppm.set_postproc_type( PPTYPE_PRICURVE );
	if( flg_usecp==0 ){
		This->ppm.flg_interpolate = 0;
	}
	This->refresh( 0 );
	This->ppm.flg_interpolate = 1;
	This->ppm.set_postproc_type( PPTYPE_UNDEF );

	This->rb_param[prm]->setonly();
	This->rb_param[prm]->do_callback();
}

void ControlPanel::cb_vs_divnum(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);
	int divnum = This->gwin->divnum = This->vs_divnum->value();

#if 1	// right: add curve, left: trim
	double crv = This->vs_divcrv->value();

	if( This->rb_fold->value() ){
		ppm->setSpherePieceBoundary( divnum, crv, CTYPE_FOLD );
	} else {
		ppm->setSpherePieceBoundary( divnum, crv, CTYPE_TRIM );
	}

	This->rb_fix[CMODE_B]->setonly();
	This->rb_param[P_FLDA]->setonly();
	This->ppm.set_postproc_type( PPTYPE_ADDCREASE );
	This->refresh( 0 );
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
#endif

#if 0	// trim
	double angle = 0.5*M_PI - 2.0*M_PI/(double)divnum;
	ppm->dccnt=0; // âº
	if( divnum>4 ){
		int cvcnt=2;
		double *cvx0 = ppm->dcurve[ppm->dccnt].cvx;
		double *cvy0 = ppm->dcurve[ppm->dccnt].cvy;
		cvx0[0] = 0.0;
		cvy0[0] = 0.0;
		cvx0[1] = sin(angle)*300;
		cvy0[1] = cos(angle)*300;
		ppm->dcurve[ppm->dccnt].cvcnt = cvcnt;
		ppm->dcurve[ppm->dccnt].ctype = CTYPE_TRIM;
		ppm->dccnt++;
	}
	//This->ppm.set_postproc_type( PPTYPE_PPEDGE );
	This->ppm.set_postproc_type( PPTYPE_ADDCREASE );
	This->refresh( 0 );
	//This->ppm.postproc();
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
#endif
	This->gwin->redraw();
	This->gwin_cp->redraw();
}

void ControlPanel::cb_vs_divcrv(Fl_Widget *wgt, void *idx)
{
	//ControlPanel *This = (ControlPanel *)idx;
	//This->vs_divnum->do_callback();
}

void ControlPanel::cb_rb_divtype(Fl_Widget *wgt, void *idx)
{
}

