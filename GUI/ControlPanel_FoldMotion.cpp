#include "ControlPanel.h"

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

	//if( 0<=f && f<fcnt ){
	if( 0<=f ){
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
