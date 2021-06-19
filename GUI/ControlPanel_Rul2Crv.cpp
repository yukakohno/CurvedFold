#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <windows.h>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

void ControlPanel::cb_btn_R2TA0(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int mode = This->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C, CMODE_R
	int prm = This->value_grpparam();	// P_CV2D, P_CV3D, P_TRSN, P_FLDA, P_RULL, P_RULR, ...
	if( mode == CMODE_R && ( prm == P_RULL || prm == P_RULR || prm == P_RULL1 || prm == P_RULR1 ) )
	{
		papermodel *ppm = &(This->ppm);
		crease *c = &(ppm->crs[0]);
		double mang = This->vs_xmang0->value()/180.0*M_PI;

		ppm->re_sidx = 0;
		ppm->re_eidx = c->Xcnt;
		int ret = c->calcR_TA( 1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, mang, 0 );

		if( ret==0 ){
			ppm->set_postproc_type( PPTYPE_PRICURVE );
			ppm->postproc();
			ppm->set_postproc_type( PPTYPE_UNDEF );
			if( This->cb_optmat->value() && ppm->tgcnt>3 ){
				ppm->optMat( CMODE_R );
			} else if (This->cb_optrot->value() && ppm->tgcnt > 3) {
				ppm->optMatRot(CMODE_R); // fixed origin
			}

			This->gwin->redraw();
			This->gwin_cp->redraw();
			This->gwin_gr->redraw();
#if 1
		} else if( ret==-1 ){
			This->acnt_inc = false;
		} else if( ret==-2 ){
			This->acnt_inc = true;
#endif
		}
		This->push_hist(c->Pcnt, CMODE_R, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);
	}
}

void ControlPanel::cb_vs_xmang(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int val = This->vs_xmang1->value();
	This->vs_xmang0->value( val );
	if( val==0 ){
		This->gwin->disp_R = This->gwin_cp->disp_R = 0; // don't display rulings in folding angle=0
	} else if( This->cb_disp[D_R]->value() != 0 ){
		This->gwin->disp_R = This->gwin_cp->disp_R = 1;
	}
	This->btn_R2TA0->do_callback();
}

void ControlPanel::cb_btn_switchRuling(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	crease* c = &(This->ppm.crs[0]);
	std::cout << This->fcnt << std::endl;
	switch (This->fcnt) {
			case 0: c->loadPb("input/rulings00.txt"); break;
			case 1: c->loadPb("input/rulings01.txt"); break;
			case 2: c->loadPb("input/rulings02.txt"); break;
			case 3: c->loadPb("input/rulings03.txt"); break;
			case 4: c->loadPb("input/rulings04.txt"); break;
			case 5: c->loadPb("input/rulings05.txt"); break;
			case 6: c->loadPb("input/rulings06.txt"); break;
		}
		This->fcnt++;
	if (This->fcnt > 6) { This->fcnt = 0; }
		//if( This->fcnt>1 ){ This->fcnt=0; }

	This->btn_R2TA0->do_callback();
}

void ControlPanel::idle(void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;

	if (This->acnt >= This->vs_xmang1->maximum()) {
		This->acnt_inc = false;
	}
	if (This->acnt <= This->vs_xmang1->minimum()) {
		This->acnt_inc = true;
	}

	if (This->acnt_inc) {
		This->acnt++;
	}
	else {
		This->acnt--;
	}
	Sleep(1);

	//printf("acnt=%d\n", This->acnt);
	This->vs_xmang1->value(This->acnt);
	This->vs_xmang1->do_callback();

	printf( "average gap = %f\n", This->ppm.avetgap );
}

//-----------------------------

void ControlPanel::cb_btn_optmat(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);
	int mode = This->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C, MODE_T, MODE_FA, MODE_TFA

	int ret = 0;
	if (This->cb_optmat->value() && ppm->tgcnt > 3) {
		ret = ppm->optMat(mode);
	}
	else if (This->cb_optrot->value() && ppm->tgcnt > 3) {
		ret = ppm->optMatRot(mode); // fixed origin
	}

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}

//-----------------------------

void ControlPanel::cb_btn_start(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	//This->acnt = (int)This->vs_xmang1->value();
	papermodel *ppm = &(This->ppm);
	crease *c = &(ppm->crs[0]);

	This->rb_fix[CMODE_R]->setonly();
	This->rb_param[P_RULL]->setonly();

	if( !This->flg_idle_active ){
		This->flg_idle_active = true;
#if 0
		{	// clear file
			FILE *fout = fopen( "output/paramresult.csv","w");
			if( fout ){
				fprintf( fout, "res_rul, crossarea,," );
				for( int i=0; i<c->Pcnt; i++ ){ fprintf( fout, "torsion%d,", i ); }
				fprintf( fout, "," );
				for( int i=0; i<c->Pcnt; i++ ){ fprintf( fout, "alpha%d, ", i ); }
				fprintf( fout, "\n" );
				fclose(fout);
			}
		}
#endif
		Fl::add_idle( ControlPanel::idle, This );
	}
}

void ControlPanel::cb_btn_stop(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	if( This->flg_idle_active ){
		This->flg_idle_active = false;
		Fl::remove_idle( ControlPanel::idle, This );
	}
}
