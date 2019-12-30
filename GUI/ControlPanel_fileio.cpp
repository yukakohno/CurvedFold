#define _USE_MATH_DEFINES
#include <math.h>
#include "ControlPanel.h"

void ControlPanel::cb_btn_load(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindow3DCF *gwin = This->gwin;
	Fl_File_Chooser *fc = This->fc;
	char wgtlbl[256];	strcpy(wgtlbl, wgt->label());
	int count;

	if(!strcmp(wgtlbl,"load"))
	{
		fc->show();

		while (fc->visible())
			Fl::wait();

		count = fc->count();

		if( count > 0 )
		{
			int sts = 0;
			char fname[1024], fname_m2m3[1024];
			strcpy(fname, fc->value());

			if( strstr( fname, ".txt") )
			{
				strcpy(fname_m2m3, fname);
				int fnamelen = strlen(fname_m2m3);
				sprintf( &(fname_m2m3[ fnamelen-4 ]), "_m2m3.txt\0" );

				// load data
				int filetype = This->ppm.crs[0].load(fname);
				int ret0 = This->ppm.crs[0].loadm2m3(fname_m2m3);
				This->ppm.set_postproc_type( PPTYPE_OPEN );
				//gwin->initObject();

				if(filetype == 1){ // xy(sCP),xyz(3D)
					This->ppm.crs[0].interpolate_spline2(0);
					This->ppm.crs[0].calcTN2d();
					This->ppm.crs[0].calcTNB();
					This->ppm.crs[0].Pcnt = 10;
					This->ppm.crs[0].setP_k(-1); // kv -> Px[_Pcnt]
					This->ppm.crs[0].setP_t(-1); // tr -> Py[_Pcnt]
					This->ppm.crs[0].setP_k2(-1); // kv2d -> Px2d[_Pcnt]
					This->rb_fix[CMODE_C]->setonly();
					if( This->rb_param[P_FLDA]->value() ){
						This->rb_param[P_CV3D]->setonly();
					}
					This->refresh(1);
				} else if( filetype == 2) { // kta
					This->rb_fix[CMODE_A]->setonly();
					if( This->rb_param[P_CV2D]->value() ){
						This->rb_param[P_CV3D]->setonly();
					}
					This->refresh(1);
				} else if( filetype == 3) { // xy(CP),a
					This->rb_fix[CMODE_B]->setonly();
					if( This->rb_param[P_CV3D]->value() || This->rb_param[P_FLDA]->value() ){
						This->rb_param[P_CV2D]->setonly();
					}
					This->refresh(1);
				} else if( filetype == 4) { // xyz(3D),a
					This->ppm.crs[0].interpolate_spline2(2); // 2:xyz‚Ì‚Ý•âŠÔ
					This->ppm.crs[0].calcTNB(); // X -> T,N,B,kv,tr
					This->ppm.crs[0].Pcnt = 10;
					This->ppm.crs[0].setP_k(-1); // kv -> Px[_Pcnt]
					This->ppm.crs[0].setP_t(-1); // tr -> Py[_Pcnt]
					This->ppm.crs[0].setA( 45.0/180.0*M_PI ); // alpha
					This->ppm.crs[0].setP_a(-1); // alpha -> Pa[_Pcnt]
					This->rb_fix[CMODE_A]->setonly();
					if( This->rb_param[P_CV2D]->value() ){ // kv2d
						This->rb_param[P_CV3D]->setonly(); // kv
					}
					This->refresh(1);
				} else if( 10 <= filetype && filetype < 10+6 ) {
					int mode=filetype-10;
					This->rb_fix[mode]->setonly();
					switch(mode){ // CMODE_A, CMODE_B, CMODE_C, CMODE_T, CMODE_FA, CMODE_TFA
						case CMODE_A: This->rb_param[P_CV3D]->setonly(); break;
						case CMODE_B: This->rb_param[P_CV2D]->setonly(); break;
						case CMODE_C: This->rb_param[P_CV3D]->setonly(); break;
							//case CMODE_TFA: This->rb_param[P_RULL]->setonly(); break;
					}
					This->refresh(1);
				}
				This->ppm.set_postproc_type( PPTYPE_UNDEF );
			}
			if( sts==0 ){
				wgt->label("clr");
			}
		}
		This->vs_ppos->bounds( 0, This->ppm.crs[0].Pcnt-1 );
	}
	else if(!strcmp(wgtlbl,"clr"))
	{
		// clear data
		gwin->ppm = NULL;
		gwin->clear();
		gwin->initObject();
		wgt->label("load");
	}
}

void ControlPanel::cb_btn_dump( Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int mode = This->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C, MODE_T, MODE_FA, MODE_TFA
	papermodel *ppm = &(This->ppm);
	crease *c0 = &(This->ppm.crs[0]);

	c0->dumpm2m3("output/m2m3.txt");
	c0->dumpP_("output/kta_current.txt", mode);
	c0->dumpP("output/P.txt", mode);
	c0->dumpX("output/X.csv");
	//c0->checkCP3D("output/checkCP3D.csv");

	ppm->dumpcv0( "output/cv0.txt" );
	ppm->dumpsvg0("output/CP0.svg");
	ppm->dumpsvg1("output/CP1.svg");
	ppm->makeObj3( 10.0 );
	ppm->dumpObj("output/CP.obj", 0.4);

	ppm->check180( "output/d180.csv" );	// ’¸“_Žü‚è‚ÌŠp“x
	ppm->checkquadplane( "output/quat.csv" );	// quad ‚Ì•½–Ê“x

	This->gwin_gr->redraw();
}
