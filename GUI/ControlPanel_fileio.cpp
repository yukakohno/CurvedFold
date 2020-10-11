#define _USE_MATH_DEFINES
#include <math.h>
#include <sys\stat.h>
#include <direct.h>
#include <time.h>

#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

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
				strcpy( This->filepath, fname );
				for( int i=strlen(This->filepath)-2; i>0; i-- ){
					if( This->filepath[i]=='/' ){
						This->filepath[i+1] = 0;
						break;
					}
				}
				sprintf( fname_m2m3, "%sm2m3.txt", This->filepath );

				// load data
				int filetype = This->ppm.crs[0].load(fname);
				int ret0 = This->ppm.crs[0].loadm2m3(fname_m2m3);
				//for( int i=0; i<16; i++ ){ This->gwin->mObject[i] = This->ppm.crs[0].m3[i]; }
				//unit_m44( This->ppm.crs[0].m3 );
				int cur_fix = This->value_grpfix();
				int cur_param = This->value_grpparam();

				This->ppm.set_postproc_type( PPTYPE_OPEN );
				//gwin->initObject();

				if(filetype == 1){ // xy(CP),xyz(3D)
					This->ppm.crs[0].interpolate_spline2(0);
					This->ppm.crs[0].calcTN2d();
					This->ppm.crs[0].calcTNB();
					This->ppm.crs[0].Pcnt = 10;
					This->ppm.crs[0].setP_k(-1); // kv -> Px[_Pcnt]
					This->ppm.crs[0].setP_t(-1); // tr -> Py[_Pcnt]
					This->ppm.crs[0].setP_k2(-1); // kv2d -> Px2d[_Pcnt]
					This->ppm.crs[0].setP_Bl(-1);
					This->ppm.crs[0].setP_Br(-1);
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
					This->ppm.crs[0].interpolate_spline2(2); // 2:xyzのみ補間
					This->ppm.crs[0].calcTNB(); // X -> T,N,B,kv,tr
					This->ppm.crs[0].Pcnt = 10;
					This->ppm.crs[0].setP_k(-1); // kv -> Px[_Pcnt]
					This->ppm.crs[0].setP_t(-1); // tr -> Py[_Pcnt]
					This->ppm.crs[0].setA( 45.0/180.0*M_PI ); // alpha
					This->ppm.crs[0].setP_a(-1); // alpha -> Pa[_Pcnt]
					This->ppm.crs[0].setP_Bl(-1);
					This->ppm.crs[0].setP_Br(-1);
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
					}
					This->refresh(1);
				}
				This->ppm.set_postproc_type( PPTYPE_UNDEF );
				This->rb_fix[cur_fix]->setonly();
				This->rb_param[cur_param]->setonly();
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

void ControlPanel::cb_btn_savelog( Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int mode = This->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C, MODE_T, MODE_FA, MODE_TFA
	papermodel *ppm = &(This->ppm);
	crease *c0 = &(This->ppm.crs[0]);
	char dirname[128], fname[128];
    struct stat statBuf;
	time_t t0 = time(NULL);
	struct tm *t1 = localtime(&t0);

	sprintf( dirname, "output/%02d%02d%02d/", t1->tm_hour, t1->tm_min, t1->tm_sec );
	if( stat(dirname, &statBuf)==0 || _mkdir( dirname )==0 ){

		sprintf(fname, "%sm2m3.txt", dirname );	c0->dumpm2m3( fname );
		//c0->dumpm2m3("output/m2m3.txt", This->gwin->mObject );
		//c0->dumpP_("output/kta_current.txt", mode);
		sprintf(fname, "%sP.txt", dirname );	c0->dumpP( fname, mode);
		sprintf(fname, "%sX.csv", dirname );	c0->dumpX( fname );
		//c0->checkCP3D("output/checkCP3D.csv");
		//sprintf(fname, "%smotion.txt", dirname );	c0->dumpMotion( fname, flg_usecp );
		sprintf(fname, "%sCP1.svg", dirname );	ppm->dumpsvg1( fname );
		ppm->makeObj3( 10.0 );
		sprintf(fname, "%sCP.obj", dirname );	ppm->dumpObj( fname, 0.4);

		//ppm->check180( "output/d180.csv" );	// 頂点周りの角度
		//ppm->checkquadplane( "output/quad.csv" );	// quad の平面度
		sprintf(fname, "%sd180.csv", dirname );	c0->check180( fname );	// 頂点周りの角度
		sprintf(fname, "%squad.csv", dirname );	c0->checkquadplane( fname );	// quad の平面度

		sprintf(fname, "%smodel.bmp", dirname );	This->gwin->exportImage( fname );
	}
	//This->gwin_gr->redraw();
}

void ControlPanel::cb_btn_savescreen( Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c0 = &(This->ppm.crs[0]);
	time_t t0 = time(NULL);
	struct tm *t1 = localtime(&t0);
	char fname[128];

	sprintf(fname, "output/model_%02d%02d%02d.bmp", t1->tm_hour, t1->tm_min, t1->tm_sec );
	This->gwin->exportImage( fname );
}

void ControlPanel::cb_btn_loadrul( Fl_Widget *wgt, void *idx )
{
	ControlPanel *This = (ControlPanel *)idx;
	crease* c = &(This->ppm.crs[0]);
	Fl_File_Chooser* fc = This->fc_rul;

	fc->show();

	while (fc->visible())
		Fl::wait();

	if(fc->count() > 0)
	{
		char fname[1024];
		strcpy(fname, fc->value());
		c->loadPb(fname);
		This->btn_R2TA0->do_callback();
		//This->gwin->redraw();
		//This->gwin_cp->redraw();
	}
}

void ControlPanel::cb_btn_saverul( Fl_Widget *wgt, void *idx )
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	c->dumpPb("output/rulings.txt");
}

void ControlPanel::cb_btn_loadtpt( Fl_Widget *wgt, void *idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	Fl_File_Chooser* fc = This->fc_tpt;

	fc->show();

	while (fc->visible())
		Fl::wait();

	if (fc->count() > 0)
	{
		char fname[1024];
		strcpy(fname, fc->value());
		strcpy(This->fname_tpt, fc->value());
		ppm->loadTgt(fname);
	This->gwin->redraw();
		This->gwin_cp->redraw();
	}
}

void ControlPanel::cb_btn_loadtptmask(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	Fl_File_Chooser* fc = This->fc_tptmask;

	fc->show();

	while (fc->visible())
		Fl::wait();

	if (fc->count() > 0)
	{
		char fname[1024];
		strcpy(fname, fc->value());
		ppm->loadTgtMask(This->fname_tpt, fname);
		This->gwin->redraw();
		This->gwin_cp->redraw();
	}
}

void ControlPanel::cb_btn_savetpt( Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);
	ppm->saveTgt("output/target.txt");
	//This->gwin->redraw();
}

void ControlPanel::cb_btn_loadall(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	GraphWindow3DCF* gwin = This->gwin;
	Fl_File_Chooser* fc = This->fc_all;
	char wgtlbl[256];	strcpy(wgtlbl, wgt->label());
	int count;

	fc->show();

	while (fc->visible())
		Fl::wait();

	count = fc->count();

	if (count > 0)
	{
		int sts = 0;
		char fname[1024], fname_m2m3[1024];
		strcpy(fname, fc->value());
		strcpy(This->filepath, fname);
		for (int i = strlen(This->filepath) - 2; i > 0; i--) {
			if (This->filepath[i] == '/') {
				This->filepath[i + 1] = 0;
				break;
			}
		}
		sprintf(fname_m2m3, "%sm2m3.txt", This->filepath);

		This->ppm.crs[0].loadAll(fname);
		This->ppm.crs[0].loadm2m3(fname_m2m3);

		This->ppm.set_postproc_type(PPTYPE_OPEN);
		This->ppm.postproc();
		This->ppm.set_postproc_type(PPTYPE_UNDEF);

		This->gwin->redraw();
		This->gwin_cp->redraw();
		This->gwin_gr->redraw();
	}
}

void ControlPanel::cb_btn_saveall(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	crease* c = &(This->ppm.crs[0]);

	//	int fileioflg[30];
	//	enum{ PX2, XX2, TX2, NX2, D2, K2, PX3, XX3, TX3, NX3, BX3, D3, K3, T3, PA, ALPHA, ALPHA0, DA, PB, BETA, BETA0, R3, R2 };
	memset(c->fileioflg, 1, sizeof(int) * 30);
	c->fileioflg[crease::PX2] = 0;
	c->fileioflg[crease::D2] = 0;
	c->fileioflg[crease::PX3] = 0;
	c->fileioflg[crease::D3] = 0;
	c->fileioflg[crease::PA] = 0;
	c->fileioflg[crease::DA] = 0;
	c->fileioflg[crease::PB] = 0;

	c->dumpAll("./output/result.csv");
}

