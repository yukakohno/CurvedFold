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

				strcpy(fname_m2m3, fname);
				for (int i = strlen(fname_m2m3) - 1; i > 0; i--) {
					if (fname_m2m3[i] == '.') {
						fname_m2m3[i] = 0;
						break;
					}
				}
				sprintf( fname_m2m3, "%s_m2m3.txt", fname_m2m3);

				// load data
				int filetype = This->ppm.crs[0].load(fname);
				int ret0 = This->ppm.crs[0].loadm2m3(fname_m2m3);
				//for( int i=0; i<16; i++ ){ This->gwin->mObject[i] = This->ppm.crs[0].m3[i]; }
				//unit_m44( This->ppm.crs[0].m3 );

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
	int divnum = This->gwin->divnum = This->vs_divnum->value();
	int frm = This->vs_fmot2->value() + c0->FM_fidx_org;
	int flg_usecp = This->cb_usecp->value();
	char dirname[128], fname[128];
    struct stat statBuf;
	time_t t0 = time(NULL);
	struct tm *t1 = localtime(&t0);

	sprintf( dirname, "output/%02d%02d%02d_f%02d/", t1->tm_hour, t1->tm_min, t1->tm_sec, frm );
	if( stat(dirname, &statBuf)==0 || _mkdir( dirname )==0 ){

		sprintf(fname, "%sm2m3.txt", dirname );	c0->dumpm2m3( fname );
		//c0->dumpm2m3("output/m2m3.txt", This->gwin->mObject );
		//c0->dumpP_("output/kta_current.txt", mode);
		sprintf(fname, "%sP.txt", dirname );	c0->dumpP( fname, mode);
		sprintf(fname, "%sX.csv", dirname );	c0->dumpX( fname );
		//c0->checkCP3D("output/checkCP3D.csv");
		//sprintf(fname, "%smotion.txt", dirname );	c0->dumpMotion( fname, flg_usecp );
		if( flg_usecp ){
			sprintf( fname, "%smotion_CP.txt", dirname );
		} else {
			sprintf( fname, "%smotion_ALL.txt", dirname );
		}
		c0->dumpMotion( fname, flg_usecp );
		sprintf(fname, "%sCP00.svg", dirname );	ppm->dumpsvg00( fname, divnum );
		sprintf(fname, "%sCP1.svg", dirname );	ppm->dumpsvg1( fname );
		ppm->makeObj3( 10.0 );
		sprintf(fname, "%sCP.obj", dirname );	ppm->dumpObj( fname, 0.4);

		//ppm->check180( "output/d180.csv" );	// 頂点周りの角度
		//ppm->checkquadplane( "output/quad.csv" );	// quad の平面度
		sprintf(fname, "%sd180_%02d.csv", dirname, frm );	c0->check180( fname );	// 頂点周りの角度
		sprintf(fname, "%squad_%02d.csv", dirname, frm );	c0->checkquadplane( fname );	// quad の平面度
		sprintf(fname, "%sgap_%02d.csv", dirname, frm );	ppm->checkGap( fname, divnum );	// ピース間の隙間

		sprintf(fname, "%smodel.bmp", dirname );	This->gwin->exportImage( fname );
	}
	//This->gwin_gr->redraw();
}

void ControlPanel::cb_btn_savescreen( Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c0 = &(This->ppm.crs[0]);
	int frm = This->vs_fmot2->value() + c0->FM_fidx_org;
	time_t t0 = time(NULL);
	struct tm *t1 = localtime(&t0);
	char fname[128];

	sprintf(fname, "output/model_%02d%02d%02d_f%02d.bmp", t1->tm_hour, t1->tm_min, t1->tm_sec, frm );
	This->gwin->exportImage( fname );
}

void ControlPanel::cb_btn_saveerr( Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);
	crease *c0 = &(This->ppm.crs[0]);
	int divnum = This->gwin->divnum = This->vs_divnum->value();
	double errdata[ MAX_STP_CNT*21 ]; // max(MAX_SPCNT,MAX_STP_CNT)
	double err180[ MAX_SPCNT*21 ];
	double errquad[ MAX_SPCNT*21*2 ];
	double errgap[ MAX_STP_CNT*21 ];
	int edsize = MAX_STP_CNT*21;
	int cntgap[21];
	int ret=0;

	for( int j=0; j<=20; j++ ){
		printf( "errors on frame %d...\n", j );
		This->vs_fmot2->value( j-c0->FM_fidx_org );
		This->vs_fmot2->do_callback();
		
		memset( errdata, 0, sizeof(double)*edsize );
		ret = c0->check180( errdata, c0->Xeidx+1, 9 );
		for( int i=c0->Xsidx; i<=c0->Xeidx; i++ ){
			err180[ MAX_SPCNT*j+i ] = errdata[i*9+8];	// atotal-2pi
		}

		memset( errdata, 0, sizeof(double)*edsize );
		ret = c0->checkquadplane( errdata, c0->Xeidx+1, 17 );
		for( int i=c0->Xsidx; i<c0->Xeidx; i++ ){
			errquad[ (MAX_SPCNT*j*2  ) +i ] = errdata[i*17+9];	// ipL
			errquad[ (MAX_SPCNT*j*2+1) +i ] = errdata[i*17+16];	// ipR
		}

		memset( errdata, 0, sizeof(double)*edsize );
		cntgap[j] = ppm->checkGap( errdata, MAX_STP_CNT, 7, divnum );
		for( int i=0; i<cntgap[j]; i++ ){
			errgap[ MAX_STP_CNT*j +i ] = errdata[i*7];	// dif
		}
	}

	FILE *fp = fopen( "output/d180_all.csv", "w" );
	if( fp ){
		int cnt=0;
		double ave=0.0, max=0.0;
		for( int j=0; j<=20; j++ ){
			for( int i=c0->Xsidx; i<=c0->Xeidx; i++ ){
				double err = fabs(err180[MAX_SPCNT*j+i]);
				ave += err;
				max = err > max ? err : max;
				cnt++;
			}
		}
		ave/=(double)cnt;
		fprintf( fp, "ave,%f\n", ave );
		fprintf( fp, "max,%f\n", max );
#if 0
		fprintf( fp, "frame,data\n" );
		for( int j=0; j<=20; j++ ){
			fprintf( fp, "%d,", j );
			for( int i=c0->Xsidx; i<=c0->Xeidx; i++ ){
				fprintf( fp, "%f,", err180[MAX_SPCNT*j+i] );
			}
			fprintf( fp, "\n" );
		}
#else
		for( int j=0; j<=20; j++ ){
			fprintf( fp, ",%d", j );
		}
		fprintf( fp, "\n" );
		for( int i=c0->Xsidx; i<=c0->Xeidx; i++ ){
			fprintf( fp, "%d,", i );
			//for( int j=0; j<=20; j++ ){
			for( int j=20; j>=0; j-- ){
				fprintf( fp, "%f,", err180[ MAX_SPCNT*j+i ] );
			}
			fprintf( fp, "\n" );
		}
#endif
		fclose( fp ); fp=NULL;
	}

	fp = fopen( "output/quad_all.csv", "w" );
	if( fp ){
		int cnt=0;
		double ave=0.0, max=0.0;
		for( int j=0; j<=20; j++ ){
			for( int i=c0->Xsidx; i<c0->Xeidx; i++ ){
				double errL = fabs(errquad[(MAX_SPCNT*j*2  )+i]);
				double errR = fabs(errquad[(MAX_SPCNT*j*2+1)+i]);
				ave += errL + errR;
				max = errL > max ? errL : max;
				max = errR > max ? errR : max;
				cnt+=2;
			}
		}
		ave/=(double)cnt;
		fprintf( fp, "ave,%f\n", ave );
		fprintf( fp, "max,%f\n", max );
#if 0
		fprintf( fp, "frame,data\n" );
		for( int j=0; j<=20; j++ ){
			fprintf( fp, "%dL,", j );
			for( int i=c0->Xsidx; i<c0->Xeidx; i++ ){
				fprintf( fp, "%f,", errquad[(MAX_SPCNT*j*2  )+i] );
			}
			fprintf( fp, "\n" );
			fprintf( fp, "%dR,", j );
			for( int i=c0->Xsidx; i<c0->Xeidx; i++ ){
				fprintf( fp, "%f,", errquad[(MAX_SPCNT*j*2+1)+i] );
			}
			fprintf( fp, "\n" );
		}
#else
		for( int j=0; j<=20; j++ ){
			fprintf( fp, ",%dL,%dR", j, j );
		}
		fprintf( fp, "\n" );
		for( int i=c0->Xsidx; i<c0->Xeidx; i++ ){
			fprintf( fp, "%d,", i );
			//for( int j=0; j<=20; j++ ){
			for( int j=20; j>=0; j-- ){
				fprintf( fp, "%f,%f,", errquad[ (MAX_SPCNT*j*2  ) +i ], errquad[ (MAX_SPCNT*j*2+1) +i ] );
			}
			fprintf( fp, "\n" );
		}
#endif
		fclose( fp ); fp=NULL;
	}

	fp = fopen( "output/gap_all.csv", "w" );
	if( fp ){
#if 1
		fprintf( fp, "frame,ave,max,data\n" );
		//for( int j=0; j<=20; j++ ){
		for( int j=20; j>=0; j-- ){
			int cnt=0;
			double ave=0.0, max=0.0;
			for( int i=0; i<cntgap[j]; i++ ){
				double gap=errgap[MAX_STP_CNT*j+i];
				if( gap<0.0 ){
					continue;
				}
				ave += gap;
				max = gap > max ? gap : max;
				cnt++;
			}
			ave/=(double)cnt;

			fprintf( fp, "%d,%f,%f,", 20-j, ave, max );
			for( int i=0; i<cntgap[j]; i++ ){
				double gap=errgap[MAX_STP_CNT*j+i];
				if( gap<0.0 ){
					fprintf( fp, "," );
				} else{
					fprintf( fp, "%f,", gap );
				}
			}
			fprintf( fp, "\n" );
		}
#else
		for( int j=0; j<=20; j++ ){
			fprintf( fp, ",%d", j );
		}
		fprintf( fp, "\n" );
		for( int i=0; i<=MAX_STP_CNT; i++ ){
			fprintf( fp, "%d,", i );
			for( int j=0; j<=20; j++ ){
				if( i<cntgap[j] ){
					fprintf( fp, "%f,", errgap[ MAX_STP_CNT*j+i ] );
				} else {
					fprintf( fp, "," );
				}
			}
			fprintf( fp, "\n" );
		}
#endif
		fclose( fp ); fp=NULL;
	}

	//This->gwin_gr->redraw();
}

void ControlPanel::cb_btn_savemotion(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int flg_usecp = This->cb_usecp->value();
	if( flg_usecp ){
		c->dumpMotion("output/motion_CP.txt", flg_usecp);
	} else {
		c->dumpMotion("output/motion_ALL.txt", flg_usecp);
	}
}

void ControlPanel::cb_btn_loadmotion(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	char fname[256];
	int flg_usecp = This->cb_usecp->value();
	if( flg_usecp ){
		sprintf( fname, "%smotion_CP.txt", This->filepath );
	} else {
		sprintf( fname, "%smotion_ALL.txt", This->filepath );
	}
	//int flg_usecp = c->loadMotion("input/motion.txt");
	flg_usecp = c->loadMotion( fname );
	This->cb_usecp->value( flg_usecp );
	c->updateFM( flg_usecp );
	int frm = This->vs_fmot2->value() + c->FM_fidx_org;
	This->update_bx_FM(c, frm);
}

void ControlPanel::cb_btn_saveframe(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int frm = This->vs_fmot2->value() + c->FM_fidx_org;
	int flg_usecp = This->cb_usecp->value();
	c->dumpMotionFrame( frm, "output/motionframe.txt", flg_usecp );
}

void ControlPanel::cb_btn_loadframe(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	crease *c = &(This->ppm.crs[0]);
	int frm = This->vs_fmot2->value() + c->FM_fidx_org;
	int flg_usecp = c->loadMotionFrame( frm, "output/motionframe.txt" );
	c->updateFM( 1 );
	This->update_bx_FM(c, frm);
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

