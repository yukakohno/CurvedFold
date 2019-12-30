#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Windows.h>
#include "ControlPanel.h"

ControlPanel::ControlPanel(int X, int Y, int W, int H, GraphWindow3DCF *_gwin )
: Fl_Window(X, Y, W, H)
{
	ControlPanel(X, Y, W, H, "", _gwin );
}

ControlPanel::ControlPanel(int X, int Y, int W, int H, const char *L, GraphWindow3DCF *_gwin )
: Fl_Window(X, Y, W, H, L)
{
	gwin = _gwin;

	kmin = -0.05;
	kmax = 0.05;
	kstep = 0.0005;
	kstepk = 0.01;

	tmin = -0.1;
	tmax = 0.1;
	tstep = 0.0005;
	tstepk = 0.01;

	amin = 0.0;
	amax = 90.0;
	astep = 0.25;
	astepk = 3.0;
	an_fcnt = 20;

	createPanel();
}

ControlPanel::~ControlPanel()
{
	gwin = NULL;
}

int ControlPanel::value_grpfix()
{
	int fix = -1;
	for( int i=0; i<N_C_MODE; i++ ){
		if( this->rb_fix[i]->value() > 0 ){
			fix = i; // CMODE_A, CMODE_B, CMODE_C
			break;
		}
	}
	return fix;
}

int ControlPanel::value_grpparam()
{
	int prm = -1;
	for( int i=0; i<N_P_IDX; i++ ){
		if( this->rb_param[i]->value() > 0 ){
			prm = i; // P_CV2D, P_CV3D, P_TRSN, P_FLDA, P_RULL, P_RULR
			break;
		}
	}
	return prm;
}

void ControlPanel::refresh(int init)
{
	// A: 3D‹Èü‚ÆÜ‚èŠp“x‚©‚çAÜ‚èü‚ð‹‚ß‚é
	// B: Ü‚èü‚ÆÜ‚èŠp“x‚©‚çA3D‹Èü‚ð‹‚ß‚é
	// C: Ü‚èü‚Æ3D‹Èü‚©‚çAÜ‚èŠp“x‚ð‹‚ß‚é
	int mode = this->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C, MODE_T, MODE_FA, MODE_TFA
	int prm = this->value_grpparam();	// P_CV2D, P_CV3D, P_TRSN, P_FLDA, P_RULL, P_RULR
	papermodel *_ppm = &(this->ppm);
	crease *_c = &(_ppm->crs[0]);

	switch( mode ){
		case CMODE_A:	// A: 3D‹Èü‚ÆÜ‚èŠp“x‚©‚çAÜ‚èü‚ð‹‚ß‚é
			if( prm==P_CV3D || prm==P_TRSN || prm==P_FLDA ){ // curvature, torsion, angle
				_c->calcXA_CP( _ppm->flg_interpolate, &_ppm->rp );
				_ppm->postproc();
			}
			break;
		case CMODE_B:	// B: Ü‚èü‚ÆÜ‚èŠp“x‚©‚çA3D‹Èü‚ð‹‚ß‚é
			if( prm==P_CV2D || prm==P_FLDA ){ // curv2d, angle
				_c->calcCPA_X( _ppm->flg_interpolate, &_ppm->rp );
				_ppm->postproc();
			}
			break;
		case CMODE_C:	// C: Ü‚èü‚Æ3D‹Èü‚©‚çAÜ‚èŠp“x‚ð‹‚ß‚é
			if( prm==P_CV2D || prm==P_CV3D || prm==P_TRSN ){ // curv2d, curvature, torsion
				_c->calcCPX_A( _ppm->flg_interpolate, &_ppm->rp );
				_ppm->postproc();
			}
			break;
	}
	this->gwin->ppm = &(this->ppm);
	this->gwin->redraw();
	this->gwin_cp->redraw();
	this->gwin_gr->redraw();
	if( init ){
		this->setorg();
	}
}

void ControlPanel::setorg()
{
	crease *c = &(this->ppm.crs[0]);
	this->vs_fmot->value(0);
	memcpy( this->k2d_org, c->k2d, sizeof(double)*c->Xcnt );
	memcpy( this->kv_org, c->kv, sizeof(double)*c->Xcnt );
	memcpy( this->tr_org, c->tr, sizeof(double)*c->Xcnt );
	memcpy( this->alpha_org, c->alpha, sizeof(double)*c->Xcnt );
	for( int i=0; i<c->Xcnt; i++ ){
		this->k2d_dif0[i] = 0;
		this->kv_dif0[i] = 0;
		this->tr_dif0[i] = c->tr[i]/(double)this->an_fcnt;
		this->alpha_dif0[i] = c->alpha[i]/(double)this->an_fcnt;
		this->k2d_dif1[i] = 0;
		this->kv_dif1[i] = 0;
		this->tr_dif1[i] = 0;
		this->alpha_dif1[i] = ( M_PI*0.5 - c->alpha[i] )/(double)this->an_fcnt;
	}
}

void ControlPanel::createPanel()
{
	int wgt_x=0, wgt_y=0;

	Fl_Tabs* tab = new Fl_Tabs(0, 0, this->w(), this->h());
	{
		Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h()-20, "IO/DSP");
		{ 
			wgt_x = 10;
			wgt_y = 20;

			// ------------------------- fILE_IO -------------------------------------------
			Fl_Box *bx_FILE_IO= new Fl_Box(0, wgt_y, g->w(), 20, "--- FILE_IO ---");

			wgt_y += 20;

			btn_load = new Fl_Button(wgt_x, wgt_y, 45, 20, "load");
			btn_load->callback(cb_btn_load, (void*)this);
			fc = new Fl_File_Chooser(RECDIR, "(*.{txt})", Fl_File_Chooser::SINGLE, "Fl_File_Chooser");
			g->add(btn_load);

			btn_loadtex = new Fl_Button(wgt_x+45, wgt_y, 45, 20, "ld_tex");
			//btn_loadtex->callback(cb_btn_loadtex, (void*)this);
			g->add(btn_loadtex);

			btn_dump = new Fl_Button(wgt_x+90, wgt_y, 45, 20, "dump");
			btn_dump->callback(cb_btn_dump, (void*)this);
			g->add(btn_dump);

			wgt_x = 10;
			wgt_y += 25;

			// ------------------------- DISPLAY -------------------------------------------
			Fl_Box *bx_DISPLAY= new Fl_Box(0, wgt_y, g->w(), 20, "--- DISPLAY ---");
			wgt_y += 20;

			btn_resetRot = new Fl_Button(wgt_x, wgt_y, 45, 20, "Y-DN");
			btn_resetRot->callback(cb_rb_gwin, (void*)this);
			btn_resetScale = new Fl_Button(wgt_x+45, wgt_y, 45, 20, "SCL");
			btn_resetScale->callback(cb_rb_gwin, (void*)this);
			btn_resetTrans = new Fl_Button(wgt_x+90, wgt_y, 45, 20, "TRNS");
			btn_resetTrans->callback(cb_rb_gwin, (void*)this);
			g->add(btn_resetRot);
			g->add(btn_resetScale);
			g->add(btn_resetTrans);

			// axis •\Ž¦^”ñ•\Ž¦
			cb_dispaxis = new Fl_Check_Button(wgt_x+135, wgt_y, 45, 20, "axis");
			cb_dispaxis->value(this->gwin->disp_axis);
			cb_dispaxis->callback(cb_cb_dispaxis, (void*)this);
			g->add(cb_dispaxis);

			wgt_x = 10;
			wgt_y += 20;

			// ’¸“_ •\Ž¦^”ñ•\Ž¦
			cb_disp[D_X] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "X");
			cb_disp[D_X]->value(this->gwin->disp_X);
			cb_disp[D_X]->callback(cb_cb_dispX, (void*)this);

			wgt_x += 45;

			// mesh •\Ž¦^”ñ•\Ž¦
			cb_disp[D_TNB] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "TNB");
			cb_disp[D_TNB]->value(this->gwin->disp_TNB);
			cb_disp[D_TNB]->callback(cb_cb_dispTNB, (void*)this);

			wgt_x += 45;

			// face •\Ž¦^”ñ•\Ž¦
			cb_disp[D_R] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "R");
			cb_disp[D_R]->value(this->gwin->disp_R);
			cb_disp[D_R]->callback(cb_cb_dispR, (void*)this);

			wgt_x += 45;

			cb_disp[D_RLEN] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "R0");
			cb_disp[D_RLEN]->value(this->gwin->disp_maxrlen);
			cb_disp[D_RLEN]->callback(cb_cb_dispRlen, (void*)this);

			g->add(cb_disp[D_X]);
			g->add(cb_disp[D_TNB]);
			g->add(cb_disp[D_R]);
			g->add(cb_disp[D_RLEN]);

			wgt_x = 10;	wgt_y += 20;

			cb_disp[D_PLY] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "PLY");
			cb_disp[D_PLY]->value(this->gwin->disp_PLY);
			cb_disp[D_PLY]->callback(cb_cb_dispPLY, (void*)this);

			wgt_x += 45;

			cb_disp[D_PTN] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "PTN");
			cb_disp[D_PTN]->value(this->gwin->disp_PTN);
			cb_disp[D_PTN]->callback(cb_cb_dispPTN, (void*)this);

			wgt_x += 45;

			cb_disp[D_CP] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "CT_P");
			cb_disp[D_CP]->value(this->gwin->disp_CP);
			cb_disp[D_CP]->callback(cb_cb_dispCP, (void*)this);

			wgt_x += 45;

			cb_CurveEnd = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "CE");
			if( ppm.flg_postproc_type == PPTYPE_X ){
				cb_CurveEnd->value(0);
			} else {
				cb_CurveEnd->value(1);
			}
			cb_CurveEnd->callback(cb_cb_CurveEnd, (void*)this);

			g->add(cb_disp[D_PLY]);
			g->add(cb_disp[D_PTN]);
			g->add(cb_disp[D_CP]);
			g->add(cb_CurveEnd);

			wgt_x = 10;	wgt_y += 25;
		}
		g->end();
	}
	{
		Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h()-20, "DSGN1");
		g->hide();
		{ 
			// ------------------------- DESIGN -------------------------------------------

			wgt_x = 10;	wgt_y = 20;
			int y_space = 25, grp_num=0;

			Fl_Box *bx_MODE= new Fl_Box(0, wgt_y, g->w(), 20, "--- CONTROL MODE ---");	wgt_y+=20;

			grp_fix = new Fl_Group(0, wgt_y, g->w(), 3*y_space);	grp_num = 0;
			rb_fix[CMODE_A] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "3D + ANGLE -> 2D");	wgt_y+=y_space; grp_num++;
			rb_fix[CMODE_B] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "2D + ANGLE -> 3D");	wgt_y+=y_space; grp_num++;
			rb_fix[CMODE_C] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "2D + 3D -> ANGLE");	wgt_y+=y_space; grp_num++;
#if 0
			rb_fix[CMODE_T] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "R-T");	wgt_y+=y_space; grp_num++;
			rb_fix[CMODE_FA] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "R-A");	wgt_y+=y_space; grp_num++;
			rb_fix[CMODE_TFA] = new Fl_Round_Button(wgt_x0, wgt_y, 60, 20, "R-TA");	wgt_y+=y_space; grp_num++;
#endif
			for( int i=0; i<grp_num; i++ ){
				rb_fix[i]->type(FL_RADIO_BUTTON);
				grp_fix->add(rb_fix[i]);
			}
			grp_fix->end();
			rb_fix[CMODE_A]->setonly();
			g->add(grp_fix);

			Fl_Box *bx_PARAMETER= new Fl_Box(0, wgt_y, g->w(), 20, "--- PARAMETER ---");	wgt_y+=20;

			grp_param = new Fl_Group(0, wgt_y, g->w(), 4*y_space);	grp_num=0;
			rb_param[P_CV2D] = new Fl_Round_Button(wgt_x, wgt_y, 45, 20, "CURVATURE (2D)");	wgt_y+=y_space; grp_num++;
			rb_param[P_CV3D] = new Fl_Round_Button(wgt_x, wgt_y, 45, 20, "CURVATURE (3D)");	wgt_y+=y_space; grp_num++;
			rb_param[P_TRSN] = new Fl_Round_Button(wgt_x, wgt_y, 45, 20, "TORSION");			wgt_y+=y_space; grp_num++;
			rb_param[P_FLDA] = new Fl_Round_Button(wgt_x, wgt_y, 45, 20, "FOLDING ANGLE");	wgt_y+=y_space; grp_num++;
			for( int i=0; i<grp_num; i++ ){
				rb_param[i]->type(FL_RADIO_BUTTON);
				rb_param[i]->callback(cb_vs_ppos, (void*)this);
				grp_param->add(rb_param[i]);
			}
			grp_param->end();
			rb_param[P_CV3D]->setonly();
			g->add(grp_param);

			Fl_Box *bx_CPINDEX= new Fl_Box(0, wgt_y, g->w(), 20, "--- CONTROL POINT INDEX ---");	wgt_y+=20;

			vs_ppos = new Fl_Value_Slider(wgt_x, wgt_y, 180, 20);
			vs_ppos->bounds(0,1);	vs_ppos->step(1);	vs_ppos->value(gwin->ppos);
			vs_ppos->align(FL_ALIGN_LEFT);
			vs_ppos->type(FL_HORIZONTAL);
			vs_ppos->callback(cb_vs_ppos, (void*)this);
			g->add( vs_ppos );

			wgt_y+=y_space;

			Fl_Box *bx_PARAMVAL= new Fl_Box(0, wgt_y, g->w(), 20, "--- PARAMETER VALUE ---");	wgt_y+=20;

			vs_pval = new Fl_Value_Slider(wgt_x, wgt_y, 180, 20);
			vs_pval->bounds(1,0);	vs_pval->step(0.01);	vs_pval->value(gwin->pprm);
			vs_pval->align(FL_ALIGN_LEFT);
			vs_pval->type(FL_HORIZONTAL);
			vs_pval->callback(cb_vs_pval, (void*)this);
			g->add( vs_pval );

			wgt_y+=y_space;

			// ------------------------- FOLD_MOTION -------------------------------------------
			Fl_Box *bx_FMOTION= new Fl_Box(0, wgt_y, g->w(), 20, "--- FOLDING_MOTION ---");	wgt_y+=20;
			vs_fmot = new Fl_Value_Slider(wgt_x, wgt_y, 180, 20);
			vs_fmot->bounds(-an_fcnt+1,an_fcnt-1);	vs_fmot->step(1);	vs_fmot->value(0);
			vs_fmot->align(FL_ALIGN_LEFT);
			vs_fmot->type(FL_HORIZONTAL);
			vs_fmot->callback(cb_vs_fmot, (void*)this);
			g->add( vs_fmot );

			wgt_y+=y_space;

			// ------------------------- ADD CREASE & TRIM ------------------------------------

			Fl_Box *bx_TRIM= new Fl_Box(0, wgt_y, g->w(), 20, "--- ADD FOLD / TRIM ---");

			wgt_y+=y_space;

			cb_foldcurve = new Fl_Check_Button(wgt_x, wgt_y, 50, 20, "FOLD");
			cb_foldcurve->value( (int)(this->gwin->flg_addcurve == 2) );
			cb_foldcurve->callback(cb_cb_addcurve, (void*)this);
			cb_foldcurve->deactivate();

			wgt_y+=y_space;

			cb_trimcurve = new Fl_Check_Button(wgt_x, wgt_y, 50, 20, "TRIM");
			cb_trimcurve->value( (int)(this->gwin->flg_addcurve == 1) );
			cb_trimcurve->callback(cb_cb_addcurve, (void*)this);

			wgt_y -= y_space;
			wgt_x = 80;

			btn_proccurve = new Fl_Button(wgt_x, wgt_y, 100, 20, "ADD CREASE");
			btn_proccurve->callback(cb_btn_proccurve, (void*)this);

			g->add(cb_trimcurve);
			g->add(cb_foldcurve);
			g->add(btn_proccurve);
		}
		g->end();
	}
	{
		Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h()-20, "DSGN2");
		g->hide();
		{
			// ------------------------- RECTIFY ------------------------------------

			wgt_x = 10;	wgt_y = 20;
			int y_space = 25, grp_num=0;
			
			Fl_Box *bx_RECTIFY= new Fl_Box(0, wgt_y, g->w(), 20, "--- RECTIFY ---");	wgt_y+=20;
			cb_rectifyA = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "FOLDING ANGLE");	wgt_y+=y_space;
			cb_rectifyA->value( this->ppm.rp.flg_rectifyA );
			cb_rectifyA->callback(cb_cb_rectifyA, (void*)this );
			cb_rectifyT = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "TORSION");	wgt_y+=y_space;
			cb_rectifyT->value( this->ppm.rp.flg_rectifyT );
			cb_rectifyT->callback(cb_cb_rectifyT, (void*)this );
			cb_rectifyR = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "Eliminate RULINGS");	wgt_y+=y_space;
			cb_rectifyR->value( this->ppm.rp.flg_rectifyR );
			cb_rectifyR->callback(cb_cb_rectifyR, (void*)this );
			g->add(cb_rectifyA);
			g->add(cb_rectifyT);
			g->add(cb_rectifyR);		
		
			btn_updateChoice = new Fl_Button(wgt_x, wgt_y, 90, 20, "update prm");
			btn_updateChoice->callback( cb_btn_updateChoice, (void*)this );
			g->add(btn_updateChoice);

			btn_setSEParam = new Fl_Button(wgt_x+90, wgt_y, 90, 20, "set prm");
			btn_setSEParam->callback( cb_btn_setSEParam, (void*)this );
			g->add(btn_setSEParam);

			wgt_y+=y_space;

			Fl_Box *bx_FANGLE= new Fl_Box(0, wgt_y, g->w(), 20, "--- FOLDING ANGLE PARAMS ---");	wgt_y+=20;

			ch_A_method = new Fl_Choice(wgt_x+50, wgt_y, 50, 20, "method");
			ch_A_method->add("Linear");
			ch_A_method->add("Bezier");
			ch_A_method->value( this->ppm.rp.rectT.method-1 );
			ch_A_method->callback( cb_ch_A_method, (void*)this );

			vi_A_kvthres = new Fl_Value_Input(wgt_x+140, wgt_y, 40, 20, "thres");
			vi_A_kvthres->value( this->ppm.rp.rectT.kvthres );
			vi_A_kvthres->callback( cb_vi_A_kvthres, (void*)this );

			wgt_y+=y_space;

			vi_A_s1mgn = new Fl_Value_Input(wgt_x+50, wgt_y, 40, 20, "s1 mgn");
			vi_A_s1mgn->value( this->ppm.rp.rectT.s1mgn );
			vi_A_s1mgn->callback( cb_vi_A_s1mgn, (void*)this );
			vi_A_s2mgn = new Fl_Value_Input(wgt_x+140, wgt_y, 40, 20, "s2 mgn"); 
			vi_A_s2mgn->value( this->ppm.rp.rectT.s2mgn );
			vi_A_s2mgn->callback( cb_vi_A_s2mgn, (void*)this );

			wgt_y+=y_space;

			ch_A_se = new Fl_Choice(wgt_x+50, wgt_y, 50, 20, "se idx");
			ch_A_se->clear();
			ch_A_se->callback( cb_ch_A_se, (void*)this );

			ch_A_src = new Fl_Choice(wgt_x+140, wgt_y, 50, 20, "src");
			ch_A_src->add("0");
			ch_A_src->add("1");
			ch_A_src->add("2");
			ch_A_src->value(0);
			ch_A_src->callback( cb_ch_A_src, (void*)this );

			wgt_y+=y_space;

			vi_A_s2 = new Fl_Value_Input(wgt_x+40, wgt_y, 35, 20, "s1e2");	vi_A_s2->callback( cb_vi_A_s2, (void*)this );
			vi_A_s1 = new Fl_Value_Input(wgt_x+75, wgt_y, 35, 20 );			vi_A_s1->callback( cb_vi_A_s1, (void*)this );
			vi_A_e1 = new Fl_Value_Input(wgt_x+110, wgt_y, 35, 20 );		vi_A_e1->callback( cb_vi_A_e1, (void*)this );
			vi_A_e2 = new Fl_Value_Input(wgt_x+145, wgt_y, 35, 20 );		vi_A_e2->callback( cb_vi_A_e2, (void*)this );

			wgt_y+=y_space;

			vi_A_sv2 = new Fl_Value_Input(wgt_x+40, wgt_y, 35, 20, "veclen");	vi_A_sv2->callback( cb_vi_A_sv2, (void*)this );
			vi_A_sv1 = new Fl_Value_Input(wgt_x+75, wgt_y, 35, 20 );		vi_A_sv1->callback( cb_vi_A_sv1, (void*)this );
			vi_A_ev1 = new Fl_Value_Input(wgt_x+110, wgt_y, 35, 20 );		vi_A_ev1->callback( cb_vi_A_ev1, (void*)this );
			vi_A_ev2 = new Fl_Value_Input(wgt_x+145, wgt_y, 35, 20 );		vi_A_ev2->callback( cb_vi_A_ev2, (void*)this );

			wgt_y+=y_space;

			g->add(ch_A_method);
			g->add(vi_A_kvthres);
			g->add(vi_A_s1mgn);	g->add(vi_A_s2mgn);
			g->add(ch_A_se);	g->add(ch_A_src);
			g->add(vi_A_s2);	g->add(vi_A_s1);	g->add(vi_A_e1);	g->add(vi_A_e2);
			g->add(vi_A_sv2);	g->add(vi_A_sv1);	g->add(vi_A_ev1);	g->add(vi_A_ev2);

			Fl_Box *bx_TORSION= new Fl_Box(0, wgt_y, g->w(), 20, "--- TORSION PARAMS ---");	wgt_y+=20;

			ch_T_method = new Fl_Choice(wgt_x+50, wgt_y, 50, 20, "method");
			ch_T_method->add("Linear");
			ch_T_method->add("Bezier");
			ch_T_method->value( this->ppm.rp.rectT.method-1 );
			ch_T_method->callback( cb_ch_T_method, (void*)this );

			vi_T_kvthres = new Fl_Value_Input(wgt_x+140, wgt_y, 40, 20, "thres");
			vi_T_kvthres->value( this->ppm.rp.rectT.kvthres );
			vi_T_kvthres->callback( cb_vi_T_kvthres, (void*)this );

			wgt_y+=y_space;

			vi_T_s1mgn = new Fl_Value_Input(wgt_x+50, wgt_y, 40, 20, "s1 mgn");
			vi_T_s1mgn->value( this->ppm.rp.rectT.s1mgn );
			vi_T_s1mgn->callback( cb_vi_T_s1mgn, (void*)this );
			vi_T_s2mgn = new Fl_Value_Input(wgt_x+140, wgt_y, 40, 20, "s2 mgn"); 
			vi_T_s2mgn->value( this->ppm.rp.rectT.s2mgn );
			vi_T_s2mgn->callback( cb_vi_T_s2mgn, (void*)this );

			wgt_y+=y_space;

			ch_T_se = new Fl_Choice(wgt_x+50, wgt_y, 50, 20, "se idx");
			ch_T_se->clear();
			ch_T_se->callback( cb_ch_T_se, (void*)this );

			ch_T_src = new Fl_Choice(wgt_x+140, wgt_y, 50, 20, "src");
			ch_T_src->add("0");
			ch_T_src->add("1");
			ch_T_src->add("2");
			ch_T_src->value(0);
			ch_T_src->callback( cb_ch_T_src, (void*)this );

			wgt_y+=y_space;

			vi_T_s2 = new Fl_Value_Input(wgt_x+40, wgt_y, 35, 20, "s1e2");	vi_T_s2->callback( cb_vi_T_s2, (void*)this );
			vi_T_s1 = new Fl_Value_Input(wgt_x+75, wgt_y, 35, 20 );			vi_T_s1->callback( cb_vi_T_s1, (void*)this );
			vi_T_e1 = new Fl_Value_Input(wgt_x+110, wgt_y, 35, 20 );		vi_T_e1->callback( cb_vi_T_e1, (void*)this );
			vi_T_e2 = new Fl_Value_Input(wgt_x+145, wgt_y, 35, 20 );		vi_T_e2->callback( cb_vi_T_e2, (void*)this );

			wgt_y+=y_space;
		
			vi_T_sv2 = new Fl_Value_Input(wgt_x+40, wgt_y, 35, 20, "veclen");	vi_T_sv2->callback( cb_vi_T_sv2, (void*)this );
			vi_T_sv1 = new Fl_Value_Input(wgt_x+75, wgt_y, 35, 20 );		vi_T_sv1->callback( cb_vi_T_sv1, (void*)this );
			vi_T_ev1 = new Fl_Value_Input(wgt_x+110, wgt_y, 35, 20 );		vi_T_ev1->callback( cb_vi_T_ev1, (void*)this );
			vi_T_ev2 = new Fl_Value_Input(wgt_x+145, wgt_y, 35, 20 );		vi_T_ev2->callback( cb_vi_T_ev2, (void*)this );

			wgt_y+=y_space;

			g->add(ch_T_method);
			g->add(vi_T_kvthres);
			g->add(vi_T_s1mgn);	g->add(vi_T_s2mgn);
			g->add(ch_T_se);	g->add(ch_T_src);
			g->add(vi_T_sv2);	g->add(vi_T_sv1);	g->add(vi_T_ev1);	g->add(vi_T_ev2);
			g->add(vi_T_s2);	g->add(vi_T_s1);	g->add(vi_T_e1);	g->add(vi_T_e2);

			Fl_Box *bx_RULING= new Fl_Box(0, wgt_y, g->w(), 20, "--- RULING PARAMS ---");	wgt_y+=20;

			vi_R_kvthres = new Fl_Value_Input(wgt_x+50, wgt_y, 40, 20, "thres");
			vi_R_kvthres->value( this->ppm.rp.rectifyR_kvthres );
			vi_R_kvthres->callback( cb_vi_R_kvthres, (void*)this );

			g->add(vi_R_kvthres);
		}
		g->end();
	}	
	tab->end();
}

int ControlPanel::handle(int event)
{
	int key, pos=-1, prm=-1;
	crease *c = &(this->ppm.crs[0]);

	switch(event){
	case FL_FOCUS:
		//return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_UNFOCUS:
		//return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_KEYDOWN:
		//printf("FL_KEYDOWN\n");
		key = Fl::event_key();
		pos = this->vs_ppos->value();
		prm = this->value_grpparam();// P_CV2D, P_CV3D, P_TRSN, P_FLDA, P_RULL, P_RULR

		if(key=='u'){
			this->ppm.set_postproc_type( PPTYPE_PRICURVE );
			switch( prm ){
	case P_CV2D: // curv2d
		c->Px2d[pos] += this->kstepk;
		this->refresh(1);
		this->vs_pval->value( c->Px2d[pos] );
		break;
	case P_CV3D: // curvature
		c->Px[pos] += this->kstepk;
		this->refresh(1);
		this->vs_pval->value( c->Px[pos] );
		break;
	case P_TRSN: // torsion
		c->Py[pos] += this->tstepk;
		this->refresh(1);
		this->vs_pval->value( c->Py[pos] );
		break;
	case P_FLDA: // angle
		c->Pa[pos] += this->astepk*M_PI/180.0;
		this->refresh(1);
		this->vs_pval->value( c->Pa[pos] *180.0/M_PI );
		break;
			}
			this->ppm.set_postproc_type( PPTYPE_UNDEF );
		} else if(key=='d'){
			this->ppm.set_postproc_type( PPTYPE_PRICURVE );
			switch( prm ){
	case P_CV2D: // curv2d
		c->Px2d[pos] -= this->kstepk;
		this->refresh(1);
		this->vs_pval->value( c->Px2d[pos] );
		break;
	case P_CV3D: // curvature
		c->Px[pos] -= this->kstepk;
		this->refresh(1);
		this->vs_pval->value( c->Px[pos] );
		break;
	case P_TRSN: // torsion
		c->Py[pos] -= this->tstepk;
		this->refresh(1);
		this->vs_pval->value( c->Py[pos] );
		break;
	case P_FLDA: // angle
		c->Pa[pos] -= this->astepk*M_PI/180.0;
		this->refresh(1);
		this->vs_pval->value( c->Pa[pos] *180.0/M_PI );
		break;
			}
			this->ppm.set_postproc_type( PPTYPE_UNDEF );
		} else if(key=='0'){
			this->ppm.set_postproc_type( PPTYPE_PRICURVE );
			switch( prm ){
	case P_CV2D: // curv2d
		c->Px2d[pos] = 0.0;
		this->refresh(1);
		this->vs_pval->value( c->Px2d[pos] );
		break;
	case P_CV3D: // curvature
		c->Px[pos] = 0.0;
		this->refresh(1);
		this->vs_pval->value( c->Px[pos] );
		break;
	case P_TRSN: // torsion
		c->Py[pos] = 0.0;
		this->refresh(1);
		this->vs_pval->value( c->Py[pos] );
		break;
	case P_FLDA: // angle
		c->Pa[pos] = M_PI/4.0;
		this->refresh(1);
		this->vs_pval->value( c->Pa[pos] );
		break;
			}
			this->ppm.set_postproc_type( PPTYPE_UNDEF );
		}
		break;
	default:
		//return Fl_Window::handle(event);
		break;
	}
	return Fl_Window::handle(event);
}
