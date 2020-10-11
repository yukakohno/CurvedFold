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

	fcnt = acnt = 0;
	acnt_inc = true;
	flg_idle_active = false;
	flg_idlerul_active = false;
	flg_idletgt_active = false;
	idlerul_idx = 0;
	idletgt_idx = 0;
	listtgcnt = 0;

	hist_head = hist_size = 0;

	sprintf( filepath, "./input/" );

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
			fix = i; // CMODE_A, CMODE_B, CMODE_C, CMODE_R
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
			prm = i; // P_CV2D, P_CV3D, P_TRSN, P_FLDA, P_RULL, P_RULR, P_TRSN1, P_FLDA1, P_RULL1, P_RULR1
			break;
		}
	}
	return prm;
}

void ControlPanel::cb_rb_fix0(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;

	if( This->rb_fix[CMODE_A]->value()>0 ){ // 3D + ANGLE -> 2D
		This->rb_param[P_CV2D]->value(0);	This->rb_param[P_CV2D]->deactivate();
		This->rb_param[P_RULL]->value(0);	This->rb_param[P_RULL]->deactivate();
		This->rb_param[P_RULR]->value(0);	This->rb_param[P_RULR]->deactivate();
		This->rb_param[P_CV3D]->activate();
		This->rb_param[P_TRSN]->activate();
		This->rb_param[P_FLDA]->activate();
		if( This->rb_param[P_CV3D]->value()==0
			&& This->rb_param[P_TRSN]->value()==0
			&& This->rb_param[P_FLDA]->value()==0 )
		{
			This->rb_param[P_CV3D]->setonly();
		}
	} else if( This->rb_fix[CMODE_B]->value()>0 ){ // 2D + ANGLE -> 3D
		This->rb_param[P_CV3D]->value(0);	This->rb_param[P_CV3D]->deactivate();
		//This->rb_param[P_TRSN]->value(0);	This->rb_param[P_TRSN]->deactivate();
		This->rb_param[P_RULL]->value(0);	This->rb_param[P_RULL]->deactivate();
		This->rb_param[P_RULR]->value(0);	This->rb_param[P_RULR]->deactivate();
		This->rb_param[P_CV2D]->activate();
		This->rb_param[P_TRSN]->activate();
		This->rb_param[P_FLDA]->activate();
		if( This->rb_param[P_CV2D]->value()==0
			&& This->rb_param[P_TRSN]->value()==0
			&& This->rb_param[P_FLDA]->value()==0 )
		{
			This->rb_param[P_CV2D]->setonly();
		}
	} else if( This->rb_fix[CMODE_C]->value()>0 ){ // 2D + 3D -> ANGLE
		This->rb_param[P_FLDA]->value(0);	This->rb_param[P_FLDA]->deactivate();
		This->rb_param[P_RULL]->value(0);	This->rb_param[P_RULL]->deactivate();
		This->rb_param[P_RULR]->value(0);	This->rb_param[P_RULR]->deactivate();
		This->rb_param[P_CV3D]->activate();
		This->rb_param[P_TRSN]->activate();
		This->rb_param[P_CV2D]->activate();
		if( This->rb_param[P_CV3D]->value()==0
			&& This->rb_param[P_TRSN]->value()==0
			&& This->rb_param[P_CV2D]->value()==0 )
		{
			This->rb_param[P_CV3D]->setonly();
		}
	} else if( This->rb_fix[CMODE_R]->value()>0 ){ // ruling
		This->rb_param[P_CV2D]->value(0);	This->rb_param[P_CV2D]->deactivate();
		This->rb_param[P_CV3D]->value(0);	This->rb_param[P_CV3D]->deactivate();
		This->rb_param[P_TRSN]->value(0);	This->rb_param[P_TRSN]->deactivate();
		This->rb_param[P_FLDA]->value(0);	This->rb_param[P_FLDA]->deactivate();
		This->rb_param[P_RULL]->activate();
		This->rb_param[P_RULR]->activate();
		if( This->rb_param[P_RULL]->value()==0 && This->rb_param[P_RULR]->value()==0 )
		{
			This->rb_param[P_RULL]->setonly();
		}
	}
	This->rb_param[P_TRSN1]->value(0);
	This->rb_param[P_FLDA1]->value(0);
	This->rb_param[P_RULL1]->value(0);
	This->rb_param[P_RULR1]->value(0);
}

void ControlPanel::cb_rb_param0(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	//This->rb_param[P_RULL]->value(0);
	//This->rb_param[P_RULR]->value(0);
	This->rb_param[P_TRSN1]->value(0);
	This->rb_param[P_FLDA1]->value(0);
	This->rb_param[P_RULL1]->value(0);
	This->rb_param[P_RULR1]->value(0);
	if( This->rb_fix[CMODE_A]->value()>0 || This->rb_fix[CMODE_B]->value()>0
		|| This->rb_fix[CMODE_C]->value()>0 || This->rb_fix[CMODE_R]->value()>0 )
	{
		This->vs_ppos->do_callback();
	}
}

void ControlPanel::cb_rb_param2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	if( This->rb_param[P_TRSN1]->value()>0 || This->rb_param[P_FLDA1]->value()>0
		|| This->rb_param[P_RULL1]->value()>0 || This->rb_param[P_RULR1]->value()>0 )
	{
		This->rb_param[P_CV2D]->value(0);
		This->rb_param[P_CV3D]->value(0);
		This->rb_param[P_TRSN]->value(0);
		This->rb_param[P_FLDA]->value(0);
		This->rb_param[P_RULL]->value(0);
		This->rb_param[P_RULR]->value(0);
		This->vs_ppos2->do_callback();
	}
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
			if( prm==P_CV2D || prm==P_FLDA || prm==P_TRSN ){ // curv2d, angle
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
		case CMODE_R:	// ruling angle -> torsion, alpha -> 3D‹ÈüAÜ‚èŠp“x
			if( prm==P_RULL || prm==P_RULR ){ // Ruling Left, Ruling Right
				_c->calcR_TA( 0/*flg_interpolate*/, &_ppm->rp, -1, -1, 2*M_PI, 2 );
				_ppm->postproc();
			}
			break;
	}
	this->push_hist(_c->Pcnt, mode, _c->Px2d, _c->Py2d, _c->Px, _c->Py, _c->Pz, _c->Pa, _c->Pbl, _c->Pbr, _c->m3);
	this->gwin->ppm = &(this->ppm);
	this->gwin->redraw();
	this->gwin_cp->redraw();
	this->gwin_gr->redraw();
	if( init ){
		this->setorg();
		_c->initFM();
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
		Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h()-20, "D0");
		g->hide();
		{ 
			wgt_x = 10;
			wgt_y = 20;

			// ------------------------- fILE_IO -------------------------------------------
			Fl_Box *bx_FILE_I = new Fl_Box(0, wgt_y, g->w(), 20, "--- LOAD ---");

			wgt_y += 20;

			btn_load = new Fl_Button(wgt_x, wgt_y, 45, 20, "load");
			btn_load->callback(cb_btn_load, (void*)this);
			fc = new Fl_File_Chooser(RECDIR, "(*.{txt})", Fl_File_Chooser::SINGLE, "Fl_File_Chooser");
			g->add(btn_load);

			btn_loadtex = new Fl_Button(wgt_x+45, wgt_y, 45, 20, "texture");
			btn_loadtex->deactivate();
			//btn_loadtex->callback(cb_btn_loadtex, (void*)this);
			g->add(btn_loadtex);

			btn_loadtpt = new Fl_Button(wgt_x+90, wgt_y, 45, 20, "tgt_pt");
			btn_loadtpt->callback(cb_btn_loadtpt, (void*)this);
			fc_tpt = new Fl_File_Chooser(RECDIR, "({target}*.{txt})", Fl_File_Chooser::SINGLE, "Fl_File_Chooser");
			g->add(btn_loadtpt);

			btn_loadrul = new Fl_Button(wgt_x+135, wgt_y, 45, 20, "rulings");
			btn_loadrul->callback(cb_btn_loadrul, (void*)this);
			fc_rul = new Fl_File_Chooser(RECDIR, "([rulings]*.{txt})", Fl_File_Chooser::SINGLE, "Fl_File_Chooser");
			g->add(btn_loadrul);
			
			wgt_x = 10;
			wgt_y += 20;

			btn_loadall = new Fl_Button(wgt_x, wgt_y, 45, 20, "all");
			btn_loadall->callback(cb_btn_loadall, (void*)this);
			fc_all = new Fl_File_Chooser(RECDIR, "(*.{csv})", Fl_File_Chooser::SINGLE, "Fl_File_Chooser");
			g->add(btn_loadall);

			btn_loadtptmask = new Fl_Button(wgt_x + 90, wgt_y, 45, 20, "tmask");
			btn_loadtptmask->callback(cb_btn_loadtptmask, (void*)this);
			fc_tptmask = new Fl_File_Chooser(RECDIR, "({tmask}*.{txt})", Fl_File_Chooser::SINGLE, "Fl_File_Chooser");
			g->add(btn_loadtptmask);

			wgt_x = 10;
			wgt_y += 20;

			Fl_Box *bx_FILE_O = new Fl_Box(0, wgt_y, g->w(), 20, "--- SAVE ---");

			wgt_y += 20;

			btn_savelog = new Fl_Button(wgt_x, wgt_y, 45, 20, "log");
			btn_savelog->callback(cb_btn_savelog, (void*)this);
			g->add(btn_savelog);

			btn_savescreen = new Fl_Button(wgt_x+45, wgt_y, 45, 20, "scshot");
			btn_savescreen->callback(cb_btn_savescreen, (void*)this);
			g->add(btn_savescreen);

			btn_savetpt = new Fl_Button(wgt_x+90, wgt_y, 45, 20, "tgt_pt");
			btn_savetpt->callback(cb_btn_savetpt, (void*)this);
			g->add(btn_savetpt);
#if 0
			btn_saveerr = new Fl_Button(wgt_x+90, wgt_y, 45, 20, "err");
			btn_saveerr->callback(cb_btn_saveerr, (void*)this);
			g->add(btn_saveerr);
#endif
			btn_saverul = new Fl_Button(wgt_x+135, wgt_y, 45, 20, "rulings");
			btn_saverul->callback(cb_btn_saverul, (void*)this);
			g->add(btn_saverul);

			wgt_x = 10;
			wgt_y += 20;

			btn_saveall = new Fl_Button(wgt_x, wgt_y, 45, 20, "all");
			btn_saveall->callback(cb_btn_saveall, (void*)this);
			g->add(btn_saveall);

			wgt_x = 10;
			wgt_y += 25;

			// ------------------------- EVALUATE -------------------------------------------
			Fl_Box *bx_EVAL= new Fl_Box(0, wgt_y, g->w(), 20, "--- EVALUATE ---");
			wgt_y += 20;

			btn_eval_gap = new Fl_Button(wgt_x, wgt_y, 60, 20, "GAP");
			btn_eval_gap->callback(cb_btn_eval_gap, (void*)this);

			btn_eval_collision = new Fl_Button(wgt_x+60, wgt_y, 60, 20, "COLLIS");
			btn_eval_collision->callback(cb_btn_eval_collision, (void*)this);

			btn_eval_rulingcross = new Fl_Button(wgt_x+120, wgt_y, 60, 20, "RULING");
			btn_eval_rulingcross->callback(cb_btn_eval_rulingcross, (void*)this);

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

			cb_dispaxis = new Fl_Check_Button(wgt_x+135, wgt_y, 45, 20, "W_AX");
			cb_dispaxis->value(this->gwin->disp_axis);
			cb_dispaxis->callback(cb_cb_dispaxis, (void*)this);
			g->add(cb_dispaxis);

			wgt_x = 10;
			wgt_y += 20;

			cb_disp[D_X] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "X");
			cb_disp[D_X]->value(this->gwin->disp_X);
			cb_disp[D_X]->callback(cb_cb_dispX, (void*)this);

			wgt_x += 45;

			cb_disp[D_TNB] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "TNB");
			cb_disp[D_TNB]->value(this->gwin->disp_TNB);
			cb_disp[D_TNB]->callback(cb_cb_dispTNB, (void*)this);

			wgt_x += 45;

			cb_disp[D_R] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "RUL");
			cb_disp[D_R]->value(this->gwin->disp_R);
			cb_disp[D_R]->callback(cb_cb_dispR, (void*)this);

			wgt_x += 45;

			cb_disp[D_AX2] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "L_AX");
			cb_disp[D_AX2]->value(this->gwin->disp_axis2);
			cb_disp[D_AX2]->callback(cb_cb_dispAx2, (void*)this);

			g->add(cb_disp[D_X]);
			g->add(cb_disp[D_TNB]);
			g->add(cb_disp[D_R]);
			g->add(cb_disp[D_AX2]);

			wgt_x = 10;	wgt_y += 20;

			cb_disp[D_PLY] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "PLY");
			cb_disp[D_PLY]->value(this->gwin->disp_PLY);
			cb_disp[D_PLY]->callback(cb_cb_dispPLY, (void*)this);

			wgt_x += 45;

			cb_disp[D_PTN] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "PTN");
			cb_disp[D_PTN]->value(this->gwin->disp_PTN);
			cb_disp[D_PTN]->callback(cb_cb_dispPTN, (void*)this);

			wgt_x += 45;

			cb_disp[D_CP] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "CPT");
			cb_disp[D_CP]->value(this->gwin->disp_CP);
			cb_disp[D_CP]->callback(cb_cb_dispCP, (void*)this);

			wgt_x += 45;

			g->add(cb_disp[D_PLY]);
			g->add(cb_disp[D_PTN]);
			g->add(cb_disp[D_CP]);

			wgt_x = 10;	wgt_y += 20;

			cb_disp[D_ONE] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "ONE");
			cb_disp[D_ONE]->value(this->gwin->disp_ONE);
			cb_disp[D_ONE]->callback(cb_cb_dispONE, (void*)this);
			g->add(cb_disp[D_ONE]);

			wgt_x += 45;

			cb_disp[D_PRI] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "PRI");
			cb_disp[D_PRI]->value(this->gwin->disp_PRI);
			cb_disp[D_PRI]->callback(cb_cb_dispPRI, (void*)this);
			g->add(cb_disp[D_PRI]);

			wgt_x += 45;

			cb_disp[D_ST] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "stch");
			cb_disp[D_ST]->value(this->gwin->disp_stitch);
			cb_disp[D_ST]->callback(cb_cb_dispST, (void*)this);
			g->add(cb_disp[D_ST]);

			wgt_x += 45;

			cb_disp[D_TGT] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "TGT");
			cb_disp[D_TGT]->value(this->gwin->disp_TGT);
			cb_disp[D_TGT]->callback(cb_cb_dispTGT, (void*)this);
			g->add(cb_disp[D_TGT]);
#if 1
			wgt_x = 10;
			wgt_y += 25;

			cb_CurveEnd = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "CE");
			if( ppm.flg_postproc_type == PPTYPE_X ){
				cb_CurveEnd->value(0);
			} else {
				cb_CurveEnd->value(1);
			}
			cb_CurveEnd->callback(cb_cb_CurveEnd, (void*)this);
			g->add(cb_CurveEnd);
#endif
			wgt_x += 45;

			cb_disp[D_LSMT] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "SMT");
			cb_disp[D_LSMT]->value(this->gwin->disp_LIN_SMOOTH);
			cb_disp[D_LSMT]->callback(cb_cb_dispLSMT, (void*)this);
			g->add(cb_disp[D_LSMT]);

			wgt_x += 45;

			cb_disp[D_OFST] = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "OFST");
			cb_disp[D_OFST]->value(this->gwin->disp_POLY_OFFSET);
			cb_disp[D_OFST]->callback(cb_cb_dispOFST, (void*)this);
			g->add(cb_disp[D_OFST]);

			wgt_x = 10;
			wgt_y += 25;

			// ------------------------- CONFIGURATION -------------------------------------------
			Fl_Box *bx_RulRes= new Fl_Box(0, wgt_y, g->w(), 20, "--- RULING RESOLUTION ---");	wgt_y += 20;
			vs_rulres = new Fl_Value_Slider(wgt_x, wgt_y, 180, 20);
			vs_rulres->bounds( 10, 100 );	vs_rulres->step(10);	vs_rulres->value(XCNT_DEF);
			vs_rulres->align(FL_ALIGN_LEFT);
			vs_rulres->type(FL_HORIZONTAL);
			vs_rulres->callback(cb_vs_rulres, (void*)this);
			g->add( vs_rulres );

		}
		g->end();
	}
	{
		Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h()-20, "D1");
		g->hide();
		{ 
			// ------------------------- DESIGN -------------------------------------------

			wgt_x = 10;	wgt_y = 20;
			int y_space = 25, grp_num=0;

			Fl_Box *bx_MODE= new Fl_Box(0, wgt_y, g->w(), 20, "--- CONTROL MODE ---");	wgt_y+=20;

			grp_fix = new Fl_Group(0, wgt_y, g->w(), 4*y_space);	grp_num = 0;
			rb_fix[CMODE_A] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "3D + ANGLE -> 2D");	wgt_y+=y_space; grp_num++;
			rb_fix[CMODE_B] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "2D + ANGLE -> 3D");	wgt_y+=y_space; grp_num++;
			rb_fix[CMODE_C] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "2D + 3D -> ANGLE");	wgt_y+=y_space; grp_num++;
			rb_fix[CMODE_R] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "Ruling Angle");	wgt_y+=y_space; grp_num++;
			for( int i=0; i<grp_num; i++ ){
				rb_fix[i]->type(FL_RADIO_BUTTON);
				rb_fix[i]->callback( cb_rb_fix0, (void*)this );
				grp_fix->add(rb_fix[i]);
			}
			grp_fix->end();
			rb_fix[CMODE_A]->setonly();
			g->add(grp_fix);

			Fl_Box *bx_PARAMETER= new Fl_Box(0, wgt_y, g->w(), 20, "--- PARAMETER ---");	wgt_y+=20;

			grp_param0 = new Fl_Group(0, wgt_y, g->w(), 6*y_space);	grp_num=0;
			rb_param[P_CV2D] = new Fl_Round_Button(wgt_x, wgt_y, 150, 20, "CURVATURE (2D)");	wgt_y+=y_space; grp_num++;
			rb_param[P_CV3D] = new Fl_Round_Button(wgt_x, wgt_y, 150, 20, "CURVATURE (3D)");	wgt_y+=y_space; grp_num++;
			rb_param[P_TRSN] = new Fl_Round_Button(wgt_x, wgt_y, 150, 20, "TORSION");		wgt_y+=y_space; grp_num++;
			rb_param[P_FLDA] = new Fl_Round_Button(wgt_x, wgt_y, 150, 20, "FOLDING ANGLE");	wgt_y+=y_space; grp_num++;
			rb_param[P_RULL] = new Fl_Round_Button(wgt_x, wgt_y, 110, 20, "RULING LEFT");	grp_num++;
			rb_param[P_RULR] = new Fl_Round_Button(wgt_x+120, wgt_y, 60, 20, "RIGHT");		wgt_y+=y_space; grp_num++;

			for( int i=0; i<grp_num; i++ ){
				rb_param[i]->type(FL_RADIO_BUTTON);
				rb_param[i]->callback( cb_rb_param0, (void*)this );
				grp_param0->add(rb_param[i]);
			}
			grp_param0->end();
			rb_param[P_CV3D]->setonly();
			g->add(grp_param0);

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
			vs_pval->bounds(0,1);	vs_pval->step(0.01);	vs_pval->value(gwin->pprm);
			vs_pval->align(FL_ALIGN_LEFT);
			vs_pval->type(FL_HORIZONTAL);
			vs_pval->callback(cb_vs_pval, (void*)this);
			g->add( vs_pval );

			wgt_y+=y_space;

			btn_RPar = new Fl_Button(wgt_x,  wgt_y, 150, 20, "Set Ruling Horizontal");
			btn_RPar->callback(cb_btn_RPar, (void*)this);
			g->add(btn_RPar);

			wgt_y+=y_space;

			btn_R2TA = new Fl_Button(wgt_x,  wgt_y, 150, 20, "Apply Ruling Angle");
			btn_R2TA->callback(cb_btn_R2TA, (void*)this);
			g->add(btn_R2TA);

			wgt_y+=y_space;

			// ------------------------- FOLD_MOTION -------------------------------------------
			Fl_Box *bx_FMOTION= new Fl_Box(0, wgt_y, g->w(), 20, "--- FOLDING_MOTION ---");	wgt_y+=20;
			vs_fmot = new Fl_Value_Slider(wgt_x, wgt_y, 140, 20);
			vs_fmot->bounds(-an_fcnt,an_fcnt*2);	vs_fmot->step(1);	vs_fmot->value(0);
			vs_fmot->align(FL_ALIGN_LEFT);
			vs_fmot->type(FL_HORIZONTAL);
			vs_fmot->callback(cb_vs_fmot, (void*)this);
			g->add( vs_fmot );

			btn_apply = new Fl_Button(wgt_x+140,  wgt_y, 40, 20, "Apply");
			btn_apply->callback(cb_btn_apply, (void*)this);
			g->add(btn_apply);

		}
		g->end();
	}
	{
		Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h()-20, "D2");
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

			// ------------------------- FOLD_MOTION -------------------------------------------
			Fl_Box *bx_FMOTION= new Fl_Box(0, wgt_y, g->w(), 20, "--- FOLDING_MOTION ---");	wgt_y+=20;
			vs_fmot = new Fl_Value_Slider(wgt_x, wgt_y, 140, 20);
			vs_fmot->bounds(-an_fcnt,an_fcnt);	vs_fmot->step(1);	vs_fmot->value(0);
			vs_fmot->align(FL_ALIGN_LEFT);
			vs_fmot->type(FL_HORIZONTAL);
			vs_fmot->callback(cb_vs_fmot, (void*)this);
			g->add( vs_fmot );

			btn_apply = new Fl_Button(wgt_x+140,  wgt_y, 40, 20, "Apply");
			btn_apply->callback(cb_btn_apply, (void*)this);
			g->add(btn_apply);

		}
		g->end();
	}
	{
		Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h()-20, "D3");
		g->hide();
		{ 
			// ------------------------- FOLDING ANGLE ------------------------------------

			wgt_x = 10;	wgt_y = 20;
			int y_space = 25, y_size = 20, grp_num=0;

			crease *c = &(this->ppm.crs[0]);

			Fl_Box *bx_DIVNUM= new Fl_Box(0, wgt_y, g->w(), y_size, "--- DIVISION NUMBER ---");	wgt_y+=y_size;
			vs_divnum = new Fl_Value_Slider(wgt_x, wgt_y, 180, y_size);
			vs_divnum->bounds(1,10);	vs_divnum->step(1);	vs_divnum->value( gwin->divnum );
			vs_divnum->align(FL_ALIGN_LEFT);
			vs_divnum->type(FL_HORIZONTAL);
			vs_divnum->callback(cb_vs_divnum, (void*)this);
			g->add( vs_divnum );

			wgt_y+=y_space;

#if 1
			Fl_Box *bx_DIVCRV= new Fl_Box(0, wgt_y, g->w(), y_size, "--- DIVISION CURVE ---");	wgt_y+=y_size;
			vs_divcrv = new Fl_Value_Slider(wgt_x, wgt_y, 180, y_size);
			vs_divcrv->bounds( kmin, kmax );	vs_divcrv->step( kstep );	vs_divcrv->value( 0 );
			vs_divcrv->align(FL_ALIGN_LEFT);
			vs_divcrv->type(FL_HORIZONTAL);
			//vs_divcrv->callback(cb_vs_divcrv, (void*)this);
			vs_divcrv->callback(cb_vs_divnum, (void*)this);
			g->add( vs_divcrv );
#else
			Fl_Box *bx_DIVTYPE= new Fl_Box(0, wgt_y, g->w(), 20, "--- DIVISION TYPE---");	wgt_y+=y_size;

			grp_divtype = new Fl_Group(0, wgt_y, g->w(), 2*y_space);	grp_num=0;

			rb_divtype[0] = new Fl_Round_Button(wgt_x, wgt_y, 45, 20, "cylinder");
			rb_divtype[0]->callback( cb_rb_divtype, (void*)this );
			rb_divtype[0]->type(FL_RADIO_BUTTON);
			rb_divtype[0]->deactivate();
			grp_divtype->add(rb_divtype[0]);

			rb_divtype[1] = new Fl_Round_Button(wgt_x+100, wgt_y, 45, 20, "disk");
			rb_divtype[1]->callback( cb_rb_divtype, (void*)this );
			rb_divtype[1]->type(FL_RADIO_BUTTON);
			rb_divtype[1]->setonly();
			grp_divtype->add(rb_divtype[1]);

			grp_divtype->end();
			g->add(grp_divtype);
#endif
			wgt_y+=y_space;

			grp_foldtrim = new Fl_Group(0, wgt_y, g->w(), y_space);
			rb_fold = new Fl_Round_Button(wgt_x,wgt_y,50,y_space,"fold");
			rb_trim = new Fl_Round_Button(wgt_x+50,wgt_y,50,y_space,"trim");
			rb_fold->type(FL_RADIO_BUTTON);
			rb_trim->type(FL_RADIO_BUTTON);
			rb_fold->callback(cb_vs_divnum, (void*)this);
			rb_trim->callback(cb_vs_divnum, (void*)this);
			rb_fold->setonly();
			grp_foldtrim->add(rb_fold);
			grp_foldtrim->add(rb_trim);
			grp_foldtrim->end();
			g->add(grp_foldtrim);

			cb_usecp = new Fl_Check_Button(wgt_x+100, wgt_y, 45, 20, "CTRL PT");
			cb_usecp->value(1);
			g->add(cb_usecp);

			wgt_y+=y_space;

			Fl_Box *bx_FMOTION= new Fl_Box(0, wgt_y, g->w(), y_size, "--- FOLDING_MOTION ---");	wgt_y+=y_size;
			vs_fmot2 = new Fl_Value_Slider(0, wgt_y, g->w(), y_size);
			vs_fmot2->bounds(0,1);	vs_fmot2->step(1);	vs_fmot2->value(0);
			vs_fmot2->align(FL_ALIGN_LEFT);
			vs_fmot2->type(FL_HORIZONTAL);
			vs_fmot2->callback(cb_vs_fmot2, (void*)this);
			g->add( vs_fmot2 );

			wgt_y+=y_size;
#if 1
			//vs_fmot2->bounds(-19,19);
			vs_fmot2->bounds(-20,20);
			int min = (int)vs_fmot2->minimum();
			int max = (int)vs_fmot2->maximum();
			int size = max-min+1;
			int dx = (int)( (float)(vs_fmot2->w() -30) / (float)size + 0.5);
			for( int j=min, x=36; j<=max; j++, x+=dx ){
				//int frm = This->vs_fmot2->value() + c->FM_fidx_org;
				int i = j + MAX_FRAME/2; //c->FM_fidx_org;
				bx_FM_Pt[i]= new Fl_Box(x, wgt_y,    dx, 15, "");
				bx_FM_Pt[i]->box(FL_FLAT_BOX); bx_FM_Pt[i]->color(FL_WHITE); g->add( bx_FM_Pt[i] );
				bx_FM_Pa[i]= new Fl_Box(x, wgt_y+15, dx, 15, "");
				bx_FM_Pa[i]->box(FL_FLAT_BOX); bx_FM_Pa[i]->color(FL_WHITE); g->add( bx_FM_Pa[i] );
				bx_FM_m3[i]= new Fl_Box(x, wgt_y+30, dx, 15, "");
				bx_FM_m3[i]->box(FL_FLAT_BOX); bx_FM_m3[i]->color(FL_WHITE); g->add( bx_FM_m3[i] );
			}
			wgt_y+=15*3;
#endif

			Fl_Box *bx_PARAMETER= new Fl_Box(0, wgt_y, g->w(), y_size, "--- PARAMETER ---");	wgt_y+=y_size;

			grp_param2 = new Fl_Group(0, wgt_y, g->w(), y_space*3);	grp_num=0;

			rb_param[P_FLDA1] = new Fl_Round_Button(wgt_x, wgt_y, 100, y_size, "Fold Angle");
			rb_param[P_FLDA1]->callback( cb_rb_param2, (void*)this );
			rb_param[P_FLDA1]->type(FL_RADIO_BUTTON);
			grp_param2->add(rb_param[P_FLDA1]);

			rb_param[P_TRSN1] = new Fl_Round_Button(wgt_x+100, wgt_y, 90, y_size, "Torsion");
			rb_param[P_TRSN1]->callback( cb_rb_param2, (void*)this );
			rb_param[P_TRSN1]->type(FL_RADIO_BUTTON);
			grp_param2->add(rb_param[P_TRSN1]);

			wgt_y+=y_space;

			btn_R2TA2 = new Fl_Button(wgt_x,  wgt_y, 60, 20, "Ruling");
			btn_R2TA2->callback(cb_btn_R2TA2, (void*)this);
			g->add(btn_R2TA2);

			rb_param[P_RULL1] = new Fl_Round_Button(wgt_x+70, wgt_y, 40, y_size, "L");
			rb_param[P_RULL1]->callback( cb_rb_param2, (void*)this );
			rb_param[P_RULL1]->type(FL_RADIO_BUTTON);
			grp_param2->add(rb_param[P_RULL1]);

			rb_param[P_RULR1] = new Fl_Round_Button(wgt_x+110, wgt_y, 40, y_size, "R");
			rb_param[P_RULR1]->callback( cb_rb_param2, (void*)this );
			rb_param[P_RULR1]->type(FL_RADIO_BUTTON);
			grp_param2->add(rb_param[P_RULR1]);

			grp_param2->end();
			g->add(grp_param2);

			wgt_y+=y_space;

			Fl_Box *bx_CPINDEX= new Fl_Box(0, wgt_y, g->w(), y_size, "--- CONTROL POINT INDEX ---");	wgt_y+=y_size;
			vs_ppos2 = new Fl_Value_Slider(wgt_x, wgt_y, 180, y_size);
			vs_ppos2->bounds(0,1);	vs_ppos2->step(1);	vs_ppos2->value(gwin->ppos);
			vs_ppos2->align(FL_ALIGN_LEFT);
			vs_ppos2->type(FL_HORIZONTAL);
			vs_ppos2->callback(cb_vs_ppos2, (void*)this);
			g->add( vs_ppos2 );

			wgt_y+=y_space;

			Fl_Box *bx_MODE= new Fl_Box(0, wgt_y, g->w(), y_size, "--- PARAMETER VALUE ---");	wgt_y+=y_size;
			vs_pval2 = new Fl_Value_Slider(wgt_x, wgt_y, 180, y_size);
			vs_pval2->bounds(0,1);	vs_pval2->step(0.01);	vs_pval2->value(gwin->pprm);
			vs_pval2->align(FL_ALIGN_LEFT);
			vs_pval2->type(FL_HORIZONTAL);
			vs_pval2->callback(cb_vs_pval2, (void*)this);
			g->add( vs_pval2 );

			wgt_y+=y_space;

			Fl_Box *bx_PAR= new Fl_Box(wgt_x, wgt_y, 60, y_size, "param");
			btn_paramopt = new Fl_Button(wgt_x+60, wgt_y, 40, y_size, "opt");
			btn_paramopt->callback(cb_btn_paramopt, (void*)this);
			btn_paramset = new Fl_Button(wgt_x+100, wgt_y, 40, y_size, "set");
			btn_paramset->callback(cb_btn_paramset, (void*)this);
			btn_paramreset = new Fl_Button(wgt_x+140, wgt_y, 40, y_size, "reset");
			btn_paramreset->callback(cb_btn_paramreset, (void*)this);

			wgt_y+=y_size;

			Fl_Box *bx_MAT= new Fl_Box(wgt_x, wgt_y, 60, y_size, "matrix");
			btn_matopt = new Fl_Button(wgt_x+60, wgt_y, 40, y_size, "opt");
			btn_matopt->callback(cb_btn_matopt, (void*)this);
			btn_matset = new Fl_Button(wgt_x+100, wgt_y, 40, y_size, "set");
			btn_matset->callback(cb_btn_matset, (void*)this);
			btn_matreset = new Fl_Button(wgt_x+140, wgt_y, 40, y_size, "reset");
			btn_matreset->callback(cb_btn_matreset, (void*)this);

			wgt_y+=y_size;

			Fl_Box *bx_MOT= new Fl_Box(wgt_x, wgt_y, 60, y_size, "motion");
			btn_savemotion = new Fl_Button(wgt_x+60, wgt_y, 60, y_size, "save");
			btn_savemotion->callback(cb_btn_savemotion, (void*)this);
			btn_loadmotion = new Fl_Button(wgt_x+120, wgt_y, 60, y_size, "load");
			btn_loadmotion->callback(cb_btn_loadmotion, (void*)this);

			wgt_y+=y_size;

			Fl_Box *bx_FRA= new Fl_Box(wgt_x, wgt_y, 60, y_size, "frame");
			btn_saveframe = new Fl_Button(wgt_x+60, wgt_y, 60, y_size, "save");
			btn_saveframe->callback(cb_btn_saveframe, (void*)this);
			btn_loadframe = new Fl_Button(wgt_x+120, wgt_y, 60, y_size, "load");
			btn_loadframe->callback(cb_btn_loadframe, (void*)this);

		}
		g->end();
	}
	{
	Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h() - 20, "D4");
	g->show();
	{
		wgt_x = 10;	wgt_y = 20;
		int y_space = 25, grp_num = 0;

		// ------------------------- RULING 2 CURVE -------------------------------------------
		Fl_Box* bx_RulCrv = new Fl_Box(0, wgt_y, g->w(), 20, "--- RULING 2 CURVE ---");	wgt_y += 20;

		vs_xmang1 = new Fl_Value_Slider(wgt_x, wgt_y, 180, 20);
		vs_xmang1->bounds(-90, 90);	vs_xmang1->step(2);	vs_xmang1->value(40);
		vs_xmang1->align(FL_ALIGN_LEFT);
		vs_xmang1->type(FL_HORIZONTAL);
		vs_xmang1->callback(cb_vs_xmang, (void*)this);
		g->add(vs_xmang1);

		wgt_y += 25;

		vs_xmang0 = new Fl_Value_Slider(wgt_x, wgt_y, 180, 20);
		vs_xmang0->bounds(-90, 90);	vs_xmang0->step(1);	vs_xmang0->value(45);
		vs_xmang0->align(FL_ALIGN_LEFT);
		vs_xmang0->type(FL_HORIZONTAL);
		g->add(vs_xmang0);

		wgt_y += 25;

		btn_switchRuling = new Fl_Button(wgt_x, wgt_y, 60, 20, "switch");
		btn_switchRuling->callback(cb_btn_switchRuling, (void*)this);
		g->add(btn_switchRuling);

		wgt_x += 60;

		btn_R2TA0 = new Fl_Button(wgt_x, wgt_y, 60, 20, "Rul->TA");
		btn_R2TA0->callback(cb_btn_R2TA0, (void*)this);
		g->add(btn_R2TA0);

		wgt_x += 60;

		btn_optmat = new Fl_Button(wgt_x, wgt_y, 60, 20, "Opt Mat");
		btn_optmat->callback(cb_btn_optmat, (void*)this);
		g->add(btn_optmat);

		wgt_x = 10;
		wgt_y += 25;

		btn_start = new Fl_Button(wgt_x, wgt_y, 60, 20, "start");
		btn_start->callback(cb_btn_start, (void*)this);
		g->add(btn_start);

		btn_stop = new Fl_Button(wgt_x + 70, wgt_y, 60, 20, "stop");
		btn_stop->callback(cb_btn_stop, (void*)this);
		g->add(btn_stop);

		cb_optmat = new Fl_Check_Button(wgt_x + 140, wgt_y, 20, 20, "mat");
		cb_optmat->value(0);
		g->add(cb_optmat);

		wgt_y += 25;

		// ------------------------- OPTIMIZATION -------------------------------------------
		wgt_x = 10;

		Fl_Box* bx_Opt = new Fl_Box(0, wgt_y, g->w(), 20, "--- OPTIMIZATION ---");	wgt_y += 20;

		btn_optfold2 = new Fl_Button(wgt_x, wgt_y, 180, 20, "Angle with Fixed Rulings");	wgt_x += 100;
		btn_optfold2->callback(cb_btn_optfold2, (void*)this);
		g->add(btn_optfold2);

		wgt_x = 10;
		wgt_y += 25;

		btn_randrul = new Fl_Button(wgt_x, wgt_y, 120, 20, "Random Rulings");	wgt_x += 120;
		btn_randrul->callback(cb_btn_randrul, (void*)this);
		g->add(btn_randrul);

		cb_optrul = new Fl_Check_Button(wgt_x+5, wgt_y, 20, 20, "optimize");
		cb_optrul->value(0);
		g->add(cb_optrul);

		wgt_x = 10;
		wgt_y += 25;

		btn_randrul2 = new Fl_Button(wgt_x, wgt_y, 180, 20, "Rulings and Angles");	wgt_x += 100;
		btn_randrul2->callback(cb_btn_randrul2, (void*)this);
		g->add(btn_randrul2);

		wgt_x = 10;
		wgt_y += 25;

		btn_optfold = new Fl_Button(wgt_x, wgt_y, 60, 20, "Angle");	wgt_x += 60;
		btn_optfold->callback(cb_btn_optfold, (void*)this);
		g->add(btn_optfold);


		btn_optrulfold = new Fl_Button(wgt_x, wgt_y, 100, 20, "Rul + Angle");	wgt_x += 100;
		btn_optrulfold->callback(cb_btn_optrulfold, (void*)this);
		g->add(btn_optrulfold);

		wgt_x = 10;
		wgt_y += 25;

		btn_opttr = new Fl_Button(wgt_x, wgt_y, 60, 20, "Trsn");	wgt_x += 60;
		btn_opttr->callback(cb_btn_opttr, (void*)this);
		btn_opttr->deactivate();
		g->add(btn_opttr);

		btn_optrul = new Fl_Button(wgt_x, wgt_y, 60, 20, "Rul");	wgt_x += 60;
		btn_optrul->callback(cb_btn_optrul, (void*)this);
		btn_optrul->deactivate();
		g->add(btn_optrul);

		btn_optcp = new Fl_Button(wgt_x, wgt_y, 60, 20, "CP");	wgt_x += 60;
		btn_optcp->callback(cb_btn_optcp, (void*)this);
		btn_optcp->deactivate();
		g->add(btn_optcp);

		wgt_x = 10;
		wgt_y += 25;

		// ------------------------- LIST RULINGS -------------------------------------------

		Fl_Box* bx_ListRul = new Fl_Box(0, wgt_y, g->w(), 20, "--- LIST RULINGS ---");	wgt_y += 20;

		btn_startrul = new Fl_Button(wgt_x, wgt_y, 60, 20, "start");	wgt_x += 70;
		btn_startrul->callback(cb_btn_startrul, (void*)this);
		g->add(btn_startrul);

		btn_stoprul = new Fl_Button(wgt_x, wgt_y, 60, 20, "stop");	wgt_x += 70;
		btn_stoprul->callback(cb_btn_stoprul, (void*)this);
		g->add(btn_stoprul);

		btn_makelist = new Fl_Button(wgt_x, wgt_y, 40, 20, "list");
		btn_makelist->callback(cb_btn_makelist, (void*)this);
		g->add(btn_makelist);

		wgt_x = 10;
		wgt_y += 25;

		grp_list = new Fl_Group(0, wgt_y, g->w(), y_space);
		rb_listleft = new Fl_Round_Button(wgt_x, wgt_y, 50, y_space, "left");
		rb_listright = new Fl_Round_Button(wgt_x + 50, wgt_y, 50, y_space, "right");
		rb_listleft->type(FL_RADIO_BUTTON);
		rb_listright->type(FL_RADIO_BUTTON);
		rb_listleft->setonly();
		grp_list->add(rb_listleft);
		grp_list->add(rb_listright);
		grp_list->end();
		g->add(grp_list);

		wgt_y += 25;

		// ------------------------- SAMPLE TARGET POINTS --------------------------------------

		Fl_Box* bx_ListTgt = new Fl_Box(0, wgt_y, g->w(), 20, "--- SAMPLE TARGET POINTS ---");	wgt_y += 20;

		btn_starttgt = new Fl_Button(wgt_x, wgt_y, 60, 20, "start");	wgt_x += 70;
		btn_starttgt->callback(cb_btn_starttgt, (void*)this);
		g->add(btn_starttgt);

		btn_stoptgt = new Fl_Button(wgt_x, wgt_y, 60, 20, "stop");	wgt_x += 70;
		btn_stoptgt->callback(cb_btn_stoptgt, (void*)this);
		g->add(btn_stoptgt);

		btn_listtgt = new Fl_Button(wgt_x, wgt_y, 40, 20, "list");
		btn_listtgt->callback(cb_btn_listtgt, (void*)this);
		g->add(btn_listtgt);

		wgt_x = 10;
		wgt_y += 25;

		// ------------------------- HISTORY -------------------------------------------

		Fl_Box* bx_Hist = new Fl_Box(0, wgt_y, g->w(), 20, "--- HISTORY ---");	wgt_y += 20;

		cb_history = new Fl_Check_Button(wgt_x, wgt_y, 20, 20, "history");
		cb_history->value(1);
		g->add(cb_history);

		wgt_y += 25;

		vs_history = new Fl_Value_Slider(wgt_x, wgt_y, 180, 20);
		vs_history->bounds(0, MAX_HISTORY-1);	vs_history->step(1);	vs_history->value(0);
		vs_history->align(FL_ALIGN_LEFT);
		vs_history->type(FL_HORIZONTAL);
		vs_history->callback(cb_vs_history, (void*)this);
		g->add(vs_history);

	}
	g->end();
	}
}

int ControlPanel::handle(int event)
{
	int i, key, pos=-1, prm=-1;
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
		} else if(key=='r'){
			this->gwin->resetObjRot();
			this->gwin->redraw();
		} else if(key=='t'){
			this->gwin->resetObjTrans();
			this->gwin->redraw();
		}
		break;
	default:
		//return Fl_Window::handle(event);
		break;
	}
	return Fl_Window::handle(event);
}
