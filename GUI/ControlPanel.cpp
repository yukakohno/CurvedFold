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

	kmin = tmin = -0.05;
	kmax = tmax = 0.05;
	kstep = tstep = 0.0025;
	kstepk = tstepk = 0.01;
	amin = 0.0;
	amax = 90.0;
	astep = 1.0;
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

int ControlPanel::value_grpopt()
{
	int opt = -1;
	for( int i=0; i<3; i++ ){
		if( this->rb_opt[i]->value() > 0 ){
			opt = i;
			break;
		}
	}
	return opt;
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
				_c->calcXA_CP( _ppm->flg_interpolate, _ppm->flg_rectifyT, _ppm->flg_rectifyA, _ppm->flg_rectifyR );
				_ppm->postproc();
			}
			break;
		case CMODE_B:	// B: Ü‚èü‚ÆÜ‚èŠp“x‚©‚çA3D‹Èü‚ð‹‚ß‚é
			if( prm==P_CV2D || prm==P_FLDA ){ // curv2d, angle
				_c->calcCPA_X( _ppm->flg_interpolate, _ppm->flg_rectifyT, _ppm->flg_rectifyA, _ppm->flg_rectifyR );
				_ppm->postproc();
			}
			break;
		case CMODE_C:	// C: Ü‚èü‚Æ3D‹Èü‚©‚çAÜ‚èŠp“x‚ð‹‚ß‚é
			if( prm==P_CV2D || prm==P_CV3D || prm==P_TRSN ){ // curv2d, curvature, torsion
				_c->calcCPX_A( _ppm->flg_interpolate, _ppm->flg_rectifyT, _ppm->flg_rectifyA, _ppm->flg_rectifyR );
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
		Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h()-20, "DESIGN1");
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

			Fl_Box *bx_RECTIFY= new Fl_Box(0, wgt_y, g->w(), 20, "--- RECTIFY ---");	wgt_y+=20;

			cb_rectifyA = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "FOLDING ANGLE");	wgt_y+=y_space;
			cb_rectifyA->value( this->ppm.flg_rectifyA );
			cb_rectifyA->callback(cb_cb_rectifyA, (void*)this);
			cb_rectifyT = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "TORSION");	wgt_y+=y_space;
			cb_rectifyT->value( this->ppm.flg_rectifyT );
			cb_rectifyT->callback(cb_cb_rectifyT, (void*)this);
			cb_rectifyR = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "Eliminate RULINGS");	wgt_y+=y_space;
			cb_rectifyR->value( this->ppm.flg_rectifyR );
			cb_rectifyR->callback(cb_cb_rectifyR, (void*)this);
			g->add(cb_rectifyA);
			g->add(cb_rectifyT);
			g->add(cb_rectifyR);

			// ------------------------- FOLD_MOTION -------------------------------------------
			Fl_Box *bx_FMOTION= new Fl_Box(0, wgt_y, g->w(), 20, "--- FOLDING_MOTION ---");	wgt_y+=20;
			vs_fmot = new Fl_Value_Slider(wgt_x, wgt_y, 180, 20);
			vs_fmot->bounds(-an_fcnt+1,an_fcnt-1);	vs_fmot->step(1);	vs_fmot->value(0);
			vs_fmot->align(FL_ALIGN_LEFT);
			vs_fmot->type(FL_HORIZONTAL);
			vs_fmot->callback(cb_vs_fmot, (void*)this);
			g->add( vs_fmot );
		}
		g->end();
	}
	{
		Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h()-20, "DESIGN2");
		g->hide();
		{ 

			// ------------------------- ADD CREASE & TRIM ------------------------------------

			wgt_x = 10;	wgt_y = 20;
			int y_space = 25, grp_num=0;

			Fl_Box *bx_TRIM= new Fl_Box(0, wgt_y, g->w(), 20, "--- ADD FOLD / TRIM ---");

			wgt_y+=25;

			cb_foldcurve = new Fl_Check_Button(wgt_x, wgt_y, 50, 20, "FOLD");
			cb_foldcurve->value( (int)(this->gwin->flg_addcurve == 2) );
			cb_foldcurve->callback(cb_cb_addcurve, (void*)this);

			wgt_y+=25;

			cb_trimcurve = new Fl_Check_Button(wgt_x, wgt_y, 50, 20, "TRIM");
			cb_trimcurve->value( (int)(this->gwin->flg_addcurve == 1) );
			cb_trimcurve->callback(cb_cb_addcurve, (void*)this);

			wgt_y -= 25;
			wgt_x = 80;

			btn_proccurve = new Fl_Button(wgt_x, wgt_y, 100, 20, "ADD CREASE");
			btn_proccurve->callback(cb_btn_proccurve, (void*)this);

			wgt_y+=25;

			btn_resetcurve = new Fl_Button(wgt_x, wgt_y, 100, 20, "REMOVE");
			btn_resetcurve->callback(cb_btn_resetcurve, (void*)this);

			g->add(cb_trimcurve);
			g->add(cb_foldcurve);
			g->add(btn_proccurve);
			g->add(btn_resetcurve);

			wgt_x = 10;	wgt_y += 25;

			Fl_Box *bx_CRSIDX= new Fl_Box(0, wgt_y, g->w(), 20, "--- CREASE INDEX ---");

			wgt_y += 25;

			vs_cidx = new Fl_Value_Slider(wgt_x, wgt_y, 180, 20 );
			vs_cidx->bounds(-MAX_CRS_CNT,MAX_CRS_CNT);
			vs_cidx->step(1);
			vs_cidx->value(0);
			vs_cidx->align(FL_ALIGN_LEFT);
			vs_cidx->type(FL_HORIZONTAL);
			vs_cidx->callback(cb_vs_cidx, (void*)this);
			g->add( vs_cidx );

			wgt_x = 10;	wgt_y += 25;

			Fl_Box *bx_COSTFUNC = new Fl_Box(0, wgt_y, g->w(), 20, "--- COST FUNCTION ---");

			wgt_y += 25;

			grp_opt = new Fl_Group(0, wgt_y, g->w(), 3*y_space);	grp_num = 0;
			rb_opt[0] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "TORSION");	wgt_y+=y_space; grp_num++;
			rb_opt[1] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "RULING ANGLE");	wgt_y+=y_space; grp_num++;
			rb_opt[2] = new Fl_Round_Button(wgt_x, wgt_y, 60, 20, "RULING CROSSING"); grp_num++;

			for( int i=0; i<grp_num; i++ ){
				rb_opt[i]->type(FL_RADIO_BUTTON);
				grp_opt->add(rb_opt[i]);
				rb_opt[i]->callback(cb_vs_otype, (void*)this);
			}
			grp_opt->end();
			rb_opt[this->ppm.hCrs_evaltype]->setonly();
			g->add(grp_opt);

			wgt_x = 10;	wgt_y += 25;

			Fl_Box *bx_OPTMODE = new Fl_Box(0, wgt_y, g->w(), 20, "--- ADJUSTMENT SUPPORT ---");

			wgt_y += 25;

			cb_rectifyC = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "RESTRICT");
			cb_rectifyC->value( this->ppm.hCrs_cridx );
			g->add(cb_rectifyC);

#if 0
			vs_otype = new Fl_Value_Slider(wgt_x+70, wgt_y, 130, 20);
			vs_otype->bounds(0,MAX_CRS_OTYPE);
			vs_otype->step(1);
			vs_otype->value(this->ppm.hCrs_evaltype);
			vs_otype->align(FL_ALIGN_LEFT);
			vs_otype->type(FL_HORIZONTAL);
			vs_otype->callback(cb_vs_otype, (void*)this);
			g->add( vs_otype );
#endif
			wgt_y += 25;

			btn_rectifyC = new Fl_Button(wgt_x, wgt_y, 100, 20, "OPTIMIZE");
			btn_rectifyC->callback(cb_btn_rectifyC, (void*)this);

			wgt_y += 25;
#if 1
			btn_resetC = new Fl_Button(wgt_x, wgt_y, 100, 20, "RESET");
			btn_resetC->callback(cb_btn_resetC, (void*)this);

			wgt_y += 25;

			btn_fixC = new Fl_Button(wgt_x, wgt_y, 100, 20, "FIX");
			btn_fixC->callback(cb_btn_fixC, (void*)this);

			g->add(btn_rectifyC);
			g->add(btn_resetC);
			g->add(btn_fixC);
#endif
			//wgt_x = 10;	wgt_y += 25;
		}
		g->end();
	}
	tab->end();
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
		}
		break;
	default:
		//return Fl_Window::handle(event);
		break;
	}
	return Fl_Window::handle(event);
}
