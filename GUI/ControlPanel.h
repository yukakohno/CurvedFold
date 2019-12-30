
/*************************************************************
	ControlPanel.h
*************************************************************/

#include <FL/Fl.h>
#include <FL/Fl_Window.h>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Group.h>
#include <FL/Fl_Int_Input.H>
#include <FL/Fl_Round_Button.h>
#include <FL/Fl_Value_Slider.h>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Choice.h>
#include <FL/Fl_Tabs.H>
#include "GraphWindow3DCF.h"
#include "GraphWindowCP.h"
#include "GraphWindowParam.h"
#include "../CurvedFoldModel/papermodel.h"
#include "common.h"

#define RECDIR "input/"

//#define TEXFNAME "texture/no_tex.jpg"
#define TEXFNAME "texture/grid_bw3.jpg"
#define TEXWIDTH  256	// texture width
#define TEXHEIGHT 256	// texture height

#ifndef CPANEL
#define CPANEL

enum DISP { D_X, D_TNB, D_R, D_RLEN, D_PLY, D_PTN, D_CP, D_7, D_8, D_9 };

class ControlPanel : public Fl_Window
{
public:
	GraphWindow3DCF *gwin;
	GraphWindowCP *gwin_cp;
	GraphWindowParam *gwin_gr;

	papermodel ppm;

	// for animation
	double kv_org[MAX_SPCNT], k2d_org[MAX_SPCNT], tr_org[MAX_SPCNT], alpha_org[MAX_SPCNT];
	double kv_dif0[MAX_SPCNT], k2d_dif0[MAX_SPCNT], tr_dif0[MAX_SPCNT], alpha_dif0[MAX_SPCNT];
	double kv_dif1[MAX_SPCNT], k2d_dif1[MAX_SPCNT], tr_dif1[MAX_SPCNT], alpha_dif1[MAX_SPCNT];
	int an_fcnt;

	double kmin,kmax,kstep,kstepk, tmin,tmax,tstep,tstepk, amin,amax,astep,astepk;

	ControlPanel(int X, int Y, int W, int H, GraphWindow3DCF *_gwin );
	ControlPanel(int X, int Y, int W, int H, const char *L, GraphWindow3DCF *_gwin );
	~ControlPanel();
	void initPanel();
	void createPanel();
	void refresh(int init);
	//void refresh2();
	void setorg();

private:
	int handle(int event);
	int value_grpfix();
	int value_grpparam();
	int value_grpopt();

	// ------------------------- fILE_IO -------------------------------------------
public:
	Fl_Button *btn_load;
	Fl_File_Chooser *fc;
	Fl_Button *btn_loadtex;
	Fl_Button *btn_dump;
	Fl_Button *btn_rlchk;

private:
	static void cb_btn_load( Fl_Widget *wgt, void *idx);
	static void cb_btn_dump( Fl_Widget *wgt, void *idx);

	// ------------------------- DISPLAY -------------------------------------------
public:
	Fl_Button *btn_resetRot;
	Fl_Button *btn_resetScale;
	Fl_Button *btn_resetTrans;
	Fl_Check_Button *cb_dispaxis;
	Fl_Check_Button *cb_disp[20]; /* enum DISP */
	Fl_Check_Button *cb_CurveEnd;
private:
	static void cb_cb_dispaxis(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispX(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispTNB(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispR(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispRlen(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispPLY(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispPTN(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispCP(Fl_Widget *wgt, void *idx);
	static void cb_cb_CurveEnd(Fl_Widget *wgt, void *idx);

	// ------------------------- DESIGN -------------------------------------------
public:
	// MODE
	Fl_Group *grp_fix;
	Fl_Round_Button *rb_fix[N_C_MODE];

	// PARAMETER
	Fl_Group *grp_param;
	Fl_Round_Button *rb_param[N_P_IDX];
	Fl_Value_Slider *vs_ppos;
	Fl_Value_Slider *vs_pval;
	Fl_Button *btn_R2TA;
	Fl_Button *btn_RPar;

private:
	static void cb_rb_gwin(Fl_Widget *wgt, void *idx);
	static void cb_vs_ppos(Fl_Widget *wgt, void *idx);
	static void cb_vs_pval(Fl_Widget *wgt, void *idx);
	static void cb_btn_R2TA(Fl_Widget *wgt, void *idx);
	static void cb_btn_RPar(Fl_Widget *wgt, void *idx);

	// ------------------------- RECTIFY -------------------------------------------
public:
	// ï‚ê≥Ç†ÇËÅ^Ç»Çµ
	Fl_Check_Button *cb_rectifyA;
	Fl_Check_Button *cb_rectifyT;
	Fl_Check_Button *cb_rectifyR;
	Fl_Button *btn_updateChoice;
	Fl_Button *btn_setSEParam;

	Fl_Choice *ch_T_method, *ch_T_se, *ch_T_src;
	Fl_Value_Input *vi_T_kvthres, *vi_T_s1mgn, *vi_T_s2mgn; 
	Fl_Value_Input *vi_T_s2, *vi_T_s1, *vi_T_e1, *vi_T_e2;
	Fl_Value_Input *vi_T_sv2, *vi_T_sv1, *vi_T_ev1, *vi_T_ev2;

	Fl_Choice *ch_A_method, *ch_A_se, *ch_A_src;
	Fl_Value_Input *vi_A_kvthres, *vi_A_s1mgn, *vi_A_s2mgn;
	Fl_Value_Input *vi_A_s2, *vi_A_s1, *vi_A_e1, *vi_A_e2;
	Fl_Value_Input *vi_A_sv2, *vi_A_sv1, *vi_A_ev1, *vi_A_ev2;

	Fl_Value_Input *vi_R_kvthres;

private:
	static void cb_cb_rectifyA(Fl_Widget *wgt, void *idx);
	static void cb_cb_rectifyT(Fl_Widget *wgt, void *idx);
	static void cb_cb_rectifyR(Fl_Widget *wgt, void *idx);

	static void cb_btn_updateChoice(Fl_Widget *wgt, void *idx);
	static void cb_btn_setSEParam(Fl_Widget *wgt, void *idx);

	static void cb_ch_T_method(Fl_Widget *wgt, void *idx);
	static void cb_ch_T_se(Fl_Widget *wgt, void *idx);
	static void cb_ch_T_src(Fl_Widget *wgt, void *idx);
	static void cb_vi_T_kvthres(Fl_Widget *wgt, void *idx);
	static void cb_vi_T_s1mgn(Fl_Widget *wgt, void *idx);
	static void cb_vi_T_s2mgn(Fl_Widget *wgt, void *idx);
	static void cb_vi_T_s2(Fl_Widget *wgt, void *idx);
	static void cb_vi_T_s1(Fl_Widget *wgt, void *idx);
	static void cb_vi_T_e1(Fl_Widget *wgt, void *idx);
	static void cb_vi_T_e2(Fl_Widget *wgt, void *idx);
	static void cb_vi_T_sv2(Fl_Widget *wgt, void *idx);
	static void cb_vi_T_sv1(Fl_Widget *wgt, void *idx);
	static void cb_vi_T_ev1(Fl_Widget *wgt, void *idx);
	static void cb_vi_T_ev2(Fl_Widget *wgt, void *idx);

	static void cb_ch_A_method(Fl_Widget *wgt, void *idx);
	static void cb_ch_A_se(Fl_Widget *wgt, void *idx);
	static void cb_ch_A_src(Fl_Widget *wgt, void *idx);
	static void cb_vi_A_kvthres(Fl_Widget *wgt, void *idx);
	static void cb_vi_A_s1mgn(Fl_Widget *wgt, void *idx);
	static void cb_vi_A_s2mgn(Fl_Widget *wgt, void *idx);
	static void cb_vi_A_s2(Fl_Widget *wgt, void *idx);
	static void cb_vi_A_s1(Fl_Widget *wgt, void *idx);
	static void cb_vi_A_e1(Fl_Widget *wgt, void *idx);
	static void cb_vi_A_e2(Fl_Widget *wgt, void *idx);
	static void cb_vi_A_sv2(Fl_Widget *wgt, void *idx);
	static void cb_vi_A_sv1(Fl_Widget *wgt, void *idx);
	static void cb_vi_A_ev1(Fl_Widget *wgt, void *idx);
	static void cb_vi_A_ev2(Fl_Widget *wgt, void *idx);

	static void cb_vi_R_kvthres(Fl_Widget *wgt, void *idx);

	// ------------------------- ADD CREASE & TRIM ------------------------------------
public:
	Fl_Check_Button *cb_trimcurve;
	Fl_Check_Button *cb_foldcurve;
	Fl_Button *btn_proccurve;

private:
	static void cb_cb_addcurve( Fl_Widget *wgt, void *idx);
	static void cb_btn_proccurve( Fl_Widget *wgt, void *idx);
	static void cb_btn_resetcurve( Fl_Widget *wgt, void *idx);
	void trimfold();

	// ------------------------- FOLD_MOTION -------------------------------------------
public:
	Fl_Value_Slider *vs_fmot;
private:
	static void cb_vs_fmot(Fl_Widget *wgt, void *idx);

};

#endif	// CPANEL
