#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Windows.h>
#include "ControlPanel.h"

void ControlPanel::cb_cb_rectifyA(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.flg_rectifyA = ((Fl_Check_Button*)wgt)->value()!=0;
	This->ppm.set_postproc_type( PPTYPE_PRICURVE );
	This->refresh(1);
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
}

void ControlPanel::cb_cb_rectifyT(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.flg_rectifyT = ((Fl_Check_Button*)wgt)->value()!=0;
	This->ppm.set_postproc_type( PPTYPE_PRICURVE );
	This->refresh(1);
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
}

void ControlPanel::cb_cb_rectifyR(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.flg_rectifyR = ((Fl_Check_Button*)wgt)->value()!=0;
	This->ppm.set_postproc_type( PPTYPE_PPEDGE );
	//This->refresh(1);
	This->ppm.postproc();
	This->ppm.set_postproc_type( PPTYPE_UNDEF );
}

void ControlPanel::cb_btn_updateChoice(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	char strnum[10][10]={"0","1","2","3","4","5","6","7","8","9"};
	This->ch_A_se->clear();
	for( int i=0; i<This->ppm.rp.rectA.scnt; i++ ){
		This->ch_A_se->add( strnum[i] );
	}
	This->ch_A_se->value(0);
	This->ch_A_se->do_callback();
	This->ch_T_se->clear();
	for( int i=0; i<This->ppm.rp.rectT.scnt; i++ ){
		This->ch_T_se->add( strnum[i] );
	}
	This->ch_T_se->value(0);
	This->ch_T_se->do_callback();
}

void ControlPanel::cb_btn_setSEParam(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	rectify_param *rpA = &(This->ppm.rp.rectA);
	rectify_param *rpT = &(This->ppm.rp.rectT);
	int idxA = This->ch_A_se->value();
	int idxT = This->ch_T_se->value();

	if( 0 <= idxA && idxA < This->ppm.rp.rectA.scnt )
	{
		rpA->s2[idxA] = (int)This->vi_A_s2->value();
		rpA->s1[idxA] = (int)This->vi_A_s1->value();
		rpA->e1[idxA] = (int)This->vi_A_e1->value();
		rpA->e2[idxA] = (int)This->vi_A_e2->value();
		rpA->svlen0[idxA] = (double)This->vi_A_sv2->value();
		rpA->svlen1[idxA] = (double)This->vi_A_sv1->value();
		rpA->evlen0[idxA] = (double)This->vi_A_ev1->value();
		rpA->evlen1[idxA] = (double)This->vi_A_ev2->value();
		printf("index_A = %d\n", idxA );
		printf("This->ppm.rp.rectA.[s2,s1,e1,e2] = %d, %d, %d, %d\n",
			rpA->s2[idxA], rpA->s1[idxA], rpA->e1[idxA], rpA->e2[idxA] );
		printf("This->ppm.rp.rectA.svlen* = %f, %f, %f, %f\n",
			rpA->svlen0[idxA], rpA->svlen1[idxA], rpA->evlen0[idxA], rpA->evlen1[idxA] );
	}

	if( 0 <= idxT && idxT < This->ppm.rp.rectT.scnt )
	{
		rpT->s2[idxT] = (int)This->vi_T_s2->value();
		rpT->s1[idxT] = (int)This->vi_T_s1->value();
		rpT->e1[idxT] = (int)This->vi_T_e1->value();
		rpT->e2[idxT] = (int)This->vi_T_e2->value();
		rpT->svlen0[idxT] = (double)This->vi_T_sv2->value();
		rpT->svlen1[idxT] = (double)This->vi_T_sv1->value();
		rpT->evlen0[idxT] = (double)This->vi_T_ev1->value();
		rpT->evlen1[idxT] = (double)This->vi_T_ev2->value();
		printf("index_T = %d\n", idxT );
		printf("This->ppm.rp.rectT.[s2,s1,e1,e2] = %d, %d, %d, %d\n",
			rpT->s2[idxT], rpT->s1[idxT], rpT->e1[idxT], rpT->e2[idxT] );
		printf("This->ppm.rp.rectT.svlen* = %f, %f, %f, %f\n",
			rpT->svlen0[idxT], rpT->svlen1[idxT], rpT->evlen0[idxT], rpT->evlen1[idxT] );
	}
}

void ControlPanel::cb_ch_A_method(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.rectA.method = ((Fl_Choice*)wgt)->value()+1;
	printf("This->ppm.rp.rectA.method = %d\n", This->ppm.rp.rectA.method);
}
void ControlPanel::cb_vi_A_kvthres(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.rectA.kvthres = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectA.kvthres = %f\n", This->ppm.rp.rectA.kvthres);
}
void ControlPanel::cb_vi_A_s1mgn(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.rectA.s1mgn = (int)((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectA.s1mgn = %d\n", This->ppm.rp.rectA.s1mgn);
}
void ControlPanel::cb_vi_A_s2mgn(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.rectA.s2mgn = (int)((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectA.s2mgn = %d\n", This->ppm.rp.rectA.s2mgn);
}
void ControlPanel::cb_ch_A_se(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int se_index = ((Fl_Choice*)wgt)->value();
	printf("se_index = %d\n", se_index);
	if( 0<=se_index && se_index<This->ppm.rp.rectA.scnt ){
		This->vi_A_s2->value( This->ppm.rp.rectA.s2[se_index] );
		This->vi_A_s1->value( This->ppm.rp.rectA.s1[se_index] );
		This->vi_A_e1->value( This->ppm.rp.rectA.e1[se_index] );
		This->vi_A_e2->value( This->ppm.rp.rectA.e2[se_index] );
		This->vi_A_sv2->value( This->ppm.rp.rectA.svlen0[se_index] );
		This->vi_A_sv1->value( This->ppm.rp.rectA.svlen1[se_index] );
		This->vi_A_ev1->value( This->ppm.rp.rectA.evlen0[se_index] );
		This->vi_A_ev2->value( This->ppm.rp.rectA.evlen1[se_index] );
	}
}
void ControlPanel::cb_ch_A_src(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.rectA.src = ((Fl_Choice*)wgt)->value();
	printf("This->ppm.rp.rectA.src = %d\n", This->ppm.rp.rectA.src);
}
void ControlPanel::cb_vi_A_s2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_A_se->value();
	This->ppm.rp.rectA.s2[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectA.s2[%d] = %d\n", index, This->ppm.rp.rectA.s2[index]);
}
void ControlPanel::cb_vi_A_s1(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_A_se->value();
	This->ppm.rp.rectA.s1[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectA.s1[%d] = %d\n", index, This->ppm.rp.rectA.s1[index]);
}
void ControlPanel::cb_vi_A_e1(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_A_se->value();
	This->ppm.rp.rectA.e1[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectA.e1[%d] = %d\n", index, This->ppm.rp.rectA.e1[index]);
}
void ControlPanel::cb_vi_A_e2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_A_se->value();
	This->ppm.rp.rectA.e2[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectA.e2[%d] = %d\n", index, This->ppm.rp.rectA.e2[index]);
}

void ControlPanel::cb_vi_A_sv2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_A_se->value();
	This->ppm.rp.rectA.svlen0[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectA.svlen0[%d] = %f\n", index, This->ppm.rp.rectA.svlen0[index]);
}
void ControlPanel::cb_vi_A_sv1(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_A_se->value();
	This->ppm.rp.rectA.svlen1[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectA.svlen1[%d] = %f\n", index, This->ppm.rp.rectA.svlen1[index]);
}
void ControlPanel::cb_vi_A_ev1(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_A_se->value();
	This->ppm.rp.rectA.evlen0[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectA.evlen0[%d] = %f\n", index, This->ppm.rp.rectA.evlen0[index]);
}
void ControlPanel::cb_vi_A_ev2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_A_se->value();
	This->ppm.rp.rectA.evlen1[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectA.evlen1[%d] = %f\n", index, This->ppm.rp.rectA.evlen1[index]);
}

void ControlPanel::cb_ch_T_method(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.rectT.method = ((Fl_Choice*)wgt)->value()+1;
	printf("This->ppm.rp.rectT.method = %d\n", This->ppm.rp.rectT.method);
}
void ControlPanel::cb_vi_T_kvthres(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.rectT.kvthres = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectT.kvthres = %f\n", This->ppm.rp.rectT.kvthres);
}
void ControlPanel::cb_vi_T_s1mgn(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.rectT.s1mgn = (int)((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectT.s1mgn = %d\n", This->ppm.rp.rectT.s1mgn);
}
void ControlPanel::cb_vi_T_s2mgn(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.rectT.s2mgn = (int)((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectT.s2mgn = %d\n", This->ppm.rp.rectT.s2mgn);
}
void ControlPanel::cb_ch_T_se(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int se_index = ((Fl_Choice*)wgt)->value();
	printf("se_index = %d\n", se_index);
	if( 0<=se_index && se_index<This->ppm.rp.rectT.scnt ){
		This->vi_T_s2->value( This->ppm.rp.rectT.s2[se_index] );
		This->vi_T_s1->value( This->ppm.rp.rectT.s1[se_index] );
		This->vi_T_e1->value( This->ppm.rp.rectT.e1[se_index] );
		This->vi_T_e2->value( This->ppm.rp.rectT.e2[se_index] );
		This->vi_T_sv2->value( This->ppm.rp.rectT.svlen0[se_index] );
		This->vi_T_sv1->value( This->ppm.rp.rectT.svlen1[se_index] );
		This->vi_T_ev1->value( This->ppm.rp.rectT.evlen0[se_index] );
		This->vi_T_ev2->value( This->ppm.rp.rectT.evlen1[se_index] );
	}
}
void ControlPanel::cb_ch_T_src(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.rectT.src = ((Fl_Choice*)wgt)->value();
	printf("This->ppm.rp.rectT.src = %d\n", This->ppm.rp.rectT.src);
}
void ControlPanel::cb_vi_T_s2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_T_se->value();
	This->ppm.rp.rectT.s2[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectT.s2[%d] = %d\n", index, This->ppm.rp.rectT.s2[index]);
}
void ControlPanel::cb_vi_T_s1(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_T_se->value();
	This->ppm.rp.rectT.s1[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectT.s1[%d] = %d\n", index, This->ppm.rp.rectT.s1[index]);
}
void ControlPanel::cb_vi_T_e1(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_T_se->value();
	This->ppm.rp.rectT.e1[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectT.e1[%d] = %d\n", index, This->ppm.rp.rectT.e1[index]);
}
void ControlPanel::cb_vi_T_e2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_T_se->value();
	This->ppm.rp.rectT.e2[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectT.e2[%d] = %d\n", index, This->ppm.rp.rectT.e2[index]);
}

void ControlPanel::cb_vi_T_sv2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_T_se->value();
	This->ppm.rp.rectT.svlen0[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectT.svlen0[%d] = %f\n", index, This->ppm.rp.rectT.svlen0[index]);
}
void ControlPanel::cb_vi_T_sv1(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_T_se->value();
	This->ppm.rp.rectT.svlen1[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectT.svlen1[%d] = %f\n", index, This->ppm.rp.rectT.svlen1[index]);
}
void ControlPanel::cb_vi_T_ev1(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_T_se->value();
	This->ppm.rp.rectT.evlen0[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectT.evlen0[%d] = %f\n", index, This->ppm.rp.rectT.evlen0[index]);
}
void ControlPanel::cb_vi_T_ev2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	int index = This->ch_T_se->value();
	This->ppm.rp.rectT.evlen1[index] = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectT.evlen1[%d] = %f\n", index, This->ppm.rp.rectT.evlen1[index]);
}

void ControlPanel::cb_vi_R_kvthres(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	This->ppm.rp.rectifyR_kvthres = ((Fl_Value_Input*)wgt)->value();
	printf("This->ppm.rp.rectifyR_kvthres = %f\n", This->ppm.rp.rectifyR_kvthres);
}
