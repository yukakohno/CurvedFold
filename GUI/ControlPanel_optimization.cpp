#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <stdio.h>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

void ControlPanel::cb_btn_optfold2(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
	int mode = This->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C, CMODE_R
	int prm = This->value_grpparam();	// P_CV2D, P_CV3D, P_TRSN, P_FLDA, P_RULL, P_RULR, ...

	if ( !(mode == CMODE_R && (prm == P_RULL || prm == P_RULR || prm == P_RULL1 || prm == P_RULR1)) )
	{
		return;
	}

	int minval=-180;
	double mintgap = 10000;
	for (int val = This->vs_xmang1->minimum(); val < This->vs_xmang1->maximum(); val++)
	{
		//This->vs_xmang0->value(val);
		//This->btn_R2TA0->do_callback();

		double mang = val / 180.0 * M_PI;
		int ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, 0, c->Xcnt, mang, 0);

		if (ret == 0) {
			ppm->set_postproc_type(PPTYPE_PRICURVE);
			ppm->postproc();
			ppm->set_postproc_type(PPTYPE_UNDEF);
			if (This->cb_optmat->value() && ppm->tgcnt > 3) {
				ppm->calcAvetgap();
				//ppm->optMat(CMODE_R);
			}
			//std::cout << "ppm->avetgap=" << ppm->avetgap << std::endl;
			if (mintgap > ppm->avetgap) {
				minval = val;
				mintgap = ppm->avetgap;
			}
		}
	}
	if (minval > -180) {
		std::cout << "gap=" << mintgap << std::endl;
		This->vs_xmang0->value(minval);
		This->vs_xmang1->value(minval);
		This->btn_R2TA0->do_callback(); // includes This->push_hist()
	}
	//This->push_hist(c->Pcnt, CMODE_R, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);
}

#define MAX_DB_SIZE 100

int set_dPb(double (&dPb)[MAX_DB_SIZE][MAX_CPCNT], int cpcnt)
{
	int size = 0;
	for (int i = 0; i < cpcnt; i++) { dPb[size][i] = 0.1; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPb[size][i] = -0.1; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPb[size][i] = 0.2 * (double)i/(double)cpcnt; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPb[size][i] = -0.2 * (double)i / (double)cpcnt; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPb[size][i] = 0.2 * (double)(cpcnt - i) / (double)cpcnt; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPb[size][i] = -0.2 * (double)(cpcnt - i) / (double)cpcnt; }	size++;
	for (int j = 0; j < cpcnt; j++) {
		for (int i = 0; i < cpcnt; i++) { dPb[size][i] = 0.0; }	dPb[size][j] = 0.1;	size++;
		for (int i = 0; i < cpcnt; i++) { dPb[size][i] = 0.0; }	dPb[size][j] = -0.1;	size++;
	}
	return size;
}

void ControlPanel::cb_btn_randrul(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
	int mode = This->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C, CMODE_R
	bool flg_optimize = This->cb_optrul->value();
	bool flg_rulOK = false;
	
	double Px2d[MAX_CPCNT], Pa[MAX_CPCNT], Pbr[MAX_CPCNT], Pbl[MAX_CPCNT];
	memcpy(Px2d, c->Px2d, sizeof(double) * c->Pcnt);
	memcpy(Pa, c->Pa, sizeof(double) * c->Pcnt);
	memcpy(Pbl, c->Pbl, sizeof(double) * c->Pcnt);
	memcpy(Pbr, c->Pbr, sizeof(double) * c->Pcnt);
	//This->push_hist(c->Pcnt, mode, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);
	double avetgap_prev = ppm->avetgap;

	int dPb_size = 10;
	double dPb[MAX_DB_SIZE][MAX_CPCNT];
	dPb_size = set_dPb( dPb, c->Pcnt );

	for (int i = 0; i < 100; i++)
	{
		// control point of the rulings
		int i0 = rand() % dPb_size;
		int i1 = rand() % dPb_size;
		for (int j = 0; j <= c->Pcnt; j++) {
			c->Pbl[j] = Pbl[j] + dPb[i0][j];
		}
		for (int j = 0; j <= c->Pcnt; j++) {
			c->Pbr[j] = Pbr[j] + dPb[i1][j];
		}

		// make new ruling
		crease::interpolate_spline_RulAngle(c->Pbl, c->Pcnt, c->betal, c->Xcnt, c->cotbl, c->cosbl, c->sinbl);
		for (int j = c->Xsidx; j <= c->Xeidx; j++) {
			c->rlx_cp[j] = c->cosbl[j] * c->Tx2d[j] - c->sinbl[j] * c->Ty2d[j];
			c->rly_cp[j] = c->sinbl[j] * c->Tx2d[j] + c->cosbl[j] * c->Ty2d[j];
			normalize_v2(&(c->rlx_cp[j]), &(c->rly_cp[j]));
			c->rllen[j] = ppm->pw * ppm->ph;
		}
		crease::interpolate_spline_RulAngle(c->Pbr, c->Pcnt, c->betar, c->Xcnt, c->cotbr, c->cosbr, c->sinbr);
		for (int j = c->Xsidx; j <= c->Xeidx; j++) {
			c->rrx_cp[j] = c->cosbr[j] * c->Tx2d[j] + c->sinbr[j] * c->Ty2d[j];
			c->rry_cp[j] = -c->sinbr[j] * c->Tx2d[j] + c->cosbr[j] * c->Ty2d[j];
			normalize_v2(&(c->rrx_cp[j]), &(c->rry_cp[j]));
			c->rrlen[j] = ppm->pw * ppm->ph;
		}
		c->calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);

		// check crossing
		int result = ppm->checkRulCross();
		if (result) { continue; }
		result = ppm->checkRulAngle();
		if (result) { continue; }
		result = ppm->checkRulCreaseCross();
		if (result) { continue; }

		flg_rulOK = true;
		break;
	}

	// if no valid ruling
	if (!flg_rulOK)
	{
		memcpy(c->Px2d, Px2d, sizeof(double) * c->Pcnt);
		memcpy(c->Pa, Pa, sizeof(double) * c->Pcnt);
		memcpy(c->Pbl, Pbl, sizeof(double) * c->Pcnt);
		memcpy(c->Pbr, Pbr, sizeof(double) * c->Pcnt);
		//int ret = This->access_hist(1, &(c->Pcnt), &mode, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);
		This->btn_R2TA0->do_callback();
		return;
	}

	// try & evaluate
	This->btn_optfold2->do_callback();
	This->push_hist(c->Pcnt, mode, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);

	if ( flg_optimize && ppm->avetgap > avetgap_prev ) {
		memcpy(c->Px2d, Px2d, sizeof(double) * c->Pcnt);
		memcpy(c->Pa, Pa, sizeof(double) * c->Pcnt);
		memcpy(c->Pbl, Pbl, sizeof(double) * c->Pcnt);
		memcpy(c->Pbr, Pbr, sizeof(double) * c->Pcnt);
		//int ret = This->access_hist(1, &(c->Pcnt), &mode, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);
		This->btn_R2TA0->do_callback();
	}
}

void ControlPanel::cb_btn_randrul2(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
	int mode = This->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C, CMODE_R
	if ( mode != CMODE_R ) {
		return;
	}
	This->cb_optrul->value(true);
	double thres_avetgap = 0.3;
	//double avetgap = 1000;

	for (int i = 0; i < 100; i++)
	{
		if (ppm->avetgap < thres_avetgap) {
			break;
		}
		This->btn_randrul->do_callback();
		//avetgap = ppm->avetgap;
		std::cout << "i = " << i << ", avetgap = " << ppm->avetgap << std::endl;
	}
}