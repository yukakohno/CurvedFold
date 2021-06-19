#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <windows.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

#define TARGET_GAP_TYPE 0 // 0: average, 1: max
//#define INPUT_PARAMETERS
//#define PRINT_OPT_PROC
#define PRINT_RAND_PROC
//#define DUMP_OPT_RESULT

void ControlPanel::cb_btn_optfold(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);

	int ret = ppm->optFold();
	if (This->cb_optmat->value() && ppm->tgcnt > 3) {
		ppm->optMat(CMODE_B); // B: Ü‚èü‚ÆÜ‚èŠp“x‚©‚çA3D‹Èü‚ð‹‚ß‚é
	} else if (This->cb_optrot->value() && ppm->tgcnt > 3) {
		ppm->optMatRot(CMODE_B);
	}
	crease* c = &(This->ppm.crs[0]);
	This->push_hist(c->Pcnt, CMODE_B, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}

void ControlPanel::cb_btn_optfold2(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
	int mode = This->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C, CMODE_R
	int prm = This->value_grpparam();	// P_CV2D, P_CV3D, P_TRSN, P_FLDA, P_RULL, P_RULR, ...

#if 1
	This->rb_fix[CMODE_R]->setonly();
	This->rb_param[P_RULL]->setonly();
#else
	if ( !(mode == CMODE_R && (prm == P_RULL || prm == P_RULR || prm == P_RULL1 || prm == P_RULR1)) )
	{
		return;
	}
#endif
#if 1
	int crcnt, lcrcnt, rcrcnt, dccnt, fccnt, tccnt;
	dccnt = ppm->dccnt; ppm->dccnt = 0;
	fccnt = ppm->fccnt; ppm->fccnt = 0;
	tccnt = ppm->tccnt; ppm->tccnt = 0;
	crcnt = ppm->crcnt; ppm->crcnt = 1;
	lcrcnt = ppm->lcrcnt; ppm->lcrcnt = 0;
	rcrcnt = ppm->rcrcnt; ppm->rcrcnt = 0;
#endif
	int minval=-180;
	double mintgap = 10000;
	for (int val = This->vs_xmang1->minimum(); val < This->vs_xmang1->maximum(); val++)
	{
		//This->vs_xmang0->value(val);
		//This->btn_R2TA0->do_callback();
		//std::cout << "val=" << val;

		double mang = val / 180.0 * M_PI;
		ppm->re_sidx = 0;
		ppm->re_eidx = c->Xcnt;
		int ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, mang, 0);
		//std::cout << ", ret=" << ret;

		if (ret == 0) {
			ppm->set_postproc_type(PPTYPE_PRICURVE);
			ppm->postproc();
			ppm->set_postproc_type(PPTYPE_UNDEF);
			// recalculate average gap
			if (This->cb_optmat->value() && ppm->tgcnt > 3) {
				ppm->calcAvetgapMat();
			} else if (This->cb_optrot->value() && ppm->tgcnt > 3) {
				ppm->calcAvetgapRot(); // calc gap with fixed origin
			}
			//std::cout << ", avetgap=" << ppm->avetgap << ", maxtgap=" << ppm->maxtgap;
#if TARGET_GAP_TYPE == 0 // 0: average, 1: max
			if (mintgap > ppm->avetgap) {
				minval = val;
				mintgap = ppm->avetgap;
				//std::cout << ", minval=" << minval << ", mintgap=" << mintgap;
			}
#else
			if (mintgap > ppm->maxtgap) {
				minval = val;
				mintgap = ppm->maxtgap;
				//std::cout << ", minval=" << minval << ", mintgap=" << mintgap;
			}
#endif
		}
		//std::cout << std::endl;
	}
#if 1
	ppm->dccnt = dccnt;
	ppm->fccnt = fccnt;
	ppm->tccnt = tccnt;
	ppm->crcnt = crcnt;
	ppm->lcrcnt = lcrcnt;
	ppm->rcrcnt = rcrcnt;
#endif
	if (minval > -180) {
		std::cout << "gap=" << mintgap << std::endl;
		This->vs_xmang0->value(minval);
		This->vs_xmang1->value(minval);
		This->btn_R2TA0->do_callback(); // includes This->push_hist()
	}
	//This->push_hist(c->Pcnt, CMODE_R, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);
}
#if 0
void ControlPanel::cb_btn_opttr(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);

	int ret = ppm->optTorsion();
	crease* c = &(This->ppm.crs[0]);
	This->push_hist(c->Pcnt, CMODE_R, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}

void ControlPanel::cb_btn_optrul(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);

	int ret = ppm->optRulings();
	crease* c = &(This->ppm.crs[0]);
	This->push_hist(c->Pcnt, CMODE_R, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}

void ControlPanel::cb_btn_optcp(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);

	int ret = ppm->optCP();
	crease* c = &(This->ppm.crs[0]);
	This->push_hist(c->Pcnt, CMODE_R, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}

void ControlPanel::cb_btn_optrulfold(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
	crease c0; memcpy(&c0, c, sizeof(crease));

	if ( This->vPbl.size()==0 || This->vPbr.size()==0 ) {
		return;
	}

	int minval = -180, minli = -1, minri = -1;
	double mintgap = 10000;
	int limax = This->vPbl.size() / c->Pcnt, rimax = This->vPbr.size() / c->Pcnt;
	for (int li = 0; li < limax; li++ ) {
	//for (int li = 300; li < limax; li+=limax) {
		for (int i = 0; i < c->Pcnt; i++) {
			c->Pbl[i] = This->vPbl[ li * c->Pcnt + i ];
		}
		for (int ri = 0; ri < rimax; ri++) {
		//for (int ri = 300; ri < rimax; ri+=rimax) {
			for (int i = 0; i < c->Pcnt; i++) {
				c->Pbr[i] = This->vPbr[ ri * c->Pcnt + i ];
			}
			int valcnt = 0, valtotal = This->vs_xmang1->maximum() - This->vs_xmang1->minimum();
			for (int val = This->vs_xmang1->minimum(); val < This->vs_xmang1->maximum(); val++) {
				double mang = val / 180.0 * M_PI;
				ppm->re_sidx = 0;
				ppm->re_eidx = c->Xcnt;
				int ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, mang, 0);

				if (ret == 0) {
					ppm->set_postproc_type(PPTYPE_PRICURVE);
					ppm->postproc();
					ppm->set_postproc_type(PPTYPE_UNDEF);
					if (This->cb_optmat->value() && ppm->tgcnt > 3) {
						ppm->calcAvetgapMat();
					} else if (This->cb_optrot->value() && ppm->tgcnt > 3) {
						ppm->calcAvetgapRot();
					}
					//std::cout << "ppm->avetgap=" << ppm->avetgap << ", ppm->maxtgap=" << ppm->maxtgap << std::endl;
#if TARGET_GAP_TYPE == 0 // 0: average, 1: max
					if (mintgap > ppm->avetgap) {
						minval = val;
						minli = li;
						minri = ri;
						mintgap = ppm->avetgap;
					}
#else
					if (mintgap > ppm->maxtgap) {
						minval = val;
						minli = li;
						minri = ri;
						mintgap = ppm->maxtgap;
					}
#endif
					valcnt++;
				}
			}
			std::cout << "li=" << li << "/" << limax << ", ri=" << ri << "/" << rimax
				<< ", valcnt=" << valcnt << "/" << valtotal
				<< ", mintgap=" << mintgap << ", minli=" << minli << ", minri=" << minri << ", minval=" << minval << std::endl;
		}
	}

	if (minval > -180 && minli > -1 && minri > -1)
	{
		for (int i = 0; i < c->Pcnt; i++) {
			c->Pbl[i] = This->vPbl[minli * c->Pcnt + i];
		}
		for (int i = 0; i < c->Pcnt; i++) {
			c->Pbr[i] = This->vPbr[minri * c->Pcnt + i];
		}
		This->vs_xmang0->value(minval);
		This->vs_xmang1->value(minval);
		This->btn_R2TA0->do_callback();
	} else {
		memcpy(c, &c0, sizeof(crease));
	}
}
#endif
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

int set_dPa(double(&dPa)[MAX_DB_SIZE][MAX_CPCNT], int cpcnt)
{
	int size = 0;
	double astep = 0.25*M_PI/180.0;
	for (int i = 0; i < cpcnt; i++) { dPa[size][i] = astep; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPa[size][i] = -astep; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPa[size][i] = astep*2.0 * (double)i / (double)cpcnt; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPa[size][i] = -astep*2.0 * (double)i / (double)cpcnt; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPa[size][i] = astep*2.0 * (double)(cpcnt - i) / (double)cpcnt; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPa[size][i] = -astep*2.0 * (double)(cpcnt - i) / (double)cpcnt; }	size++;
	for (int j = 0; j < cpcnt; j++) {
		for (int i = 0; i < cpcnt; i++) { dPa[size][i] = 0.0; }	dPa[size][j] = astep;	size++;
		for (int i = 0; i < cpcnt; i++) { dPa[size][i] = 0.0; }	dPa[size][j] = -astep;	size++;
	}
	return size;
}

int set_dPy(double(&dPy)[MAX_DB_SIZE][MAX_CPCNT], int cpcnt)
{
	int size = 0;
	double tstep = 0.0005;
	for (int i = 0; i < cpcnt; i++) { dPy[size][i] = tstep; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPy[size][i] = -tstep; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPy[size][i] = tstep*2.0 * (double)i / (double)cpcnt; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPy[size][i] = -tstep*2.0 * (double)i / (double)cpcnt; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPy[size][i] = tstep*2.0 * (double)(cpcnt - i) / (double)cpcnt; }	size++;
	for (int i = 0; i < cpcnt; i++) { dPy[size][i] = -tstep*2.0 * (double)(cpcnt - i) / (double)cpcnt; }	size++;
	for (int j = 0; j < cpcnt; j++) {
		for (int i = 0; i < cpcnt; i++) { dPy[size][i] = 0.0; }	dPy[size][j] = tstep;	size++;
		for (int i = 0; i < cpcnt; i++) { dPy[size][i] = 0.0; }	dPy[size][j] = -tstep;	size++;
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
	double maxtgap_prev = ppm->maxtgap;

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
		int result = c->checkRulingCross();
		if (result) { continue; }
		result = c->checkRulingAngle();
		if (result) { continue; }
		result = c->checkRulCreaseCross();
		if (result) { continue; }

		flg_rulOK = true;
		This->optlog_trial_til_validrul[This->optlog_itr] = i;
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
		This->optlog_trial_til_validrul[This->optlog_itr] = 100;
		This->optlog_avetgap[This->optlog_itr] = -1.0;
		This->optlog_maxtgap[This->optlog_itr] = -1.0;
		This->optlog_min_avetgap[This->optlog_itr] = -1.0;
		This->optlog_min_maxtgap[This->optlog_itr] = -1.0;
		return;
	}

	// try & evaluate
	This->btn_optfold2->do_callback();
	This->push_hist(c->Pcnt, mode, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);

	This->optlog_avetgap[This->optlog_itr] = ppm->avetgap;
	This->optlog_maxtgap[This->optlog_itr] = ppm->maxtgap;

#if TARGET_GAP_TYPE == 0 // 0: average, 1: max
	if ( flg_optimize && ppm->avetgap > avetgap_prev )
#else
	if (flg_optimize && ppm->maxtgap > maxtgap_prev)
#endif
	{
		memcpy(c->Px2d, Px2d, sizeof(double) * c->Pcnt);
		memcpy(c->Pa, Pa, sizeof(double) * c->Pcnt);
		memcpy(c->Pbl, Pbl, sizeof(double) * c->Pcnt);
		memcpy(c->Pbr, Pbr, sizeof(double) * c->Pcnt);
		//int ret = This->access_hist(1, &(c->Pcnt), &mode, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);
		ppm->avetgap = avetgap_prev;
		ppm->maxtgap = maxtgap_prev;
		if (This->optlog_itr < 0 || This->optlog_itr == This->optlog_itrmax) {
			This->btn_R2TA0->do_callback();
		}
	}
	This->optlog_min_avetgap[This->optlog_itr] = ppm->avetgap;
	This->optlog_min_maxtgap[This->optlog_itr] = ppm->maxtgap;
}

void ControlPanel::cb_btn_randrul2(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
#if 1
	This->rb_fix[CMODE_R]->setonly();
	This->rb_param[P_RULL]->setonly();
#else
	int mode = This->value_grpfix();	// CMODE_A, CMODE_B, CMODE_C, CMODE_R
	if ( mode != CMODE_R ) {
		return;
	}
#endif
	This->cb_optrul->value(true);
	double thres_avetgap = 0.3;
	double thres_maxtgap = 1.0;
	//double avetgap = 1000;

	This->optlog_itr = 0;
	This->optlog_itrmax = 100;
	clock_t start_clock, end_clock;
	start_clock = clock();

	for (int i = 0; i < This->optlog_itrmax; i++)
	{
#if TARGET_GAP_TYPE == 0 // 0: average, 1: max
		if (ppm->avetgap < thres_avetgap) {
#else
		if (ppm->maxtgap < thres_maxtgap) {
#endif
			break;
		}
		This->btn_randrul->do_callback();
		//std::cout << "i = " << i << ", avetgap = " << ppm->avetgap <<  ", maxtgap = " << ppm->maxtgap << std::endl;
		std::cout << "iter = " << This->optlog_itr
			<< ", ave gap = " << This->optlog_avetgap[This->optlog_itr]
			<< ", max gap = " << This->optlog_maxtgap[This->optlog_itr]
			<< ", min ave gap = " << This->optlog_min_avetgap[This->optlog_itr]
			<< ", min max gap = " << This->optlog_min_maxtgap[This->optlog_itr]
			<< std::endl;
		This->optlog_itr++;
	}
	end_clock = clock();
	std::cout << (double)(end_clock - start_clock) / CLOCKS_PER_SEC << std::endl;
#if 1
	{
		std::ofstream ofs("./output/time0.csv", std::ios_base::app);
		ofs << "clock," << (double)(end_clock - start_clock) / CLOCKS_PER_SEC << std::endl;
		ofs.close();
	}
	{
		std::ofstream ofs("./output/err_seq0_target.csv");
		ofs << "iter,trial_validrul,err,minerr" << std::endl;
		for (int i = 0; i < This->optlog_itr; i++)
		{
			ofs << i << "," << This->optlog_trial_til_validrul[i] << ","
				<< This->optlog_avetgap[i] << "," << This->optlog_maxtgap[i] << ","
				<< This->optlog_min_avetgap[i] << "," << This->optlog_min_maxtgap[i] << std::endl;
		}
		ofs.close();
	}
#endif
	This->optlog_cnt = This->optlog_itr;
	This->optlog_itr = -1;
}

void ControlPanel::cb_btn_randrul3(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	This->rb_fix[CMODE_R]->setonly();
	This->rb_param[P_RULL]->setonly();

	double thres_avetgap = 0.3;
	double thres_maxtgap = 1.0;

	double dPb[MAX_DB_SIZE][MAX_CPCNT], dPb_val[MAX_CPCNT * 2];
	double dPb_max = 5.0 / 180.0 * M_PI; // random in [-5 degrees, 5 degrees]
	int dPb_size = set_dPb(dPb, c->Pcnt);
#ifdef PRINT_OPT_PROC
	for (int i = 0; i < dPb_size; i++) {
		for (int j = 0; j < c->Pcnt; j++) {
			std::cout << dPb[i][j] << ", ";
		}
		std::cout << std::endl;
	}
#endif

	// get parameters
	int opt_type = 2;
	int max_itr = 500;
	int max_trial_til_validrul = 100;
	bool double_check_rulings = true;
#ifdef INPUT_PARAMETERS
	std::cout << "type of parameter shift, 1: random, 2: table(default)" << std::endl;
	std::cin >> opt_type;
	if (0 < opt_type && opt_type << 3) {
		opt_type = 2;
	}
	std::cout << "opt_type = " << (opt_type == 1 ? "1: random" : "2: table") << std::endl;

	std::cout << "max iteration (default: 500)" << std::endl;
	std::cin >> max_itr;
	if (max_itr == 0) {
		max_itr = 500;
	}
	std::cout << "max_itr = " << max_itr << std::endl;

	std::cout << "max iteration for rulings (default: 100)" << std::endl;
	std::cin >> max_trial_til_validrul;
	if (max_trial_til_validrul == 0) {
		max_trial_til_validrul = 100;
	}
	std::cout << "max_trial_til_validrul = " << max_trial_til_validrul << std::endl;

	std::cout << "double check rulings, y(default) / n" << std::endl;
	char yn;
	std::cin >> yn;
	double_check_rulings = (yn == 'n' ? false : true);
#endif

	clock_t start_clock, end_clock;
	start_clock = clock();

	int minminval = c->alpha[c->Xcnt/2];	// angle with min gap in all rulings
	for (int i = 0; i < max_itr; i++, This->optlog_cnt = i)
	{
#if TARGET_GAP_TYPE == 0 // 0: average, 1: max
		if (ppm->avetgap < thres_avetgap) {
#else
		if (ppm->maxtgap < thres_maxtgap) {
#endif
			break; // i
		}
#if 0
		This->btn_randrul->do_callback();
#else
		double prevPx2d[MAX_CPCNT], prevPa[MAX_CPCNT], prevPbr[MAX_CPCNT], prevPbl[MAX_CPCNT], prev_avetgap, prev_maxtgap;
		This->optlog_avetgap[i] = -1;
		This->optlog_maxtgap[i] = -1;
		prev_avetgap = This->optlog_min_avetgap[i] = ppm->avetgap;
		prev_maxtgap = This->optlog_min_maxtgap[i] = ppm->maxtgap;
		memcpy(prevPx2d, c->Px2d, sizeof(double) * c->Pcnt);
		memcpy(prevPa, c->Pa, sizeof(double) * c->Pcnt);
		memcpy(prevPbl, c->Pbl, sizeof(double) * c->Pcnt);
		memcpy(prevPbr, c->Pbr, sizeof(double) * c->Pcnt);

		bool flg_rulOK = false;
		for (int j = 0; j < max_trial_til_validrul; j++)
		{
			// control point of the rulings
			int i0 = rand() % dPb_size;
			int i1 = rand() % dPb_size;
			for (int k = 0; k <= c->Pcnt; k++) {
				if (opt_type == 1) {
					dPb_val[k] = ((double)rand() / (double)RAND_MAX - 0.5) * 2.0 * dPb_max; // [-dPb_max, dPb_max]
				} else {
					dPb_val[k] = dPb[i0][k];
				}
				c->Pbl[k] = prevPbl[k] + dPb_val[k];
			}
			for (int k = 0; k <= c->Pcnt; k++) {
				if (opt_type == 1) {
					dPb_val[MAX_CPCNT + k] = ((double)rand() / (double)RAND_MAX - 0.5) * 2.0 * dPb_max; // [-dPb_max, dPb_max]
				} else {
					dPb_val[MAX_CPCNT + k] = dPb[i1][k];

				}
				c->Pbr[k] = prevPbr[k] + dPb_val[MAX_CPCNT + k];
			}

			// make new ruling
			crease::interpolate_spline_RulAngle(c->Pbl, c->Pcnt, c->betal, c->Xcnt, c->cotbl, c->cosbl, c->sinbl);
			for (int k = c->Xsidx; k <= c->Xeidx; k++) {
				c->rlx_cp[k] = c->cosbl[k] * c->Tx2d[k] - c->sinbl[k] * c->Ty2d[k];
				c->rly_cp[k] = c->sinbl[k] * c->Tx2d[k] + c->cosbl[k] * c->Ty2d[k];
				normalize_v2(&(c->rlx_cp[k]), &(c->rly_cp[k]));
				c->rllen[k] = ppm->pw * ppm->ph;
			}
			crease::interpolate_spline_RulAngle(c->Pbr, c->Pcnt, c->betar, c->Xcnt, c->cotbr, c->cosbr, c->sinbr);
			for (int k = c->Xsidx; k <= c->Xeidx; k++) {
				c->rrx_cp[k] = c->cosbr[k] * c->Tx2d[k] + c->sinbr[k] * c->Ty2d[k];
				c->rry_cp[k] = -c->sinbr[k] * c->Tx2d[k] + c->cosbr[k] * c->Ty2d[k];
				normalize_v2(&(c->rrx_cp[k]), &(c->rry_cp[k]));
				c->rrlen[k] = ppm->pw * ppm->ph;
			}
			c->calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);

			// check crossing
			int result = c->checkRulingCross();
			if (result) { continue; }
			result =c->checkRulingAngle();
			if (result) { continue; }
			result = c->checkRulCreaseCross();
			if (result) { continue; }
#if 0
			std::cout << "i = " << i << ", j = " << j << ": ";
			for (int k = 0; k < c->Pcnt; k++) {
				std::cout << std::setprecision(2) << dPb_val[k] << ", ";
			}
			std::cout << "/ ";
			for (int k = 0; k < c->Pcnt; k++) {
				std::cout << std::setprecision(2) << dPb_val[MAX_CPCNT + k] << ", ";
			}
			std::cout << std::endl;
#endif
			flg_rulOK = true;
			This->optlog_trial_til_validrul[i] = j;
			break; // j
		}

		// if no valid ruling
		if (!flg_rulOK)
		{
			std::cout << "i = " << i << ", ! flg_rulOK" << std::endl;
			This->optlog_trial_til_validrul[i] = -1;
			memcpy(c->Px2d, prevPx2d, sizeof(double) * c->Pcnt);
			memcpy(c->Pa, prevPa, sizeof(double) * c->Pcnt);
			memcpy(c->Pbl, prevPbl, sizeof(double) * c->Pcnt);
			memcpy(c->Pbr, prevPbr, sizeof(double) * c->Pcnt);
			This->optlog_avetgap[i] = -1;
			This->optlog_maxtgap[i] = -1;
			This->optlog_min_avetgap[i] = -1;
			This->optlog_min_maxtgap[i] = -1;
			This->optlog_cnt = i;
			break; // i
		}

		// try all folding angles & evaluate
#if 0
		This->btn_optfold2->do_callback();
#else
		int crcnt, lcrcnt, rcrcnt, dccnt, fccnt, tccnt;
		dccnt = ppm->dccnt; ppm->dccnt = 0;
		fccnt = ppm->fccnt; ppm->fccnt = 0;
		tccnt = ppm->tccnt; ppm->tccnt = 0;
		crcnt = ppm->crcnt; ppm->crcnt = 1;
		lcrcnt = ppm->lcrcnt; ppm->lcrcnt = 0;
		rcrcnt = ppm->rcrcnt; ppm->rcrcnt = 0;

		int minval = -180;	// angle with min gap within the rulings
		double min_avetgap = 10000;
		double min_maxtgap = 10000;
		for (int val = This->vs_xmang1->minimum(); val < This->vs_xmang1->maximum(); val++)
		{
			double mang = (double)val / 180.0 * M_PI;
			ppm->re_sidx = 0;
			ppm->re_eidx = c->Xcnt;
			int ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, mang, 0);

			if (ret == 0) {
				ppm->set_postproc_type(PPTYPE_PRICURVE);
				ppm->postproc();
				ppm->set_postproc_type(PPTYPE_UNDEF);

				// check crossing
				int result = 0;
				if (double_check_rulings) {
					result = c->checkRulingCross();
				}
				if (result) {
#ifdef PRINT_OPT_PROC
					std::cout << "i = " << i << ", ang = " << val << ", c->checkRulingCross() = " << result << std::endl;
#endif
				}
				else {
					if (double_check_rulings) {
						result = c->checkRulingAngle();
					}
					if (result) {
#ifdef PRINT_OPT_PROC
						std::cout << "i = " << i << ", ang = " << val << ", c->checkRulingAngle() = " << result << std::endl;
#endif
					}
					else {
						if (double_check_rulings) {
							result = c->checkRulCreaseCross();
						}
						if (result) {
#ifdef PRINT_OPT_PROC
							std::cout << "i = " << i << ", ang = " << val << ",  c->checkRulCreaseCross() = " << result << std::endl;
#endif
						}
						else {
							// recalculate average gap
							if (This->cb_optmat->value() && ppm->tgcnt > 3) {
								ppm->calcAvetgapMat();	// calculate ppm->avetgap, maxtgap;
							}
							else if (This->cb_optrot->value() && ppm->tgcnt > 3) {
								ppm->calcAvetgapRot();
							}
#if TARGET_GAP_TYPE == 0 // 0: average, 1: max
							if (min_avetgap > ppm->avetgap) {
#else
							if (min_maxtgap > ppm->maxtgap) {
#endif
								minval = val;
								min_avetgap = ppm->avetgap;
								min_maxtgap = ppm->maxtgap;
							}
						}
					}
				}
			} else {
				// no valid folding angle

			}
		}
		ppm->dccnt = dccnt;
		ppm->fccnt = fccnt;
		ppm->tccnt = tccnt;
		ppm->crcnt = crcnt;
		ppm->lcrcnt = lcrcnt;
		ppm->rcrcnt = rcrcnt;
#endif
		This->optlog_avetgap[i] = min_avetgap == 10000 ? -1 : min_avetgap;
		This->optlog_maxtgap[i] = min_maxtgap == 10000 ? -1 : min_maxtgap;
#if TARGET_GAP_TYPE == 0 // 0: average, 1: max
		if (min_avetgap > prev_avetgap)
#else
		if (min_maxtgap > prev_maxtgap)
#endif
		{
			memcpy(c->Px2d, prevPx2d, sizeof(double) * c->Pcnt);
			memcpy(c->Pa, prevPa, sizeof(double) * c->Pcnt);
			memcpy(c->Pbl, prevPbl, sizeof(double) * c->Pcnt);
			memcpy(c->Pbr, prevPbr, sizeof(double) * c->Pcnt);
			This->optlog_min_avetgap[i] = ppm->avetgap = prev_avetgap;
			This->optlog_min_maxtgap[i] = ppm->maxtgap = prev_maxtgap;
		} else {
			This->optlog_min_avetgap[i] = ppm->avetgap = min_avetgap;
			This->optlog_min_maxtgap[i] = ppm->maxtgap = min_maxtgap;
			minminval = minval;
		}
#endif

#ifdef PRINT_OPT_PROC
		std::cout << "i = " << i
			<< ", ave gap = " << This->optlog_avetgap[i]
			<< ", max gap = " << This->optlog_maxtgap[i]
			<< ", angle = " << minval
			<< ", min ave gap = " << This->optlog_min_avetgap[i]
			<< ", min max gap = " << This->optlog_min_maxtgap[i]
			<< ", minangle = " << minminval
			<< std::endl;
#endif
	} // i

#if 0
	This->btn_optfold2->do_callback();
#else
	std::cout << "min ave gap = " << ppm->avetgap << ", min max gap = " << ppm->maxtgap << std::endl;
	This->vs_xmang0->value(minminval);
	This->vs_xmang1->value(minminval);
	This->btn_R2TA0->do_callback(); // includes This->push_hist()
	ppm->calcAvetgapMat();	// calculate ppm->avetgap, maxtgap;
	std::cout << "ppm->avetgap = " << ppm->avetgap << ", maxtgap = " << ppm->maxtgap << std::endl;
#endif

	end_clock = clock();
	std::cout << "process time: " << (double)(end_clock - start_clock) / CLOCKS_PER_SEC << std::endl;

#ifdef DUMP_OPT_RESULT
	{
		std::ofstream ofs("./output/time0.csv", std::ios_base::app);
		ofs << (double)(end_clock - start_clock) / CLOCKS_PER_SEC << std::endl;
		ofs.close();
	}
	{
		std::ofstream ofs("./output/err_seq0_target.csv");
		ofs << "iter,trial_validrul,ave gap,max gap,min ave gap,min max gap" << std::endl;
		for (int i = 0; i < This->optlog_cnt; i++)
		{
			ofs << i << "," << This->optlog_trial_til_validrul[i] << ","
				<< This->optlog_avetgap[i] << "," << This->optlog_maxtgap[i] << ","
				<< This->optlog_min_avetgap[i] << "," << This->optlog_min_maxtgap[i] << std::endl;
		}
		ofs.close();
	}
#endif
}

void ControlPanel::cb_btn_randtrfold(Fl_Widget* wgt, void* idx)

{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	This->rb_fix[CMODE_B]->setonly();
	This->rb_param[P_TRSN]->setonly();
	This->cb_optmat->value(1);

	double prevPx2d[MAX_CPCNT], prevPy[MAX_CPCNT], prevPa[MAX_CPCNT], prev_avetgap, prev_maxtgap;
	memcpy(prevPx2d, c->Px2d, sizeof(double) * c->Pcnt);
	memcpy(prevPa, c->Pa, sizeof(double) * c->Pcnt);
	memcpy(prevPy, c->Py, sizeof(double) * c->Pcnt);

	for (int i = 0; i < 10000; i++)
	{
		// control point of the rulings
		for (int j = 0; j <= c->Pcnt; j++) {
			c->Pa[j] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI; // [-pi/2, pi/2];
		}
		for (int j = 0; j <= c->Pcnt; j++) {
			c->Py[j] = ((double)rand() / (double)RAND_MAX - 0.5) * 2.0 * 0.1; // [-0.1, 0.1];
		}

		// make new ruling
		c->calcCPA_X(1/*ppm->flg_interpolate*/, &ppm->rp);
		ppm->set_postproc_type(PPTYPE_PRICURVE);
		ppm->postproc();
		ppm->set_postproc_type(PPTYPE_UNDEF);

		// check crossing
		int result = c->checkRulingCross();
		if (result) {
#ifdef PRINT_RAND_PROC
			std::cout << "i = " << i << ", c->checkRulingCross() = " << result << std::endl;
#endif
		}
		else {
			result = c->checkRulingAngle();
			if (result) {
#ifdef PRINT_RAND_PROC
				std::cout << "i = " << i << ", c->checkRulingAngle() = " << result << std::endl;
#endif
			}
			else {
				result = c->checkRulCreaseCross();
				if (result) {
#ifdef PRINT_RAND_PROC
					std::cout << "i = " << i << ", c->checkRulCreaseCross() = " << result << std::endl;
#endif
				}
				else {
#ifdef PRINT_RAND_PROC
					std::cout << "i = " << i << ", OK" << result << std::endl;
#endif
				}
			}
		}
	} // i

#if 0
	{
		std::ofstream ofs("./output/time1.csv", std::ios_base::app);
		ofs << (double)(end_clock - start_clock) / CLOCKS_PER_SEC << std::endl;
		ofs.close();
	}
	{
		std::ofstream ofs("./output/err_seq1_target.csv");
		ofs << "iter,,ave gap,max gap,min ave gap,min max gap" << std::endl;
		for (int i = 0; i < This->optlog_cnt; i++)
		{
			if (This->optlog_avetgap[i] > 0.0) {
				ofs << i << ",," << This->optlog_avetgap[i] << "," << This->optlog_maxtgap[i]
					<< "," << This->optlog_min_avetgap[i] << "," << This->optlog_min_maxtgap[i] << std::endl;
			}
		}
		ofs.close();
	}
#endif

	memcpy(c->Px2d, prevPx2d, sizeof(double) * c->Pcnt);
	memcpy(c->Pa, prevPa, sizeof(double) * c->Pcnt);
	memcpy(c->Py, prevPy, sizeof(double) * c->Pcnt);
	c->calcCPA_X(1/*ppm->flg_interpolate*/, &ppm->rp);
	ppm->set_postproc_type(PPTYPE_PRICURVE);
	ppm->postproc();
	ppm->set_postproc_type(PPTYPE_UNDEF);
	This->btn_optmat->do_callback();

	This->gwin->redraw();
	This->gwin_cp->redraw();
}

void ControlPanel::cb_btn_opttrfold(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	This->rb_fix[CMODE_B]->setonly();
	This->rb_param[P_TRSN]->setonly();
	This->cb_optmat->value(1);

	//double thres_avetgap = 0.3;
	double thres_avetgap = 0.1;
	double thres_maxtgap = 1.0;
	double avetgap_org = ppm->avetgap;

	double dPa[MAX_DB_SIZE][MAX_CPCNT], dPa_val[MAX_CPCNT * 2], dPa_max = 5.0 / 180.0 * M_PI; // [-5 degrees, 5 degrees]
	double dPy[MAX_DB_SIZE][MAX_CPCNT], dPy_val[MAX_CPCNT * 2], dPy_max = 0.001; // [-0.001, 0.001]
	int dPa_size = set_dPa(dPa, c->Pcnt);
	int dPy_size = set_dPy(dPy, c->Pcnt);
#ifdef PRINT_OPT_PROC
	std::cout << "dPalpha :" << std::endl;
	for (int i = 0; i < dPa_size; i++) {
		for (int j = 0; j < c->Pcnt; j++) {
			std::cout << dPa[i][j] << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << "dPtau :" << std::endl;
	for (int i = 0; i < dPy_size; i++) {
		for (int j = 0; j < c->Pcnt; j++) {
			std::cout << dPy[i][j] << ", ";
		}
		std::cout << std::endl;
	}
#endif

	// get parameters
	int opt_type = 1;
	int max_itr = MAX_OPT_ITR;
#ifdef INPUT_PARAMETERS
	std::cout << "type of parameter shift, 1: random, 2: table(default)" << std::endl;
	std::cin >> opt_type;
	if (0 < opt_type && opt_type << 3) {
		opt_type = 2;
	}
	std::cout << "opt_type = " << (opt_type == 1 ? "1: random" : "2: table") << std::endl;

	std::cout << "max iteration (default: " << MAX_OPT_ITR << ")" << std::endl;
	std::cin >> max_itr;
	if (max_itr == 0) {
		max_itr = MAX_OPT_ITR;
	}
	std::cout << "max_itr = " << max_itr << std::endl;
#endif
	double prevPbl[MAX_CPCNT], prevPbr[MAX_CPCNT], minminval;
	memcpy(prevPbl, c->Pbl, sizeof(double)* c->Pcnt);
	memcpy(prevPbr, c->Pbr, sizeof(double)* c->Pcnt);
	minminval = c->alpha[c->Xcnt/2]*180.0/M_PI;

	clock_t start_clock, end_clock;
	start_clock = clock();

	int miniter = -1;
	for (int i = 0; i < max_itr; i++)
	{
		This->optlog_cnt = i;

#if TARGET_GAP_TYPE == 0 // 0: average, 1: max
		if (ppm->avetgap < thres_avetgap) {
#else
		if (ppm->maxtgap < thres_maxtgap) {
#endif
			break;
		}

		double prevPx2d[MAX_CPCNT], prevPy[MAX_CPCNT], prevPa[MAX_CPCNT], prev_avetgap, prev_maxtgap;
		This->optlog_avetgap[i] = -1;
		This->optlog_maxtgap[i] = -1;
		prev_avetgap = This->optlog_min_avetgap[i] = ppm->avetgap;
		prev_maxtgap = This->optlog_min_maxtgap[i] = ppm->maxtgap;
		memcpy(prevPx2d, c->Px2d, sizeof(double) * c->Pcnt);
		memcpy(prevPa, c->Pa, sizeof(double) * c->Pcnt);
		memcpy(prevPy, c->Py, sizeof(double) * c->Pcnt);

		// control point of the rulings
		int i0 = rand() % dPa_size;
		int i1 = rand() % dPy_size;
		for (int j = 0; j <= c->Pcnt; j++) {
			if (opt_type == 1) {
				dPa_val[j] = ((double)rand() / (double)RAND_MAX - 0.5) * 2.0 * dPa_max; // [-dPa_max, dPa_max]
			} else {
				dPa_val[j] = dPa[i0][j];
			}
			c->Pa[j] = prevPa[j] + dPa_val[j];
		}
		for (int j = 0; j <= c->Pcnt; j++) {
			if (opt_type == 1) {
				dPy_val[j] = ((double)rand() / (double)RAND_MAX - 0.5) * 2.0 * dPy_max; // [-dPy_max, dPy_max]
			}
			else {
				dPy_val[j] = dPy[i0][j];
			}
			c->Py[j] = prevPy[j] + dPy_val[j];
		}

		// make new ruling
		c->calcCPA_X(1/*ppm->flg_interpolate*/, &ppm->rp);
		ppm->set_postproc_type(PPTYPE_PRICURVE);
		ppm->postproc();
		ppm->set_postproc_type(PPTYPE_UNDEF);

		// check crossing
		int result = c->checkRulingCross();
		if (result) {
#ifdef PRINT_OPT_PROC
			std::cout << "i = " << i << ", c->checkRulingCross() = " << result << std::endl;
#endif
		}
		else {
			result = c->checkRulingAngle();
			if (result) {
#ifdef PRINT_OPT_PROC
				std::cout << "i = " << i << ", c->checkRulingAngle() = " << result << std::endl;
#endif
			}
			else {
				result = c->checkRulCreaseCross();
				if (result) {
#ifdef PRINT_OPT_PROC
					std::cout << "i = " << i << ", c->checkRulCreaseCross() = " << result << std::endl;
#endif
				}
				else {
					// recalculate average gap
					if (This->cb_optmat->value() && ppm->tgcnt > 3) {
						ppm->calcAvetgapMat();	// calculate ppm->avetgap, maxtgap;
					}
					else if (This->cb_optrot->value() && ppm->tgcnt > 3) {
						ppm->calcAvetgapRot();
					}
					This->optlog_avetgap[i] = ppm->avetgap;
					This->optlog_maxtgap[i] = ppm->maxtgap;

#if TARGET_GAP_TYPE == 0 // 0: average, 1: max
					if (ppm->avetgap < prev_avetgap)
#else
					if (ppm->maxtgap < prev_maxtgap)
#endif
					{
						miniter = i;
						This->optlog_min_avetgap[i] = ppm->avetgap;
						This->optlog_min_maxtgap[i] = ppm->maxtgap;
					}

				}
			}
		}

		if(miniter != i) {
			memcpy(c->Px2d, prevPx2d, sizeof(double) * c->Pcnt);
			memcpy(c->Pa, prevPa, sizeof(double) * c->Pcnt);
			memcpy(c->Py, prevPy, sizeof(double) * c->Pcnt);
			ppm->avetgap = prev_avetgap;
			ppm->maxtgap = prev_maxtgap;
		}
#ifdef PRINT_OPT_PROC
		std::cout << "i = " << i
			<< ", ave gap = " << This->optlog_avetgap[i]
			<< ", max gap = " << This->optlog_maxtgap[i]
			<< ", min ave gap = " << This->optlog_min_avetgap[i]
			<< ", min max gap = " << This->optlog_min_maxtgap[i]
			<< std::endl;
#endif
	} // i

#if 1
	std::cout << "iter_adopted = " << miniter << std::endl;
	std::cout << "min ave gap = " << ppm->avetgap << ", min max gap = " << ppm->maxtgap << std::endl;
	if ( ppm->avetgap < avetgap_org ) {
		c->calcCPA_X(1/*ppm->flg_interpolate*/, &ppm->rp);
		ppm->set_postproc_type(PPTYPE_PRICURVE);
		ppm->postproc();
		ppm->set_postproc_type(PPTYPE_UNDEF);
		This->btn_optmat->do_callback();
		ppm->calcAvetgapMat();	// calculate ppm->avetgap, maxtgap;
	}
	else {
		This->rb_fix[CMODE_R]->setonly();
		This->rb_param[P_RULL]->setonly();
		memcpy(c->Pbl, prevPbl, sizeof(double) * c->Pcnt);
		memcpy(c->Pbr, prevPbr, sizeof(double) * c->Pcnt);
		This->vs_xmang0->value(minminval);
		This->vs_xmang1->value(minminval);
		This->btn_R2TA0->do_callback(); // includes This->push_hist()
		ppm->calcAvetgapMat();	// calculate ppm->avetgap, maxtgap;
	}
	std::cout << "ppm->avetgap = " << ppm->avetgap << ", maxtgap = " << ppm->maxtgap << std::endl;
#endif

	end_clock = clock();
	std::cout << "process time: " << (double)(end_clock - start_clock) / CLOCKS_PER_SEC << std::endl;

#ifdef DUMP_OPT_RESULT
	{
		std::ofstream ofs("./output/time1.csv", std::ios_base::app);
		ofs << (double)(end_clock - start_clock) / CLOCKS_PER_SEC << std::endl;
		ofs.close();
	}
	{
		std::ofstream ofs("./output/err_seq1_target.csv");
		ofs << "iter,,ave gap,max gap,min ave gap,min max gap" << std::endl;
		for (int i = 0; i < This->optlog_cnt; i++)
		{
			if (This->optlog_avetgap[i] > 0.0) {
				ofs << i << ",," << This->optlog_avetgap[i] << "," << This->optlog_maxtgap[i]
					<< "," << This->optlog_min_avetgap[i] << "," << This->optlog_min_maxtgap[i] << std::endl;
			}
		}
		ofs.close();
	}
#endif
	This->gwin->redraw();
	This->gwin_cp->redraw();
}
