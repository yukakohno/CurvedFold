#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <windows.h>
#include <stdio.h>
#include <time.h>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

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
				ppm->calcAvetgap();
			} else if (This->cb_optrot->value() && ppm->tgcnt > 3) {
				ppm->calcAvetgapRot(); // calc gap with fixed origin
			}
			//std::cout << ", avetgap=" << ppm->avetgap;
			if (mintgap > ppm->avetgap) {
				minval = val;
				mintgap = ppm->avetgap;
				//std::cout << ", minval=" << minval << ", mintgap=" << mintgap;
			}
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
						ppm->calcAvetgap();
					} else if (This->cb_optrot->value() && ppm->tgcnt > 3) {
						ppm->calcAvetgapRot();
					}
					//std::cout << "ppm->avetgap=" << ppm->avetgap << std::endl;
					if (mintgap > ppm->avetgap) {
						minval = val;
						minli = li;
						minri = ri;
						mintgap = ppm->avetgap;
					}
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

// ------------------------- LIST RULINGS -------------------------------------------

#define IPP 10
int checkRulCross_recursive(papermodel *ppm, crease* c, int pidx, int rl, int *cnt, std::vector<double>& vPb)
{
	int idx_crossing = -1;

	if ( pidx == c->Pcnt-1 ) {

		// TODO: check starting ruling angle
		int i_start = 1;

		for (int i = i_start; i < 180; i+=IPP) { 
			if (rl < 0) {
				c->Pbl[pidx] = (double)i / 180.0 * M_PI;
				crease::interpolate_spline_RulAngle(c->Pbl, c->Pcnt, c->betal, c->Xcnt, c->cotbl, c->cosbl, c->sinbl);
				for (int j = c->Xsidx; j <= c->Xeidx; j++) {
					//c->betal[j] = atan(1.0 / c->cotbl[j]);
					//if (c->betal[j] < 0.0) {
					//	c->betal[j] += M_PI;
					//}
					c->rlx_cp[j] = c->cosbl[j] * c->Tx2d[j] - c->sinbl[j] * c->Ty2d[j];
					c->rly_cp[j] = c->sinbl[j] * c->Tx2d[j] + c->cosbl[j] * c->Ty2d[j];
					normalize_v2(&(c->rlx_cp[j]), &(c->rly_cp[j]));
					c->rllen[j] = ppm->pw * ppm->ph;
				}
			}else {
				c->Pbr[pidx] = (double)i / 180.0 * M_PI;
				crease::interpolate_spline_RulAngle(c->Pbr, c->Pcnt, c->betar, c->Xcnt, c->cotbr, c->cosbr, c->sinbr);
				for (int j = c->Xsidx; j <= c->Xeidx; j++) {
					//c->betar[j] = atan(1.0 / c->cotbr[j]);
					//if (c->betar[j] < 0.0) {
					//	c->betar[j] += M_PI;
					//}
					c->rrx_cp[j] = c->cosbr[j] * c->Tx2d[j] + c->sinbr[j] * c->Ty2d[j];
					c->rry_cp[j] = -c->sinbr[j] * c->Tx2d[j] + c->cosbr[j] * c->Ty2d[j];
					normalize_v2(&(c->rrx_cp[j]), &(c->rry_cp[j]));
					c->rrlen[j] = ppm->pw * ppm->ph;
				}
			}
			c->calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);		// ruling’·‚³C³i˜güA‹Èü‚Ü‚Å‚Ì’·‚³—Dæj

			// return index of rulings crossing, -1 if no crossing
			idx_crossing = c->checkRulingCross(rl); // -1:left, 1:right
			if ( idx_crossing > -1 ) {
#if 0
				std::cout << "idx_crossing=" << idx_crossing;
				if (rl<0) {
					for (int k = 0; k < c->Pcnt; k++) { std::cout << ", " << c->Pbl[k]; }
				} else {
					for (int k = 0; k < c->Pcnt; k++) { std::cout << ", " << c->Pbr[k]; }
				}
				std::cout << std::endl;
#endif
				break;
			}

			(*cnt)++;
#if 1
			{
				//std::cout << "OK, cnt="<< (*cnt);
				if (rl < 0) {
					for (int k = 0; k < c->Pcnt; k++) {
						//std::cout << ", " << c->Pbl[k];
						vPb.push_back(c->Pbl[k]);
					}
				}
				else {
					for (int k = 0; k < c->Pcnt; k++) {
						//std::cout << ", " << c->Pbr[k];
						vPb.push_back(c->Pbr[k]);
					}
				}
				//std::cout << std::endl;
			}
#endif
		}
	} else {
		double* Pb = NULL;
		if (rl < 0) {
			Pb = c->Pbl;
		}
		else {
			Pb = c->Pbr;
		}

		// TODO: check starting ruling angle
		int i_start = 1;

		for (int i = i_start; i < 180; i+=IPP) {
			Pb[pidx] = (double)i / 180.0 * M_PI;

			idx_crossing = checkRulCross_recursive(ppm, c, pidx + 1, rl, cnt, vPb);
			if ( idx_crossing > -1 ) {
				continue;
			}
		}
	}
	return idx_crossing;
}

void ControlPanel::cb_btn_makelist(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* _c = &(ppm->crs[0]);
	crease c; memcpy(&c, _c, sizeof(crease));
	int idx_crossing = -1, lcnt=0, rcnt=0;

	if (This->rb_listleft->value() != 0) {
		This->vPbl.clear();
		idx_crossing = checkRulCross_recursive(ppm, &c, 0, -1, &lcnt, This->vPbl);
		std::cout << "vPbl->size = " << This->vPbl.size() << ", " << This->vPbl.size() / c.Pcnt << std::endl;
	} else {
		This->vPbr.clear();
		idx_crossing = checkRulCross_recursive(ppm, &c, 0, 1, &rcnt, This->vPbr);
		std::cout << "vPbr->size = " << This->vPbr.size() << ", " << This->vPbr.size() / c.Pcnt << std::endl;
	}
}

void ControlPanel::idlerul(void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	if (This->rb_listleft->value() != 0) {

		if (This->idlerul_idx >= This->vPbl.size() / c->Pcnt) {
		This->idlerul_idx = 0;
	}
#if 0
	std::cout << "idlerul_idx = " << This->idlerul_idx;
	for (int i = 0; i < c->Pcnt; i++) {
		std::cout << ", " << This->vPbl[This->idlerul_idx * c->Pcnt + i];
	}
	std::cout << std::endl;
#endif
	for (int i = 0; i < c->Pcnt; i++) {
		c->Pbl[i] = This->vPbl[This->idlerul_idx * c->Pcnt + i];
	}
	crease::interpolate_spline_RulAngle(c->Pbl, c->Pcnt, c->betal, c->Xcnt, c->cotbl, c->cosbl, c->sinbl);
	for (int j = c->Xsidx; j <= c->Xeidx; j++) {
		c->rlx_cp[j] = c->cosbl[j] * c->Tx2d[j] - c->sinbl[j] * c->Ty2d[j];
		c->rly_cp[j] = c->sinbl[j] * c->Tx2d[j] + c->cosbl[j] * c->Ty2d[j];
		normalize_v2(&(c->rlx_cp[j]), &(c->rly_cp[j]));
		c->rllen[j] = ppm->pw * ppm->ph;
	}
	} else {
	if (This->idlerul_idx >= This->vPbr.size() / c->Pcnt) {
		This->idlerul_idx = 0;
	}
#if 0
	std::cout << "idlerul_idx = " << This->idlerul_idx;
	for (int i = 0; i < c->Pcnt; i++) {
		std::cout << ", " << This->vPbr[This->idlerul_idx * c->Pcnt + i];
	}
	std::cout << std::endl;
#endif
	for (int i = 0; i < c->Pcnt; i++) {
		c->Pbr[i] = This->vPbr[This->idlerul_idx * c->Pcnt + i];
	}
	crease::interpolate_spline_RulAngle(c->Pbr, c->Pcnt, c->betar, c->Xcnt, c->cotbr, c->cosbr, c->sinbr);
	for (int j = c->Xsidx; j <= c->Xeidx; j++) {
		c->rrx_cp[j] = c->cosbr[j] * c->Tx2d[j] + c->sinbr[j] * c->Ty2d[j];
		c->rry_cp[j] = -c->sinbr[j] * c->Tx2d[j] + c->cosbr[j] * c->Ty2d[j];
		normalize_v2(&(c->rrx_cp[j]), &(c->rry_cp[j]));
		c->rrlen[j] = ppm->pw * ppm->ph;
	}
	}
	c->calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);		// ruling’·‚³C³i˜güA‹Èü‚Ü‚Å‚Ì’·‚³—Dæj

	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
	Sleep(1);

	This->idlerul_idx++;
}
void ControlPanel::cb_btn_startrul(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	//This->acnt = (int)This->vs_xmang1->value();
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	if (!This->flg_idlerul_active) {
		if (This->rb_listleft->value() != 0 && This->vPbl.size() == 0) {
			return;
		}
		if (This->rb_listright->value() != 0 && This->vPbr.size() == 0) {
			return;
		}
		This->flg_idlerul_active = true;
		Fl::add_idle(ControlPanel::idlerul, This);
	}
}

void ControlPanel::cb_btn_stoprul(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	if (This->flg_idlerul_active) {
		This->flg_idlerul_active = false;
		Fl::remove_idle(ControlPanel::idlerul, This);
	}
}

// ------------------------- SAMPLE TARGET POINTS --------------------------------------


void ControlPanel::cb_btn_listtgt(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
	crease c0; memcpy(&c0, c, sizeof(crease));

	srand((unsigned int)time(NULL));
	int targetcnt = 100;

	This->vPbl.clear();
	This->vPbr.clear();

	// random rulings left -> This->vPbl
	for (int i = 0; i < 1000000; i++)
	{
		for (int j = 0; j < c->Pcnt; j++) {
			c->Pbl[j] = ((double)rand() / (double)RAND_MAX) * M_PI; // [0, PI]
		}
		crease::interpolate_spline_RulAngle(c->Pbl, c->Pcnt, c->betal, c->Xcnt, c->cotbl, c->cosbl, c->sinbl);
		for (int j = c->Xsidx; j <= c->Xeidx; j++) {
			c->rlx_cp[j] = c->cosbl[j] * c->Tx2d[j] - c->sinbl[j] * c->Ty2d[j];
			c->rly_cp[j] = c->sinbl[j] * c->Tx2d[j] + c->cosbl[j] * c->Ty2d[j];
			normalize_v2(&(c->rlx_cp[j]), &(c->rly_cp[j]));
			c->rllen[j] = ppm->pw * ppm->ph;
		}
		c->calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);		// ruling’·‚³C³i˜güA‹Èü‚Ü‚Å‚Ì’·‚³—Dæj
		int idx_crossing = c->checkRulingCross(-1); // -1:left, 1:right
		if (idx_crossing > -1) {
			continue;
		}
		for (int k = 0; k < c->Pcnt; k++) {
			This->vPbl.push_back(c->Pbl[k]);
		}
		if (This->vPbl.size() >= targetcnt * c->Pcnt) {
			break;
		}
	}

	// random rulings right -> This->vPbr
	for (int i = 0; i < 1000000; i++)
	{
		for (int j = 0; j < c->Pcnt; j++) {
			c->Pbr[j] = ((double)rand() / (double)RAND_MAX) * M_PI; // [0, PI]
		}
		crease::interpolate_spline_RulAngle(c->Pbr, c->Pcnt, c->betar, c->Xcnt, c->cotbr, c->cosbr, c->sinbr);
		for (int j = c->Xsidx; j <= c->Xeidx; j++) {
			c->rrx_cp[j] = c->cosbr[j] * c->Tx2d[j] + c->sinbr[j] * c->Ty2d[j];
			c->rry_cp[j] = -c->sinbr[j] * c->Tx2d[j] + c->cosbr[j] * c->Ty2d[j];
			normalize_v2(&(c->rrx_cp[j]), &(c->rry_cp[j]));
			c->rrlen[j] = ppm->pw * ppm->ph;
		}
		c->calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);		// ruling’·‚³C³i˜güA‹Èü‚Ü‚Å‚Ì’·‚³—Dæj
		int idx_crossing = c->checkRulingCross(1); // -1:left, 1:right
		if (idx_crossing > -1) {
			continue;
		}
		for (int k = 0; k < c->Pcnt; k++) {
			This->vPbr.push_back(c->Pbr[k]);
		}
		if (This->vPbr.size() >= targetcnt * c->Pcnt) {
			break;
		}
	}

	This->listtgcnt = 0;
	int cntl = (This->vPbl.size() / c->Pcnt);
	int cntr = (This->vPbr.size() / c->Pcnt);
	int targetcnt1 = cntl < cntr ? cntl : cntr;

	// This->vPbl + This->vPbl + random folding angle -> This->tgt*
	for (int j = 0; j < targetcnt1; j++) {

		for (int i = 0; i < c->Pcnt; i++) {
			c->Pbl[i] = This->vPbl[j * c->Pcnt + i];
		}
		for (int i = 0; i < c->Pcnt; i++) {
			c->Pbr[i] = This->vPbr[j * c->Pcnt + i];
		}
		bool flg = false;
		for (int i = 0; i < 1000000; i++)
		{
			double mang = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI; // [-PI/2, PI/2]
			ppm->re_sidx = 0;
			ppm->re_eidx = c->Xcnt;
			int ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, mang, 0);
			if (ret == 0) {
				flg = true;
				break;
			}
		}
		if (flg) {
			This->tgPcnt[This->listtgcnt] = c->Pcnt;
			memcpy(This->tgtPbl[This->listtgcnt], c->Pbl, sizeof(double) * MAX_CPCNT);
			memcpy(This->tgtPbr[This->listtgcnt], c->Pbr, sizeof(double) * MAX_CPCNT);
			memcpy(This->tgtPa[This->listtgcnt], c->Pa, sizeof(double) * MAX_CPCNT);

			ppm->set_postproc_type(PPTYPE_PRICURVE);
			ppm->postproc();
			ppm->set_postproc_type(PPTYPE_UNDEF);

			This->tgcnt[This->listtgcnt] = ppm->tgcnt;
			memcpy(This->tgx[This->listtgcnt], ppm->ogx, sizeof(double) * MAX_TGT_CNT);
			memcpy(This->tgy[This->listtgcnt], ppm->ogy, sizeof(double) * MAX_TGT_CNT);
			memcpy(This->tgz[This->listtgcnt], ppm->ogz, sizeof(double) * MAX_TGT_CNT);
			memcpy(This->ogx_cp[This->listtgcnt], ppm->ogx_cp, sizeof(double) * MAX_TGT_CNT);
			memcpy(This->ogy_cp[This->listtgcnt], ppm->ogy_cp, sizeof(double) * MAX_TGT_CNT);

			This->listtgcnt++;
			std::cout << "listtgcnt=" << This->listtgcnt << std::endl;

			if (This->listtgcnt >= MAX_TGTLST) {
				break;
			}
		}
	}
#if 1
	// file output
	for (int i = 0; i < This->listtgcnt; i++) {
		std::stringstream ss; ss << "./output/target" << std::setw(3) << std::setfill('0') << i << ".txt";
		std::ofstream ofs( ss.str() );
		for (int j = 0; j < This->tgcnt[i]; j++) {
			ofs << This->tgx[i][j] << "\t" << This->tgy[i][j] << "\t" << This->tgz[i][j]
				<< "\t" << This->ogx_cp[i][j] << "\t" << This->ogy_cp[i][j] << std::endl;
		}
		ofs.close();
		ss.str(""); ss << "./output/rulings" << std::setw(3) << std::setfill('0') << i << ".txt";
		ofs.open(ss.str());
		for (int j = 0; j < This->tgPcnt[i]; j++) { ofs << This->tgtPbl[i][j] << "\t"; }
		ofs << std::endl;
		for (int j = 0; j < This->tgPcnt[i]; j++) { ofs << This->tgtPbr[i][j] << "\t"; }
		ofs << std::endl;
		ofs.close();
	}
#endif
	memcpy(c, &c0, sizeof(crease));
	ppm->set_postproc_type(PPTYPE_PRICURVE);
	ppm->postproc();
	ppm->set_postproc_type(PPTYPE_UNDEF);
}

void ControlPanel::idletgt(void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	if (This->listtgcnt==0) {
		return;
	}

	if (This->idletgt_idx >= This->listtgcnt) {
		This->idletgt_idx = 0;
	}

	ppm->tgcnt = This->tgcnt[This->idletgt_idx];
	memcpy(ppm->tgx, This->tgx[This->idletgt_idx], sizeof(double)* MAX_TGT_CNT);
	memcpy(ppm->tgy, This->tgy[This->idletgt_idx], sizeof(double) * MAX_TGT_CNT);
	memcpy(ppm->tgz, This->tgz[This->idletgt_idx], sizeof(double) * MAX_TGT_CNT);
	memcpy(ppm->ogx_cp, This->ogx_cp[This->idletgt_idx], sizeof(double) * MAX_TGT_CNT);
	memcpy(ppm->ogy_cp, This->ogy_cp[This->idletgt_idx], sizeof(double) * MAX_TGT_CNT);
	ppm->getTgt2D3D();
#if 0
	for (int i = 0; i < ppm->tgcnt; i++) {
		std::cout << i << ":" << ppm->tgx[i] << "," << ppm->tgy[i] << "," << ppm->tgz[i]
			<< "," << ppm->ogx[i] << "," << ppm->ogy[i] << "," << ppm->ogz[i]
			<< "," << ppm->ogx_cp[i] << "," << ppm->ogy_cp[i] << std::endl;
	}
	std::cout << std::endl;
#endif
	This->btn_optfold2->do_callback();

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
	Sleep(1);

	This->idletgt_idx++;
}

void ControlPanel::cb_btn_starttgt(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;

	if (This->listtgcnt == 0) {
		return;
	}

	if (!This->flg_idletgt_active) {
		This->flg_idletgt_active = true;
		Fl::add_idle(ControlPanel::idletgt, This);
	}
}

void ControlPanel::cb_btn_stoptgt(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	if (This->flg_idletgt_active) {
		This->flg_idletgt_active = false;
		Fl::remove_idle(ControlPanel::idletgt, This);
	}
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
		This->optlog_err[This->optlog_itr] = -1.0;
		This->optlog_minerr[This->optlog_itr] = -1.0;
		return;
	}

	// try & evaluate
	This->btn_optfold2->do_callback();
	This->push_hist(c->Pcnt, mode, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);


	This->optlog_err[This->optlog_itr] = ppm->avetgap;

	if ( flg_optimize && ppm->avetgap > avetgap_prev ) {
		memcpy(c->Px2d, Px2d, sizeof(double) * c->Pcnt);
		memcpy(c->Pa, Pa, sizeof(double) * c->Pcnt);
		memcpy(c->Pbl, Pbl, sizeof(double) * c->Pcnt);
		memcpy(c->Pbr, Pbr, sizeof(double) * c->Pcnt);
		//int ret = This->access_hist(1, &(c->Pcnt), &mode, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);
		ppm->avetgap = avetgap_prev;
		if (This->optlog_itr < 0 || This->optlog_itr == This->optlog_itrmax) {
			This->btn_R2TA0->do_callback();
		}
	}
	This->optlog_minerr[This->optlog_itr] = ppm->avetgap;
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
	//double avetgap = 1000;

	This->optlog_itr = 0;
	This->optlog_itrmax = 100;
	clock_t start_clock, end_clock;
	start_clock = clock();
	for (int i = 0; i < 100; i++)
	{
		if (ppm->avetgap < thres_avetgap) {
			break;
		}
		This->btn_randrul->do_callback();
		//avetgap = ppm->avetgap;
		//std::cout << "i = " << i << ", avetgap = " << ppm->avetgap << std::endl;
		std::cout << "iter = " << This->optlog_itr
			<< ", err = " << This->optlog_err[This->optlog_itr]
			<< ", minerr = " << This->optlog_minerr[This->optlog_itr] << std::endl;
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
				<< This->optlog_err[i] << "," << This->optlog_minerr[i] << std::endl;
		}
		ofs.close();
	}
#endif
	This->optlog_itr = -1;
}
	}
}