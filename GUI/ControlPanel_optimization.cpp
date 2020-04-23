#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <windows.h>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

void ControlPanel::cb_btn_optfold(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);

	int ret = ppm->optFold();
	if (This->cb_optmat->value() && ppm->tgcnt > 3) {
		ppm->optMat(CMODE_B); // B: Ü‚èü‚ÆÜ‚èŠp“x‚©‚çA3D‹Èü‚ð‹‚ß‚é
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
		ppm->re_sidx = 0;
		ppm->re_eidx = c->Xcnt;
		int ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, mang);

		if (ret == 0) {
			ppm->set_postproc_type(PPTYPE_PRICURVE);
			ppm->postproc();
			ppm->set_postproc_type(PPTYPE_UNDEF);
			if (This->cb_optmat->value() && ppm->tgcnt > 3) {
				ppm->optMat(CMODE_R);
			}
			if (mintgap > ppm->avetgap) {
				minval = val;
				mintgap = ppm->avetgap;
			}
		}
	}
	if (minval > -180) {
		This->vs_xmang0->value(minval);
		This->vs_xmang1->value(minval);
		This->btn_R2TA0->do_callback();
	}
	This->push_hist(c->Pcnt, CMODE_R, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);
}

void ControlPanel::cb_btn_opttr(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);

	int ret = ppm->optTorsion();

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}
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

			if ( idx_crossing < (double)pidx/(double)c->Pcnt * (double)c->Xcnt ) {
				break;
			}
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

	idx_crossing = checkRulCross_recursive( ppm, &c, 0, -1, &lcnt, This->vPbl);
	std::cout << "vPbl->size = " << This->vPbl.size() << ", " << This->vPbl.size()/c.Pcnt << std::endl;

	//idx_crossing = checkRulCross_recursive( ppm, &c, 0, 1, &rcnt, This->vPbr);
	//std::cout << "vPbr->size = " << This->vPbr.size() << ", " << This->vPbr.size()/c.Pcnt << std::endl;

}

void ControlPanel::idlerul(void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
#if 1
	if (This->idlerul_idx >= This->vPbl.size()/ c->Pcnt) {
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
#else
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
#endif
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

