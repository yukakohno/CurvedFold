#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <float.h>

#include "papermodel.h"
#include "util.h"

double papermodel::calcAveDiff( double *tr0, double *tr1, double max, int cvcnt )
{
	double diff=0.0, df;
	memset( cverr, 0.0, sizeof(double)*MAX_SPCNT );
	for( int i=3; i<cvcnt-3; i++ ){
		df = abs(tr0[i]-tr1[i]);
		diff += df;
		cverr[i] = df/max;
	}
end:
	return diff/(double)(cvcnt-6);
}

double papermodel::calcRulCross( double *xx, double *xy, double *rx, double *ry, double *rlen, int cvcnt ) // Xx2d, Xy2d, rlx_cp, rly_cp, rllen, cnt
{
	double ex[MAX_SPCNT],ey[MAX_SPCNT], diff=0.0;
	double constlen=50.0;
	memset( cverr, 0.0, sizeof(double)*MAX_SPCNT );

	for( int i=2;  i<cvcnt-2; i++ ){
		ex[i] = xx[i]+rlen[i]*rx[i];
		ey[i] = xy[i]+rlen[i]*ry[i];
		//ex[i] = xx[i]+constlen*rx[i];
		//ey[i] = xy[i]+constlen*ry[i];
	}
#if 0
	{
		FILE *fp = fopen("output/calcRulCross.csv","w");
		if(fp){
			for( int i=4;  i<cvcnt-4; i++ ){
				fprintf(fp, "%f,%f\n%f,%f\n\n", xx[i],xy[i],ex[i],ey[i]);
			}
			fclose(fp);
		}
	}
#endif
	for( int i=4;  i<cvcnt-4; i++ ){
		double ix,iy, l0,l1, l0_,l1_, a,b,c,d, area;
		int ret = intersectionOfLine( xx[i],xy[i],ex[i],ey[i], xx[i+1],xy[i+1],ex[i+1],ey[i+1], &ix, &iy, &l0, &l1 );
		if( ret<0 ){
			continue;
		}
		// ƒTƒ‰ƒX‚ÌŒöŽ® O(0,0),A(a,b),B(c,d) -> area of OAB = abs(a*d-b*c)*0.5
		//l0_ = rlen[i]  -l0;
		//l1_ = rlen[i+1]-l1;
		l0_ = constlen-l0;
		l1_ = constlen-l1;
		a = l0_*rx[i];
		b = l0_*ry[i];
		c = l1_*rx[i+1];
		d = l1_*ry[i+1];
		area = cverr[i] = sqrt(abs(a*d-b*c))+0.01;
		diff += area;
		cverr[i] /= hCrs_cvmaxerr[hCrs_evaltype];
	}
end:
	return diff*0.1 / (double)(cvcnt-6);
}

#define ACNT 12
int papermodel::optRulLen_SplineCP( crease *_c0, crease *_c1 )
{
	int ret=0, min_cpidx=-1, min_aidx=-1;
	double  dx[ACNT], dy[ACNT], dd=0.5/*mm*/, mindiff, diff;
	crease c0, c1;

	for( int i=0; i<12; i++ ){
		double a = (double)i * 2.0*M_PI / (double)ACNT;
		dx[i]= sin(a);
		dy[i]= cos(a);
	}
	switch( hCrs_evaltype ){
	case EVTYPE_TORSION: 
		mindiff = calcAveDiff( &(_c0->tr[_c1->Xsidx]), &(_c1->tr[_c1->Xsidx]), hCrs_cvmaxerr[hCrs_evaltype], _c1->Xeidx-_c1->Xsidx+1 );
		break;
	case EVTYPE_BETALR:
		mindiff = calcAveDiff( &(_c1->betal[_c1->Xsidx]), &(_c1->betar[_c1->Xsidx]), hCrs_cvmaxerr[hCrs_evaltype], _c1->Xeidx-_c1->Xsidx+1 );
		break;
	case EVTYPE_RULCROSS:
		if( _c1->rl<0 ){
			mindiff = calcRulCross( &(_c1->Xx2d[_c1->Xsidx]), &(_c1->Xy2d[_c1->Xsidx]), &(_c1->rlx_cp[_c1->Xsidx]), &(_c1->rly_cp[_c1->Xsidx]), &(_c1->rllen[_c1->Xsidx]), _c1->Xeidx-_c1->Xsidx+1 );
		} else if( _c1->rl>0 ){
			mindiff = calcRulCross( &(_c1->Xx2d[_c1->Xsidx]), &(_c1->Xy2d[_c1->Xsidx]),&( _c1->rrx_cp[_c1->Xsidx]), &(_c1->rry_cp[_c1->Xsidx]), &(_c1->rrlen[_c1->Xsidx]), _c1->Xeidx-_c1->Xsidx+1 );
		}
		break;
	default:
		mindiff = DBL_MAX;
		break;
	}

	for( int cnt=0; cnt<1; cnt++ ){
		for( int cpidx=0; cpidx<CCNT; cpidx++ ){
			for( int aidx=0; aidx<ACNT; aidx++ ){
				c0.copy( _c0 );
				c1.copy( _c1 );

				if( c1.rl<0 ){
					for( int i=0; i<MAX_SPCNT; i++ ){ c0.rlflg[i] = RETYPE_UNDEF; }
				} else {
					for( int i=0; i<MAX_SPCNT; i++ ){ c0.rrflg[i] = RETYPE_UNDEF; }
				}
				for( int i=0; i<MAX_SPCNT; i++ ){
					c1.rlflg[i] = c1.rrflg[i] = RETYPE_UNDEF;
				}
				c1.CPx[cpidx] += dx[aidx];
				c1.CPy[cpidx] += dy[aidx];

				if( c1.rl<0 ){
					ret = setRulLen_SplineCP( c0.Xx2d, c0.Xy2d, c0.rlx_cp, c0.rly_cp, c0.Xsidx, c0.Xeidx,
						c1.CPx, c1.CPy, c0.rllen, c0.rlflg );
					c1.setL2R( &c0 ); c1.Xsidx = _c1->Xsidx; c1.Xeidx = _c1->Xeidx;
					c1.calcTNB( flg_rectifyT );
					c1.calcLeft();
					crease::calcRLenP_( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]),
						&(c1.rlx_cp[c1.Xsidx]), &(c1.rly_cp[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1, psx,psy,pex,pey,
						&(c1.rllen[c1.Xsidx]), &(c1.rlflg[c1.Xsidx]), NULL /*"rrlen.csv"*/ );
					for( int j=0; j<tccnt; j++ ){
						crease::calcRLenC_( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]),
							&(c1.rlx_cp[c1.Xsidx]), &(c1.rly_cp[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1,
							tcurve[j].cvx, tcurve[j].cvy, tcurve[j].cvcnt, &(c1.rllen[c1.Xsidx]), &(c1.rlflg[c1.Xsidx]) );
					} // j
				} else if( c1.rl>0 ){
					ret = setRulLen_SplineCP( c0.Xx2d, c0.Xy2d, c0.rrx_cp, c0.rry_cp, c0.Xsidx, c0.Xeidx,
						c1.CPx, c1.CPy, c0.rrlen, c0.rrflg );
					c1.setR2L( &c0 ); c1.Xsidx = _c1->Xsidx; c1.Xeidx = _c1->Xeidx;
					c1.calcTNB( flg_rectifyT );
					c1.calcRight();
					crease::calcRLenP_( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]),
						&(c1.rrx_cp[c1.Xsidx]), &(c1.rry_cp[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1, psx,psy,pex,pey,
						&(c1.rrlen[c1.Xsidx]), &(c1.rrflg[c1.Xsidx]), NULL /*"rrlen.csv"*/ );
					for( int j=0; j<tccnt; j++ ){
						crease::calcRLenC_( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]),
							&(c1.rrx_cp[c1.Xsidx]), &(c1.rry_cp[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1,
							tcurve[j].cvx, tcurve[j].cvy, tcurve[j].cvcnt, &(c1.rrlen[c1.Xsidx]), &(c1.rrflg[c1.Xsidx]) );
					} // j
				}

				switch( hCrs_evaltype ){
		case EVTYPE_TORSION:
			diff = calcAveDiff( &(c0.tr[c1.Xsidx]), &(c1.tr[c1.Xsidx]), hCrs_cvmaxerr[hCrs_evaltype], c1.Xeidx-c1.Xsidx+1 );
			break;
		case EVTYPE_BETALR:
			diff = calcAveDiff( &(c1.betal[c1.Xsidx]), &(c1.betar[c1.Xsidx]), hCrs_cvmaxerr[hCrs_evaltype], c1.Xeidx-c1.Xsidx+1 );
			break;
		case EVTYPE_RULCROSS:
			if( c1.rl<0 ){
				diff = calcRulCross( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]), &(c1.rlx_cp[c1.Xsidx]), &(c1.rly_cp[c1.Xsidx]), &(c1.rllen[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1 );
			} else if( c1.rl>0 ){
				diff = calcRulCross( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]),&( c1.rrx_cp[c1.Xsidx]), &(c1.rry_cp[c1.Xsidx]), &(c1.rrlen[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1 );
			}
			break;
		default:
			//diff = 1.0;
			break;
				}

				if( mindiff > diff ){
					mindiff = diff;
					min_cpidx = cpidx;
					min_aidx = aidx;
				}

				printf( "cp_i:%d, a:%.0f, diff=%.2lf, mindiff=%.2lf: (%.0f,%.0f) (%.0f,%.0f) (%.0f,%.0f) (%.0f,%.0f) (%.0f,%.0f) (%.0f,%.0f) \n",
					cpidx, (double)aidx * 360.0 / (double)ACNT, diff, mindiff,
					c1.CPx[0], c1.CPy[0], c1.CPx[1], c1.CPy[1], c1.CPx[2], c1.CPy[2],
					c1.CPx[3], c1.CPy[3], c1.CPx[4], c1.CPy[4], c1.CPx[5], c1.CPy[5] );

			} // aidx
		} // cpidx

		if( min_cpidx>-1 && min_aidx>-1 )
		{
			if( _c1->rl<0 ){
				_c0->init_rlflg();
			} else {
				_c0->init_rrflg();
			}
			_c1->init_rlflg();
			_c1->init_rrflg();
			_c1->CPx[min_cpidx] += dx[min_aidx];
			_c1->CPy[min_cpidx] += dy[min_aidx];

			if( _c1->rl<0 ){
				ret = setRulLen_SplineCP( _c0->Xx2d, _c0->Xy2d, _c0->rlx_cp, _c0->rly_cp, _c0->Xsidx, _c0->Xeidx,
					_c1->CPx, _c1->CPy, _c0->rllen, _c0->rlflg );
				c1.setL2R( _c0 );
				c1.calcTNB( flg_rectifyT );
				c1.calcLeft();
			} else if( _c1->rl>0 ){
				ret = setRulLen_SplineCP( _c0->Xx2d, _c0->Xy2d, _c0->rrx_cp, _c0->rry_cp, _c0->Xsidx, _c0->Xeidx,
					_c1->CPx, _c1->CPy, _c0->rrlen, _c0->rrflg );
				_c1->setR2L( _c0 );
				_c1->calcTNB( flg_rectifyT );
				_c1->calcRight();
			}
		}
		hCrs_mindiff[hCrs_evaltype] = mindiff;

	} // cnt

end:
	return ret;
}

int papermodel::set_rectifyCrs12( int cridx, int cpidx, double dx, double dy )
{
	int ret = -1;
	hCrs_cridx = 0;
	hCrs_cpidx = -1;
	//hCrs_mindiff = 1000.0;
	hCrs_cpdx = dx;
	hCrs_cpdy = dy;

	for( int i=1; i<rcrcnt; i++ ){
		if( rcrs[i] == &(crs[cridx]) ){
			hCrs_cridx = i;
			hCrs_cpidx = cpidx;
			ret = 0; goto end;
		}
	}
	for( int i=1; i<lcrcnt; i++ ){
		if( lcrs[i] == &(crs[cridx]) ){
			hCrs_cridx = -i;
			hCrs_cpidx = cpidx;
			ret = 0; goto end;
		}
	}
end:
	return ret;
}

int papermodel::set_resetCrs12()
{
	hCrs_cridx = 0;
	hCrs_cpidx = -1;
	hCrs_cpdx = hCrs_cpdy = 0;
	return 0;
}

int papermodel::optRulLen_SplineCP1( crease *_c0, crease *_c1, int cpidx, double cpdx, double cpdy )
{
	int ret=0;
	double diff_org=0, diff_new=0;
	crease c0; c0.copy( _c0 );
	crease c1; c1.copy( _c1 );

	if( c1.rl<0 ){
		c0.init_rlflg();
	} else {
		c0.init_rrflg();
	}
	c1.init_rlflg();
	c1.init_rrflg();
	c1.CPx[cpidx] = cpdx;
	c1.CPy[cpidx] = cpdy;

	if( c1.rl<0 ){
		ret = setRulLen_SplineCP( c0.Xx2d, c0.Xy2d, c0.rlx_cp, c0.rly_cp, c0.Xsidx, c0.Xeidx,
			c1.CPx, c1.CPy, c0.rllen, c0.rlflg );
		crease::calcRLenP_( &(c0.Xx2d[c0.Xsidx]), &(c0.Xy2d[c0.Xsidx]),
			&(c0.rlx_cp[c0.Xsidx]), &(c0.rly_cp[c0.Xsidx]), c0.Xeidx-c0.Xsidx+1, psx,psy,pex,pey,
			&(c0.rllen[c0.Xsidx]), &(c0.rlflg[c0.Xsidx]), NULL /*"rrlen.csv"*/ );
		c1.setL2R( &c0 ); c1.Xsidx = _c1->Xsidx; c1.Xeidx = _c1->Xeidx;
		c1.calcTNB( flg_rectifyT );
		c1.calcLeft();
		crease::calcRLenP_( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]),
			&(c1.rlx_cp[c1.Xsidx]), &(c1.rly_cp[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1, psx,psy,pex,pey,
			&(c1.rllen[c1.Xsidx]), &(c1.rlflg[c1.Xsidx]), NULL /*"rrlen.csv"*/ );
		for( int j=0; j<tccnt; j++ ){
			crease::calcRLenC_( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]),
				&(c1.rlx_cp[c1.Xsidx]), &(c1.rly_cp[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1,
				tcurve[j].cvx, tcurve[j].cvy, tcurve[j].cvcnt, &(c1.rllen[c1.Xsidx]), &(c1.rlflg[c1.Xsidx]) );
		} // j
	} else if( c1.rl>0 ){
		ret = setRulLen_SplineCP( c0.Xx2d, c0.Xy2d, c0.rrx_cp, c0.rry_cp, c0.Xsidx, c0.Xeidx,
			c1.CPx, c1.CPy, c0.rrlen, c0.rrflg );
		crease::calcRLenP_( &(c0.Xx2d[c0.Xsidx]), &(c0.Xy2d[c0.Xsidx]),
			&(c0.rrx_cp[c0.Xsidx]), &(c0.rry_cp[c0.Xsidx]), c0.Xeidx-c0.Xsidx+1, psx,psy,pex,pey,
			&(c0.rrlen[c0.Xsidx]), &(c0.rrflg[c0.Xsidx]), NULL /*"rrlen.csv"*/ );
		c1.setR2L( &c0 ); c1.Xsidx = _c1->Xsidx; c1.Xeidx = _c1->Xeidx;
		c1.calcTNB( flg_rectifyT );
		c1.calcRight();
		crease::calcRLenP_( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]),
			&(c1.rrx_cp[c1.Xsidx]), &(c1.rry_cp[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1, psx,psy,pex,pey,
			&(c1.rrlen[c1.Xsidx]), &(c1.rrflg[c1.Xsidx]), NULL /*"rrlen.csv"*/ );
		for( int j=0; j<tccnt; j++ ){
			crease::calcRLenC_( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]),
				&(c1.rrx_cp[c1.Xsidx]), &(c1.rry_cp[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1,
				tcurve[j].cvx, tcurve[j].cvy, tcurve[j].cvcnt, &(c1.rrlen[c1.Xsidx]), &(c1.rrflg[c1.Xsidx]) );
		} // j
	}

	switch( hCrs_evaltype ){
		case EVTYPE_TORSION:
			diff_org = calcAveDiff( &(_c0->tr[_c1->Xsidx]), &(_c1->tr[_c1->Xsidx]), hCrs_cvmaxerr[hCrs_evaltype], _c1->Xeidx-_c1->Xsidx+1 );
			diff_new = calcAveDiff( &(c0.tr[c1.Xsidx]), &(c1.tr[c1.Xsidx]), hCrs_cvmaxerr[hCrs_evaltype], c1.Xeidx-c1.Xsidx+1 );
			break;
		case EVTYPE_BETALR:
			diff_org = calcAveDiff( &(_c1->betal[_c1->Xsidx]), &(_c1->betar[_c1->Xsidx]), hCrs_cvmaxerr[hCrs_evaltype], _c1->Xeidx-_c1->Xsidx+1 );
			diff_new = calcAveDiff( &(c1.betal[c1.Xsidx]), &(c1.betar[c1.Xsidx]), hCrs_cvmaxerr[hCrs_evaltype], c1.Xeidx-c1.Xsidx+1 );
			break;
		case EVTYPE_RULCROSS:
			if( c1.rl<0 ){
				diff_org = calcRulCross( &(_c1->Xx2d[_c1->Xsidx]), &(_c1->Xy2d[_c1->Xsidx]), &(_c1->rlx_cp[_c1->Xsidx]), &(_c1->rly_cp[_c1->Xsidx]), &(_c1->rllen[_c1->Xsidx]), _c1->Xeidx-_c1->Xsidx+1 );
				diff_new = calcRulCross( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]), &(c1.rlx_cp[c1.Xsidx]), &(c1.rly_cp[c1.Xsidx]), &(c1.rllen[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1 );
			} else if( c1.rl>0 ){
				diff_org = calcRulCross( &(_c1->Xx2d[_c1->Xsidx]), &(_c1->Xy2d[_c1->Xsidx]),&( _c1->rrx_cp[_c1->Xsidx]), &(_c1->rry_cp[_c1->Xsidx]), &(_c1->rrlen[_c1->Xsidx]), _c1->Xeidx-_c1->Xsidx+1 );
				diff_new = calcRulCross( &(c1.Xx2d[c1.Xsidx]), &(c1.Xy2d[c1.Xsidx]),&( c1.rrx_cp[c1.Xsidx]), &(c1.rry_cp[c1.Xsidx]), &(c1.rrlen[c1.Xsidx]), c1.Xeidx-c1.Xsidx+1 );
			}
			break;
		default:
			diff_org = diff_new = 1.0;
			break;
	}

	if( diff_org > diff_new - hCrs_mgndiff[hCrs_evaltype] ){
		hCrs_mindiff[hCrs_evaltype] = diff_new;
		_c0->copy( &c0 );
		_c1->copy( &c1 );	c1.Xsidx = 0; c1.Xeidx = XCNT;
	} else {
		hCrs_mindiff[hCrs_evaltype] = diff_org;
	}

	printf( "cr_i:%d, cp_i:%d, diff=%f, mindiff=%f: (%.0f,%.0f) (%.0f,%.0f) (%.0f,%.0f) (%.0f,%.0f) (%.0f,%.0f) (%.0f,%.0f)\n",
		hCrs_cridx, hCrs_cpidx, diff_new, hCrs_mindiff[hCrs_evaltype],
		c1.CPx[0], c1.CPy[0], c1.CPx[1], c1.CPy[1], c1.CPx[2], c1.CPy[2],
		c1.CPx[3], c1.CPy[3], c1.CPx[4], c1.CPy[4], c1.CPx[5], c1.CPy[5] );

end:
	return ret;
}

