#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "crease.h"
#include "Spline.h"

int crease::interpolate_spline( double *_P, int _Pcnt, double *_X, int _Xcnt )
{
	int i, ret=0;
	Spline sp;
	sp.init(_P, _Pcnt);
	for( i=0; i<_Xcnt; i++ ){
		double t = (double)i/(double)_Xcnt * (double)(_Pcnt-1);
		_X[i] = sp.culc(t);
	}
end:
	return ret;
}
int crease::interpolate_spline_RulAngle( double *_P, int _Pcnt,
										double *_B, int _Xcnt, double *cotB, double *cosB, double *sinB )
{
	int i, ret=0;
	ret = interpolate_spline( _P, _Pcnt, _B, _Xcnt );

	for( i=0; i<_Xcnt; i++ ){
		if( _B[i] < 0.0 ){
			_B[i] += M_PI;
		}
		if( _B[i] > M_PI ){
			_B[i] -= M_PI;
		}
		cotB[i] = 1.0/tan(_B[i]);
		cosB[i] = cos(_B[i]);
		sinB[i] = sin(_B[i]);
	}
end:
	return ret;
}

int crease::setP_k( int Pidx=-1 )
{
	int i, ret=0, idxk[MAX_CPCNT];

	for( i=0; i<Pcnt; i++ ){
		idxk[i] = (int)((double)i*(double)Xcnt/(double)(Pcnt-1) + 0.5);
	}
	if( idxk[0] < CEMGN ){ idxk[0] = CEMGN; }
	if( idxk[Pcnt-1] > Xcnt-CEMGN-1 ){ idxk[Pcnt-1] = Xcnt-CEMGN-1; }

	if( Pidx<0 || Pcnt<=Pidx ){
		for( i=0; i<Pcnt; i++ ){
			Px[i] = kv[idxk[i]];
		}
	} else {
		Px[Pidx] = kv[idxk[Pidx]];
	}

end:
	return ret;
}

int crease::setP_t( int Pidx=-1 )
{
	int i, ret=0, idxt[MAX_CPCNT];

	for( i=0; i<Pcnt; i++ ){
		idxt[i] = (int)((double)i*(double)Xcnt/(double)(Pcnt-1) + 0.5);
	}
	if( idxt[0] < CEMGN+1 ){ idxt[0] = CEMGN+1; }
	if( idxt[Pcnt-1] > Xcnt-CEMGN-2 ){ idxt[Pcnt-1] = Xcnt-CEMGN-2; }

	if( Pidx<0 || Pcnt<=Pidx ){
		for( i=0; i<Pcnt; i++ ){
			Py[i] = tr[idxt[i]];
		}
	} else {
		Py[Pidx] = tr[idxt[Pidx]];
	}

end:
	return ret;
}

int crease::setP_a( int Pidx=-1 )
{
	int i, ret=0, idxa[MAX_CPCNT];

	for( i=0; i<Pcnt; i++ ){
		idxa[i] = (int)((double)i*(double)Xcnt/(double)(Pcnt-1) + 0.5);
	}
	if( idxa[0] < CEMGN ){ idxa[0] = CEMGN; }
	if( idxa[Pcnt-1] > Xcnt-CEMGN-1 ){ idxa[Pcnt-1] = Xcnt-CEMGN-1; }

	if( Pidx<0 || Pcnt<=Pidx ){
		for( i=0; i<Pcnt; i++ ){
			Pa[i] = alpha[idxa[i]];
		}
	} else {
		Pa[Pidx] = alpha[idxa[Pidx]];
	}

end:
	return ret;
}

int crease::setP_k2( int Pidx=-1 )
{
	int i, ret=0, idxk[MAX_CPCNT];

	for( i=0; i<Pcnt; i++ ){
		idxk[i] = (int)((double)i*(double)Xcnt/(double)(Pcnt-1) + 0.5);
	}
	if( idxk[0] < CEMGN ){ idxk[0] = CEMGN; }
	if( idxk[Pcnt-1] > Xcnt-CEMGN-1 ){ idxk[Pcnt-1] = Xcnt-CEMGN-1; }

	if( Pidx<0 || Pcnt<=Pidx ){
		for( i=0; i<Pcnt; i++ ){
			Px2d[i] = k2d[idxk[i]];
		}
	} else {
		Px2d[Pidx] = k2d[idxk[Pidx]];
	}
end:
	return ret;
}

int crease::setP_Bl( int Pidx=-1 )
{
	int i, ret=0, idxt[MAX_CPCNT];

	for( i=0; i<Pcnt; i++ ){
		idxt[i] = (int)((double)i*(double)Xcnt/(double)(Pcnt-1) + 0.5);
	}
	if( idxt[0] < CEMGN+1 ){ idxt[0] = CEMGN+1; }
	if( idxt[Pcnt-1] > Xcnt-CEMGN-2 ){ idxt[Pcnt-1] = Xcnt-CEMGN-2; }

	if( Pidx<0 || Pcnt<=Pidx ){
		for( i=0; i<Pcnt; i++ ){
			Pbl[i] = betal[idxt[i]];
		}
	} else {
		Pbl[Pidx] = betal[idxt[Pidx]];
	}
end:
	return ret;
}

int crease::setP_Br( int Pidx=-1 )
{
	int i, ret=0, idxt[MAX_CPCNT];

	for( i=0; i<Pcnt; i++ ){
		idxt[i] = (int)((double)i*(double)Xcnt/(double)(Pcnt-1) + 0.5);
	}
	if( idxt[0] < CEMGN+1 ){ idxt[0] = CEMGN+1; }
	if( idxt[Pcnt-1] > Xcnt-CEMGN-2 ){ idxt[Pcnt-1] = Xcnt-CEMGN-2; }

	if( Pidx<0 || Pcnt<=Pidx ){
		for( i=0; i<Pcnt; i++ ){
			Pbr[i] = betar[idxt[i]];
		}
	} else {
		Pbr[Pidx] = betar[idxt[Pidx]];
	}
end:
	return ret;
}

// P*をスプライン補間（sに関して等間隔）
int crease::interpolate_spline2( int mode )	// mode=0:CP&3D 1:CP 2:3D
{
	int i, j, ret=0, Tcnt = 10000;
	double *Tx2d=NULL, *Ty2d=NULL, *Tx=NULL, *Ty=NULL, *Tz=NULL, *L2d=NULL, *L=NULL;
	Spline spx2, spy2, spx, spy, spz;

	if( mode==0 || mode==1 ){
		Tx2d = new double[Tcnt];
		Ty2d = new double[Tcnt];
		L2d = new double[Tcnt];
		spx2.init(Px2d, Pcnt);
		spy2.init(Py2d, Pcnt);
	}
	if( mode==0 || mode==2 ){
		Tx = new double[Tcnt];
		Ty = new double[Tcnt];
		Tz = new double[Tcnt];
		L = new double[Tcnt];
		spx.init(Px, Pcnt);
		spy.init(Py, Pcnt);
		spz.init(Pz, Pcnt);
	}

	for( i=0; i<Tcnt; i++ ){
		double t = (double)i/(double)Tcnt * (double)(Pcnt-1);
		if( mode==0 || mode==1 ){

			Tx2d[i] = spx2.culc(t);
			Ty2d[i] = spy2.culc(t);
		}
		if( mode==0 || mode==2 ){
			Tx[i] = spx.culc(t);
			Ty[i] = spy.culc(t);
			Tz[i] = spz.culc(t);
		}
		if( i==0 ){
			if( mode==0 || mode==1 ){
				L2d[i] = 0.0;
			}
			if( mode==0 || mode==2 ){
				L[i] = 0.0;
			}
			continue;
		}
		if( mode==0 || mode==1 ){
			double dx2 = Tx2d[i] - Tx2d[i-1];
			double dy2 = Ty2d[i] - Ty2d[i-1];
			L2d[i] = sqrt( dx2*dx2 + dy2*dy2 );
			L2d[i] += L2d[i-1];
		}
		if( mode==0 || mode==2 ){
			double dx = Tx[i] - Tx[i-1];
			double dy = Ty[i] - Ty[i-1];
			double dz = Tz[i] - Tz[i-1];
			L[i] = sqrt( dx*dx + dy*dy + dz*dz );
			L[i] += L[i-1];
		}
	}

	// 媒介変数(0～num-1の実数）に対する値を計算
	Xcnt = XcntOrg;
	for( i=0; i<Xcnt; i++ ){
		double t0 = 0.0;
		double t1 = (double)i/(double)Xcnt;
		if( mode==0 || mode==1 ){
			for( j=0; j<Tcnt; j++ ){
				if( t1 < L2d[j]/L2d[Tcnt-1] ){
					double dj = (double)j;
					// jとj-1 の間で線形補間
					if(j>0){
						dj = (double)(j-1) + (t1*L2d[Tcnt-1]-L2d[j-1])/(L2d[j]-L2d[j-1]);
					}
					t0 = dj/(double)Tcnt * (double)(Pcnt-1); // 0 <= t0 <= Pcnt-1
					break;
				}
			}
		} else if( mode==2 ){
			for( j=0; j<Tcnt; j++ ){
				if( t1 < L[j]/L[Tcnt-1] ){
					double dj = (double)j;
					// jとj-1 の間で線形補間
					if(j>0){
						dj = (double)(j-1) + (t1*L[Tcnt-1]-L[j-1])/(L[j]-L[j-1]);
					}
					t0 = dj/(double)Tcnt * (double)(Pcnt-1); // 0 <= t0 <= Pcnt-1
					break;
				}
			}
		}
		if( mode==0 || mode==1 ){
			Xx2d[i] = spx2.culc(t0);
			Xy2d[i] = spy2.culc(t0);
		}
		if( mode==0 || mode==2 ){
			Xx[i] = spx.culc(t0);
			Xy[i] = spy.culc(t0);
			Xz[i] = spz.culc(t0);
		}
	}

#if 0
	{
		FILE *fp=fopen("output/dump_interpolate_spline2.csv","w");
		if( fp ){
			for( i=0; i<Xcnt; i++ ){
				fprintf(fp, "%f,%f,%f,%f,%f\n", Xx2d[i], Xy2d[i], Xx[i], Xy[i], Xz[i]);
			}
			fclose(fp); fp=NULL;
		}
	}
#endif
end:
	if( Tx2d ){ delete[] Tx2d; Tx2d=NULL; } 
	if( Ty2d ){ delete[] Ty2d; Ty2d=NULL; } 
	if( Tx ){ delete[] Tx; Tx=NULL; } 
	if( Ty ){ delete[] Ty; Ty=NULL; } 
	if( Tz ){ delete[] Tz; Tz=NULL; } 
	if( L2d ){ delete[] L2d; L2d=NULL; } 
	if( L ){ delete[] L; L=NULL; } 
	return ret;
}
