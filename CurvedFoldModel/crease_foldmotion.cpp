#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "crease.h"
#include "util.h"

#define INTP_TYPE 1 // 1:LINEAR, 2:SPLINE

void crease::initFM()
{
	memcpy( Px2d_org, Px2d, sizeof(double)*Pcnt );
	memcpy( Px_org, Px, sizeof(double)*Pcnt );
	memcpy( Py_org, Py, sizeof(double)*Pcnt );
	memcpy( Pa_org, Pa, sizeof(double)*Pcnt );
	memcpy( m3_org, m3, sizeof(double)*16 );
	FM_fidx_min = 0;
	FM_fidx_org = MAX_FRAME/2; // 41/2=20
	FM_fidx_max = MAX_FRAME-1; // 40

	memset( flg_FM_Pt, 0, sizeof(int)*MAX_FRAME );
	memset( flg_FM_Pa, 0, sizeof(int)*MAX_FRAME );
	memset( flg_FM_m3, 0, sizeof(int)*MAX_FRAME );
	flg_FM_Pt[FM_fidx_org] = flg_FM_Pa[FM_fidx_org] = flg_FM_m3[FM_fidx_org] = 1;
	memcpy( &(FM_Pt[FM_fidx_org][0]), Py, sizeof(double)*Pcnt );
	memcpy( &(FM_Pa[FM_fidx_org][0]), Pa, sizeof(double)*Pcnt );
	memcpy( &(FM_m3[FM_fidx_org][0]), m3, sizeof(double)*16 );
#if 1
	flg_FM_Pt[FM_fidx_min] = flg_FM_Pa[FM_fidx_min] = flg_FM_m3[FM_fidx_min] = 1;
	for( int j=0; j<Pcnt; j++ ){
		FM_Pt[FM_fidx_min][j] = 0.0;
		FM_Pa[FM_fidx_min][j] = 0.0;
	}
	unit_m44( &(FM_m3[FM_fidx_min][0]) );

	flg_FM_Pt[FM_fidx_max] = flg_FM_Pa[FM_fidx_max] = flg_FM_m3[FM_fidx_max] = 1;
	for( int j=0; j<Pcnt; j++ ){
		FM_Pt[FM_fidx_max][j] = FM_Pt[FM_fidx_org][j];
		FM_Pa[FM_fidx_max][j] = M_PI*0.5;
	}
	memcpy( &(FM_m3[FM_fidx_max][0]), m3, sizeof(double)*16 );
	updateFM( 1 );
#else
	double dPt0[MAX_CPCNT], dPa0[MAX_CPCNT], dPt1[MAX_CPCNT], dPa1[MAX_CPCNT];
	for( int j=0; j<Pcnt; j++ ){
		dPt0[j] = FM_Pt[FM_fidx_org][j]/(double)(FM_fidx_org-FM_fidx_min);
		dPa0[j] = FM_Pa[FM_fidx_org][j]/(double)(FM_fidx_org-FM_fidx_min);
		dPt1[j] = 0.0;
		dPa1[j] = (M_PI*0.5 - FM_Pa[FM_fidx_org][j])/(double)(FM_fidx_max-FM_fidx_org);
	}
	for( int i=FM_fidx_min; i<FM_fidx_org; i++ ){
		int ii = FM_fidx_org-i;
		if( i==FM_fidx_min ){
			flg_FM_Pt[i] = flg_FM_Pa[i] = 1;
		} else {
			flg_FM_Pt[i] = flg_FM_Pa[i] = 0;
		}
		for( int j=0; j<Pcnt; j++ ){
			FM_Pt[i][j] = FM_Pt[FM_fidx_org][j] - dPt0[j]*ii;
			FM_Pa[i][j] = FM_Pa[FM_fidx_org][j] - dPa0[j]*ii;
		}
	}
	for( int i=FM_fidx_org+1; i<=FM_fidx_max; i++ ){
		int ii = i-FM_fidx_org;
		if( i==FM_fidx_max ){
			flg_FM_Pt[i] = flg_FM_Pa[i] = 1;
		} else {
			flg_FM_Pt[i] = flg_FM_Pa[i] = 0;
		}
		for( int j=0; j<Pcnt; j++ ){
			FM_Pt[i][j] = FM_Pt[FM_fidx_org][j] + dPt1[j]*ii;
			FM_Pa[i][j] = FM_Pa[FM_fidx_org][j] + dPa1[j]*ii;
		}
	}
#endif

#if 0
	{
		FILE *fp=fopen("initFM.csv","w");
		if( fp ){
			for( int i=FM_fidx_min; i<=FM_fidx_max; i++ ){
				for( int j=0; j<Pcnt; j++ ){
					fprintf(fp, ",%f", FM_Pt[i][j]);
				}
				fprintf(fp, ",");
				for( int j=0; j<Pcnt; j++ ){
					fprintf(fp, ",%f", FM_Pa[i][j]*180.0/M_PI);
				}
				fprintf(fp, "\n");
			}
			fclose(fp);
		}
	}
#endif
}

int crease::updateFM( int flg_interpolateCP )
{
	int ret=0;

#if INTP_TYPE==1 // LINEAR

	for( int j=0; j<Pcnt; j++ ){
		double a[MAX_FRAME], t[MAX_FRAME]; // for debug
		for( int i=FM_fidx_min, ii=i; i<=FM_fidx_max; i++ ){
			if( flg_FM_Pa[i]>0 && i!=ii){
				double da = (FM_Pa[i][j] - FM_Pa[ii][j])/(double)(i-ii);
				for( int k=ii+1; k<i; k++ ){
					FM_Pa[k][j] = a[k] = FM_Pa[ii][j] + da*(k-ii);
				}
				ii=i;
			}
		}
		for( int i=FM_fidx_min, ii=i; i<=FM_fidx_max; i++ ){
			if( flg_FM_Pt[i]>0 && i!=ii){
				double dt = (FM_Pt[i][j] - FM_Pt[ii][j])/(double)(i-ii);
				for( int k=ii+1; k<i; k++ ){
					FM_Pt[k][j] = t[k] = FM_Pt[ii][j] + dt*(k-ii);
				}
				ii=i;
			}
		}
	}

#if FM_MODE==2

	if( flg_interpolateCP ){
		for( int i=FM_fidx_min; i<=FM_fidx_max; i++ ){
			if( flg_FM_Pa[i]>0 ){
				interpolate_spline( FM_Pa[i], Pcnt, FM_alpha[i], Xcnt );
			}
			if( flg_FM_Pt[i]>0 ){
				interpolate_spline( FM_Pt[i], Pcnt, FM_tr[i], Xcnt );
			}
		}
	}
	for( int j=0; j<Xcnt; j++ ){
		double a[MAX_FRAME], t[MAX_FRAME]; // for debug
		for( int i=FM_fidx_min, ii=i; i<=FM_fidx_max; i++ ){
			if( flg_FM_Pa[i]>0 && i!=ii){
				double da = (FM_alpha[i][j] - FM_alpha[ii][j])/(double)(i-ii);
				for( int k=ii+1; k<i; k++ ){
					FM_alpha[k][j] = a[k] = FM_alpha[ii][j] + da*(k-ii);
				}
				ii=i;
			}
		}
		for( int i=FM_fidx_min, ii=i; i<=FM_fidx_max; i++ ){
			if( flg_FM_Pt[i]>0 && i!=ii){
				double dt = (FM_tr[i][j] - FM_tr[ii][j])/(double)(i-ii);
				for( int k=ii+1; k<i; k++ ){
					FM_tr[k][j] = t[k] = FM_tr[ii][j] + dt*(k-ii);
				}
				ii=i;
			}
		}
	}
	if( flg_interpolateCP==0 ){
		setP_t(-1);
		setP_a(-1);
	}
#endif

	// interpolate m3
	for( int i=FM_fidx_min, ii=i; i<=FM_fidx_max; i++ ){
		if( flg_FM_m3[i]>0 && i!=ii)
		{
			double inv_m3[16], dm[16], dq[4], dax[3], dang, da, dt[3];

			// FM_m3[i] * inv(FM_m3[ii]) -> dax[3], dang
			memcpy( dm, FM_m3[i], sizeof(double)*16 );
			memcpy( inv_m3, FM_m3[ii], sizeof(double)*16 );
			inv_m44( inv_m3 );
			mult_m44_n44( dm, inv_m3, 1);
			mat_quat( dm, dq );
			quat_axis_ang( dq, dax, &dang );

			// Å’ZŒo˜H‚É
			if( dang > M_PI ){
				dang -= 2.0*M_PI;
			}
			if( dang < -M_PI ){
				dang += 2.0*M_PI;
			}
			da = dang/(double)(i-ii);

			dt[0] = (FM_m3[i][12] - FM_m3[ii][12])/(double)(i-ii);
			dt[1] = (FM_m3[i][13] - FM_m3[ii][13])/(double)(i-ii);
			dt[2] = (FM_m3[i][14] - FM_m3[ii][14])/(double)(i-ii);

			for( int k=ii+1; k<i; k++ )
			{
				double ang = da*(k-ii);
				// dax[3], ang -> dq[4] -> dm[16] : FM_m3[k] = FM_m3[ii] * dm
				axis_ang_quat( dax, ang, dq);
				quat_mat( dq, dm );
				memcpy( FM_m3[k], FM_m3[ii], sizeof(double)*16 );
				mult_m44_n44( dm, FM_m3[k], 0 );

				FM_m3[k][12] = FM_m3[ii][12] + dt[0]*(k-ii);
				FM_m3[k][13] = FM_m3[ii][13] + dt[1]*(k-ii);
				FM_m3[k][14] = FM_m3[ii][14] + dt[2]*(k-ii);
			}
			ii=i;
		}
	}

#else if INTP_TYPE==2 // SPLINE

	Spline spa, spt;

	for( int j=0; j<Pcnt; j++ ){
		double a0[MAX_FRAME], t[MAX_FRAME];
		int acnt=0;
		for( int i=FM_fidx_min, ii=i; i<=FM_fidx_max; i++ ){
			if( flg_FM_Pa[i]>0 ){
				a0[acnt] = FM_Pa[i][j];
				for( int k=ii; k<=i; k++ ){
					t[k] = (double)(acnt-1) + (double)(k-ii)/(double)(i-ii);
				}
				ii=i;
				acnt++;
			}
		}
		spa.init( a0, acnt );
		for( int i=FM_fidx_min; i<=FM_fidx_max; i++ ){
			FM_Pa[i][j] = spa.culc( t[i] );
		}
	}

	for( int j=0; j<Pcnt; j++ ){
		double t0[MAX_FRAME], t[MAX_FRAME];
		int tcnt=0;
		for( int i=FM_fidx_min, ii=i; i<=FM_fidx_max; i++ ){
			if( flg_FM_Pt[i]>0 ){
				t0[tcnt] = FM_Pt[i][j];
				for( int k=ii; k<=i; k++ ){
					t[k] = (double)(tcnt-1) + (double)(k-ii)/(double)(i-ii);
				}
				ii=i;
				tcnt++;
			}
		}
		spt.init( t0, tcnt );
		for( int i=FM_fidx_min; i<=FM_fidx_max; i++ ){
			FM_Pt[i][j] = spt.culc( t[i] );
		}
	}

#endif
	return ret;
}

