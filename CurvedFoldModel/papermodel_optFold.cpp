#include <stdio.h>
#include <iostream>
#include <fstream>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include "papermodel.h"
#include "util.h"
#include "polynomial.h"

int papermodel::saveTgt( char *fname )
{
	int ret=0;
	FILE *fp=NULL;
	char buf[256];

	if( fname==NULL ){
		ret = -1; return ret;
	}
	fp = fopen( fname, "w" );
	if( !fp ){
		ret = -1; return ret;
	}

	if (tgcnt == 0) {
		for (int j = 20; j < ph; j += 20) {
			for (int i = 20; i < pw; i += 20) {
				ogx_cp[tgcnt] = i;
				ogy_cp[tgcnt] = j;
				tgcnt++;
			}
		}
		getTgt2D3D();
	}

#if 1
	for( int i=0; i<tgcnt; i++ ){
		fprintf( fp, "%f %f %f %f %f\n", ogx[i], ogy[i], ogz[i], ogx_cp[i], ogy_cp[i] );
	}
#else // data for debug
	double a = 10.0/180.0*M_PI;
	double m[16] = {cos(a),sin(a),0,0,
					-sin(a),cos(a),0,0,
					0,0,1,10,
					0,0,0,1};
	for( int i=0; i<tgcnt; i++ ){
		//fprintf( fp, "%f %f %f %f %f\n", ogx[i], ogy[i], ogz[i], ogx_cp[i], ogy_cp[i] );
		double x = m[0]*ogx[i] + m[1]*ogy[i] + m[2]*ogz[i] + m[3];
		double y = m[4]*ogx[i] + m[5]*ogy[i] + m[6]*ogz[i] + m[7];
		double z = m[8]*ogx[i] + m[9]*ogy[i] + m[10]*ogz[i] + m[11];
		fprintf( fp, "%f %f %f %f %f\n", x, y, z, ogx_cp[i], ogy_cp[i] );
	}
#endif
	fclose( fp );

	tgcnt = 0;

	return ret;
}

int papermodel::loadTgt( char *fname )
{
	int ret=0;
	FILE *fp=NULL;
	char buf[256];

	if( fname==NULL ){
		ret = -1; return ret;
	}
	fp = fopen( fname, "r" );
	if( !fp ){
		ret = -1; return ret;
	}

	tgcnt = 0;
	while( fgets( buf, 1024, fp) ){
		if( buf[0] == '#' ){
			continue;
		}
		int ret = sscanf( buf, "%lf %lf %lf %lf %lf", &tgx[tgcnt], &tgy[tgcnt], &tgz[tgcnt], &ogx_cp[tgcnt], &ogy_cp[tgcnt] );
		if(ret<5){
			continue;
		}
		//printf( "%d: %f, %f, %f, %f, %f\n", tgcnt, tgx[tgcnt], tgy[tgcnt], tgz[tgcnt], ogx_cp[tgcnt], ogy_cp[tgcnt] );
		tgcnt++;
	}

	fclose( fp );

	getTgt2D3D();

	return ret;
}

int papermodel::loadTgtMask(char* fname_tpt, char* fname_mask)
{
	int ret = 0;
	FILE* fp = NULL;
	char buf[1024];
	int tmask[MAX_TGT_CNT]; memset(tmask, 1, sizeof(MAX_TGT_CNT));
	int tmaskcnt = 0;

	if (fname_tpt == NULL || fname_mask == NULL) {
		ret = -1; return ret;
	}
	fp = fopen(fname_mask, "r");
	if (!fp) {
		ret = -1; return ret;
	}
	while (fgets(buf, 1024, fp)) {
		char *pbuf = buf, retbuf[16];
		int ret = 1;
		while (ret) {
			ret = csvread(&pbuf, retbuf, ' ', 1);
			if (strlen(retbuf) == 0) {
				break;
			}
			tmask[tmaskcnt] = atoi(retbuf);	tmaskcnt++;
			if (tmaskcnt >= MAX_TGT_CNT) {
				break;
			}
		}
		if (tmaskcnt >= MAX_TGT_CNT) {
			break;
		}
	}
	fclose(fp);

	fp = fopen(fname_tpt, "r");
	if (!fp) {
		ret = -1; return ret;
	}

	tgcnt = 0;
	tmaskcnt = 0;
	while (fgets(buf, 1024, fp)) {
		if (buf[0] == '#') {
			continue;
		}
		int ret = sscanf(buf, "%lf %lf %lf %lf %lf", &tgx[tgcnt], &tgy[tgcnt], &tgz[tgcnt], &ogx_cp[tgcnt], &ogy_cp[tgcnt]);
		if (ret < 5) {
			continue;
		}
		if (tmask[tmaskcnt++] == 0) {
			continue;
		}
		//printf("%d: %f, %f, %f, %f, %f\n", tgcnt, tgx[tgcnt], tgy[tgcnt], tgz[tgcnt], ogx_cp[tgcnt], ogy_cp[tgcnt]);
		tgcnt++;
	}

	fclose(fp);

	getTgt2D3D();

	return ret;
}

int papermodel::getTgt2D3D()
{
	int ret=0, avecnt=0;
	avetgap = 0.0;
	maxtgap = 0.0;

	for( int i=0; i<tgcnt; i++ ){
		bool flg = false;
		int minplidx = -1;
		double minpldist = pw + ph;
		for( int j=0; j<plcnt; j++){

			if (pl_cridx[j] != 0) {
				continue;
			}

			int left[4], right[4], lcnt=0, rcnt=0;
			memset( left, 0, sizeof(int)*4 );
			memset( right, 0, sizeof(int)*4 );
			for( int k=0; k<plvcnt[j]; k++ ){
				double vx0, vy0, vx1, vy1, lr;
				vx0 = plx_cp[j*4+(k+1)%plvcnt[j]] - plx_cp[j*4+k];
				vy0 = ply_cp[j*4+(k+1)%plvcnt[j]] - ply_cp[j*4+k];
				vx1 = ogx_cp[i] - plx_cp[j*4+k];
				vy1 = ogy_cp[i] - ply_cp[j*4+k];
				lr = vx0*vy1 - vx1*vy0;
				if( lr >= 0.0 ){ left[k] = 1; lcnt++; }
				if( lr <= 0.0 ){ right[k] = 1; rcnt++; }

				double pldist = line_point_dist_2d(0.0, 0.0, vx0, vy0, vx1, vy1);
				if ( minpldist > pldist ) {
					minpldist = pldist;
					minplidx = j;
				}
			}
			if( lcnt>=plvcnt[j] || rcnt>=plvcnt[j] ){
				double *pmat = &plmat[j*16];
				ogx[i] = pmat[0]*ogx_cp[i] + pmat[4]*ogy_cp[i] + pmat[12];
				ogy[i] = pmat[1]*ogx_cp[i] + pmat[5]*ogy_cp[i] + pmat[13];
				ogz[i] = pmat[2]*ogx_cp[i] + pmat[6]*ogy_cp[i] + pmat[14];

				double dx = tgx[i] - ogx[i];
				double dy = tgy[i] - ogy[i];
				double dz = tgz[i] - ogz[i];
				tgap[i] = sqrt( dx*dx + dy*dy + dz*dz );
				avetgap += tgap[i]; avecnt++;
				maxtgap = maxtgap > tgap[i] ? maxtgap : tgap[i];

				if (dx * pmat[8] + dy * pmat[9] + dz * pmat[10] < 0) {
					tgap[i] = -tgap[i];
				}

				//std::cout << "target: " << i << ", ply" << j << std::endl;
				flg = true;
				break; // j
			}
		}
		if( !flg ){
			if (minplidx > -1 ) {
				// plot on the nearest plane with max margin X mm
				double* pmat = &plmat[minplidx * 16];
				ogx[i] = pmat[0] * ogx_cp[i] + pmat[4] * ogy_cp[i] + pmat[12];
				ogy[i] = pmat[1] * ogx_cp[i] + pmat[5] * ogy_cp[i] + pmat[13];
				ogz[i] = pmat[2] * ogx_cp[i] + pmat[6] * ogy_cp[i] + pmat[14];

				double dx = tgx[i] - ogx[i];
				double dy = tgy[i] - ogy[i];
				double dz = tgz[i] - ogz[i];
				tgap[i] = sqrt(dx * dx + dy * dy + dz * dz);
				avetgap += tgap[i]; avecnt++;
				maxtgap = maxtgap > tgap[i] ? maxtgap : tgap[i];

				if (dx * pmat[8] + dy * pmat[9] + dz * pmat[10] < 0) {
					tgap[i] = -tgap[i];
				}

			} else {
				ogx[i] = ogy[i] = ogz[i] = 0.0;
			}
		}
		//std::cout << "tgap[" << i << "] = " << tgap[i] << std::endl;
	}
	if( avecnt>0 ){
		avetgap /= (double)avecnt;
	}

	return ret;
}

int papermodel::checkCollision()
{
	int ret=0;
	double mgn = 1.0; // 判定のマージン
	double minx[MAX_PL_FCNT], maxx[MAX_PL_FCNT], miny[MAX_PL_FCNT], maxy[MAX_PL_FCNT], minz[MAX_PL_FCNT], maxz[MAX_PL_FCNT];
	double pix[4], piy[4], piz[4], pi0x[4], pi0y[4], pjx[4], pjy[4], pjz[4], pj0x[4], pj0y[4];

	for( int i=0; i<plcnt; i++ ){
		minx[i] = DBL_MAX;		maxx[i] = -DBL_MAX;
		miny[i] = DBL_MAX;		maxy[i] = -DBL_MAX;
		minz[i] = DBL_MAX;		maxz[i] = -DBL_MAX;
		for( int j=0; j<plvcnt[i]; j++){
			minx[i] = minx[i] < plx[i*4+j] ? minx[i] : plx[i*4+j];
			maxx[i] = maxx[i] > plx[i*4+j] ? maxx[i] : plx[i*4+j];
			miny[i] = miny[i] < ply[i*4+j] ? miny[i] : ply[i*4+j];
			maxy[i] = maxy[i] > ply[i*4+j] ? maxy[i] : ply[i*4+j];
			minz[i] = minz[i] < plz[i*4+j] ? minz[i] : plz[i*4+j];
			maxz[i] = maxz[i] > plz[i*4+j] ? maxz[i] : plz[i*4+j];
		}
	}

	for( int i=0; i<plcnt; i++ ){
		for( int j=0; j<i; j++){
			// check collision
			if( maxx[i] < minx[j] ){
				continue;
			}
			if( minx[i] > maxx[j] ){
				continue;
			}
			if( maxy[i] < miny[j] ){
				continue;
			}
			if( miny[i] > maxy[j] ){
				continue;
			}
			if( maxz[i] < minz[j] ){
				continue;
			}
			if( minz[i] > maxz[j] ){
				continue;
			}
			bool left[4]={false}, right[4]={false};
			int lcnt=0, rcnt=0;
			for( int k=0; k<plvcnt[i]; k++){
				double dx = plx[i*4+k] - plx[j*4];
				double dy = ply[i*4+k] - ply[j*4];
				double dz = plz[i*4+k] - plz[j*4];
				double ip = dx*plmat[j*16+8] + dy*plmat[j*16+9] + dz*plmat[j*16+10];
				if( ip<=0 ){ left[k]=true; lcnt++; }
				if( ip>=0 ){ right[k]=true; rcnt++; }
			}
			if( lcnt >= plvcnt[i] || rcnt >= plvcnt[i] ){
				continue;
			}
			memset(left,0,sizeof(bool)*4); memset(right,0,sizeof(bool)*4); lcnt=0; rcnt=0;
			for( int k=0; k<plvcnt[j]; k++){
				double dx = plx[j*4+k] - plx[i*4];
				double dy = ply[j*4+k] - ply[i*4];
				double dz = plz[j*4+k] - plz[i*4];
				double ip = dx*plmat[i*16+8] + dy*plmat[i*16+9] + dz*plmat[i*16+10];
				if( ip<=0 ){ left[k]=true; lcnt++; }
				if( ip>=0 ){ right[k]=true; rcnt++; }
			}
			if( lcnt >= plvcnt[i] || rcnt >= plvcnt[i] ){
				continue;
			}

			//double pix[4], piy[4], piz[4], pi0x[4], pi0y[4], pjx[4], pjy[4], pjz[4], pj0x[4], pj0y[4];
			double *pmi = &plmat[i*16];
			double *pmj = &plmat[j*16];
 			for( int k=0; k<plvcnt[i]; k++){
				double dx = plx[i*4+k] - pmj[12];
				double dy = ply[i*4+k] - pmj[13];
				double dz = plz[i*4+k] - pmj[14];
				pix[k] = dx*pmj[0] + dy*pmj[1] + dz*pmj[2];
				piy[k] = dx*pmj[4] + dy*pmj[5] + dz*pmj[6];
				piz[k] = dx*pmj[8] + dy*pmj[9] + dz*pmj[10];
			}
			int iscnti=0;
			for( int k=0; k<plvcnt[i]; k++ ){
				int k1 = (k+1)%plvcnt[i];
				if( piz[k]>=0.0 && piz[k1]>=0.0 || piz[k]<=0.0 && piz[k1]<=0.0 ){
					continue;
				}
				double denom = fabs(piz[k]) + fabs(piz[k1]);
				double z0 = fabs(piz[k]) / denom;
				double z1 = fabs(piz[k1]) / denom;
				pi0x[iscnti] = z0 * pix[k1] + z1 * pix[k];
				pi0y[iscnti] = z0 * piy[k1] + z1 * piy[k];
				iscnti++;
			}
			for( int ii=0; ii<iscnti; ii++ ){
				int left[4], right[4], lcnt=0, rcnt=0;
				memset( left, 0, sizeof(int)*4 );
				memset( right, 0, sizeof(int)*4 );
				for( int k=0; k<plvcnt[j]; k++ ){
					int k0 = j*4+k;
					int k1 = j*4+(k+1)%plvcnt[j];
					double vx0, vy0, vx1, vy1, lr;
					vx0 = plx_cp[k1] - plx_cp[k0];
					vy0 = ply_cp[k1] - ply_cp[k0];
					vx1 = pi0x[ii] - plx_cp[k0];
					vy1 = pi0y[ii] - ply_cp[k0];
					lr = vx0*vy1 - vx1*vy0;
					if( lr > 0.0 ){
						lr /= sqrt(vx0*vx0+vy0*vy0);
						if( lr > mgn ){
							left[k] = 1; lcnt++;
						}
					}
					if( lr < 0.0 ){
						lr /= sqrt(vx0*vx0+vy0*vy0);
						if( lr < -mgn ){
							right[k] = 1; rcnt++; 
						}
					}
				}
				if( lcnt>=plvcnt[j] || rcnt>=plvcnt[j] ){
					printf( "collision between polygon[%d] and polygon[%d]\n", i, j );
					ret = 1;
					//goto end;
					break;
				}
			}

 			for( int k=0; k<plvcnt[j]; k++){
				double dx = plx[j*4+k] - pmi[12];
				double dy = ply[j*4+k] - pmi[13];
				double dz = plz[j*4+k] - pmi[14];
				pjx[k] = dx*pmi[0] + dy*pmi[1] + dz*pmi[2];
				pjy[k] = dx*pmi[4] + dy*pmi[5] + dz*pmi[6];
				pjz[k] = dx*pmi[8] + dy*pmi[9] + dz*pmi[10];
			}
			int iscntj=0;
 			for( int k=0; k<plvcnt[j]; k++){
				int k1 = (k+1)%plvcnt[j];
				if( pjz[k]>=0.0 && pjz[k1]>=0.0 || pjz[k]<=0.0 && pjz[k1]<=0.0 ){
					continue;
				}
				double denom = fabs(pjz[k]) + fabs(pjz[k1]);
				double z0 = fabs(pjz[k]) / denom;
				double z1 = fabs(pjz[k1]) / denom;
				pj0x[iscntj] = z0 * pjx[k1] + z1 * pjx[k];
				pj0y[iscntj] = z0 * pjy[k1] + z1 * pjy[k];
				iscntj++;
			}
			for( int jj=0; jj<iscntj; jj++ ){
				//printf("jj=%d\n",jj);
				int left[4], right[4], lcnt=0, rcnt=0;
				memset( left, 0, sizeof(int)*4 );
				memset( right, 0, sizeof(int)*4 );
				for( int k=0; k<plvcnt[i]; k++ ){
					int k0 = i*4+k;
					int k1 = i*4+(k+1)%plvcnt[i];
					double vx0, vy0, vx1, vy1, lr;
					vx0 = plx_cp[k1] - plx_cp[k0];
					vy0 = ply_cp[k1] - ply_cp[k0];
					vx1 = pj0x[jj] - plx_cp[k0];
					vy1 = pj0y[jj] - ply_cp[k0];
					lr = vx0*vy1 - vx1*vy0;
					if( lr > 0.0 ){
						lr /= sqrt(vx0*vx0+vy0*vy0);
						if( lr > mgn ){
							left[k] = 1; lcnt++;
						}
					}
					if( lr < 0.0 ){
						lr /= sqrt(vx0*vx0+vy0*vy0);
						if( lr < -mgn ){
							right[k] = 1; rcnt++; 
						}
					}
				}
				if( lcnt>=plvcnt[i] || rcnt>=plvcnt[i] ){
					printf( "collision between polygon[%d] and polygon[%d]\n", j, i );
#if 0
					{
						FILE *fp=fopen("output/debug.obj","w");
#if 1
						fprintf(fp, "v %f %f %f\n", plx_cp[i*4+0], ply_cp[i*4+0], 0.0);
						fprintf(fp, "v %f %f %f\n", plx_cp[i*4+1], ply_cp[i*4+1], 0.0);
						fprintf(fp, "v %f %f %f\n", plx_cp[i*4+2], ply_cp[i*4+2], 0.0);
						fprintf(fp, "v %f %f %f\n", plx_cp[i*4+3], ply_cp[i*4+3], 0.0);
						fprintf(fp, "v %f %f %f\n", pjx[0], pjy[0], pjz[0]);
						fprintf(fp, "v %f %f %f\n", pjx[1], pjy[1], pjz[1]);
						fprintf(fp, "v %f %f %f\n", pjx[2], pjy[2], pjz[2]);
						fprintf(fp, "v %f %f %f\n", pjx[3], pjy[3], pjz[3]);
						fprintf(fp, "v %f %f %f\n", pj0x[jj], pj0y[jj], 0.0);
#else
						fprintf(fp, "v %f %f %f %f\n", plx[i*4+0], ply[i*4+0], plz[i*4+0]);
						fprintf(fp, "v %f %f %f %f\n", plx[i*4+1], ply[i*4+1], plz[i*4+1]);
						fprintf(fp, "v %f %f %f %f\n", plx[i*4+2], ply[i*4+2], plz[i*4+2]);
						fprintf(fp, "v %f %f %f %f\n", plx[i*4+3], ply[i*4+3], plz[i*4+3]);
						fprintf(fp, "v %f %f %f %f\n", plx[j*4+0], ply[j*4+0], plz[j*4+0]);
						fprintf(fp, "v %f %f %f %f\n", plx[j*4+1], ply[j*4+1], plz[j*4+1]);
						fprintf(fp, "v %f %f %f %f\n", plx[j*4+2], ply[j*4+2], plz[j*4+2]);
						fprintf(fp, "v %f %f %f %f\n", plx[j*4+3], ply[j*4+3], plz[j*4+3]);
#endif
						fprintf(fp, "f 1// 2// 3// 4//\n");
						fprintf(fp, "f 5// 6// 7// 8//\n");
						//fprintf(fp, "f 9// 10// 11//\n");
						fclose(fp);
					}
#endif
					ret = 1;
					//goto end;
					break;
				}
			}
		}
	}

end:
	return ret;
}

int papermodel::optMat( int mode )
{
	int ret=0;

	if (tgcnt == 0) {
		std::cout << "no target point." << std::endl;
		return -1;
	}

	double tgx0[MAX_TGT_CNT], tgy0[MAX_TGT_CNT], tgz0[MAX_TGT_CNT];
	double ogx0[MAX_TGT_CNT], ogy0[MAX_TGT_CNT], ogz0[MAX_TGT_CNT];
	memcpy( tgx0, tgx, sizeof(double)*tgcnt );
	memcpy( tgy0, tgy, sizeof(double)*tgcnt );
	memcpy( tgz0, tgz, sizeof(double)*tgcnt );
	memcpy( ogx0, ogx, sizeof(double)*tgcnt );
	memcpy( ogy0, ogy, sizeof(double)*tgcnt );
	memcpy( ogz0, ogz, sizeof(double)*tgcnt );

	double mat[16];
	getMat( tgcnt, ogx0,ogy0,ogz0, tgx0,tgy0,tgz0, mat);
	transpose_m44( mat );

	// dst=1 : n44[16] = n44[16] * m44[16]
	mult_m44_n44( crs[0].m3, mat, 1);

	switch( mode ){
		case 0/*CMODE_A*/:	crs[0].calcXA_CP( 1/*flg_interpolate*/, &rp );	break;
		case 1/*CMODE_B*/:	crs[0].calcCPA_X( 1/*flg_interpolate*/, &rp );	break;
		case 2/*CMODE_C*/:	crs[0].calcCPX_A( 1/*flg_interpolate*/, &rp );	break;
		case 3/*CMODE_R*/:	crs[0].calcR_TA( 0/*flg_interpolate*/, &rp, -1, -1, 2*M_PI, 0 );	break;
	}
	set_postproc_type( PPTYPE_PRICURVE );
	postproc();
	set_postproc_type( PPTYPE_UNDEF );

	return ret;
}

int papermodel::optMatRot(int mode)
{
	int ret = 0;

	if (tgcnt == 0) {
		std::cout << "no target point." << std::endl;
		return -1;
	}

	double tgx0[MAX_TGT_CNT], tgy0[MAX_TGT_CNT], tgz0[MAX_TGT_CNT];
	double ogx0[MAX_TGT_CNT], ogy0[MAX_TGT_CNT], ogz0[MAX_TGT_CNT];
	memcpy(tgx0, tgx, sizeof(double) * tgcnt);
	memcpy(tgy0, tgy, sizeof(double) * tgcnt);
	memcpy(tgz0, tgz, sizeof(double) * tgcnt);
	memcpy(ogx0, ogx, sizeof(double) * tgcnt);
	memcpy(ogy0, ogy, sizeof(double) * tgcnt);
	memcpy(ogz0, ogz, sizeof(double) * tgcnt);
	for (int i = 0; i < tgcnt; i++ ) {
		tgx0[i] -= crs[0].m3[12];
		tgy0[i] -= crs[0].m3[13];
		tgz0[i] -= crs[0].m3[14];
		ogx0[i] -= crs[0].m3[12];
		ogy0[i] -= crs[0].m3[13];
		ogz0[i] -= crs[0].m3[14];
	}

	double mat[16], tmat0[16], tmat1[16];
	getMatRot(tgcnt, ogx0, ogy0, ogz0, tgx0, tgy0, tgz0, mat);
	transpose_m44(mat);
	unit_m44(tmat0); tmat0[12] = -crs[0].m3[12]; tmat0[13] = -crs[0].m3[13]; tmat0[14] = -crs[0].m3[14];
	unit_m44(tmat1); tmat1[12] =  crs[0].m3[12]; tmat1[13] =  crs[0].m3[13]; tmat1[14] =  crs[0].m3[14];
	//memcpy(tmat1, tmat0, sizeof(double) * 16); inv_m44(tmat1);
	// dst=1 : n44[16] = n44[16] * m44[16]
	mult_m44_n44(crs[0].m3, tmat0, 1);
	mult_m44_n44(crs[0].m3, mat, 1);
	mult_m44_n44(crs[0].m3, tmat1, 1);

	switch (mode) {
	case 0/*CMODE_A*/:	crs[0].calcXA_CP(1/*flg_interpolate*/, &rp);	break;
	case 1/*CMODE_B*/:	crs[0].calcCPA_X(1/*flg_interpolate*/, &rp);	break;
	case 2/*CMODE_C*/:	crs[0].calcCPX_A(1/*flg_interpolate*/, &rp);	break;
	case 3/*CMODE_R*/:	crs[0].calcR_TA(0/*flg_interpolate*/, &rp, -1, -1, 2 * M_PI, 0);	break;
	}
	set_postproc_type(PPTYPE_PRICURVE);
	postproc();
	set_postproc_type(PPTYPE_UNDEF);

	return ret;
}

int papermodel::calcAvetgap()
{
	int ret = 0;

	if (tgcnt == 0) {
		std::cout << "no target point." << std::endl;
		return -1;
	}

	double dx, dy, dz, diff, diff0 = 0.0, diff1 = 0.0, maxdiff = 0.0;
	for (int i = 0; i < tgcnt; i++ ) {
		dx = ogx[i] - tgx[i];
		dy = ogy[i] - tgy[i];
		dz = ogz[i] - tgz[i];
		diff = sqrt(dx * dx + dy * dy + dz * dz);
		diff0 += diff;
		maxdiff = maxdiff > diff ? maxdiff : diff;
	}
	//std::cout << "diff0=" << diff0/(double)tgcnt << ", diff1=" << diff1/(double)tgcnt << std::endl;
	this->avetgap = diff0 / (double)tgcnt;
	this->maxtgap = maxdiff;

	return ret;
}

int papermodel::calcAvetgapMat()
{
	int ret = 0;

	if (tgcnt == 0) {
		std::cout << "no target point." << std::endl;
		return -1;
	}

	double tgx0[MAX_TGT_CNT], tgy0[MAX_TGT_CNT], tgz0[MAX_TGT_CNT];
	double ogx0[MAX_TGT_CNT], ogy0[MAX_TGT_CNT], ogz0[MAX_TGT_CNT];
	memcpy(tgx0, tgx, sizeof(double) * tgcnt);
	memcpy(tgy0, tgy, sizeof(double) * tgcnt);
	memcpy(tgz0, tgz, sizeof(double) * tgcnt);
	memcpy(ogx0, ogx, sizeof(double) * tgcnt);
	memcpy(ogy0, ogy, sizeof(double) * tgcnt);
	memcpy(ogz0, ogz, sizeof(double) * tgcnt);

	double mat[16];
	getMat(tgcnt, ogx0, ogy0, ogz0, tgx0, tgy0, tgz0, mat);

	double dx, dy, dz, diff, diff0 = 0.0, diff1 = 0.0, maxdiff = 0.0;
#if 1
	for (int i = 0; i < tgcnt; i++ ) {
		ogx0[i] = mat[0] * ogx[i] + mat[1] * ogy[i] + mat[2] * ogz[i] + mat[3];
		ogy0[i] = mat[4] * ogx[i] + mat[5] * ogy[i] + mat[6] * ogz[i] + mat[7];
		ogz0[i] = mat[8] * ogx[i] + mat[9] * ogy[i] + mat[10] * ogz[i] + mat[11];
		dx = ogx0[i] - tgx[i];
		dy = ogy0[i] - tgy[i];
		dz = ogz0[i] - tgz[i];
		diff = sqrt(dx * dx + dy * dy + dz * dz);
		diff0 += diff;
		maxdiff = maxdiff > diff ? maxdiff : diff;
	}
#else
	for (int i = 0; i < tgcnt; i++) {
		tgx0[i] = mat[0] * tgx[i] + mat[1] * tgy[i] + mat[2] * tgz[i] + mat[3];
		tgy0[i] = mat[4] * tgx[i] + mat[5] * tgy[i] + mat[6] * tgz[i] + mat[7];
		tgz0[i] = mat[8] * tgx[i] + mat[9] * tgy[i] + mat[10] * tgz[i] + mat[11];
		dx = tgx0[i] - ogx[i];
		dy = tgy0[i] - ogy[i];
		dz = tgz0[i] - ogz[i];
		diff = sqrt(dx * dx + dy * dy + dz * dz);
		diff1 += diff;
		maxdiff = maxdiff > diff ? maxdiff : diff;
	}
#endif
	//std::cout << "diff0=" << diff0/(double)tgcnt << ", diff1=" << diff1/(double)tgcnt << std::endl;
	this->avetgap = diff0 / (double)tgcnt;
	this->maxtgap = maxdiff;

	return ret;
}

int papermodel::calcAvetgapRot()
{
	int ret = 0;

	if (tgcnt == 0) {
		std::cout << "no target point." << std::endl;
		return -1;
	}

	double tgx0[MAX_TGT_CNT], tgy0[MAX_TGT_CNT], tgz0[MAX_TGT_CNT];
	double ogx0[MAX_TGT_CNT], ogy0[MAX_TGT_CNT], ogz0[MAX_TGT_CNT];
	memcpy(tgx0, tgx, sizeof(double) * tgcnt);
	memcpy(tgy0, tgy, sizeof(double) * tgcnt);
	memcpy(tgz0, tgz, sizeof(double) * tgcnt);
	memcpy(ogx0, ogx, sizeof(double) * tgcnt);
	memcpy(ogy0, ogy, sizeof(double) * tgcnt);
	memcpy(ogz0, ogz, sizeof(double) * tgcnt);
	for (int i = 0; i < tgcnt; i++) {
		tgx0[i] -= crs[0].m3[12];
		tgy0[i] -= crs[0].m3[13];
		tgz0[i] -= crs[0].m3[14];
		ogx0[i] -= crs[0].m3[12];
		ogy0[i] -= crs[0].m3[13];
		ogz0[i] -= crs[0].m3[14];
	}

	double mat[16];
	getMatRot(tgcnt, ogx0, ogy0, ogz0, tgx0, tgy0, tgz0, mat);

	double dx, dy, dz, diff, diff0 = 0.0, diff1 = 0.0, maxdiff = 0.0, tmpx, tmpy, tmpz;
#if 1
	for (int i = 0; i < tgcnt; i++) {
		tmpx = mat[0] * ogx0[i] + mat[1] * ogy0[i] + mat[2] * ogz0[i];
		tmpy = mat[4] * ogx0[i] + mat[5] * ogy0[i] + mat[6] * ogz0[i];
		tmpz = mat[8] * ogx0[i] + mat[9] * ogy0[i] + mat[10] * ogz0[i];
		dx = tmpx - tgx0[i];
		dy = tmpy - tgy0[i];
		dz = tmpz - tgz0[i];
		diff = sqrt(dx * dx + dy * dy + dz * dz);
		diff0 += diff;
		maxdiff = maxdiff > diff ? maxdiff : diff;
	}
#else
	for (int i = 0; i < tgcnt; i++) {
		tgx0[i] = mat[0] * tgx[i] + mat[1] * tgy[i] + mat[2] * tgz[i] + mat[3];
		tgy0[i] = mat[4] * tgx[i] + mat[5] * tgy[i] + mat[6] * tgz[i] + mat[7];
		tgz0[i] = mat[8] * tgx[i] + mat[9] * tgy[i] + mat[10] * tgz[i] + mat[11];
		dx = tgx0[i] - ogx[i];
		dy = tgy0[i] - ogy[i];
		dz = tgz0[i] - ogz[i];
		diff = sqrt(dx * dx + dy * dy + dz * dz);
		diff1 += diff;
		maxdiff = maxdiff > diff ? maxdiff : diff;
	}
#endif
	//std::cout << "diff0=" << diff0/(double)tgcnt << ", diff1=" << diff1/(double)tgcnt << std::endl;
	this->avetgap = diff0 / (double)tgcnt;
	this->maxtgap = maxdiff;

	return ret;
}

int onStrip(double Xx0, double Xy0, double Xx1, double Xy1,
	double rx0, double ry0, double rx1, double ry1,
	double px, double py )
{
#if 0
	X1--->rx1, ry1;
	| px, py;
	X0--->rx0, ry0;
#endif
	double px0 = px - Xx0, py0 = py - Xy0, px1 = px - Xx1, py1 = py - Xy1;
	double tx = Xx1 - Xx0, ty = Xy1 - Xy0;
	double lrx = tx * py0 - py0 * ty; // >0:right
	double lr0 = rx0 * py0 - px0 * ry0; // <0:left
	double lr1 = rx1 * py1 - px1 * ry1; // >0:right

	if ( lrx>=0 && lr0 <= 0 && lr1 >=0 ) {
		return 1;
	}
	return 0;
}

int papermodel::optFold()
{
	int ret=0;
	crease* c = &(this->crs[0]);

	double wgt = 0.9;
	int polydegree = 3;

	int stcnt = 0;
	double sttgx_[MAX_TGT_CNT], sttgy_[MAX_TGT_CNT], sttgz_[MAX_TGT_CNT];
	double stogx_[MAX_TGT_CNT], stogy_[MAX_TGT_CNT], stogz_[MAX_TGT_CNT];

	int stlcnt[MAX_SPCNT], strcnt[MAX_SPCNT];
	double *stltgx[MAX_SPCNT], *stltgy[MAX_SPCNT], *stltgz[MAX_SPCNT];
	double *stlogx[MAX_SPCNT], *stlogy[MAX_SPCNT], *stlogz[MAX_SPCNT];
	double *strtgx[MAX_SPCNT], *strtgy[MAX_SPCNT], *strtgz[MAX_SPCNT];
	double *strogx[MAX_SPCNT], *strogy[MAX_SPCNT], *strogz[MAX_SPCNT];
	memset(stlcnt, 0, sizeof(int) * MAX_SPCNT);
	memset(strcnt, 0, sizeof(int) * MAX_SPCNT);

	int pcnt=0;
	double p_x[MAX_SPCNT], p_y[MAX_SPCNT];

	// for file out
	std::ofstream fout("./output/optFold.csv");
	double tgtAlphaL0_[MAX_TGT_CNT];	memset(tgtAlphaL0_, 0, sizeof(double) * MAX_TGT_CNT);
	double tgtAlphaR0_[MAX_TGT_CNT];	memset(tgtAlphaR0_, 0, sizeof(double) * MAX_TGT_CNT);
	double *tgtAlphaL0[MAX_SPCNT], *tgtAlphaR0[MAX_SPCNT];
	double orgAlpha[MAX_SPCNT];	memset(orgAlpha, 0, sizeof(double)*MAX_SPCNT);
	double tgtAlpha[MAX_SPCNT];	memset(tgtAlpha, 0, sizeof(double) * MAX_SPCNT);
	double revAlpha[MAX_SPCNT];	memset(revAlpha, 0, sizeof(double) * MAX_SPCNT);
	double orgAlpha2[MAX_SPCNT];	memset(orgAlpha2, 0, sizeof(double) * MAX_SPCNT);
	double revAlpha2[MAX_SPCNT];	memset(revAlpha2, 0, sizeof(double) * MAX_SPCNT);

	for (int i = c->Xsidx; i < c->Xeidx; i++ ) {
		stltgx[i] = &(sttgx_[stcnt]);
		stltgy[i] = &(sttgy_[stcnt]);
		stltgz[i] = &(sttgz_[stcnt]);
		stlogx[i] = &(stogx_[stcnt]);
		stlogy[i] = &(stogy_[stcnt]);
		stlogz[i] = &(stogz_[stcnt]);
		tgtAlphaL0[i] = &(tgtAlphaL0_[stcnt]);
		for (int j = 0; j < tgcnt; j++) {
			if (onStrip(c->Xx2d[i], c->Xy2d[i], c->Xx2d[i + 1], c->Xy2d[i + 1],
				c->rlx_cp[i], c->rly_cp[i], c->rlx_cp[i + 1], c->rly_cp[i + 1], ogx_cp[j], ogy_cp[j]))
			{
				double tmpx, tmpy, tmpz;
				tmpx = tgx[j] - c->Xx[i];
				tmpy = tgy[j] - c->Xy[i];
				tmpz = tgz[j] - c->Xz[i];
				sttgx_[stcnt] = tmpx * c->Tx[i] + tmpy * c->Ty[i] + tmpz * c->Tz[i];
				sttgy_[stcnt] = tmpx * c->Nx[i] + tmpy * c->Ny[i] + tmpz * c->Nz[i];
				sttgz_[stcnt] = tmpx * c->Bx[i] + tmpy * c->By[i] + tmpz * c->Bz[i];
				tmpx = ogx[j] - c->Xx[i];
				tmpy = ogy[j] - c->Xy[i];
				tmpz = ogz[j] - c->Xz[i];
				stogx_[stcnt] = tmpx * c->Tx[i] + tmpy * c->Ty[i] + tmpz * c->Tz[i];
				stogy_[stcnt] = tmpx * c->Nx[i] + tmpy * c->Ny[i] + tmpz * c->Nz[i];
				stogz_[stcnt] = tmpx * c->Bx[i] + tmpy * c->By[i] + tmpz * c->Bz[i];
				stcnt++;
				stlcnt[i]++;
#if 0
				std::cout << i << "_left: (" << ogx_cp[j] << "," << ogy_cp[j] << ")" << std::endl;
#endif
			}
		}
		strtgx[i] = &(sttgx_[stcnt]);
		strtgy[i] = &(sttgy_[stcnt]);
		strtgz[i] = &(sttgz_[stcnt]);
		strogx[i] = &(stogx_[stcnt]);
		strogy[i] = &(stogy_[stcnt]);
		strogz[i] = &(stogz_[stcnt]);
		tgtAlphaR0[i] = &(tgtAlphaR0_[stcnt]);
		for (int j = 0; j < tgcnt; j++) {
			if (onStrip(c->Xx2d[i+1], c->Xy2d[i+1], c->Xx2d[i], c->Xy2d[i],
				c->rrx_cp[i+1], c->rry_cp[i+1], c->rrx_cp[i], c->rry_cp[i], ogx_cp[j], ogy_cp[j]))
			{
				double tmpx, tmpy, tmpz;
				tmpx = tgx[j] - c->Xx[i];
				tmpy = tgy[j] - c->Xy[i];
				tmpz = tgz[j] - c->Xz[i];
				sttgx_[stcnt] = tmpx * c->Tx[i] + tmpy * c->Ty[i] + tmpz * c->Tz[i];
				sttgy_[stcnt] = tmpx * c->Nx[i] + tmpy * c->Ny[i] + tmpz * c->Nz[i];
				sttgz_[stcnt] = tmpx * c->Bx[i] + tmpy * c->By[i] + tmpz * c->Bz[i];
				tmpx = ogx[j] - c->Xx[i];
				tmpy = ogy[j] - c->Xy[i];
				tmpz = ogz[j] - c->Xz[i];
				stogx_[stcnt] = tmpx * c->Tx[i] + tmpy * c->Ty[i] + tmpz * c->Tz[i];
				stogy_[stcnt] = tmpx * c->Nx[i] + tmpy * c->Ny[i] + tmpz * c->Nz[i];
				stogz_[stcnt] = tmpx * c->Bx[i] + tmpy * c->By[i] + tmpz * c->Bz[i];
				stcnt++;
				strcnt[i]++;
#if 0
				std::cout << i << "_right: (" << ogx_cp[j] << "," << ogy_cp[j] << ")" << std::endl;
#endif
			}
		}
	}

#if 0
	for (int i = c->Xsidx; i < c->Xeidx; i++) {
		for (int j = 0; j < stlcnt[i]; j++) {
			std::cout << i << "-L: (" << stltgx[i][j] << "," << stltgy[i][j] << "," << stltgz[i][j] << "), ("
				<< stlogx[i][j] << "," << stlogy[i][j] << "," << stlogz[i][j] << ")" << std::endl;
		}
	}
	for (int i = c->Xsidx; i < c->Xeidx; i++) {
		for (int j = 0; j < strcnt[i]; j++) {
			std::cout << i << "-R: (" << strtgx[i][j] << "," << strtgy[i][j] << "," << strtgz[i][j] << "), ("
				<< stlogx[i][j] << "," << strogy[i][j] << "," << strogz[i][j] << ")" << std::endl;
		}
	}
#endif

	// calculate new folding angle for each vertex
	pcnt = 0;
	for (int i = c->Xsidx; i < c->Xeidx; i++) {
		int acnt = 0;
		double tgt_alpha = 0.0, rev_alpha=0.0;

		for (int j = 0; j < stlcnt[i]; j++) {
			double alpha = atan2(stltgz[i][j], stltgy[i][j]);
			//std::cout << "  L(" << i << "," << j << ")=" << alpha * 180.0 / M_PI << std::endl;
			tgtAlphaL0[i][j] = alpha * 180.0/M_PI;
			tgt_alpha += alpha;
			acnt++;
		}
		for (int j = 0; j < strcnt[i]; j++) {
			double alpha = atan2(strtgz[i][j], -strtgy[i][j]);
			//std::cout << "  R(" << i << "," << j << ")=" << alpha * 180.0 / M_PI << std::endl;
			tgtAlphaR0[i][j] = alpha * 180.0/M_PI;
			tgt_alpha += alpha;
			acnt++;
		}
		if (acnt == 0) {
			continue;
		}
		tgt_alpha /= (double)acnt;
		rev_alpha = wgt * tgt_alpha + (1 - wgt) * c->alpha[i];
#if 0
		std::cout << i << ": org alpha = " << c->alpha[i] * 180.0 / M_PI
			<< ", target alpha = " << tgt_alpha * 180.0 / M_PI
			<< ", revised alpha = " << rev_alpha * 180.0 / M_PI << std::endl;
#endif
		orgAlpha[i] = c->alpha[i] * 180.0 / M_PI;
		tgtAlpha[i] = tgt_alpha * 180.0 / M_PI;
		revAlpha[i] = rev_alpha * 180.0 / M_PI;

		p_x[pcnt] = i;
		p_y[pcnt] = rev_alpha;
		pcnt++;
	}

	// interpolate & update

	double * coeff = least(polydegree, pcnt, p_x, p_y);

	if (coeff != NULL) {
		for (int i = 0; i < polydegree + 1; i++) {
			std::cout << "coeffi of degree " << polydegree-i << " = " << coeff[i] << std::endl;
		}
		for (int i = c->Xsidx; i <= c->Xeidx; i++) {
			double x = (double)i;
			double xx = x*x;
			double xxx = xx*x;
			double alpha = coeff[0]*xxx + coeff[1]*xx + coeff[2]*x + coeff[3];
#if 0
			std::cout << i << ": org alpha = " << c->alpha[i] * 180.0 / M_PI
				<< ", revised alpha = " << alpha * 180.0 / M_PI << std::endl;
#endif
			orgAlpha2[i] = c->alpha[i] * 180.0 / M_PI;
			revAlpha2[i] = alpha * 180.0 / M_PI;

			c->alpha[i] = alpha;
		}
		delete[] coeff;
	} else {
		ret = -1;
	}
	c->setP_a(-1);
	c->calcCPA_X(1/*flg_interpolate*/, &rp);

	set_postproc_type(PPTYPE_PRICURVE);
	postproc();
	set_postproc_type(PPTYPE_UNDEF);

	fout << "i,orgAlpha,tgtAlpha,revAlpha,orgAlpha2,revAlpha2,c->alpha,,tgtLeft,,tgtRight" << std::endl;
	for (int i = c->Xsidx; i <= c->Xeidx; i++) {
		fout << i << "," << orgAlpha[i]
			<< "," << tgtAlpha[i]
			<< "," << revAlpha[i]
			<< "," << orgAlpha2[i]
			<< "," << revAlpha2[i]
			<< "," << c->alpha[i]*180.0/M_PI << std::endl;
	}
	for (int i = c->Xsidx; i <= c->Xeidx; i++) {
		for (int j = 0; j < stlcnt[i]; j++) {
			fout << ",,,,,,," << i << "," << tgtAlphaL0[i][j] << std::endl;
		}
		for (int j = 0; j < strcnt[i]; j++) {
			fout << ",,,,,,,,," << i << "," << tgtAlphaR0[i][j] << std::endl;

		}
	}
	return ret;
}

int papermodel::optTorsion()
{
	int ret = 0;
	crease* c = &(this->crs[0]);

	return ret;
}

int papermodel::optRulings()
{
	int ret = 0;

	return ret;
}

int papermodel::optCP()
{
	int ret = 0;

	return ret;
}
