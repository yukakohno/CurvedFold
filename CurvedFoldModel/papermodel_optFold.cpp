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

	tgcnt = 0;
	for( int j=20; j<ph; j+=20 ){
		for( int i=20; i<pw; i+=20 ){
			ogx_cp[tgcnt] = i;
			ogy_cp[tgcnt] = j;
			tgcnt++;
		}
	}

	getTgt2D3D();
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
		printf( "%d: %f, %f, %f, %f, %f\n", tgcnt, tgx[tgcnt], tgy[tgcnt], tgz[tgcnt], ogx_cp[tgcnt], ogy_cp[tgcnt] );
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
		printf("%d: %f, %f, %f, %f, %f\n", tgcnt, tgx[tgcnt], tgy[tgcnt], tgz[tgcnt], ogx_cp[tgcnt], ogy_cp[tgcnt]);
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

	for( int i=0; i<tgcnt; i++ ){
		bool sign = true;
		bool flg = false;
		for( int j=0; j<plcnt; j++){
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
			}
			if( lcnt>=plvcnt[j] || rcnt>=plvcnt[j] ){
				double *pmat = &plmat[j*16];
				ogx[i] = pmat[0]*ogx_cp[i] + pmat[4]*ogy_cp[i] + pmat[12];
				ogy[i] = pmat[1]*ogx_cp[i] + pmat[5]*ogy_cp[i] + pmat[13];
				ogz[i] = pmat[2]*ogx_cp[i] + pmat[6]*ogy_cp[i] + pmat[14];
				flg = true;

				double dx = tgx[i] - ogx[i];
				double dy = tgy[i] - ogy[i];
				double dz = tgz[i] - ogz[i];
				if( dx*pmat[8] + dy*pmat[9] + dz*pmat[10] <0 ){
					sign = false;
				}

				tgap[i] = sqrt( dx*dx + dy*dy + dz*dz );
				if( sign==false ){
					tgap[i] = -tgap[i];
				}
				avetgap += fabs(tgap[i]); avecnt++;

				break; // j
			}
		}
		if( !flg ){
			//printf("ogx[%d] = ogy[%d] = ogz[%d] = 0.0\n", i,i,i );
			ogx[i] = ogy[i] = ogz[i] = 0.0;
		} else {
//			double dx = tgx[i] - ogx[i];
//			double dy = tgy[i] - ogy[i];
//			double dz = tgz[i] - ogz[i];
//			tgap[i] = sqrt( dx*dx + dy*dy + dz*dz );
//			if( sign==false ){
//				tgap[i] = -tgap[i];
//			}
//			avetgap += fabs(tgap[i]); avecnt++;
		}
	}
	if( avecnt>0 ){
		avetgap /= (double)avecnt;
	}

	return ret;
}

double papermodel::calcRulCross( double *xx, double *xy, double *rx, double *ry, double *rlen, int cvcnt ) // Xx2d, Xy2d, rlx_cp, rly_cp, rllen, cnt
{
	double ex[MAX_SPCNT],ey[MAX_SPCNT], diff=0.0;
	double constlen=50.0;

	//for( int i=2;  i<cvcnt-2; i++ ){
	for( int i=0;  i<cvcnt; i++ ){
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
	//for( int i=4;  i<cvcnt-4; i++ ){
	for( int i=0;  i<cvcnt-1; i++ ){
		double ix,iy, l0,l1, l0_,l1_, a,b,c,d, area;
		int ret = intersectionOfLine( xx[i],xy[i],ex[i],ey[i], xx[i+1],xy[i+1],ex[i+1],ey[i+1], &ix, &iy, &l0, &l1 );
		if( ret<0 ){
			continue;
		}
		// サラスの公式 O(0,0),A(a,b),B(c,d) -> area of OAB = abs(a*d-b*c)*0.5
		l0_ = rlen[i]  -l0;
		l1_ = rlen[i+1]-l1;
		//l0_ = constlen-l0;
		//l1_ = constlen-l1;
		a = l0_*rx[i];
		b = l0_*ry[i];
		c = l1_*rx[i+1];
		d = l1_*ry[i+1];
		area = sqrt(fabs(a*d-b*c))+0.01;
		diff += area;
	}
end:
	//return diff*0.1 / (double)cvcnt;
	return diff / (double)cvcnt;
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

int papermodel::checkRulAngle()
{
	int ret=0;
	crease *c = &crs[0];

	for( int i=c->Xsidx; i<=c->Xeidx; i++ ){
		double dx0=0, dy0=0, dx1=0, dy1=0;
		if( c->Xsidx<i ){
			dx0 = c->Xx2d[i] - c->Xx2d[i-1];
			dy0 = c->Xy2d[i] - c->Xy2d[i-1];
		}
		if( i<c->Xeidx ){
			dx1 = c->Xx2d[i+1] - c->Xx2d[i];
			dy1 = c->Xy2d[i+1] - c->Xy2d[i];
		}
		double ip = c->rlx_cp[i] * c->Tx2d[i] + c->rly_cp[i] * c->Ty2d[i];
		double rl=0;
		if( ip>0.0 && i<c->Xeidx ){
			rl = c->rlx_cp[i] * dy1 - c->rly_cp[i] * dx1;
		} else if( ip<0.0 && c->Xsidx<i ){
			rl = c->rlx_cp[i] * dy0 - c->rly_cp[i] * dx0;
		} else {
			rl = c->rlx_cp[i] * c->Ty2d[i] - c->rly_cp[i] * c->Tx2d[i];
		}
		if( rl >= 0.0 ){
			ret = -1;
			break;
		}

		ip = c->rrx_cp[i] * c->Tx2d[i] + c->rry_cp[i] * c->Ty2d[i];
		if( ip>0.0 && i<c->Xeidx ){
			rl = c->rrx_cp[i] * dy1 - c->rry_cp[i] * dx1;
		} else if( ip<0.0 && c->Xsidx<i ){
			rl = c->rrx_cp[i] * dy0 - c->rry_cp[i] * dx0;
		} else {
			rl = c->rrx_cp[i] * c->Ty2d[i-1] - c->rry_cp[i] * c->Tx2d[i-1];
		}
		if( rl < 0.0 ){
			ret = 1;
			break;
		}
	}

	return ret;
}

int papermodel::checkRulCross()
{
	double area;
	int result = checkRulCross( area );
	//printf("ruling crossing = %f\n", area);
	return result;
}
int papermodel::checkRulCross( double &area )
{
	int ret=0;
	crease *c = &crs[0];
	//int ii = c->Xsidx;
	//area = calcRulCross( &(c->Xx2d[ii]), &(c->Xy2d[ii]), &(c->rlx_cp[ii]), &(c->rly_cp[ii]), &(c->rllen[ii]), c->Xeidx-ii+1 );
	//area = calcRulCross( c->Xx2d, c->Xy2d, c->rlx_cp, c->rly_cp, c->rllen, c->Xcnt );
	double area_rl = calcRulCross( &(c->Xx2d[c->Xsidx]), &(c->Xy2d[c->Xsidx]),
		&(c->rlx_cp[c->Xsidx]), &(c->rly_cp[c->Xsidx]),
		&(c->rllen[c->Xsidx]), c->Xeidx-c->Xsidx+1 );
	double area_rr = calcRulCross( &(c->Xx2d[c->Xsidx]), &(c->Xy2d[c->Xsidx]),
		&(c->rrx_cp[c->Xsidx]), &(c->rry_cp[c->Xsidx]),
		&(c->rrlen[c->Xsidx]), c->Xeidx-c->Xsidx+1 );
	area = area_rl + area_rr;
	if( area>0.0 ){ ret = 1; }
	return ret;
}

int papermodel::checkRulCreaseCross()
{
	int ret=0;
	crease *c = &crs[0];

	for( int i=c->Xsidx; i<=c->Xeidx; i++ ){
		double lex, ley, rex, rey;
		lex = c->Xx2d[i] + c->rllen[i] * c->rlx_cp[i];
		ley = c->Xy2d[i] + c->rllen[i] * c->rly_cp[i];
		rex = c->Xx2d[i] + c->rrlen[i] * c->rrx_cp[i];
		rey = c->Xy2d[i] + c->rrlen[i] * c->rry_cp[i];
		for( int j=c->Xsidx; j<c->Xeidx; j++ ){
			double ix, iy, l0, l1;
			int res = intersectionOfLine( c->Xx2d[i], c->Xy2d[i], lex, ley,
				c->Xx2d[j], c->Xy2d[j],c->Xx2d[j+1], c->Xy2d[j+1], &ix, &iy, &l0, &l1 );
			if( res==0 && l0>0.0 && l1>0.0 ){
				ret = -1;
				break;
			}
		}
		if( ret<0 ){
			break;
		}
		for( int j=c->Xsidx; j<c->Xeidx; j++ ){
			double ix, iy, l0, l1;
			int res = intersectionOfLine( c->Xx2d[i], c->Xy2d[i], rex, rey,
				c->Xx2d[j], c->Xy2d[j],c->Xx2d[j+1], c->Xy2d[j+1], &ix, &iy, &l0, &l1 );
			if( res==0 && l0>0.0 && l1>0.0 ){
				ret = -1;
				break;
			}
		}
		if( ret<0 ){
			break;
		}
	}

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

int papermodel::calcAvetgap()
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

	double dx, dy, dz, diff, diff0 = 0.0, diff1 = 0.0;
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
	}
#endif
	//std::cout << "diff0=" << diff0/(double)tgcnt << ", diff1=" << diff1/(double)tgcnt << std::endl;
	this->avetgap = diff0 / (double)tgcnt;

	return ret;
}

