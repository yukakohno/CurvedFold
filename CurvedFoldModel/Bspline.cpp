#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "Bspline.h"
#include "util.h"

int Bspline::calcT(int vcnt, double tmax, int toption, double *t, double *vx, double *vy, double *vz)
{
	int i,j,k,ret=0;

	switch( toption ){
	case 0: // parameter t（等間隔）
		for( i=0; i<vcnt-1; i++ ){
			t[i] = (double)i / (double)(vcnt-1) * (double)tmax;
		}
		t[vcnt-1] = (double)tmax - 0.0001;
		break;
	case 1: // 点間距離に応じて
		if( vx==NULL && vy==NULL && vz==NULL ){
			ret = -1; goto end;
		}
		t[0] = 0.0;
		for( i=1; i<vcnt; i++ ){
			double dx=0, dy=0, dz=0;
			if(vx){ dx = vx[i]-vx[i-1]; }
			if(vy){ dy = vy[i]-vy[i-1]; }
			if(vz){ dz = vz[i]-vz[i-1]; }
			t[i] = t[i-1] + sqrt(dx*dx + dy*dy + dz*dz);
		}
		for( i=0; i<vcnt; i++ ){
			t[i] = t[i] / t[vcnt-1] * (double)tmax;
		}
		t[vcnt-1] -= 0.0001;
		break;
	}
end:
	return ret;
}

double Bspline::baseN(int i,int k,double t,int *nv)
{
	double w1=0.0,w2=0.0;
	if(k==1){
		if(t >=nv[i] && t<nv[i+1]) return 1.0;
		//if(t >=nv[i] && t<=nv[i+1]) return 1.0;
		else return 0.0;
	}
	else {
		if((nv[i+k]-nv[i+1])!=0){
			w1=((nv[i+k]-t)/(nv[i+k]-nv[i+1])) * baseN(i+1,k-1,t,nv);
		}
		if((nv[i+k-1]-nv[i]) !=0 ){
			w2=((t-nv[i])/(nv[i+k-1]-nv[i])) * baseN(i,k-1,t,nv);
		}
		return (w1+w2);
	}
}

int Bspline::calcBSpline(double *cx, double *cy, double *cz, int ccnt,
						 int *nv, int nvcnt, int deg, int pcnt, double *t,
						 double *vx, double *vy, double *vz )
{
	int i,j;
	double *bn = new double[ccnt*pcnt];

	// r(t)=ΣNi,k(t）Pi
	for( i=0; i<pcnt; i++ ){
		if( vx ){ vx[i] = 0.0; }
		if( vy ){ vy[i] = 0.0; }
		if( vz ){ vz[i] = 0.0; }
		for( j=0; j<ccnt; j++ ){
			bn[i*ccnt+j] = baseN( j, deg, t[i], nv);
		}
		if( t[i]==nv[nvcnt-1] ){ // last
			for( j=0; j<ccnt-1; j++ ){
				bn[i*ccnt+j] = 0.0;
			}
			bn[i*ccnt+j] = 1.0;
		}
		for( j=0; j<ccnt; j++ ){
			//vx[i] += baseN( j, deg, t[i], nv) * cx[j];
			//vy[i] += baseN( j, deg, t[i], nv) * cy[j];
			if( cx && vx ){ vx[i] += bn[i*ccnt+j] * cx[j]; }
			if( cy && vy ){ vy[i] += bn[i*ccnt+j] * cy[j]; }
			if( cz && vz ){ vz[i] += bn[i*ccnt+j] * cz[j]; }
		}
	}
#if 0
	FILE *fp = fopen("output/debug_calcBSpline.csv", "w");
	if(fp){
		for( i=0; i<pcnt; i++ ){
			fprintf(fp, "%f,", t[i]);
			for( j=0; j<ccnt; j++ ){
				fprintf(fp, "%f,", bn[i*ccnt+j]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		for( i=0; i<ccnt; i++ ){
			if( cx ){ fprintf(fp, "%f,", cx[i]); }
			if( cy ){ fprintf(fp, "%f,", cy[i]); }
			if( cz ){ fprintf(fp, "%f,", cz[i]); }
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		for( i=0; i<pcnt; i++ ){
			if( vx ){ fprintf(fp, "%f,", vx[i]); }
			if( vy ){ fprintf(fp, "%f,", vy[i]); }
			if( vz ){ fprintf(fp, "%f,", vz[i]); }
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
#endif

end:
	if(bn){ delete bn; bn=NULL; }
	return 0;
}


int Bspline::calcCP( double *vx, double *vy, double *vz, int vcnt,
					int *nv, int nvcnt, int deg, double *t, int ccnt,
					double *cx, double *cy, double *cz )
{
	int i,j,k,tmax,ret=0;
	double *bn=NULL, *qx=NULL, *qy=NULL, *qz=NULL, *q=NULL, *m=NULL, *mi=NULL, *mtmp=NULL;
	bn = new double[ccnt*vcnt];
	if( vx && cx ){ qx = new double[vcnt]; }
	if( vy && cy ){ qy = new double[vcnt]; }
	if( vz && cz ){ qz = new double[vcnt]; }
	q = new double[ccnt*3];
	m = new double[ccnt*ccnt];
	mi = new double[ccnt*ccnt];
	mtmp = new double[(ccnt-2)*(ccnt-2)];

	// matrix N (base)
	for( i=0; i<vcnt; i++ ){ // data points -> row i
		for( j=0; j<ccnt; j++ ){ // control points -> col j
			bn[i*ccnt+j] = baseN( j, deg, t[i], nv);
		}
	}
	// matrix M=NtN
	for( i=0; i<ccnt*ccnt; i++ ){ m[i] = 0.0; }
	for( i=1; i<ccnt-1; i++ ){
		for( j=1; j<ccnt-1; j++ ){
			for( k=1; k<vcnt-1; k++ ){
				m[i*ccnt+j] += bn[k*ccnt+i]*bn[k*ccnt+j];
			}
		}
	}

	// matrix M-1
	for( i=1; i<ccnt-1; i++ ){
		for( j=1; j<ccnt-1; j++ ){
			mtmp[(i-1)*(ccnt-2)+(j-1)] = m[i*ccnt+j];
		}
	}
	ret = m_inverse( mtmp, ccnt-2 );
	if( ret < 0 ){
		goto end;
	}
	for( i=0; i<ccnt*ccnt; i++ ){ mi[i] = 0.0; }
	for( i=1; i<ccnt-1; i++ ){
		for( j=1; j<ccnt-1; j++ ){
			mi[i*ccnt+j] = mtmp[(i-1)*(ccnt-2)+(j-1)];
		}
	}

	if( cx && vx ){ cx[0] = vx[0]; }
	if( cy && vy ){ cy[0] = vy[0]; }
	if( cz && vz ){ cz[0] = vz[0]; }
	if( cx && vx ){ cx[ccnt-1] = vx[vcnt-1]; }
	if( cy && vy ){ cy[ccnt-1] = vy[vcnt-1]; }
	if( cz && vz ){ cz[ccnt-1] = vz[vcnt-1]; }

	for( i=1; i<vcnt-1; i++ ){
		if( qx && vx ){ qx[i] = vx[i] - bn[i*ccnt+0]*vx[0] - bn[i*ccnt+(ccnt-1)]*vx[vcnt-1]; }
		if( qy && vy ){ qy[i] = vy[i] - bn[i*ccnt+0]*vy[0] - bn[i*ccnt+(ccnt-1)]*vy[vcnt-1]; }
		if( qz && vz ){ qz[i] = vz[i] - bn[i*ccnt+0]*vz[0] - bn[i*ccnt+(ccnt-1)]*vz[vcnt-1]; }
	}

	// matrix Q
	for( i=0; i<ccnt*3; i++ ){ q[i] = 0.0; }
	for( j=1; j<ccnt-1; j++ ){
		for( i=1; i<vcnt-1; i++ ){
			if( qx ){ q[j*3  ] += bn[i*ccnt+j] * qx[i]; }
			if( qy ){ q[j*3+1] += bn[i*ccnt+j] * qy[i]; }
			if( qz ){ q[j*3+2] += bn[i*ccnt+j] * qz[i]; }
		}
	}

	// MP=Q -> P=M-1 Q
	for( j=1; j<ccnt-1; j++ ){
		if( cx ){ cx[j] = 0.0; }
		if( cy ){ cy[j] = 0.0; }
		if( cz ){ cz[j] = 0.0; }
		for( i=1; i<ccnt-1; i++ ){
			if( cx ){ cx[j] += mi[i*ccnt+j] * q[i*3  ]; }
			if( cy ){ cy[j] += mi[i*ccnt+j] * q[i*3+1]; }
			if( cz ){ cz[j] += mi[i*ccnt+j] * q[i*3+2]; }
		}
	}

end:
#if 0
	FILE *fp = fopen("output/debug_calcCP.csv", "w");
	if(fp){
		for( i=0; i<ccnt; i++ ){
			for( j=0; j<ccnt; j++ ){
				fprintf(fp, "%f,", m[i*ccnt+j]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		for( i=0; i<ccnt; i++ ){
			for( j=0; j<ccnt; j++ ){
				fprintf(fp, "%f,", mi[i*ccnt+j]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		for( j=0; j<ccnt; j++ ){
			fprintf(fp, "%f,%f,%f\n", q[j*3  ], q[j*3+1], q[j*3+2]);
		}
		fprintf(fp, "\n");
		for( i=0; i<ccnt; i++ ){
			if( cx ){ fprintf(fp, "%f,", cx[i]); }
			if( cy ){ fprintf(fp, "%f,", cy[i]); }
			if( cz ){ fprintf(fp, "%f,", cz[i]); }
			fprintf(fp, "\n");
		}
		fclose(fp); fp=NULL;
	}
#endif
	if(bn){ delete bn; bn=NULL; }
	if(qx){ delete qx; qx=NULL; }
	if(qy){ delete qy; qy=NULL; }
	if(qz){ delete qz; qz=NULL; }
	if(q){ delete q; q=NULL; }
	if(m){ delete m; m=NULL; }
	if(mi){ delete mi; mi=NULL; }
	if(mtmp){ delete mtmp; mtmp=NULL; }

	return ret;
}

double baseN3( int j, double t )
{
	if( t<0 || 3<t || j<0 || 4<j){
		return -1.0;
	}
	if( 0<= t && t <= 1 ){
		switch( j ){
			case 0:	return (t-1.0)*(t-1.0);	break;
			case 1:	return t*(1.0-t) + t*(2.0-t)/2.0;	break;
			case 2:	return t*t/2.0;	break;
			case 3:	return 0.0;	break;
			case 4:	return 0.0;	break;
		}
	} else if( t <= 2 ){
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return (2.0-t)*(2.0-t)/2.0;	break;
			case 2:	return t*(2.0-t)/2.0 + (3.0-t)*(t-1.0)/2.0;	break;
			case 3:	return (t-1.0)*(t-1.0)/2.0;	break;
			case 4:	return 0.0;	break;
		}
	} else {
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return 0.0;	break;
			case 2:	return (3.0-t)*(3.0-t)/2.0;	break;
			case 3:	return (t-1.0)*(3.0-t)/2.0 + (3.0-t)*(t-2.0);	break;
			case 4:	return (t-2.0)*(t-2.0);	break;
		}
	}
	return -1.0;
}

double baseN3Diff( int j, double t )
{
	if( t<0 || 3<t || j<0 || 4<j){
		return -1.0;
	}
	if( 0<= t && t <= 1 ){
		switch( j ){
			case 0:	return 2.0*t-2.0;	break;
			case 1:	return 2.0-3.0*t;	break;
			case 2:	return t;	break;
			case 3:	return 0.0;	break;
			case 4:	return 0.0;	break;
		}
	} else if( t <= 2 ){
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return t-2.0;	break;
			case 2:	return 3.0-2.0*t;	break;
			case 3:	return t-1.0;	break;
			case 4:	return 0.0;	break;
		}
	} else {
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return 0.0;	break;
			case 2:	return t-3.0;	break;
			case 3:	return 7.0-3.0*t;	break;
			case 4:	return 2.0*t-4.0;	break;
		}
	}
	return -1.0;
}

double baseN3Diff2( int j, double t )
{
	if( t<0 || 3<t || j<0 || 4<j){
		return -1.0;
	}
	if( 0<= t && t <= 1 ){
		switch( j ){
			case 0:	return 2.0;	break;
			case 1:	return -3.0;	break;
			case 2:	return 1.0;	break;
			case 3:	return 0.0;	break;
			case 4:	return 0.0;	break;
		}
	} else if( t <= 2 ){
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return 1.0;	break;
			case 2:	return -2.0;	break;
			case 3:	return 1.0;	break;
			case 4:	return 0.0;	break;
		}
	} else {
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return 0.0;	break;
			case 2:	return 1.0;	break;
			case 3:	return -3.0;	break;
			case 4:	return 2.0;	break;
		}
	}
	return -1.0;
}

double baseN4( int j, double t )
{
	if( t<0 || 3<t || j<0 || 5<j){
		return -1.0;
	}
	if( 0<= t && t <= 1 ){
		switch( j ){
			case 0:	return (1.0-t)*(1.0-t)*(1.0-t);	break;
			case 1:	return t*(1.0-t)*(1.0-t) + t*(2.0-t)*(1.0-t)/2.0 + t*(2.0-t)*(2.0-t)/4.0;	break;
			case 2:	return t*t*(1.0-t)/2.0 + t*t*(2.0-t)/4.0 + t*t*(3.0-t)/6.0;	break;
			case 3:	return t*t*t/6.0;	break;
			case 4:	return 0.0;	break;
			case 5:	return 0.0;	break;
		}
	} else if( t <= 2 ){
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return (2.0-t)*(2.0-t)*(2.0-t)/4.0;	break;
			case 2:	return t*(2.0-t)*(2.0-t)/4.0 + t*(3.0-t)*(2.0-t)/6.0 + (3.0-t)*(3.0-t)*(t-1.0)/6.0;	break;
			case 3:	return t*t*(2.0-t)/6.0 + t*(3.0-t)*(t-1.0)/6.0 + (3.0-t)*(t-1.0)*(t-1.0)/4.0;	break;
			case 4:	return (t-1.0)*(t-1.0)*(t-1.0)/4.0;	break;
			case 5:	return 0.0;	break;
		}
	} else {
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return 0.0;	break;
			case 2:	return (3.0-t)*(3.0-t)*(3.0-t)/6.0;	break;
			case 3:	return t*(3.0-t)*(3.0-t)/6.0 + (3.0-t)*(t-1.0)*(3.0-t)/4.0 + (3.0-t)*(3.0-t)*(t-2.0)/2.0;	break;
			case 4:	return (t-1.0)*(t-1.0)*(3.0-t)/4.0 + (t-1.0)*(3.0-t)*(t-2.0)/2.0 + (3.0-t)*(t-2.0)*(t-2.0);	break;
			case 5:	return (t-2.0)*(t-2.0)*(t-2.0);	break;
		}
	}
	return -1.0;
}

double baseN4Diff( int j, double t )
{
	if( t<0 || 3<t || j<0 || 5<j){
		return -1.0;
	}
	if( 0<= t && t <= 1 ){
		switch( j ){
			case 0:	return -3.0*(t-1)*(t-1);	break;
			case 1:	return t*(2.0*t-2.0)+(t*(2.0*t-4.0))/4.0+((t-1.0)*(t-2.0))/2.0+(t-1.0)*(t-1.0)+(t-2.0)*(t-2.0)/4.0+(t*(t-1.0))/2.0+(t*(t-2.0))/2.0;	break;
			case 2:	return -t*(t-1.0)-(t*(t-2.0))/2.0-(t*(t-3.0))/3.0-(11.0*t*t)/12.0;	break;
			case 3:	return t*t/2;	break;
			case 4:	return 0.0;	break;
			case 5:	return 0.0;	break;
		}
	} else if( t <= 2 ){
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return -(3.0*(t-2)*(t-2))/4.0;	break;
			case 2:	return (t*(2.0*t-4.0))/4.0+((t-2.0)*(t-3.0))/6.0+(t-2.0)*(t-2.0)/4.0+(t-3.0)*(t-3.0)/6.0+((2.0*t-6.0)*(t-1.0))/6.0+(t*(t-2.0))/6.0+(t*(t-3.0))/6.0;	break;
			case 3:	return -((t-1.0)*(t-3.0))/6.0-(t-1.0)*(t-1.0)/4.0-((2.0*t-2.0)*(t-3.0))/4.0-(t*(t-1))/6.0-(t*(t-2.0))/3.0-(t*(t-3.0))/6.0-t*t/6.0;	break;
			case 4:	return (3*(t-1)*(t-1))/4;	break;
			case 5:	return 0.0;	break;
		}
	} else {
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return 0.0;	break;
			case 2:	return -(t-3.0)*(t-3.0)/2.0;	break;
			case 3:	return (t*(2.0*t-6.0))/6.0+(11.0*(t-3.0)*(t-3.0))/12.0+((2.0*t-6.0)*(t-1.0))/4.0+((2.0*t-6.0)*(t-2.0))/2.0;	break;
			case 4:	return -((t-1.0)*(t-2.0))/2.0-((t-1.0)*(t-3.0))/2.0-((t-2.0)*(t-3.0))/2.0-(t-1.0)*(t-1.0)/4.0-(t-2.0)*(t-2.0)-((2.0*t-2.0)*(t-3.0))/4.0-(2.0*t-4.0)*(t-3.0);	break;
			case 5:	return 3.0*(t-2.0)*(t-2.0);	break;
		}
	}
	return -1.0;
}

double baseN4Diff2( int j, double t )
{
	if( t<0 || 3<t || j<0 || 5<j){
		return -1.0;
	}
	if( 0<= t && t <= 1 ){
		switch( j ){
			case 0:	return 6.0-6.0*t;	break;
			case 1:	return (21.0*t)/2.0-9.0;	break;
			case 2:	return 3.0-(11.0*t)/2.0;	break;
			case 3:	return t;	break;
			case 4:	return 0.0;	break;
			case 5:	return 0.0;	break;
		}
	} else if( t <= 2 ){
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return 3.0-(3.0*t)/2.0;	break;
			case 2:	return (7.0*t)/2.0-6.0;	break;
			case 3:	return 9.0/2.0-(7.0*t)/2.0;	break;
			case 4:	return (3.0*t)/2.0-3.0/2.0;	break;
			case 5:	return 0.0;	break;
		}
	} else {
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return 0.0;	break;
			case 2:	return 3-t;	break;
			case 3:	return (11*t)/2-27/2;	break;
			case 4:	return 45/2-(21*t)/2;	break;
			case 5:	return 6*t-12;	break;
		}
	}
	return -1.0;
}

double baseN4Diff3( int j, double t )
{
	if( t<0 || 3<t || j<0 || 5<j){
		return -1.0;
	}
	if( 0<= t && t <= 1 ){
		switch( j ){
			case 0:	return -6.0;	break;
			case 1:	return 21.0/2.0;	break;
			case 2:	return -11.0/2.0;	break;
			case 3:	return 1.0;	break;
			case 4:	return 0.0;	break;
			case 5:	return 0.0;	break;
		}
	} else if( t <= 2 ){
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return -3.0/2.0;	break;
			case 2:	return 7.0/2.0;	break;
			case 3:	return -7.0/2.0;	break;
			case 4:	return 3.0/2.0;	break;
			case 5:	return 0.0;	break;
		}
	} else {
		switch( j ){
			case 0:	return 0.0;	break;
			case 1:	return 0.0;	break;
			case 2:	return -1.0;	break;
			case 3:	return 11.0/2.0;	break;
			case 4:	return -21.0/2.0;	break;
			case 5:	return 6.0;	break;
		}
	}
	return -1.0;
}
int Bspline::calcBSDiff(double *cx, double *cy, double *cz, int ccnt,
						int *nv, int nvcnt, int deg, int vcnt, double *t,
						double *vx, double *vy, double *vz, int diffn )
{
	int i,j, ret=0;

	if(deg==3){
		if( nvcnt!=8 || nv[0]!=0 || nv[1]!=0 || nv[2]!=0 || nv[3]!=1 || nv[4]!=2 || nv[5]!=3 || nv[6]!=3 || nv[7]!=3 ){
			ret = -1; goto end;
		}
	} else if(deg==4){

	} else{
		ret = -1; goto end;
	}

	double *bn = new double[ccnt*vcnt];

	for( i=0; i<vcnt; i++ ){
		if( vx ){ vx[i] = 0.0; }
		if( vy ){ vy[i] = 0.0; }
		if( vz ){ vz[i] = 0.0; }
		for( j=0; j<ccnt; j++ ){
			if( deg==3 ){
				switch( diffn ){
					case 0: bn[i*ccnt+j] = baseN3( j, t[i]); break;
					case 1: bn[i*ccnt+j] = baseN3Diff( j, t[i]); break;
					case 2: bn[i*ccnt+j] = baseN3Diff2( j, t[i]); break;
				}
			} else if( deg==4 ){
				switch( diffn ){
					case 0: bn[i*ccnt+j] = baseN4( j, t[i]); break;
					case 1: bn[i*ccnt+j] = baseN4Diff( j, t[i]); break;
					case 2: bn[i*ccnt+j] = baseN4Diff2( j, t[i]); break;
					case 3: bn[i*ccnt+j] = baseN4Diff3( j, t[i]); break;
				}
			}
		}
		for( j=0; j<ccnt; j++ ){
			if( cx && vx ){ vx[i] += bn[i*ccnt+j] * cx[j]; }
			if( cy && vy ){ vy[i] += bn[i*ccnt+j] * cy[j]; }
			if( cz && vz ){ vz[i] += bn[i*ccnt+j] * cz[j]; }
		}
	}
#if 0
	FILE *fp = fopen("output/debug_calcBSplineDiff.csv", "w");
	if(fp){
		for( i=0; i<vcnt; i++ ){
			fprintf(fp, "%f,", t[i]);
			for( j=0; j<ccnt; j++ ){
				fprintf(fp, "%f,", bn[i*ccnt+j]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		for( i=0; i<ccnt; i++ ){
			if( cx ){ fprintf(fp, "%f,", cx[i]); }
			if( cy ){ fprintf(fp, "%f,", cy[i]); }
			if( cz ){ fprintf(fp, "%f,", cz[i]); }
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		for( i=0; i<vcnt; i++ ){
			if( vx ){ fprintf(fp, "%f,", vx[i]); }
			if( vy ){ fprintf(fp, "%f,", vy[i]); }
			if( vz ){ fprintf(fp, "%f,", vz[i]); }
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
#endif

end:
	if(bn){ delete bn; bn=NULL; }
	return 0;
}
