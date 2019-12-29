#include <stdio.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <memory.h>
#include "util.h"
#include "Quaternion.h"

// return value 0: end, 1: continue
int csvread(char **_in,		// input line 
			char *ret,		// item
			char sep,		// separator charactor
			int flg_con)	// 1: egnore consecutive separator
{
	int cnt=0;
	char c, *in = *_in;

	while (1){
		c = in[cnt];

		if(!c || c == 10) // 10: new line
		{
			ret[cnt]='\0';
			*_in = &in[cnt];
			return 0;
		}

		// separater
		if (c == sep)
		{
			int i = 1;
			ret[cnt]='\0';
			// continue until NOT separater
			if(flg_con){
				while((c = in[cnt+i]) == sep){
					i++;
				}
				if(!c || c == 10) // NULL / new line
				{
					*_in = &in[cnt];
					return 0;
				}
			}
			*_in = &in[cnt+i];
			return 1;
		}

		ret[cnt]=(char) c;
		cnt++;
	}
}

void unit_m44(float *m44)
{
	m44[0] = 1.0;	m44[1] = 0.0;	m44[2] = 0.0;	m44[3] = 0.0;
	m44[4] = 0.0;	m44[5] = 1.0;	m44[6] = 0.0;	m44[7] = 0.0;
	m44[8] = 0.0;	m44[9] = 0.0;	m44[10] = 1.0;	m44[11] = 0.0;
	m44[12] = 0.0;	m44[13] = 0.0;	m44[14] = 0.0;	m44[15] = 1.0;
}
void unit_m44(double *m44)
{
	m44[0] = 1.0;	m44[1] = 0.0;	m44[2] = 0.0;	m44[3] = 0.0;
	m44[4] = 0.0;	m44[5] = 1.0;	m44[6] = 0.0;	m44[7] = 0.0;
	m44[8] = 0.0;	m44[9] = 0.0;	m44[10] = 1.0;	m44[11] = 0.0;
	m44[12] = 0.0;	m44[13] = 0.0;	m44[14] = 0.0;	m44[15] = 1.0;
}

void unit_m33(float *m33)
{
	m33[0] = 1.0;	m33[1] = 0.0;	m33[2] = 0.0;
	m33[3] = 0.0;	m33[4] = 1.0;	m33[5] = 0.0;
	m33[6] = 0.0;	m33[7] = 0.0;	m33[8] = 1.0;
}
void unit_m33(double *m33)
{
	m33[0] = 1.0;	m33[1] = 0.0;	m33[2] = 0.0;
	m33[3] = 0.0;	m33[4] = 1.0;	m33[5] = 0.0;
	m33[6] = 0.0;	m33[7] = 0.0;	m33[8] = 1.0;
}

void mult_m44_v4(float *_m44, float *_v4)
{
	int i,j;
	float v0[4], m[4][4], v1[4];

	memcpy(v1, _v4, sizeof(float)*4);
	memcpy(m, _m44, sizeof(float)*16);
	memset(v0, 0, sizeof(float)*4);

	for(i = 0; i < 4; i++){
		for(j = 0; j < 4; j++){
			v0[i] += m[i][j] * v1[j];
		}
	}
	memcpy(_v4, v0, sizeof(float)*4);
}

void mult_m44_v4(double *_m44, double *_v4)
{
	int i,j;
	double v0[4], m[4][4], v1[4];

	memcpy(v1, _v4, sizeof(double)*4);
	memcpy(m, _m44, sizeof(double)*16);
	memset(v0, 0, sizeof(double)*4);

	for(i = 0; i < 4; i++){
		for(j = 0; j < 4; j++){
			v0[i] += m[i][j] * v1[j];
		}
	}
	memcpy(_v4, v0, sizeof(double)*4);
}

void mult_v4_m44(float *_v4, float *_m44)
{
	int i,j;
	float v0[4], v1[4], m[4][4];

	memcpy(v1, _v4, sizeof(float)*4);
	memcpy(m, _m44, sizeof(float)*16);
	memset(v0, 0, sizeof(float)*4);

	for(j = 0; j < 4; j++){
		for(i = 0; i < 4; i++){
			v0[j] += m[i][j] * v1[i];
		}
	}
	memcpy(_v4, v0, sizeof(float)*4);
}

void mult_v4_m44(double *_v4, double *_m44)
{
	int i,j;
	double v0[4], v1[4], m[4][4];

	memcpy(v1, _v4, sizeof(double)*4);
	memcpy(m, _m44, sizeof(double)*16);
	memset(v0, 0, sizeof(double)*4);

	for(j = 0; j < 4; j++){
		for(i = 0; i < 4; i++){
			v0[j] += m[i][j] * v1[i];
		}
	}
	memcpy(_v4, v0, sizeof(double)*4);
}

void mult_m44_n44(float *_n44, float *_m44, int dst)
{
	int i,j,k;
	float m0[16], m1[16], m2[16];

	memcpy(m1, _n44, sizeof(float)*16);
	memcpy(m2, _m44, sizeof(float)*16);
	memset(m0, 0, sizeof(float)*16);

	for(i = 0; i < 4; i++){
		for(j = 0; j < 4; j++){
			for(k = 0; k < 4; k++){
				m0[i*4+j] += m1[i*4+k] * m2[k*4+j];
			}
		}
	}
	if(dst == 0){
		memcpy(_m44, m0, sizeof(float)*16);
	} else{
		memcpy(_n44, m0, sizeof(float)*16);
	}
}

void mult_m44_n44(double *_n44, double *_m44, int dst)
{
	int i,j,k;
	double m0[16], m1[16], m2[16];

	memcpy(m1, _n44, sizeof(double)*16);
	memcpy(m2, _m44, sizeof(double)*16);
	memset(m0, 0, sizeof(double)*16);

	for(i = 0; i < 4; i++){
		for(j = 0; j < 4; j++){
			for(k = 0; k < 4; k++){
				m0[i*4+j] += m1[i*4+k] * m2[k*4+j];
			}
		}
	}
	if(dst == 0){
		memcpy(_m44, m0, sizeof(double)*16);
	} else{
		memcpy(_n44, m0, sizeof(double)*16);
	}
}

void mult_m33_n33(float *_n33, float *_m33, int dst)
{
	int i,j,k;
	float m0[9], m1[9], m2[9];

	memcpy(m1, _n33, sizeof(float)*9);
	memcpy(m2, _m33, sizeof(float)*9);
	memset(m0, 0, sizeof(float)*9);

	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			for(k = 0; k < 3; k++){
				m0[i*3+j] += m1[i*3+k] * m2[k*3+j];
			}
		}
	}
	if(dst == 0){
		memcpy(_m33, m0, sizeof(float)*9);
	} else{
		memcpy(_n33, m0, sizeof(float)*9);
	}
}

void mult_m33_n33(double *_n33, double *_m33, int dst)
{
	int i,j,k;
	double m0[9], m1[9], m2[9];

	memcpy(m1, _n33, sizeof(double)*9);
	memcpy(m2, _m33, sizeof(double)*9);
	memset(m0, 0, sizeof(double)*9);

	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			for(k = 0; k < 3; k++){
				m0[i*3+j] += m1[i*3+k] * m2[k*3+j];
			}
		}
	}
	if(dst == 0){
		memcpy(_m33, m0, sizeof(double)*9);
	} else{
		memcpy(_n33, m0, sizeof(double)*9);
	}
}

void mult_m44_n44_rot(float *_n44, float *_m44, int dst)
{
	int i,j,k;
	float m0[16], m1[16], m2[16];

	memcpy(m1, _n44, sizeof(float)*16);
	memcpy(m2, _m44, sizeof(float)*16);
	memset(m0, 0, sizeof(float)*16);

	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			for(k = 0; k < 3; k++){
				m0[i*4+j] += m1[i*4+k] * m2[k*4+j];
			}
		}
	}

	if(dst == 0){
		m0[3]=_m44[3]; m0[7]=_m44[7]; m0[11]=_m44[11];
		memcpy(_m44, m0, sizeof(float)*12);
	} else{
		m0[3]=_m44[3]; m0[7]=_m44[7]; m0[11]=_m44[11];
		memcpy(_n44, m0, sizeof(float)*12);
	}
}

void mult_m44_n44_rot(double *_n44, double *_m44, int dst)
{
	int i,j,k;
	double m0[16], m1[16], m2[16];

	memcpy(m1, _n44, sizeof(double)*16);
	memcpy(m2, _m44, sizeof(double)*16);
	memset(m0, 0, sizeof(double)*16);

	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			for(k = 0; k < 3; k++){
				m0[i*4+j] += m1[i*4+k] * m2[k*4+j];
			}
		}
	}

	if(dst == 0){
		m0[3]=_m44[3]; m0[7]=_m44[7]; m0[11]=_m44[11];
		memcpy(_m44, m0, sizeof(double)*12);
	} else{
		m0[3]=_n44[3]; m0[7]=_n44[7]; m0[11]=_n44[11];
		memcpy(_n44, m0, sizeof(double)*12);
	}
}

int inv_m44(float *_m)
{
	int i, ret;
	double m[16];

	for(i = 0; i < 16; i++){ m[i] = (double)_m[i]; }

	ret = m_inverse( m, 4 );

	for(i = 0; i < 16; i++){ _m[i] = (float)m[i]; }

	return ret;
}

int inv_m44(double *_m)
{
	return m_inverse( _m, 4 );
}

int inv_m33(float *m)
{
	float det, mat[9];
	det = m[0]*m[4]*m[8] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7] - m[0]*m[5]*m[7] - m[1]*m[3]*m[8] - m[2]*m[4]*m[6];
	if(fabs(det) < DBL_MIN)
		return -1;
	mat[0] = (m[4]*m[8]-m[5]*m[7])/det;
	mat[1] = (m[2]*m[7]-m[1]*m[8])/det;
	mat[2] = (m[1]*m[5]-m[2]*m[4])/det;
	mat[3] = (m[5]*m[6]-m[3]*m[8])/det;
	mat[4] = (m[0]*m[8]-m[2]*m[6])/det;
	mat[5] = (m[2]*m[3]-m[0]*m[5])/det;
	mat[6] = (m[3]*m[7]-m[4]*m[6])/det;
	mat[7] = (m[1]*m[6]-m[0]*m[7])/det;
	mat[8] = (m[0]*m[4]-m[1]*m[3])/det;
	memcpy(m,mat,sizeof(float)*9);
	return 0;
}

int inv_m33(double *m)
{
	double det, mat[9];
	det = m[0]*m[4]*m[8] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7] - m[0]*m[5]*m[7] - m[1]*m[3]*m[8] - m[2]*m[4]*m[6];
	if(fabs(det) < DBL_MIN)
		return -1;
	mat[0] = (m[4]*m[8]-m[5]*m[7])/det;
	mat[1] = (m[2]*m[7]-m[1]*m[8])/det;
	mat[2] = (m[1]*m[5]-m[2]*m[4])/det;
	mat[3] = (m[5]*m[6]-m[3]*m[8])/det;
	mat[4] = (m[0]*m[8]-m[2]*m[6])/det;
	mat[5] = (m[2]*m[3]-m[0]*m[5])/det;
	mat[6] = (m[3]*m[7]-m[4]*m[6])/det;
	mat[7] = (m[1]*m[6]-m[0]*m[7])/det;
	mat[8] = (m[0]*m[4]-m[1]*m[3])/det;
	memcpy(m,mat,sizeof(double)*9);
	return 0;
}

/* original code written by Y.Fukui
m[][30]: s—ñ‚ª“ü‚é”z—ñi‚R‚O~‚R‚OŽŸŒ³‚Ì‚à‚Ì‚Ü‚Å‰Â”\j
‚±‚ÌŠÖ”‚Ìˆ—‚ÌŒ‹‰ÊA‹ts—ñ‚ª“ü‚é”z—ñi“üo—Íˆø”j
l:@‹ts—ñ‚ð‹‚ß‚é‚×‚«‚Ž~‚Ž(n<30)‚Ì³•ûs—ñ‚ÌŽŸ”i“ü—Íˆø”j*/
int m_inverse( double *m, int l )
{
	double	eps,p,aa,w,a[30][30];
	int	i,j,nn,ip,nw,noseq[30];

	eps=1.e-10;
	for( i=0; i<l; i++ ) {
		for( j=0; j<l; j++ ) a[i][j] = m[i*l+j]; // revised by Y.Kohno
	}
	for( i=0; i< l; i++ ) noseq[i] = i;
	for( nn=0; nn< l; nn++ ) {
		p = 0.0;
		for( i=nn; i < l; i++ ) {
			aa = fabs( a[i][0] );
			if( p < aa ) { p = aa; ip = i;	}
		}
		if( p < eps ) return( -1 );
		nw = noseq[ip];
		noseq[ip] = noseq[nn];
		noseq[nn] = nw;
		for( j=0; j< l; j++ ) {
			w=a[ip][j]; a[ip][j]=a[nn][j]; a[nn][j]=w;
		}
		w = a[nn][0];
		for( j=1; j< l; j++ ) a[nn][j-1] = a[nn][j] / w;
		a[nn][l-1] = 1.0 / w;

		for( i=0; i< l; i++ ){
			if( i != nn ) {
				w = a[i][0];
				for(j=1; j< l; j++) a[i][j-1] = a[i][j] - w*a[nn][j-1];
				a[i][l-1] = -w*a[nn][l-1];
			}
		}
	}
	for( nn=0; nn< l; nn++ ) {
		for( j=nn; j < l; j++ ) {
			if( noseq[j] == nn ) break;
		}
		noseq[j] = noseq[nn];
		for( i=0; i<l; i++ ) {
			w=a[i][j]; a[i][j]=a[i][nn]; a[i][nn]=w;
		}
	}
	for( i=0; i<l; i++ ) {
		for( j=0; j<l; j++ ) m[i*l+j] = a[i][j]; // revised by Y.Kohno
	}
	return( 0 );
}

int mat_quat(float *mat, float *quat, int matsize)
{
	float m11,m12,m13, m21,m22,m23, m31,m32,m33, qx,qy,qz,qw;
	if(matsize == 4){
		m11 = mat[0];	m12 = mat[1];	m13 = mat[2];
		m21 = mat[4];	m22 = mat[5];	m23 = mat[6];
		m31 = mat[8];	m32 = mat[9];	m33 = mat[10];
	} else if(matsize == 3){
		m11 = mat[0];	m12 = mat[1];	m13 = mat[2];
		m21 = mat[3];	m22 = mat[4];	m23 = mat[5];
		m31 = mat[6];	m32 = mat[7];	m33 = mat[8];
	} else{
		printf("invalid input @mat_quat(), matsize = %d\n", matsize);
	}

	float elem[ 4 ]; // 0:x, 1:y, 2:z, 3:w
	elem[ 0 ] = m11 - m22 - m33 + 1.0f;
	elem[ 1 ] = -m11 + m22 - m33 + 1.0f;
	elem[ 2 ] = -m11 - m22 + m33 + 1.0f;
	elem[ 3 ] = m11 + m22 + m33 + 1.0f;

	unsigned biggestIndex = 0;
	for ( int i = 1; i < 4; i++ ) {
		if ( elem[i] > elem[biggestIndex] )
			biggestIndex = i;
	}

	if ( elem[biggestIndex] < 0.0f )
		return -1;

	float *q[4] = {&qx, &qy, &qz, &qw};
	float v = sqrtf( elem[biggestIndex] ) * 0.5f;
	*q[biggestIndex] = v;
	float mult = 0.25f / v;

	switch ( biggestIndex ) {
	case 0: // x
		*q[1] = (m12 + m21) * mult;
		*q[2] = (m31 + m13) * mult;
		*q[3] = (m23 - m32) * mult;
		break;
	case 1: // y
		*q[0] = (m12 + m21) * mult;
		*q[2] = (m23 + m32) * mult;
		*q[3] = (m31 - m13) * mult;
		break;
	case 2: // z
		*q[0] = (m31 + m13) * mult;
		*q[1] = (m23 + m32) * mult;
		*q[3] = (m12 - m21) * mult;
		break;
	case 3: // w
		*q[0] = (m23 - m32) * mult;
		*q[1] = (m31 - m13) * mult;
		*q[2] = (m12 - m21) * mult;
		break;
	}
	quat[0] = qx;
	quat[1] = qy;
	quat[2] = qz;
	quat[3] = qw;
	return 0;
}

int mat_quat(double *mat, double *quat, int matsize)
{
	double m11,m12,m13, m21,m22,m23, m31,m32,m33, qx,qy,qz,qw;
	if(matsize == 4){
		m11 = mat[0];	m12 = mat[1];	m13 = mat[2];
		m21 = mat[4];	m22 = mat[5];	m23 = mat[6];
		m31 = mat[8];	m32 = mat[9];	m33 = mat[10];
	} else if(matsize == 3){
		m11 = mat[0];	m12 = mat[1];	m13 = mat[2];
		m21 = mat[3];	m22 = mat[4];	m23 = mat[5];
		m31 = mat[6];	m32 = mat[7];	m33 = mat[8];
	} else{
		printf("invalid input @mat_quat(), matsize = %d\n", matsize);
	}

	double elem[ 4 ]; // 0:x, 1:y, 2:z, 3:w
	elem[ 0 ] = m11 - m22 - m33 + 1.0f;
	elem[ 1 ] = -m11 + m22 - m33 + 1.0f;
	elem[ 2 ] = -m11 - m22 + m33 + 1.0f;
	elem[ 3 ] = m11 + m22 + m33 + 1.0f;

	unsigned biggestIndex = 0;
	for ( int i = 1; i < 4; i++ ) {
		if ( elem[i] > elem[biggestIndex] )
			biggestIndex = i;
	}

	if ( elem[biggestIndex] < 0.0f )
		return -1;

	double *q[4] = {&qx, &qy, &qz, &qw};
	double v = sqrt( elem[biggestIndex] ) * 0.5;
	*q[biggestIndex] = v;
	double mult = 0.25f / v;

	switch ( biggestIndex ) {
	case 0: // x
		*q[1] = (m12 + m21) * mult;
		*q[2] = (m31 + m13) * mult;
		*q[3] = (m23 - m32) * mult;
		break;
	case 1: // y
		*q[0] = (m12 + m21) * mult;
		*q[2] = (m23 + m32) * mult;
		*q[3] = (m31 - m13) * mult;
		break;
	case 2: // z
		*q[0] = (m31 + m13) * mult;
		*q[1] = (m23 + m32) * mult;
		*q[3] = (m12 - m21) * mult;
		break;
	case 3: // w
		*q[0] = (m23 - m32) * mult;
		*q[1] = (m31 - m13) * mult;
		*q[2] = (m12 - m21) * mult;
		break;
	}
	quat[0] = qx;
	quat[1] = qy;
	quat[2] = qz;
	quat[3] = qw;
	return 0;
}

void quat_mat(float *quat, float *mat, int matsize)
{
	float r11,r12,r13, r21,r22,r23, r31,r32,r33, qx,qy,qz,qw;

	qx = quat[0];
	qy = quat[1];
	qz = quat[2];
	qw = quat[3];
	r11 = 1.0f - 2.0f * qy * qy - 2.0f * qz * qz;
	r12 = 2.0f * qx * qy + 2.0f * qw * qz;
	r13 = 2.0f * qx * qz - 2.0f * qw * qy;
	r21 = 2.0f * qx * qy - 2.0f * qw * qz;
	r22 = 1.0f - 2.0f * qx * qx - 2.0f * qz * qz;
	r23 = 2.0f * qy * qz + 2.0f * qw * qx;
	r31 = 2.0f * qx * qz + 2.0f * qw * qy;
	r32 = 2.0f * qy * qz - 2.0f * qw * qx;
	r33 = 1.0f - 2.0f * qx * qx - 2.0f * qy * qy;

	if(matsize == 4){
		mat[0] = r11;	mat[1] = r12;	mat[2] = r13;	mat[3] = 0;
		mat[4] = r21;	mat[5] = r22;	mat[6] = r23;	mat[7] = 0;
		mat[8] = r31;	mat[9] = r32;	mat[10] = r33;	mat[11] = 0;
		mat[12] = 0;	mat[13] = 0;	mat[14] = 0;	mat[15] = 1;
	} else if(matsize == 3){
		mat[0] = r11;	mat[1] = r12;	mat[2] = r13;
		mat[3] = r21;	mat[4] = r22;	mat[5] = r23;
		mat[6] = r31;	mat[7] = r32;	mat[8] = r33;
	} else{
		printf("invalid input @mat_quat(), matsize = %d\n", matsize);
	}
}

void quat_mat(double *quat, double *mat, int matsize)
{
	double r11,r12,r13, r21,r22,r23, r31,r32,r33, qx,qy,qz,qw;

	qx = quat[0];
	qy = quat[1];
	qz = quat[2];
	qw = quat[3];
	r11 = 1.0f - 2.0f * qy * qy - 2.0f * qz * qz;
	r12 = 2.0f * qx * qy + 2.0f * qw * qz;
	r13 = 2.0f * qx * qz - 2.0f * qw * qy;
	r21 = 2.0f * qx * qy - 2.0f * qw * qz;
	r22 = 1.0f - 2.0f * qx * qx - 2.0f * qz * qz;
	r23 = 2.0f * qy * qz + 2.0f * qw * qx;
	r31 = 2.0f * qx * qz + 2.0f * qw * qy;
	r32 = 2.0f * qy * qz - 2.0f * qw * qx;
	r33 = 1.0f - 2.0f * qx * qx - 2.0f * qy * qy;

	if(matsize == 4){
		mat[0] = r11;	mat[1] = r12;	mat[2] = r13;	mat[3] = 0;
		mat[4] = r21;	mat[5] = r22;	mat[6] = r23;	mat[7] = 0;
		mat[8] = r31;	mat[9] = r32;	mat[10] = r33;	mat[11] = 0;
		mat[12] = 0;	mat[13] = 0;	mat[14] = 0;	mat[15] = 1;
	} else if(matsize == 3){
		mat[0] = r11;	mat[1] = r12;	mat[2] = r13;
		mat[3] = r21;	mat[4] = r22;	mat[5] = r23;
		mat[6] = r31;	mat[7] = r32;	mat[8] = r33;
	} else{
		printf("invalid input @mat_quat(), matsize = %d\n", matsize);
	}
}

double normalize_v2( double *x, double *y )
{
	double len = sqrt( (*x)*(*x) + (*y)*(*y) );
	*x = (*x)/len;
	*y = (*y)/len;
	return len;
}

double normalize_v3( double *x, double *y, double *z )
{
	double len = sqrt( (*x)*(*x) + (*y)*(*y) + (*z)*(*z) );
	*x = (*x)/len;
	*y = (*y)/len;
	*z = (*z)/len;
	return len;
}

int getMat(int vcnt,
		   double *vx0, double *vy0, double *vz0,
		   double *vx1, double *vy1, double *vz1,
		   double *_mat)
{
	int i, cnt, ret=0;
	double mat[16], rm[16], tm0[16], tm1[16], eps = 0.0000001;
	double rm00, rm01, rm02, rm10, rm11, rm12, rm20, rm21, rm22;
	double tv00, tv01, tv02, tv10, tv11, tv12;
	double w, sw, cw, p, sp, cp, k, ck, sk, Fx, Fy, Fz;
	double dfdp00, dfdp01, dfdp02, dfdp10, dfdp11, dfdp12, dfdp20, dfdp21, dfdp22;
	double m[9], b0,b1,b2, dw,dp,dk;

	if(vcnt<3){
		ret = -1;
		goto end;
	}

	// translation
	tv00 = tv01 = tv02 = tv10 = tv11 = tv12 = 0;
	for( i=0; i<vcnt; i++ ){
		tv00+=vx0[i];	tv01+=vy0[i];	tv02+=vz0[i];
		tv10+=vx1[i];	tv11+=vy1[i];	tv12+=vz1[i];
	}
	tv00/=vcnt;	tv01/=vcnt;	tv02/=vcnt;
	tv10/=vcnt;	tv11/=vcnt;	tv12/=vcnt;
	for( i=0; i<vcnt; i++ ){
		vx0[i]-=tv00;	vy0[i]-=tv01;	vz0[i]-=tv02;
		vx1[i]-=tv10;	vy1[i]-=tv11;	vz1[i]-=tv12;
	}

	// initial value for rotation matrix
	w=0;	cw=cos(w);	sw=sin(w);
	p=0;	cp=cos(p);	sp=sin(p);
	k=0;	ck=cos(k);	sk=sin(k);
	cnt = 0;
	dw = dp = dk = DBL_MAX;
	while( (fabs(dw)>eps || fabs(dp)>eps || fabs(dk)>eps) && cnt<100 ){
		memset( m, 0, sizeof(double)*9 );
		b0 = b1 = b2 = 0;
		for( i=0; i<vcnt; i++){
			rm00 = cp*ck;			rm01 = -cp*sk;			rm02 = sp;
			rm10 = cw*sk+sw*sp*ck;	rm11 = cw*ck-sw*sp*sk;	rm12 = -sw*cp;
			rm20 = sw*sk-cw*sp*ck;	rm21 = sw*ck+cw*sp*sk;	rm22 = cw*cp;

			Fx = rm00*vx0[i] + rm01*vy0[i] + rm02*vz0[i] - vx1[i];
			Fy = rm10*vx0[i] + rm11*vy0[i] + rm12*vz0[i] - vy1[i];
			Fz = rm20*vx0[i] + rm21*vy0[i] + rm22*vz0[i] - vz1[i];

			dfdp00 = 0;
			dfdp01 = (-sp*ck)			*vx0[i] + (sp*sk)			*vy0[i]	+ (cp)		*vz0[i];
			dfdp02 = (-cp*sk)			*vx0[i]	+ (-cp*ck)			*vy0[i];
			dfdp10 = (-sw*sk+cw*sp*ck)	*vx0[i]	+ (-sw*ck-cw*sp*sk)	*vy0[i]	+ (-cw*cp)	*vz0[i];
			dfdp11 = (sw*cp*ck)			*vx0[i]	+ (-sw*cp*sk)		*vy0[i]	+ (sw*sp)	*vz0[i];
			dfdp12 = (cw*ck-sw*sp*sk)	*vx0[i]	+ (-cw*sk-sw*sp*ck)	*vy0[i];
			dfdp20 = (cw*sk+sw*sp*ck)	*vx0[i]	+ (cw*ck-sw*sp*sk)	*vy0[i]	+ (-sw*cp)	*vz0[i];
			dfdp21 = (-cw*cp*ck)		*vx0[i]	+ (cw*cp*sk)		*vy0[i]	+ (-cw*sp)	*vz0[i];
			dfdp22 = (sw*ck+cw*sp*sk)	*vx0[i]	+ (-sw*sk+cw*sp*ck)	*vy0[i];

			m[0] += dfdp00 * dfdp00 + dfdp10 * dfdp10 + dfdp20 * dfdp20;
			m[1] += dfdp00 * dfdp01 + dfdp10 * dfdp11 + dfdp20 * dfdp21;
			m[2] += dfdp00 * dfdp02 + dfdp10 * dfdp12 + dfdp20 * dfdp22;
			m[3] += dfdp01 * dfdp00 + dfdp11 * dfdp10 + dfdp21 * dfdp20;
			m[4] += dfdp01 * dfdp01 + dfdp11 * dfdp11 + dfdp21 * dfdp21;
			m[5] += dfdp01 * dfdp02 + dfdp11 * dfdp12 + dfdp21 * dfdp22;
			m[6] += dfdp02 * dfdp00 + dfdp12 * dfdp10 + dfdp22 * dfdp20;
			m[7] += dfdp02 * dfdp01 + dfdp12 * dfdp11 + dfdp22 * dfdp21;
			m[8] += dfdp02 * dfdp02 + dfdp12 * dfdp12 + dfdp22 * dfdp22;
			b0 += dfdp00 * Fx + dfdp10 * Fy + dfdp20 * Fz;
			b1 += dfdp01 * Fx + dfdp11 * Fy + dfdp21 * Fz;
			b2 += dfdp02 * Fx + dfdp12 * Fy + dfdp22 * Fz;
		}
		inv_m33(m);
		dw = m[0]*b0 + m[1]*b1 + m[2]*b2;
		dp = m[3]*b0 + m[4]*b1 + m[5]*b2;
		dk = m[6]*b0 + m[7]*b1 + m[8]*b2;
		w=w-dw;	sw=sin(w);	cw=cos(w);
		p=p-dp;	sp=sin(p);	cp=cos(p);
		k=k-dk;	sk=sin(k);	ck=cos(k);

		cnt++;
	}
	rm00 = cp*ck;			rm01 = -cp*sk;			rm02 = sp;
	rm10 = cw*sk+sw*sp*ck;	rm11 = cw*ck-sw*sp*sk;	rm12 = -sw*cp;
	rm20 = sw*sk-cw*sp*ck;	rm21 = sw*ck+cw*sp*sk;	rm22 = cw*cp;

	unit_m44(tm0);	tm0[3] = -tv00;	tm0[7] =-tv01;	tm0[11] = -tv02;
	unit_m44(tm1);	tm1[3] = tv10;	tm1[7] = tv11;	tm1[11] = tv12;
	rm[0] = rm00;	rm[1] = rm01;	rm[2] = rm02;	rm[3] = 0;
	rm[4] = rm10;	rm[5] = rm11;	rm[6] = rm12;	rm[7] = 0;
	rm[8] = rm20;	rm[9] = rm21;	rm[10] = rm22;	rm[11] = 0;
	rm[12] = 0;		rm[13] = 0;		rm[14] = 0;		rm[15] = 1;

	memcpy(mat,tm1,sizeof(double)*16);
	mult_m44_n44(mat, rm, 1);	// dst=1 : n44[16] = n44[16] * m44[16] 
	mult_m44_n44(mat, tm0, 1);
	if(_mat) memcpy(_mat,mat,sizeof(double)*16);
end:
	return ret;
}

// return -1: no intersection
int intersectionOfLine( double sx0, double sy0, double ex0, double ey0,
					   double sx1, double sy1, double ex1, double ey1,
					   double *ix, double *iy, double *l0, double *l1 )
{
	int ret=0;
	double dx0=ex0-sx0, dy0=ey0-sy0, dx1=ex1-sx1, dy1=ey1-sy1;
	double len0 = sqrt(dx0*dx0+dy0*dy0), len1 = sqrt(dx1*dx1+dy1*dy1);
	double d = dx0*dy1-dy0*dx1;
	if (d == 0){
		ret=-1; goto end;
	}
	double u = ((sx1-sx0)*dy1 - (sy1-sy0)*dx1) /d;
	double v = ((sx1-sx0)*dy0 - (sy1-sy0)*dx0) /d;
	if (u < 0.0 || u > 1.0){
		ret=-1; goto end; // intersection point not between p1 and p2
	}
	if (v < 0.0 || v > 1.0){
		ret=-1; goto end; // intersection point not between p3 and p4
	}
	if(ix){ *ix = sx0 + u * dx0; }
	if(iy){ *iy = sy0 + u * dy0; }
	if(l0){ *l0 = u * len0; }
	if(l1){ *l1 = v * len1; }
end:
	return ret;
}
