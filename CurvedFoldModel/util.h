#pragma once

// return value 0: end, 1: continue
int csvread(char **_in,			// input line 
			char *ret,			// item
			char sep = ',',		// separator charactor
			int flg_con = 0);	// 1: ignore consecutive separator

extern void unit_m44(float *m44);
extern void unit_m44(double *m44);
extern void unit_m33(float *m33);
extern void unit_m33(double *m33);

// v4[4] = m44[16] * v4[4]
extern void mult_m44_v4(float *_m44, float *_v4);
extern void mult_m44_v4(double *_m44, double *_v4);

// v4[4] = v4[4] * m44[16]
extern void mult_v4_m44(float *_v4, float *_m44);
extern void mult_v4_m44(double *_v4, double *_m44);

// dst=0 : m44[16] = n44[16] * m44[16]
// dst=1 : n44[16] = n44[16] * m44[16]
extern void mult_m44_n44(float *_n44, float *_m44, int dst=0);
extern void mult_m44_n44(double *_n44, double *_m44, int dst=0);
extern void mult_m33_n33(float *_n33, float *_m33, int dst=0);
extern void mult_m33_n33(double *_n33, double *_m33, int dst=0);

// dst=0 : m44[16] = n44[16] * m44[16]
// dst=1 : n44[16] = n44[16] * m44[16] 
extern void mult_m44_n44_rot(float *_n44, float *_m44, int dst=0);
extern void mult_m44_n44_rot(double *_n44, double *_m44, int dst=0);

// return value  0:OK, -1:NG
extern int inv_m44(float *m);
extern int inv_m44(double *m);
extern int inv_m33(float *m);
extern int inv_m33(double *m);
extern int m_inverse( double *m, int l );

// dst=0 : q2[4] = q1[4] * q2[4]
// dst=1 : q1[4] = q1[4] * q2[4]
extern void mult_quat(float *q1, float *q2, int dst = 0);

// rotmat[9] or [16] -> quat[4], matsize: 3 or 4
// reference URL: http://marupeke296.com/DXG_No58_RotQuaternionTrans.html
extern int mat_quat( float *mat, float *quat, int matsize=4 );
extern int mat_quat( double *mat, double *quat, int matsize=4 );

// quat[4] -> rotmat[9] or [16], matsize: 3 or 4
// reference URL: http://marupeke296.com/DXG_No58_RotQuaternionTrans.html
extern void quat_mat( float *quat, float *mat, int matsize=4 );
extern void quat_mat( double *quat, double *mat, int matsize=4 );

extern double normalize_v2( double *x, double *y );
extern double normalize_v3( double *x, double *y, double *z );

extern int getMat( int vcnt,
		   double *vx0, double *vy0, double *vz0,
		   double *vx1, double *vy1, double *vz1,
				  double *_mat );

// return -1: no intersection
extern int intersectionOfLine( double sx0, double sy0, double ex0, double ey0,
					   double sx1, double sy1, double ex1, double ey1,
					   double *ix, double *iy, double *l0, double *l1 );


