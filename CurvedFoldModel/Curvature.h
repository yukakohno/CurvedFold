/*************************************************************
	Curvature.h
	PCA�ɂ��厲�����̐ݒ�
	2���Ȗ�a*x^2+b*y^2+h*x*y=0�ɋߎ�
	��NPAL���_�v���O�������痬�p
*************************************************************/

#ifndef CURVATURE_H
#define CURVATURE_H

//#include "PC3D.h"

typedef struct{
	double m0[16];		// �ŗL�x�N�g���{�ʒu�i��]�O�j
	double ev[3];		// �ŗL�l
	double a,b,h;		// �񎟋Ȗʂ̌W���Cz = a*x*x + b*y*y + 2*h*x*y
	double sint,cost;	// pivot�p�x
	double m[16];		// �ŗL�x�N�g���{�ʒu
	double A,B;			// �񎟋Ȗʂ̌W���CZ = A*X*X + B*Y*Y
	double err;			// �񎟋Ȗʋߎ��덷 090514�ǉ�
	double rho0,rho1;	// �ȗ�
	double shapeindex;	
} CurvInfo;

/*		�R�����̓_�Q �� �听������
����	double cen[3];			// ���S���W
		double x[], y[], z[];	// �_�Q���W
		int n;					// �_��
�o��	double eval[3];			// �ŗL�l
		double evec[4*4];		// �ŗL�x�N�g���C���S���W���܂ޕϊ��s��
�Ԃ�l	0:����I���C1:�ُ�I��
*/
extern int PCA3(double *x, double *y, double *z, int n, double *_cen,
				double *eval, double *evec);


/* �񎟋Ȗ� a*x^2 + b*y^2 + 2h*x*y - z = 0 �ɋߎ�����
����	double cen[3];			// ���S���W
		double x[], y[], z[];	// �_�Q���W
		int n;					// �_��
�o��	double a,b,h;			// �񎟋ȖʌW��
�Ԃ�l	0:����I���C1:�ُ�I��
*/
extern int QuadricSurface(double *px, double *py, double *pz, int n,
						  double *a, double *b, double *h);

extern double calcError(double *px, double *py, double *pz, int n, double a, double b, double h);

extern int Curvature(double *cen, double *nrm, double *_px, double *_py, double *_pz, int n,
			  CurvInfo *cinfo);

extern void dumpCurvInfo(CurvInfo *info);
extern void dumpCurvInfo(FILE *fp, CurvInfo *info);
extern void readCurvInfo(FILE *fp, CurvInfo *info);


extern int PCA2(double *x, double *y, int n, double *cen, double *eval, double *evec);


#endif // CURVATURE_H
