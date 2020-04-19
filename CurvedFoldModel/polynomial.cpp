// source:
// https://www.sist.ac.jp/~suganuma/kougi/other_lecture/SE/predict/program/multi/least/cpp.htm

/****************************/
/* �ŏ����@�i�������ߎ��j */
/*      coded by Y.Suganuma */
/****************************/

#include <stdio.h>
#include <math.h>
#include "polynomial.h"

#if 0
int main()
{
	double* x, * y, * z;
	int i1, m, n;

	scanf("%d %d", &m, &n);   // �������̎����ƃf�[�^�̐�

	x = new double[n];
	y = new double[n];

	for (i1 = 0; i1 < n; i1++)   // �f�[�^
		scanf("%lf %lf", &x[i1], &y[i1]);

	z = least(m, n, x, y);

	if (z != NULL) {
		printf("����\n");
		for (i1 = 0; i1 < m + 1; i1++)
			printf("   %d ���̌W�� %f\n", m - i1, z[i1]);
		delete[] z;
	}
	else
		printf("***error  �t�s������߂邱�Ƃ��ł��܂���ł���\n");

	delete[] x;
	delete[] y;

	return 0;
}
#endif

/******************************************/
/* �ŏ��Q��@                             */
/*      m : �������̎���                  */
/*      n : �f�[�^�̐�                    */
/*      x,y : �f�[�^                      */
/*      return : �������̌W���i��������j */
/*               �G���[�̏ꍇ��NULL��Ԃ� */
/******************************************/
double* least(int m, int n, double* x, double* y)
{
	double** A, ** w, * z, x1, x2;
	int i1, i2, i3, sw;

	m++;
	z = new double[m];
	w = new double* [m];
	for (i1 = 0; i1 < m; i1++)
		w[i1] = new double[m + 1];
	A = new double* [n];

	for (i1 = 0; i1 < n; i1++) {
		A[i1] = new double[m];
		A[i1][m - 2] = x[i1];
		A[i1][m - 1] = 1.0;
		x1 = A[i1][m - 2];
		x2 = x1;
		for (i2 = m - 3; i2 >= 0; i2--) {
			x2 *= x1;
			A[i1][i2] = x2;
		}
	}

	for (i1 = 0; i1 < m; i1++) {
		for (i2 = 0; i2 < m; i2++) {
			w[i1][i2] = 0.0;
			for (i3 = 0; i3 < n; i3++)
				w[i1][i2] += A[i3][i1] * A[i3][i2];
		}
	}

	for (i1 = 0; i1 < m; i1++) {
		w[i1][m] = 0.0;
		for (i2 = 0; i2 < n; i2++)
			w[i1][m] += A[i2][i1] * y[i2];
	}

	sw = gauss(w, m, 1, 1.0e-10);

	if (sw == 0) {
		for (i1 = 0; i1 < m; i1++)
			z[i1] = w[i1][m];
	}
	else
		z = NULL;

	for (i1 = 0; i1 < n; i1++)
		delete[] A[i1];
	for (i1 = 0; i1 < m; i1++)
		delete[] w[i1];
	delete[] A;
	delete[] w;

	return z;
}

/*******************************************************/
/* ���`�A���������������i�t�s������߂�j              */
/*      w : �������̍��Ӌy�щE��                       */
/*      n : �������̐�                                 */
/*      m : �������̉E�ӂ̗�̐�                       */
/*      eps : �������𔻒肷��K��                     */
/*      return : =0 : ����                             */
/*               =1 : �t�s�񂪑��݂��Ȃ�               */
/*******************************************************/
int gauss(double** w, int n, int m, double eps)
{
	double y1, y2;
	int ind = 0, nm, m1, m2, i1, i2, i3;

	nm = n + m;

	for (i1 = 0; i1 < n && ind == 0; i1++) {

		y1 = .0;
		m1 = i1 + 1;
		m2 = 0;

		for (i2 = i1; i2 < n; i2++) {
			y2 = fabs(w[i2][i1]);
			if (y1 < y2) {
				y1 = y2;
				m2 = i2;
			}
		}

		if (y1 < eps)
			ind = 1;

		else {

			for (i2 = i1; i2 < nm; i2++) {
				y1 = w[i1][i2];
				w[i1][i2] = w[m2][i2];
				w[m2][i2] = y1;
			}

			y1 = 1.0 / w[i1][i1];

			for (i2 = m1; i2 < nm; i2++)
				w[i1][i2] *= y1;

			for (i2 = 0; i2 < n; i2++) {
				if (i2 != i1) {
					for (i3 = m1; i3 < nm; i3++)
						w[i2][i3] -= w[i2][i1] * w[i1][i3];
				}
			}
		}
	}

	return(ind);
}
