#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

#include "crease.h"
#include "Bezier.h"

void crease::setDefault( rectify_param *rp )
{
	rp->flg_src_s1e1 = 0;	// 0:k2d, 1:kv*cos(alpha)
	rp->kvthres = MIN_CURV;	// �ȗ�0�Ƃ݂Ȃ�臒l
	rp->s1mgn = 2;			// �ȗ�0��Ԃ��L���镝 >=0
	rp->s2mgn = -1;			// ��ԋ��  -1: �ɒl�܂�, 0: ��ԂȂ�, >0: ��ԕ�
	rp->method = 2;			// ��ԕ��@  1: linear, 2:Bezier
	rp->src = 0;			// 0: ���ׂĎZ�o  1: s1e1�͊���, 2: s2e2�͊���
	rp->scnt = 0;
	memset( rp->s1, 0, sizeof(int)*RECT_ASIZE );
	memset( rp->e1, 0, sizeof(int)*RECT_ASIZE );
	memset( rp->s2, 0, sizeof(int)*RECT_ASIZE );
	memset( rp->e2, 0, sizeof(int)*RECT_ASIZE );
	for( int i=0; i<RECT_ASIZE; i++ ){
		rp->svlen0[i] = 10.0;
		rp->svlen1[i] = 10.0;
		rp->evlen0[i] = 10.0;
		rp->evlen1[i] = 10.0;
	}
}

int crease::gets1e1( int Xcnt, int Xs, int Xe, double *val, double thres_val, int mgn,
					int *_s1, int *_e1, int semax )
{
	int i, j, s0=-1,e0=-1, scnt=0, scnt0=0;
	int flg[MAX_SPCNT]; memset( flg, 0, sizeof(int)*MAX_SPCNT );

	// curv=0 �ɋ߂��͈͂����o -> s0,e0
	for( i=Xs; i<=Xe; i++ ){
		double aval = fabs(val[i]);
		if( s0>-1 ){
			flg[i] = 1;
		}
		if( s0<0 && aval<thres_val ){
			_s1[scnt] = i-1;
			flg[i-1] = flg[i] = 1;
			s0 = 1;
		} else if( s0>0 && aval>thres_val ){
			_e1[scnt] = i;
			scnt++;
			s0 = -1;
			if( scnt >= semax){
				printf("scnt overflow");
				scnt = -1; goto end;
			}
		}
	}
	if( s0>0 ){
		_e1[scnt] = i;
		scnt++;
		s0 = -1;
	}

	// + �� -�A- �� + ��ǉ��Ō��o�i���z���傫���ꍇ�j -> s0,e0
	for( i=Xs; i<Xe; i++ ){
		if( flg[i] || flg[i+1] ){
			continue;
		}
		if( val[i] < 0 && val[i+1] > 0 || val[i] > 0 && val[i+1] < 0 ){

			if( scnt==0 ){
				_s1[scnt] = i;
				_e1[scnt] = i+1;
				scnt++;
				continue;
			}

			// insert
			scnt++;
			if( scnt > semax){
				printf("scnt overflow");
				scnt = -1; goto end;
			}
			int ins=scnt-1;
			for( j=0; j<scnt-1; j++ ){
				if( i < _s1[j] ){
					ins = j;
					break;
				}
			}
			for( j=scnt-2; j>=ins; j-- ){
				_s1[j+1] = _s1[j];
				_e1[j+1] = _e1[j];
			}
			_s1[ins] = i;
			_e1[ins] = i+1;
		}
	}
	if( scnt==0 ){
		goto end;
	}

	//	for( i=0; i<scnt; i++ ){
	//		printf( "%d: s0=%d, e0=%d\n", i, _s1[i], _e1[i] );
	//	}

	// curv=0 �O��̕␳�͈͂�ݒ� -> s1,e1
	for( i=0; i<scnt; i++ ){
		_s1[i] = _s1[i]-mgn;
		if( _s1[i]<0 ){ _s1[i]=0; }
		_e1[i] = _e1[i]+mgn;
		if( _e1[i]>=Xcnt ){ _e1[i]=Xcnt-1; }
	}

	// �ׂƂ��������ꍇ�͓���
	for( i=1; i<scnt; i++ ){
		if( _s1[i] < _e1[i-1] ){
			_s1[i] = _s1[i-1];
			_s1[i-1] = _e1[i-1] = -1;
		}
	}
	for( i=0; i<scnt; i++ ){
		if( _s1[i] > -1 && _e1[i] > -1 ){
			_s1[scnt0] = _s1[i];
			_e1[scnt0] = _e1[i];
			scnt0++;
		}
	}
	scnt = scnt0;

	//	for( i=0; i<scnt; i++ ){
	//		printf( "%d: s1=%d, e1=%d\n", i, _s1[i], _e1[i] );
	//	}
end:
	return scnt;
}

// abs(k2d[i])>min_kv �͈̔́A�������ς��_�A��mgn�����L�����͈͂����߂�
int crease::gets1e1( int flg_src_s1e1, double thres, int mgn, int *_s1, int *_e1, int semax )
{
	if( flg_src_s1e1==0 ){
		return gets1e1( Xcnt, CEMGN, Xcnt-CEMGN-1, k2d, thres, mgn, _s1, _e1, semax );
	} else {
		double ak2d[MAX_SPCNT] = {0.0};
		for( int i=CEMGN; i<Xcnt-CEMGN; i++ ){
			ak2d[i] = kv[i]*cos(alpha[i]);
		}
		return gets1e1( Xcnt, CEMGN, Xcnt-CEMGN-1, ak2d, thres, mgn, _s1, _e1, semax );
	}
}

int crease::gets2e2( int Xcnt, int *s1, int *e1, int scnt, double *value, int *s2, int *e2 )
{
	int m[10], ret=0;
	for( int i=0; i<scnt-1; i++ ){
		m[i] = (int)( ( e1[i] + s1[i+1] ) /2.0 );
	}
	for( int i=0; i<scnt; i++ ){
		// �܂��ő�܂ōL����
		if( i==0 ){
			s2[i] = CEMGN;
		} else {
			s2[i] = m[i-1]+1;
		}
		if( i==scnt-1 ){
			e2[i] = Xcnt-CEMGN-1;
		} else {
			e2[i] = m[i]-1;
		}
		// �ł��߂��ɑ�^�ɏ��܂ŋ��߂�
		double sgn = value[ s1[i]-1 ] - value[ s1[i] ];
		for( int j=s1[i]; j>s2[i]; j-- ){
			if( sgn>=0 && value[j-1]<value[j] || sgn<0 && value[j-1]>value[j] ){
				s2[i] = j;
				break;
			}
		}
		sgn = value[ s1[i]+1 ] - value[ s1[i] ];
		for( int j=e1[i]; j<e2[i]; j++ ){
			if( sgn>=0 && value[j+1]<value[j] || sgn<0 && value[j+1]>value[j] ){
				e2[i] = j;
				break;
			}
		}
	}
	return ret;
}
