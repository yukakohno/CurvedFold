#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "curvedraw.h"
#include "util.h"

curvedraw::curvedraw()
{
	init();
}

curvedraw::~curvedraw()
{
	init();
}

void curvedraw::init()
{
	ctype = CTYPE_UNDEF;
	cvcnt = 0;
	memset( cvx, 0, sizeof(double)*MAXDCX );
	memset( cvy, 0, sizeof(double)*MAXDCX );
}

void curvedraw::closecv( double th_d /* 5.0 mm */ )
{
	if( ctype != CTYPE_TRIM ){
		return;
	}
	double th_d2 = th_d*th_d; // mm
	double sx = cvx[0];
	double sy = cvy[0];
	double ex = cvx[cvcnt-1];
	double ey = cvy[cvcnt-1];
	double dx=ex-sx, dy=ey-sy, d2=dx*dx+dy*dy;
	if( d2<th_d2 ){
		if( cvcnt<MAXDCX-1 ){
			cvx[cvcnt] = sx;
			cvy[cvcnt] = sy;
			cvcnt++;
		} else {
			cvx[cvcnt-1] = sx;
			cvy[cvcnt-1] = sy;
		}
	}
}


int *val_cidx;
int compCidx( const void *a, const void *b )
{
	int ia = *(int*)a;
	int ib = *(int*)b;
	return val_cidx[ia] - val_cidx[ib];
}

int curvedraw::IS_drawcurve( curvedraw *dc, crease *c, int psx, int psy, int pex, int pey,
							double *_ISx, double *_ISy, int ISsize, int &_IScnt, 
							cvertex_type *_Xtype, int *_Xidx, double *_Xlen, int *_Cidx )
{
	int ret=0;

	int IScnt = 0;
	double *ISx = new double[ISsize];
	double *ISy = new double[ISsize]; 
	cvertex_type *Xtype = new cvertex_type[ISsize];
	int *Xidx = new int[ISsize];
	int *Cidx = new int[ISsize];
	double *Xlen = new double[ISsize];
	int *ISodr = new int[ISsize];

	for( int i=0; i<ISsize; i++ ){
		ISx[i] = ISy[i] = 0.0;
		Xtype[i] = CVTYPE_UNDEF;
		Xidx[i] = Cidx[i] = -1;
		Xlen[i] = 0.0;
	}

	//
	// intersection of curvedraw and rulings
	//
	for( int rl=0; rl<2; rl++ )
	{
		double *rx=NULL, *ry=NULL;
		switch( rl ){
			case 0:	rx=c->rlx_cp;	ry=c->rly_cp;	break;	// ¶‘¤ ruling|‹ÈüŒð“_
			case 1:	rx=c->rrx_cp;	ry=c->rry_cp;	break;	// ‰E‘¤ ruling|‹ÈüŒð“_
		}
		for( int i=c->Xsidx; i<=c->Xeidx; i++ ){
			int intersect=0;
			double maxlen=100000;
			for( int j=0; j<dc->cvcnt-1; j++ ){
				double ix,iy,l0,l1;
				int ret = intersectionOfLine( c->Xx2d[i], c->Xy2d[i],
					c->Xx2d[i]+rx[i]*maxlen, c->Xy2d[i]+ry[i]*maxlen,
					dc->cvx[j], dc->cvy[j], dc->cvx[j+1], dc->cvy[j+1],
					&ix, &iy, &l0, &l1 );
				if( ret<0 ){
					continue;
				}
				// ruling‚Æ‹Èü‚ª•¡”‰ñŒð·‚·‚éê‡‚ÍÅ‚à’Z‚¢‚à‚Ì‚ðŽg‚¤
				if( maxlen > l0 ){
					intersect = 1;
					Cidx[IScnt] = j;
					Xlen[IScnt] = maxlen = l0;
					ISx[IScnt] = ix;
					ISy[IScnt] = iy;
				}
			}
			if( intersect ){
				switch( rl ){
					case 0:	Xtype[IScnt]=CVTYPE_RUL_LEFT;	break;
					case 1:	Xtype[IScnt]=CVTYPE_RUL_RIGHT;	break;
				}
				Xidx[IScnt] = i;
				IScnt++;
			}
		}
	} // rl

	//
	// intersection of curvedraw and paper edge
	//
	for( int j=0; j<dc->cvcnt-1; j++ ){
		int ret[4];
		double ix[4],iy[4],l0[4],l1[4];
		double cvx0=dc->cvx[j], cvy0=dc->cvy[j], cvx1=dc->cvx[j+1], cvy1=dc->cvy[j+1];
		ret[0] = intersectionOfLine( psx, psy, pex, psy, cvx0, cvy0, cvx1, cvy1, &ix[0], &iy[0], &l0[0], &l1[0] ); // TOP
		ret[1] = intersectionOfLine( pex, psy, pex, pey, cvx0, cvy0, cvx1, cvy1, &ix[1], &iy[1], &l0[1], &l1[1] ); // RIGHT
		ret[2] = intersectionOfLine( pex, pey, psx, pey, cvx0, cvy0, cvx1, cvy1, &ix[2], &iy[2], &l0[2], &l1[2] ); // BOTTOM
		ret[3] = intersectionOfLine( psx, pey, psx, psy, cvx0, cvy0, cvx1, cvy1, &ix[3], &iy[3], &l0[3], &l1[3] ); // LEFT
		for( int i=0; i<4; i++ ){
			if( ret[i]<0 || l1[i]==0.0 ){
				continue;
			}
			Cidx[IScnt] = j;
			switch(i){
				case 0:	Xtype[IScnt] = CVTYPE_PAPER_TOP; break;
				case 1:	Xtype[IScnt] = CVTYPE_PAPER_RIGHT; break;
				case 2:	Xtype[IScnt] = CVTYPE_PAPER_BOTTOM; break;
				case 3:	Xtype[IScnt] = CVTYPE_PAPER_LEFT; break;
			}
			//Xidx[IScnt] = 0;
			Xlen[IScnt] = l0[i];
			ISx[IScnt] = ix[i];
			ISy[IScnt] = iy[i];
			IScnt++;
		}
	}

	//
	// intersection of curvedraw and crease
	//
	for( int i=c->Xsidx; i<c->Xeidx; i++ ){
		for( int j=0; j<dc->cvcnt-1; j++ ){
			double ix,iy,l0,l1;
			int ret = intersectionOfLine( c->Xx2d[i], c->Xy2d[i], c->Xx2d[i+1], c->Xy2d[i+1],
				dc->cvx[j], dc->cvy[j], dc->cvx[j+1], dc->cvy[j+1], &ix, &iy, &l0, &l1 );
			if( ret<0 ){
				continue;
			}
			Cidx[IScnt] = j;
			Xtype[IScnt] = CVTYPE_CREASE; // START/END undefined
			Xlen[IScnt] = l0;
			ISx[IScnt] = ix;
			ISy[IScnt] = iy;
			IScnt++;
			break;
		}
	}

	// sort ISodr[] by Cidx
	for( int i=0; i<IScnt; i++ ){
		ISodr[i] = i;
	}
	val_cidx = Cidx;
	qsort( ISodr, IScnt, sizeof(int), compCidx );

	_IScnt = IScnt;
	for( int i=0; i<IScnt; i++ ){
		int j = ISodr[i];
		_ISx[i] = ISx[j];
		_ISy[i] = ISy[j];
		_Xtype[i] = Xtype[j];
		_Xidx[i] = Xidx[j];
		_Xlen[i] = Xlen[j];
		_Cidx[i] = Cidx[j];
		//printf("%d ", Cidx[j]); // check
	}
#if 0
	{
		FILE *fp = fopen( "output/curve.csv", "w" );
		if( fp ){
			for( int i=0; i<IScnt; i++ ){
				fprintf( fp, "%f,%f,%d,%d,%f,%d\n", _ISx[i],_ISy[i],_Xtype[i],_Xidx[i],_Xlen[i], _Cidx[i] );
			}
			fprintf( fp, "\n" );
			for( int i=c->Xsidx; i<c->Xeidx+1; i++ ){
				fprintf( fp, "%f,%f\n", c->Xx2d[i], c->Xy2d[i] );
			}
			fclose(fp); fp=NULL;
		}
	}
#endif

end:
	if( ISx ){ delete[] ISx; }
	if( ISy ){ delete[] ISy; }
	if( Xtype ){ delete[] Xtype; }
	if( Xidx ){ delete[] Xidx; }
	if( Cidx ){ delete[] Cidx; }
	if( Xlen ){ delete[] Xlen; }
	if( ISodr ){ delete[] ISodr; }
	return ret;
}

int curvedraw::IS_check_cw( curvedraw *dc, crease *c, int psx, int psy, int pex, int pey,
						   double *ISx, double *ISy, int ISsize, int &IScnt, 
						   cvertex_type *Xtype, int *Xidx, double *Xlen, int *Cidx )
{
	int i, cw=1;
	for( i=0; i<IScnt-1; i++ ){
		if( Xtype[i]==CVTYPE_CREASE_START ){
			double vx0,vy0,vx1,vy1;
			int xi = Xidx[i];
			vx0 = c->Xx2d[xi] - c->Xx2d[xi+1];
			vy0 = c->Xy2d[xi] - c->Xy2d[xi+1];
			vx1 = ISx[i+1] - ISx[i];
			vy1 = ISy[i+1] - ISy[i];
			if( vx0*vy1-vx1*vy0 < 0 ){ // ¶‘¤
				//inverse_c(k);
				cw=-1;
			}
			break;
		}
		if( Xtype[i]==CVTYPE_CREASE_END ){
			double vx0,vy0,vx1,vy1;
			int xi = Xidx[i];
			vx0 = c->Xx2d[xi] - c->Xx2d[xi-1];
			vy0 = c->Xy2d[xi] - c->Xy2d[xi-1];
			vx1 = ISx[i+1] - ISx[i];
			vy1 = ISy[i+1] - ISy[i];
			if( vx0*vy1-vx1*vy0 < 0 ){
				//inverse_c(k);
				cw=-1;
			}
			break;
		}
	}
	if( i == IScnt-1 ){ // break‚µ‚Ä‚¢‚È‚¢Ü‚èü‚ÆŒð‚í‚ç‚È‚¢
		for( i=0; i<IScnt-1; i++ ){
			if( Xtype[i]==CVTYPE_RUL_LEFT || Xtype[i]==CVTYPE_RUL_RIGHT ){
				double vx0,vy0,vx1,vy1;
				vx0 = ISx[i] - c->Xx2d[Xidx[i]];
				vy0 = ISy[i] - c->Xy2d[Xidx[i]];
				vx1 = ISx[i+1] - ISx[i];
				vy1 = ISy[i+1] - ISy[i];
				if( vx0*vy1-vx1*vy0 < 0 ){
					//inverse_c(k);
					cw=-1;
				}
				break;
			}
		}
	}

	return cw;
}
