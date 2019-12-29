#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "papermodel.h"
#include "Bspline.h"
#include "util.h"

#define MIN(x,y) x<y?x:y

int papermodel::addPly( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp,
					   double x0_3d, double y0_3d, double z0_3d, double x1_3d, double y1_3d, double z1_3d,
					   double x2_3d, double y2_3d, double z2_3d )
{
	int ret=0;
	double x2[4], y2[4], z2[4];
	double x3[4], y3[4], z3[4];
	double mat[16];

	if( plcnt >= MAX_PL_FCNT-1 ){
		return -1;
	}

	plvcnt[plcnt] = 3;
	x2[0] = plx_cp[plcnt*4  ] = x0_cp;	y2[0] = ply_cp[plcnt*4  ] = y0_cp;
	x2[1] = plx_cp[plcnt*4+1] = x1_cp;	y2[1] = ply_cp[plcnt*4+1] = y1_cp;
	x2[2] = plx_cp[plcnt*4+2] = x2_cp;	y2[2] = ply_cp[plcnt*4+2] = y2_cp;
	//x2[3] = plx_cp[plcnt*4+3] = x3_cp;	y2[3] = ply_cp[plcnt*4+3] = y3_cp;
	z2[0] = z2[1] = z2[2] = z2[3] = 0.0;
	x3[0] = plx[plcnt*4  ] = x0_3d;	y3[0] = ply[plcnt*4  ] = y0_3d;	z3[0] = plz[plcnt*4  ] = z0_3d;
	x3[1] = plx[plcnt*4+1] = x1_3d;	y3[1] = ply[plcnt*4+1] = y1_3d;	z3[1] = plz[plcnt*4+1] = z1_3d;
	x3[2] = plx[plcnt*4+2] = x2_3d;	y3[2] = ply[plcnt*4+2] = y2_3d;	z3[2] = plz[plcnt*4+2] = z2_3d;
	//x3[3] = plx[plcnt*4+3] = x3_3d;	y3[3] = ply[plcnt*4+3] = y3_3d;	z3[3] = plz[plcnt*4+3] = z3_3d;
	getMat( plvcnt[plcnt], x2,y2,z2, x3,y3,z3, mat);
	plmat[plcnt*16+ 0]=mat[0]; plmat[plcnt*16+ 1]=mat[4]; plmat[plcnt*16+ 2]=mat[ 8]; plmat[plcnt*16+ 3]=mat[12];
	plmat[plcnt*16+ 4]=mat[1]; plmat[plcnt*16+ 5]=mat[5]; plmat[plcnt*16+ 6]=mat[ 9]; plmat[plcnt*16+ 7]=mat[13];
	plmat[plcnt*16+ 8]=mat[2]; plmat[plcnt*16+ 9]=mat[6]; plmat[plcnt*16+10]=mat[10]; plmat[plcnt*16+11]=mat[14];
	plmat[plcnt*16+12]=mat[3]; plmat[plcnt*16+13]=mat[7]; plmat[plcnt*16+14]=mat[11]; plmat[plcnt*16+15]=mat[15];
	plcnt++;

	return ret;
}

int papermodel::addPly( double x0_cp, double y0_cp, double x1_cp, double y1_cp,
					   double x2_cp, double y2_cp, double x3_cp, double y3_cp,
					   double x0_3d, double y0_3d, double z0_3d, double x1_3d, double y1_3d, double z1_3d,
					   double x2_3d, double y2_3d, double z2_3d, double x3_3d, double y3_3d, double z3_3d )
{
	int ret=0;
	double x2[4], y2[4], z2[4];
	double x3[4], y3[4], z3[4];
	double mat[16];

	if( plcnt >= MAX_PL_FCNT-1 ){
		return -1;
	}

	if( x3_3d==0.0 && y3_3d==0.0 && z3_3d==0.0 ){ // 3�_�����Ƃ�4�p�`�����

		ret = addPly( x0_cp, y0_cp, x1_cp, y1_cp, x2_cp, y2_cp, x0_3d, y0_3d, z0_3d, x1_3d, y1_3d, z1_3d, x2_3d, y2_3d, z2_3d );

		int pi = plcnt-1;
		plvcnt[pi] = 4;
		plx_cp[pi*4+3] = x3_cp;
		ply_cp[pi*4+3] = y3_cp;
		plx[pi*4+3] = plmat[pi*16+ 0]*x3_cp + plmat[pi*16+ 4]*y3_cp + plmat[pi*16+ 8]*0.0 + plmat[pi*16+12];	
		ply[pi*4+3] = plmat[pi*16+ 1]*x3_cp + plmat[pi*16+ 5]*y3_cp + plmat[pi*16+ 9]*0.0 + plmat[pi*16+13];	
		plz[pi*4+3] = plmat[pi*16+ 2]*x3_cp + plmat[pi*16+ 6]*y3_cp + plmat[pi*16+10]*0.0 + plmat[pi*16+14];

	} else {
		plvcnt[plcnt] = 4;
		x2[0] = plx_cp[plcnt*4  ] = x0_cp;	y2[0] = ply_cp[plcnt*4  ] = y0_cp;
		x2[1] = plx_cp[plcnt*4+1] = x1_cp;	y2[1] = ply_cp[plcnt*4+1] = y1_cp;
		x2[2] = plx_cp[plcnt*4+2] = x2_cp;	y2[2] = ply_cp[plcnt*4+2] = y2_cp;
		x2[3] = plx_cp[plcnt*4+3] = x3_cp;	y2[3] = ply_cp[plcnt*4+3] = y3_cp;
		z2[0] = z2[1] = z2[2] = z2[3] = 0.0;
		x3[0] = plx[plcnt*4  ] = x0_3d;	y3[0] = ply[plcnt*4  ] = y0_3d;	z3[0] = plz[plcnt*4  ] = z0_3d;
		x3[1] = plx[plcnt*4+1] = x1_3d;	y3[1] = ply[plcnt*4+1] = y1_3d;	z3[1] = plz[plcnt*4+1] = z1_3d;
		x3[2] = plx[plcnt*4+2] = x2_3d;	y3[2] = ply[plcnt*4+2] = y2_3d;	z3[2] = plz[plcnt*4+2] = z2_3d;
		x3[3] = plx[plcnt*4+3] = x3_3d;	y3[3] = ply[plcnt*4+3] = y3_3d;	z3[3] = plz[plcnt*4+3] = z3_3d;
		getMat( plvcnt[plcnt], x2,y2,z2, x3,y3,z3, mat);
		plmat[plcnt*16+ 0]=mat[0]; plmat[plcnt*16+ 1]=mat[4]; plmat[plcnt*16+ 2]=mat[ 8]; plmat[plcnt*16+ 3]=mat[12];
		plmat[plcnt*16+ 4]=mat[1]; plmat[plcnt*16+ 5]=mat[5]; plmat[plcnt*16+ 6]=mat[ 9]; plmat[plcnt*16+ 7]=mat[13];
		plmat[plcnt*16+ 8]=mat[2]; plmat[plcnt*16+ 9]=mat[6]; plmat[plcnt*16+10]=mat[10]; plmat[plcnt*16+11]=mat[14];
		plmat[plcnt*16+12]=mat[3]; plmat[plcnt*16+13]=mat[7]; plmat[plcnt*16+14]=mat[11]; plmat[plcnt*16+15]=mat[15];
		plcnt++;
	}

	return ret;
}

int papermodel::addPlymat( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp, double *mat )
{
	int ret=0;

	if( plcnt >= MAX_PL_FCNT-1 ){
		return -1;
	}

	plvcnt[plcnt] = 3;
	plx_cp[plcnt*4  ] = x0_cp;	ply_cp[plcnt*4  ] = y0_cp;
	plx_cp[plcnt*4+1] = x1_cp;	ply_cp[plcnt*4+1] = y1_cp;
	plx_cp[plcnt*4+2] = x2_cp;	ply_cp[plcnt*4+2] = y2_cp;

	double *m = &(plmat[plcnt*16]);
	memcpy( m, mat, sizeof(double)*16 );

	plx[plcnt*4+0] = m[0]*x0_cp + m[4]*y0_cp + m[ 8]*0.0 + m[12];	
	ply[plcnt*4+0] = m[1]*x0_cp + m[5]*y0_cp + m[ 9]*0.0 + m[13];	
	plz[plcnt*4+0] = m[2]*x0_cp + m[6]*y0_cp + m[10]*0.0 + m[14];
	plx[plcnt*4+1] = m[0]*x1_cp + m[4]*y1_cp + m[ 8]*0.1 + m[12];	
	ply[plcnt*4+1] = m[1]*x1_cp + m[5]*y1_cp + m[ 9]*0.1 + m[13];	
	plz[plcnt*4+1] = m[2]*x1_cp + m[6]*y1_cp + m[10]*0.1 + m[14];
	plx[plcnt*4+2] = m[0]*x2_cp + m[4]*y2_cp + m[ 8]*0.2 + m[12];	
	ply[plcnt*4+2] = m[1]*x2_cp + m[5]*y2_cp + m[ 9]*0.2 + m[13];	
	plz[plcnt*4+2] = m[2]*x2_cp + m[6]*y2_cp + m[10]*0.2 + m[14];

	plcnt++;
	return ret;
}

int papermodel::addEdge( double x0, double y0, double z0, double x1, double y1, double z1 )
{
	int ret=0;

	if( plecnt >= MAX_PL_ECNT-1 ){
		return -1;
	}

	plex[plecnt*2  ] = x0;
	pley[plecnt*2  ] = y0;
	plez[plecnt*2  ] = z0;
	plex[plecnt*2+1] = x1;
	pley[plecnt*2+1] = y1;
	plez[plecnt*2+1] = z1;
	plecnt++;
	return ret;
}

int papermodel::addEdge( double x0, double y0, double z0, double x1, double y1, double z1,
						double x2, double y2, double z2 )
{
	int ret=0;

	plex[plecnt*2  ] = x0;
	if( plecnt >= MAX_PL_ECNT-1 ){
		return -1;
	}

	pley[plecnt*2  ] = y0;
	plez[plecnt*2  ] = z0;
	plex[plecnt*2+1] = x1;
	pley[plecnt*2+1] = y1;
	plez[plecnt*2+1] = z1;
	plecnt++;
	plex[plecnt*2  ] = x1;
	pley[plecnt*2  ] = y1;
	plez[plecnt*2  ] = z1;
	plex[plecnt*2+1] = x2;
	pley[plecnt*2+1] = y2;
	plez[plecnt*2+1] = z2;
	plecnt++;
	return ret;
}

int papermodel::calcRulingPly()
{
	int i, li, ret=0;
	plcnt = plecnt = 0;

	crease *c = &(crs[0]);
	for( i=c->Xsidx+1, li=c->Xsidx; i<=c->Xeidx; i++ ){
		if( c->rllen[i]==0.0 ){
			continue;
		}
		addPly( c->Xx2d[ i], c->Xy2d[ i], c->Xx2d[li], c->Xy2d[li],
			c->Xx2d[li] + c->rlx_cp[li]*c->rllen[li], c->Xy2d[li] + c->rly_cp[li]*c->rllen[li],
			c->Xx2d[ i] + c->rlx_cp[ i]*c->rllen[ i], c->Xy2d[ i] + c->rly_cp[ i]*c->rllen[ i],
			c->Xx[ i], c->Xy[ i], c->Xz[ i], c->Xx[li],  c->Xy[li], c->Xz[li],
			c->Xx[li] + c->rlx[li]*c->rllen[li], c->Xy[li] + c->rly[li]*c->rllen[li], c->Xz[li] + c->rlz[li]*c->rllen[li],
			c->Xx[ i] + c->rlx[ i]*c->rllen[ i], c->Xy[ i] + c->rly[ i]*c->rllen[ i], c->Xz[ i] + c->rlz[ i]*c->rllen[ i] );
#if 1 // add corner triangle
		int pi = plcnt-1;
		if( /*c->rlflg[li] == RETYPE_EGLEFT && c->rlflg[i] == RETYPE_EGTOP ||*/
			c->rlflg[li] == RETYPE_EGTOP && c->rlflg[i] == RETYPE_EGLEFT ){
				addPlymat(
					c->Xx2d[ i] + c->rlx_cp[ i]*c->rllen[ i], c->Xy2d[ i] + c->rly_cp[ i]*c->rllen[ i],
					c->Xx2d[li] + c->rlx_cp[li]*c->rllen[li], c->Xy2d[li] + c->rly_cp[li]*c->rllen[li],
					psx,psy,&(plmat[pi*16]) );
		} else if( /*c->rlflg[li] == RETYPE_EGTOP && c->rlflg[i] == RETYPE_EGRIGHT ||*/
			c->rlflg[li] == RETYPE_EGRIGHT && c->rlflg[i] == RETYPE_EGTOP ){
				addPlymat(
					c->Xx2d[ i] + c->rlx_cp[ i]*c->rllen[ i], c->Xy2d[ i] + c->rly_cp[ i]*c->rllen[ i],
					c->Xx2d[li] + c->rlx_cp[li]*c->rllen[li], c->Xy2d[li] + c->rly_cp[li]*c->rllen[li],
					pex,psy,&(plmat[pi*16]) );
		} else if( /*c->rlflg[li] == RETYPE_EGRIGHT && c->rlflg[i] == RETYPE_EGBOTTOM ||*/
			c->rlflg[li] == RETYPE_EGBOTTOM && c->rlflg[i] == RETYPE_EGRIGHT ){
				addPlymat(
					c->Xx2d[ i] + c->rlx_cp[ i]*c->rllen[ i], c->Xy2d[ i] + c->rly_cp[ i]*c->rllen[ i],
					c->Xx2d[li] + c->rlx_cp[li]*c->rllen[li], c->Xy2d[li] + c->rly_cp[li]*c->rllen[li],
					pex,pey,&(plmat[pi*16]) );
		} else if( /*c->rlflg[li] == RETYPE_EGBOTTOM && c->rlflg[i] == RETYPE_EGLEFT ||*/
			c->rlflg[li] == RETYPE_EGLEFT && c->rlflg[i] == RETYPE_EGBOTTOM ){
				addPlymat(
					c->Xx2d[ i] + c->rlx_cp[ i]*c->rllen[ i], c->Xy2d[ i] + c->rly_cp[ i]*c->rllen[ i],
					c->Xx2d[li] + c->rlx_cp[li]*c->rllen[li], c->Xy2d[li] + c->rly_cp[li]*c->rllen[li],
					psx,pey,&(plmat[pi*16]) );
		}
#endif
#if 0
		if( c->rlflg[li] == RETYPE_EGLEFT && c->rlflg[i] == RETYPE_EGLEFT )
		{
			addEdge1( c->Xx[ i] + c->rlx[ i]*c->rllen[ i],
				c->Xy[ i] + c->rly[ i]*c->rllen[ i],
				c->Xz[ i] + c->rlz[ i]*c->rllen[ i],
				c->Xx[li] + c->rlx[li]*c->rllen[li],
				c->Xy[li] + c->rly[li]*c->rllen[li],
				c->Xz[li] + c->rlz[li]*c->rllen[li] );

		} else
#endif
			if( !(c->rlflg[li] == RETYPE_UNDEF || c->rlflg[li] == RETYPE_CURVE)
				&& !(c->rlflg[i] == RETYPE_UNDEF || c->rlflg[i] == RETYPE_CURVE) )
			{
				addEdge( c->Xx[ i] + c->rlx[ i]*c->rllen[ i],
					c->Xy[ i] + c->rly[ i]*c->rllen[ i],
					c->Xz[ i] + c->rlz[ i]*c->rllen[ i],
					c->Xx[li] + c->rlx[li]*c->rllen[li],
					c->Xy[li] + c->rly[li]*c->rllen[li],
					c->Xz[li] + c->rlz[li]*c->rllen[li] );
			}
			li=i;
	}

	for( i=c->Xsidx+1, li=c->Xsidx; i<=c->Xeidx; i++ ){
		if( c->rrlen[i]==0.0 ){
			continue;
		}
		addPly( c->Xx2d[li], c->Xy2d[li], c->Xx2d[ i], c->Xy2d[ i],
			c->Xx2d[ i] + c->rrx_cp[ i]*c->rrlen[ i], c->Xy2d[ i] + c->rry_cp[ i]*c->rrlen[ i],
			c->Xx2d[li] + c->rrx_cp[li]*c->rrlen[li], c->Xy2d[li] + c->rry_cp[li]*c->rrlen[li],
			c->Xx[li], c->Xy[li], c->Xz[li], c->Xx[ i],  c->Xy[ i], c->Xz[ i],
			c->Xx[ i] + c->rrx[ i]*c->rrlen[ i], c->Xy[ i] + c->rry[ i]*c->rrlen[ i], c->Xz[ i] + c->rrz[ i]*c->rrlen[ i],
			c->Xx[li] + c->rrx[li]*c->rrlen[li], c->Xy[li] + c->rry[li]*c->rrlen[li], c->Xz[li] + c->rrz[li]*c->rrlen[li] );
#if 1 // add corner triangle
		int pi = plcnt-1;
		if( c->rrflg[li] == RETYPE_EGLEFT && c->rrflg[i] == RETYPE_EGTOP
			/*|| c->rrflg[li] == RETYPE_EGTOP && c->rrflg[i] == RETYPE_EGLEFT*/ ){
				addPlymat(
					c->Xx2d[li] + c->rrx_cp[li]*c->rrlen[li], c->Xy2d[li] + c->rry_cp[li]*c->rrlen[li],
					c->Xx2d[ i] + c->rrx_cp[ i]*c->rrlen[ i], c->Xy2d[ i] + c->rry_cp[ i]*c->rrlen[ i],
					psx,psy,&(plmat[pi*16]) );

		} else if( c->rrflg[li] == RETYPE_EGTOP && c->rrflg[i] == RETYPE_EGRIGHT
			/*|| c->rrflg[li] == RETYPE_EGRIGHT && c->rrflg[i] == RETYPE_EGTOP*/ ){
				addPlymat(
					c->Xx2d[li] + c->rrx_cp[li]*c->rrlen[li], c->Xy2d[li] + c->rry_cp[li]*c->rrlen[li],
					c->Xx2d[ i] + c->rrx_cp[ i]*c->rrlen[ i], c->Xy2d[ i] + c->rry_cp[ i]*c->rrlen[ i],
					pex,psy,&(plmat[pi*16]) );

		} else if( c->rrflg[li] == RETYPE_EGRIGHT && c->rrflg[i] == RETYPE_EGBOTTOM
			/*|| c->rrflg[li] == RETYPE_EGBOTTOM && c->rrflg[i] == RETYPE_EGRIGHT*/ ){
				addPlymat(
					c->Xx2d[li] + c->rrx_cp[li]*c->rrlen[li], c->Xy2d[li] + c->rry_cp[li]*c->rrlen[li],
					c->Xx2d[ i] + c->rrx_cp[ i]*c->rrlen[ i], c->Xy2d[ i] + c->rry_cp[ i]*c->rrlen[ i],
					pex,pey,&(plmat[pi*16]) );

		} else if( c->rrflg[li] == RETYPE_EGBOTTOM && c->rrflg[i] == RETYPE_EGLEFT
			/*|| c->rrflg[li] == RETYPE_EGLEFT && c->rrflg[i] == RETYPE_EGBOTTOM*/ ){
				addPlymat(
					c->Xx2d[li] + c->rrx_cp[li]*c->rrlen[li], c->Xy2d[li] + c->rry_cp[li]*c->rrlen[li],
					c->Xx2d[ i] + c->rrx_cp[ i]*c->rrlen[ i], c->Xy2d[ i] + c->rry_cp[ i]*c->rrlen[ i],
					psx,pey,&(plmat[pi*16]) );
		}
#endif
#if 0
		if( c->rlflg[li] == RETYPE_EGLEFT && c->rlflg[i] == RETYPE_EGLEFT )
		{
			addEdge1( c->Xx[li] + c->rrx[li]*c->rrlen[li],
				c->Xy[li] + c->rry[li]*c->rrlen[li],
				c->Xz[li] + c->rrz[li]*c->rrlen[li],
				c->Xx[ i] + c->rrx[ i]*c->rrlen[ i],
				c->Xy[ i] + c->rry[ i]*c->rrlen[ i],
				c->Xz[ i] + c->rrz[ i]*c->rrlen[ i] );

		} else
#endif
			if( !(c->rrflg[li] == RETYPE_UNDEF || c->rrflg[li] == RETYPE_CURVE)
				&& !(c->rrflg[i] == RETYPE_UNDEF || c->rrflg[i] == RETYPE_CURVE) )
			{
				addEdge( c->Xx[li] + c->rrx[li]*c->rrlen[li],
					c->Xy[li] + c->rry[li]*c->rrlen[li],
					c->Xz[li] + c->rrz[li]*c->rrlen[li],
					c->Xx[ i] + c->rrx[ i]*c->rrlen[ i],
					c->Xy[ i] + c->rry[ i]*c->rrlen[ i],
					c->Xz[ i] + c->rrz[ i]*c->rrlen[ i] );
			}
			li=i;
	}

#if 1	// �[�_�i�܂���n�_�E�����j
	crease *c0 = &(crs[0]);

	if( (c0->Xstype==CETYPE_EGTOP0 || c0->Xstype==CETYPE_EGTOP1) && c0->rlflg[c0->Xsidx]==RETYPE_EGTOP
		|| (c0->Xstype==CETYPE_EGRIGHT0 || c0->Xstype==CETYPE_EGRIGHT1) && c0->rlflg[c0->Xsidx]==RETYPE_EGRIGHT
		|| (c0->Xstype==CETYPE_EGBOTTOM0 || c0->Xstype==CETYPE_EGBOTTOM1) && c0->rlflg[c0->Xsidx]==RETYPE_EGBOTTOM
		|| (c0->Xstype==CETYPE_EGLEFT0 || c0->Xstype==CETYPE_EGLEFT1) && c0->rlflg[c0->Xsidx]==RETYPE_EGLEFT
		|| c0->rlflg[c0->Xsidx]==RETYPE_TRIM )
	{
		// �܂���[�_��ruling�[�_�������Ӂ^�Ȑ� -> �O�p�`��\��
		//if( c0->Xxs!=0.0 || c0->Xys!=0.0 || c0->Xzs!=0.0 )
		{
			addPly( c0->Xx2d[c0->Xsidx], c0->Xy2d[c0->Xsidx],
				c0->Xxs2d, c0->Xys2d,
				c0->Xx2d[c0->Xsidx] + c0->rlx_cp[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy2d[c0->Xsidx] + c0->rly_cp[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xx[c0->Xsidx], c0->Xy[c0->Xsidx], c0->Xz[c0->Xsidx],
				c0->Xxs, c0->Xys, c0->Xzs,
				c0->Xx[c0->Xsidx] + c0->rlx[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rly[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rlz[c0->Xsidx]*c0->rllen[c0->Xsidx] );
			addEdge( c0->Xx[c0->Xsidx] + c0->rlx[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rly[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rlz[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xxs, c0->Xys, c0->Xzs );
		}
	}
	else if( (c0->Xstype==CETYPE_EGTOP0 || c0->Xstype==CETYPE_EGTOP1) && c0->rlflg[c0->Xsidx]==RETYPE_EGLEFT )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`
		//if( c0->Xxs!=0.0 || c0->Xys!=0.0 || c0->Xzs!=0.0 )
		{
			addPly(
				c0->Xx2d[c0->Xsidx] + c0->rlx_cp[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy2d[c0->Xsidx] + c0->rly_cp[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xx2d[c0->Xsidx], c0->Xy2d[c0->Xsidx],
				c0->Xxs2d, c0->Xys2d,
				psx, psy, // top-left
				c0->Xx[c0->Xsidx] + c0->rlx[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rly[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rlz[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xx[c0->Xsidx], c0->Xy[c0->Xsidx], c0->Xz[c0->Xsidx],
				c0->Xxs, c0->Xys, c0->Xzs,
				0,0,0); // undefined
			addEdge( c0->Xx[c0->Xsidx] + c0->rlx[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rly[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rlz[c0->Xsidx]*c0->rllen[c0->Xsidx],
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xxs, c0->Xys, c0->Xzs );
		}
	}
	else if( (c0->Xstype==CETYPE_EGRIGHT0 || c0->Xstype==CETYPE_EGRIGHT1) && c0->rlflg[c0->Xsidx]==RETYPE_EGTOP )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`
		//if( c0->Xxs!=0.0 || c0->Xys!=0.0 || c0->Xzs!=0.0 )
		{
			addPly(
				c0->Xx2d[c0->Xsidx] + c0->rlx_cp[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy2d[c0->Xsidx] + c0->rly_cp[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xx2d[c0->Xsidx], c0->Xy2d[c0->Xsidx],
				c0->Xxs2d, c0->Xys2d,
				pex, psy, // top-right
				c0->Xx[c0->Xsidx] + c0->rlx[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rly[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rlz[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xx[c0->Xsidx], c0->Xy[c0->Xsidx], c0->Xz[c0->Xsidx],
				c0->Xxs, c0->Xys, c0->Xzs,
				0,0,0); // undefined
			addEdge( c0->Xx[c0->Xsidx] + c0->rlx[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rly[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rlz[c0->Xsidx]*c0->rllen[c0->Xsidx],
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xxs, c0->Xys, c0->Xzs );
		}
	}
	else if( (c0->Xstype==CETYPE_EGBOTTOM0 || c0->Xstype==CETYPE_EGBOTTOM1) && c0->rlflg[c0->Xsidx]==RETYPE_EGRIGHT )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`
		//if( c0->Xxs!=0.0 || c0->Xys!=0.0 || c0->Xzs!=0.0 )
		{
			addPly(
				c0->Xx2d[c0->Xsidx] + c0->rlx_cp[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy2d[c0->Xsidx] + c0->rly_cp[c0->Xsidx]*c0->rllen[c0->Xsidx],					
				c0->Xx2d[c0->Xsidx], c0->Xy2d[c0->Xsidx],
				c0->Xxs2d, c0->Xys2d,
				pex, pey, // bottom-right
				c0->Xx[c0->Xsidx] + c0->rlx[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rly[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rlz[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xx[c0->Xsidx], c0->Xy[c0->Xsidx], c0->Xz[c0->Xsidx],
				c0->Xxs, c0->Xys, c0->Xzs,
				0,0,0); // undefined
			addEdge( c0->Xx[c0->Xsidx] + c0->rlx[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rly[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rlz[c0->Xsidx]*c0->rllen[c0->Xsidx],
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xxs, c0->Xys, c0->Xzs );
		}
	}
	else if( (c0->Xstype==CETYPE_EGLEFT0 || c0->Xstype==CETYPE_EGLEFT1) && c0->rlflg[c0->Xsidx]==RETYPE_EGBOTTOM )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`
		//if( c0->Xxs!=0.0 || c0->Xys!=0.0 || c0->Xzs!=0.0 )
		{
			addPly(
				c0->Xx2d[c0->Xsidx] + c0->rlx_cp[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy2d[c0->Xsidx] + c0->rly_cp[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xx2d[c0->Xsidx], c0->Xy2d[c0->Xsidx],
				c0->Xxs2d, c0->Xys2d,
				psx, pey, // bottom-left
				c0->Xx[c0->Xsidx] + c0->rlx[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rly[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rlz[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xx[c0->Xsidx], c0->Xy[c0->Xsidx], c0->Xz[c0->Xsidx],
				c0->Xxs, c0->Xys, c0->Xzs,
				0,0,0); // undefined
			addEdge( c0->Xx[c0->Xsidx] + c0->rlx[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rly[c0->Xsidx]*c0->rllen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rlz[c0->Xsidx]*c0->rllen[c0->Xsidx],
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xxs, c0->Xys, c0->Xzs );
		}
	}
#endif	// �[�_�i�܂���n�_�E�����j

#if 1	// �[�_�i�܂���n�_�E�E���j
	if( (c0->Xstype==CETYPE_EGTOP0 || c0->Xstype==CETYPE_EGTOP1) && c0->rrflg[c0->Xsidx]==RETYPE_EGTOP
		|| (c0->Xstype==CETYPE_EGRIGHT0 || c0->Xstype==CETYPE_EGRIGHT1) && c0->rrflg[c0->Xsidx]==RETYPE_EGRIGHT
		|| (c0->Xstype==CETYPE_EGBOTTOM0 || c0->Xstype==CETYPE_EGBOTTOM1) && c0->rrflg[c0->Xsidx]==RETYPE_EGBOTTOM
		|| (c0->Xstype==CETYPE_EGLEFT0 || c0->Xstype==CETYPE_EGLEFT1) && c0->rrflg[c0->Xsidx]==RETYPE_EGLEFT
		|| c0->rrflg[c0->Xsidx]==RETYPE_TRIM )
	{
		// �܂���[�_��ruling�[�_�������Ӂ^�Ȑ� -> �O�p�`��\��
		//if( c0->Xxs!=0.0 || c0->Xys!=0.0 || c0->Xzs!=0.0 )
		{
			addPly( c0->Xxs2d, c0->Xys2d,
				c0->Xx2d[c0->Xsidx], c0->Xy2d[c0->Xsidx],
				c0->Xx2d[c0->Xsidx] + c0->rrx_cp[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy2d[c0->Xsidx] + c0->rry_cp[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xxs, c0->Xys, c0->Xzs,
				c0->Xx[c0->Xsidx], c0->Xy[c0->Xsidx], c0->Xz[c0->Xsidx],
				c0->Xx[c0->Xsidx] + c0->rrx[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rry[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rrz[c0->Xsidx]*c0->rrlen[c0->Xsidx] );
			addEdge( c0->Xxs, c0->Xys, c0->Xzs,
				c0->Xx[c0->Xsidx] + c0->rrx[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rry[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rrz[c0->Xsidx]*c0->rrlen[c0->Xsidx]	);
		}
	}
	else if((c0->Xstype==CETYPE_EGTOP0 || c0->Xstype==CETYPE_EGTOP1) && c0->rrflg[c0->Xsidx]==RETYPE_EGRIGHT )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxs!=0.0 || c0->Xys!=0.0 || c0->Xzs!=0.0 )
		{
			addPly(	c0->Xxs2d, c0->Xys2d,
				c0->Xx2d[c0->Xsidx], c0->Xy2d[c0->Xsidx],
				c0->Xx2d[c0->Xsidx] + c0->rrx_cp[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy2d[c0->Xsidx] + c0->rry_cp[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				pex, psy, // top-right
				c0->Xxs, c0->Xys, c0->Xzs,
				c0->Xx[c0->Xsidx], c0->Xy[c0->Xsidx], c0->Xz[c0->Xsidx],
				c0->Xx[c0->Xsidx] + c0->rrx[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rry[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rrz[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xxs, c0->Xys, c0->Xzs,
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xx[c0->Xsidx] + c0->rrx[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rry[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rrz[c0->Xsidx]*c0->rrlen[c0->Xsidx] );
		}
	}
	else if((c0->Xstype==CETYPE_EGRIGHT0 || c0->Xstype==CETYPE_EGRIGHT1) && c0->rrflg[c0->Xsidx]==RETYPE_EGBOTTOM )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxs!=0.0 || c0->Xys!=0.0 || c0->Xzs!=0.0 )
		{
			addPly(	c0->Xxs2d, c0->Xys2d,
				c0->Xx2d[c0->Xsidx], c0->Xy2d[c0->Xsidx],
				c0->Xx2d[c0->Xsidx] + c0->rrx_cp[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy2d[c0->Xsidx] + c0->rry_cp[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				pex, pey, // bottom-right
				c0->Xxs, c0->Xys, c0->Xzs,
				c0->Xx[c0->Xsidx], c0->Xy[c0->Xsidx], c0->Xz[c0->Xsidx],
				c0->Xx[c0->Xsidx] + c0->rrx[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rry[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rrz[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xxs, c0->Xys, c0->Xzs,
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xx[c0->Xsidx] + c0->rrx[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rry[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rrz[c0->Xsidx]*c0->rrlen[c0->Xsidx] );
		}
	}
	else if((c0->Xstype==CETYPE_EGBOTTOM0 || c0->Xstype==CETYPE_EGBOTTOM1) && c0->rrflg[c0->Xsidx]==RETYPE_EGLEFT )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxs!=0.0 || c0->Xys!=0.0 || c0->Xzs!=0.0 )
		{
			addPly(	c0->Xxs2d, c0->Xys2d,
				c0->Xx2d[c0->Xsidx], c0->Xy2d[c0->Xsidx],
				c0->Xx2d[c0->Xsidx] + c0->rrx_cp[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy2d[c0->Xsidx] + c0->rry_cp[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				psx, pey, // bottom-left
				c0->Xxs, c0->Xys, c0->Xzs,
				c0->Xx[c0->Xsidx], c0->Xy[c0->Xsidx], c0->Xz[c0->Xsidx],
				c0->Xx[c0->Xsidx] + c0->rrx[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rry[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rrz[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xxs, c0->Xys, c0->Xzs,
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xx[c0->Xsidx] + c0->rrx[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rry[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rrz[c0->Xsidx]*c0->rrlen[c0->Xsidx] );
		}
	}
	else if((c0->Xstype==CETYPE_EGLEFT0 || c0->Xstype==CETYPE_EGLEFT1) && c0->rrflg[c0->Xsidx]==RETYPE_EGTOP )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxs!=0.0 || c0->Xys!=0.0 || c0->Xzs!=0.0 )
		{
			addPly(	c0->Xxs2d, c0->Xys2d,
				c0->Xx2d[c0->Xsidx], c0->Xy2d[c0->Xsidx],
				c0->Xx2d[c0->Xsidx] + c0->rrx_cp[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy2d[c0->Xsidx] + c0->rry_cp[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				psx, psy, // top-left
				c0->Xxs, c0->Xys, c0->Xzs,
				c0->Xx[c0->Xsidx], c0->Xy[c0->Xsidx], c0->Xz[c0->Xsidx],
				c0->Xx[c0->Xsidx] + c0->rrx[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rry[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rrz[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xxs, c0->Xys, c0->Xzs,
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xx[c0->Xsidx] + c0->rrx[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xy[c0->Xsidx] + c0->rry[c0->Xsidx]*c0->rrlen[c0->Xsidx],
				c0->Xz[c0->Xsidx] + c0->rrz[c0->Xsidx]*c0->rrlen[c0->Xsidx] );
		}
	}
#endif	// �[�_�i�܂���n�_�E�E���j

#if 1	// �[�_�i�܂���I�_�E�����j

	if( (c0->Xetype==CETYPE_EGTOP0 || c0->Xetype==CETYPE_EGTOP1) && c0->rlflg[c0->Xeidx]==RETYPE_EGTOP
		|| (c0->Xetype==CETYPE_EGRIGHT0 || c0->Xetype==CETYPE_EGRIGHT1) && c0->rlflg[c0->Xeidx]==RETYPE_EGRIGHT
		|| (c0->Xetype==CETYPE_EGBOTTOM0 || c0->Xetype==CETYPE_EGBOTTOM1) && c0->rlflg[c0->Xeidx]==RETYPE_EGBOTTOM
		|| (c0->Xetype==CETYPE_EGLEFT0 || c0->Xetype==CETYPE_EGLEFT1) && c0->rlflg[c0->Xeidx]==RETYPE_EGLEFT
		|| c0->rlflg[c0->Xeidx]==RETYPE_TRIM )
	{
		// �܂���[�_��ruling�[�_�������Ӂ^�Ȑ� -> �O�p�`��\��
		//if( c0->Xxe!=0.0 || c0->Xye!=0.0 || c0->Xze!=0.0 )
		{
			addPly( c0->Xxe2d, c0->Xye2d,
				c0->Xx2d[c0->Xeidx], c0->Xy2d[c0->Xeidx],
				c0->Xx2d[c0->Xeidx] + c0->rlx_cp[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy2d[c0->Xeidx] + c0->rly_cp[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xxe, c0->Xye, c0->Xze,
				c0->Xx[c0->Xeidx], c0->Xy[c0->Xeidx], c0->Xz[c0->Xeidx],
				c0->Xx[c0->Xeidx] + c0->rlx[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rly[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rlz[c0->Xeidx]*c0->rllen[c0->Xeidx] );
			addEdge( c0->Xxe, c0->Xye, c0->Xze,
				c0->Xx[c0->Xeidx] + c0->rlx[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rly[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rlz[c0->Xeidx]*c0->rllen[c0->Xeidx] );
		}
	}
	else if((c0->Xetype==CETYPE_EGTOP0 || c0->Xetype==CETYPE_EGTOP1) && c0->rlflg[c0->Xeidx]==RETYPE_EGRIGHT )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxe!=0.0 || c0->Xye!=0.0 || c0->Xze!=0.0 )
		{
			addPly( c0->Xxe2d, c0->Xye2d,
				c0->Xx2d[c0->Xeidx], c0->Xy2d[c0->Xeidx],
				c0->Xx2d[c0->Xeidx] + c0->rlx_cp[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy2d[c0->Xeidx] + c0->rly_cp[c0->Xeidx]*c0->rllen[c0->Xeidx],
				pex, psy, // top-right
				c0->Xxe, c0->Xye, c0->Xze,
				c0->Xx[c0->Xeidx], c0->Xy[c0->Xeidx], c0->Xz[c0->Xeidx],
				c0->Xx[c0->Xeidx] + c0->rlx[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rly[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rlz[c0->Xeidx]*c0->rllen[c0->Xeidx],
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xxe, c0->Xye, c0->Xze,
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xx[c0->Xeidx] + c0->rlx[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rly[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rlz[c0->Xeidx]*c0->rllen[c0->Xeidx] );
		}
	}
	else if((c0->Xetype==CETYPE_EGRIGHT0 || c0->Xetype==CETYPE_EGRIGHT1) && c0->rlflg[c0->Xeidx]==RETYPE_EGBOTTOM )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxe!=0.0 || c0->Xye!=0.0 || c0->Xze!=0.0 )
		{
			addPly( c0->Xxe2d, c0->Xye2d,
				c0->Xx2d[c0->Xeidx], c0->Xy2d[c0->Xeidx],
				c0->Xx2d[c0->Xeidx] + c0->rlx_cp[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy2d[c0->Xeidx] + c0->rly_cp[c0->Xeidx]*c0->rllen[c0->Xeidx],
				pex, pey, // bottom-right
				c0->Xxe, c0->Xye, c0->Xze,
				c0->Xx[c0->Xeidx], c0->Xy[c0->Xeidx], c0->Xz[c0->Xeidx],
				c0->Xx[c0->Xeidx] + c0->rlx[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rly[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rlz[c0->Xeidx]*c0->rllen[c0->Xeidx],
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xxe, c0->Xye, c0->Xze,
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xx[c0->Xeidx] + c0->rlx[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rly[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rlz[c0->Xeidx]*c0->rllen[c0->Xeidx] );
		}
	}
	else if((c0->Xetype==CETYPE_EGBOTTOM0 || c0->Xetype==CETYPE_EGBOTTOM1) && c0->rlflg[c0->Xeidx]==RETYPE_EGLEFT )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxe!=0.0 || c0->Xye!=0.0 || c0->Xze!=0.0 )
		{
			addPly( c0->Xxe2d, c0->Xye2d,
				c0->Xx2d[c0->Xeidx], c0->Xy2d[c0->Xeidx],
				c0->Xx2d[c0->Xeidx] + c0->rlx_cp[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy2d[c0->Xeidx] + c0->rly_cp[c0->Xeidx]*c0->rllen[c0->Xeidx],
				psx, pey, // bottom-left
				c0->Xxe, c0->Xye, c0->Xze,
				c0->Xx[c0->Xeidx], c0->Xy[c0->Xeidx], c0->Xz[c0->Xeidx],
				c0->Xx[c0->Xeidx] + c0->rlx[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rly[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rlz[c0->Xeidx]*c0->rllen[c0->Xeidx],
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xxe, c0->Xye, c0->Xze,
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xx[c0->Xeidx] + c0->rlx[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rly[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rlz[c0->Xeidx]*c0->rllen[c0->Xeidx] );
		}
	}
	else if((c0->Xetype==CETYPE_EGLEFT0 || c0->Xetype==CETYPE_EGLEFT1) && c0->rlflg[c0->Xeidx]==RETYPE_EGTOP )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxe!=0.0 || c0->Xye!=0.0 || c0->Xze!=0.0 )
		{
			addPly( c0->Xxe2d, c0->Xye2d,
				c0->Xx2d[c0->Xeidx], c0->Xy2d[c0->Xeidx],
				c0->Xx2d[c0->Xeidx] + c0->rlx_cp[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy2d[c0->Xeidx] + c0->rly_cp[c0->Xeidx]*c0->rllen[c0->Xeidx],
				psx, psy, // top-left
				c0->Xxe, c0->Xye, c0->Xze,
				c0->Xx[c0->Xeidx], c0->Xy[c0->Xeidx], c0->Xz[c0->Xeidx],
				c0->Xx[c0->Xeidx] + c0->rlx[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rly[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rlz[c0->Xeidx]*c0->rllen[c0->Xeidx],
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xxe, c0->Xye, c0->Xze,
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xx[c0->Xeidx] + c0->rlx[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rly[c0->Xeidx]*c0->rllen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rlz[c0->Xeidx]*c0->rllen[c0->Xeidx] );
		}
	}
#endif	// �[�_�i�܂���I�_�E�����j

#if 1	// �[�_�i�܂���I�_�E�E���j
	if( (c0->Xetype==CETYPE_EGTOP0 || c0->Xetype==CETYPE_EGTOP1) && c0->rrflg[c0->Xeidx]==RETYPE_EGTOP
		|| (c0->Xetype==CETYPE_EGRIGHT0 || c0->Xetype==CETYPE_EGRIGHT1) && c0->rrflg[c0->Xeidx]==RETYPE_EGRIGHT
		|| (c0->Xetype==CETYPE_EGBOTTOM0 || c0->Xetype==CETYPE_EGBOTTOM1) && c0->rrflg[c0->Xeidx]==RETYPE_EGBOTTOM
		|| (c0->Xetype==CETYPE_EGLEFT0 || c0->Xetype==CETYPE_EGLEFT1) && c0->rrflg[c0->Xeidx]==RETYPE_EGLEFT
		|| c0->rrflg[c0->Xeidx]==RETYPE_TRIM )
	{
		// �܂���[�_��ruling�[�_�������Ӂ^�Ȑ� -> �O�p�`��\��
		//if( c0->Xxe!=0.0 || c0->Xye!=0.0 || c0->Xze!=0.0 )
		{

			addPly( c0->Xx2d[c0->Xeidx], c0->Xy2d[c0->Xeidx],
				c0->Xxe2d, c0->Xye2d,
				c0->Xx2d[c0->Xeidx] + c0->rrx_cp[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy2d[c0->Xeidx] + c0->rry_cp[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xx[c0->Xeidx], c0->Xy[c0->Xeidx], c0->Xz[c0->Xeidx],
				c0->Xxe, c0->Xye, c0->Xze,
				c0->Xx[c0->Xeidx] + c0->rrx[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rry[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rrz[c0->Xeidx]*c0->rrlen[c0->Xeidx] );
			addEdge( c0->Xx[c0->Xeidx] + c0->rrx[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rry[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rrz[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xxe, c0->Xye, c0->Xze );
		}
	}
	else if((c0->Xetype==CETYPE_EGTOP0 || c0->Xetype==CETYPE_EGTOP1) && c0->rrflg[c0->Xeidx]==RETYPE_EGLEFT )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxe!=0.0 || c0->Xye!=0.0 || c0->Xze!=0.0 )
		{
			addPly(
				c0->Xx2d[c0->Xeidx] + c0->rrx_cp[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy2d[c0->Xeidx] + c0->rry_cp[c0->Xeidx]*c0->rrlen[c0->Xeidx],				
				c0->Xx2d[c0->Xeidx], c0->Xy2d[c0->Xeidx],
				c0->Xxe2d, c0->Xye2d,
				psx, psy, // top-left
				c0->Xx[c0->Xeidx] + c0->rrx[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rry[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rrz[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xx[c0->Xeidx], c0->Xy[c0->Xeidx], c0->Xz[c0->Xeidx],
				c0->Xxe, c0->Xye, c0->Xze,
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xx[c0->Xeidx] + c0->rrx[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rry[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rrz[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xxe, c0->Xye, c0->Xze );
		}
	}
	else if((c0->Xetype==CETYPE_EGRIGHT0 || c0->Xetype==CETYPE_EGRIGHT1) && c0->rrflg[c0->Xeidx]==RETYPE_EGTOP )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxe!=0.0 || c0->Xye!=0.0 || c0->Xze!=0.0 )
		{
			addPly(
				c0->Xx2d[c0->Xeidx] + c0->rrx_cp[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy2d[c0->Xeidx] + c0->rry_cp[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xx2d[c0->Xeidx], c0->Xy2d[c0->Xeidx],
				c0->Xxe2d, c0->Xye2d,
				pex, psy, // top-right
				c0->Xx[c0->Xeidx] + c0->rrx[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rry[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rrz[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xx[c0->Xeidx], c0->Xy[c0->Xeidx], c0->Xz[c0->Xeidx],
				c0->Xxe, c0->Xye, c0->Xze,
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xx[c0->Xeidx] + c0->rrx[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rry[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rrz[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xxe, c0->Xye, c0->Xze );
		}
	}
	else if((c0->Xetype==CETYPE_EGBOTTOM0 || c0->Xetype==CETYPE_EGBOTTOM1) && c0->rrflg[c0->Xeidx]==RETYPE_EGRIGHT )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxe!=0.0 || c0->Xye!=0.0 || c0->Xze!=0.0 )
		{
			addPly(
				c0->Xx2d[c0->Xeidx] + c0->rrx_cp[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy2d[c0->Xeidx] + c0->rry_cp[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xx2d[c0->Xeidx], c0->Xy2d[c0->Xeidx],
				c0->Xxe2d, c0->Xye2d,
				pex, pey, // bottom-right
				c0->Xx[c0->Xeidx] + c0->rrx[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rry[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rrz[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xx[c0->Xeidx], c0->Xy[c0->Xeidx], c0->Xz[c0->Xeidx],
				c0->Xxe, c0->Xye, c0->Xze,
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xx[c0->Xeidx] + c0->rrx[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rry[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rrz[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xxe, c0->Xye, c0->Xze );
		}
	}
	else if((c0->Xetype==CETYPE_EGLEFT0 || c0->Xetype==CETYPE_EGLEFT1) && c0->rrflg[c0->Xeidx]==RETYPE_EGBOTTOM )
	{
		// �O�p�`�ŕϊ����߂Ďl�p�`	
		//if( c0->Xxe!=0.0 || c0->Xye!=0.0 || c0->Xze!=0.0 )
		{
			addPly(
				c0->Xx2d[c0->Xeidx] + c0->rrx_cp[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy2d[c0->Xeidx] + c0->rry_cp[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xx2d[c0->Xeidx], c0->Xy2d[c0->Xeidx],
				c0->Xxe2d, c0->Xye2d,
				psx, pey, // bottom-left
				c0->Xx[c0->Xeidx] + c0->rrx[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rry[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rrz[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xx[c0->Xeidx], c0->Xy[c0->Xeidx], c0->Xz[c0->Xeidx],
				c0->Xxe, c0->Xye, c0->Xze,
				0.0, 0.0, 0.0 ); // undefined
			addEdge( c0->Xx[c0->Xeidx] + c0->rrx[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xy[c0->Xeidx] + c0->rry[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				c0->Xz[c0->Xeidx] + c0->rrz[c0->Xeidx]*c0->rrlen[c0->Xeidx],
				plx[(plcnt-1)*4+3], ply[(plcnt-1)*4+3], plz[(plcnt-1)*4+3],
				c0->Xxe, c0->Xye, c0->Xze );
		}
	}
#endif	// �[�_�i�܂���I�_�E�E���j

end:
	return ret;
}

int papermodel::makeObj3( double Dthk )
{
	int ret=0;
	double Swth = -1.0;
	double Swth_Dthk = Swth - Dthk;
	double Dx=0.0, Dy=0.0, Dz=0.0, Dlen;	// �����o������
	for( int i=crs[0].Xsidx; i<crs[0].Xeidx+1; i++ ){
		Dx += crs[0].Bx[i];
		Dy += crs[0].By[i];
		Dz += crs[0].Bz[i];
	}
	Dlen = sqrt( Dx*Dx + Dy*Dy + Dz*Dz );
	Dx/=Dlen;	Dy/=Dlen;	Dz/=Dlen;

	o_vcnt = o_fcnt = 0;
	for( int i=0; i<plcnt; i++ ){
		o_fvcnt[o_fcnt] = MIN( plvcnt[i], 4 );
		for( int j=0; j<o_fvcnt[o_fcnt]; j++ ){
			if( o_vcnt+1 < MAX_PL_FCNT*4 && o_fcnt+1 < MAX_PL_FCNT ){
				o_vx[o_vcnt] = plx[i*4 + j];	
				o_vy[o_vcnt] = ply[i*4 + j];
				o_vz[o_vcnt] = plz[i*4 + j];
				o_fi[o_fcnt][j] = o_vcnt;
				o_vcnt++;
			}
		}
		if( o_fcnt+1 < MAX_PL_FCNT ){
			o_fcnt++;
		}
#if 0 // backside
		o_fvcnt[o_fcnt] = MIN( plvcnt[i], 4 );
		for( int j=0; j<o_fvcnt[o_fcnt]; j++ ){
			o_vx[o_vcnt] = plx[i*4 + j] +Dx*Dthk;	
			o_vy[o_vcnt] = ply[i*4 + j] +Dy*Dthk;
			o_vz[o_vcnt] = plz[i*4 + j] +Dz*Dthk;
			o_fi[o_fcnt][o_fvcnt[o_fcnt]-1-j] = o_vcnt;
			o_vcnt++;
		}
		o_fcnt++;
#endif
	}
#if 0 // edge
	for( int i=0; i<plecnt; i++ ){
		o_fvcnt[o_fcnt] = 4;
		o_vx[o_vcnt] = plex[i*2  ];	o_vy[o_vcnt] = pley[i*2  ];	o_vz[o_vcnt] = plez[i*2  ];	o_fi[o_fcnt][0] = o_vcnt;	o_vcnt++;
		o_vx[o_vcnt] = plex[i*2+1];	o_vy[o_vcnt] = pley[i*2+1];	o_vz[o_vcnt] = plez[i*2+1];	o_fi[o_fcnt][1] = o_vcnt;	o_vcnt++;
		o_vx[o_vcnt] = plex[i*2+1]+Dx*Dthk;
		o_vy[o_vcnt] = pley[i*2+1]+Dy*Dthk;
		o_vz[o_vcnt] = plez[i*2+1]+Dz*Dthk;	o_fi[o_fcnt][2] = o_vcnt;	o_vcnt++;
		o_vx[o_vcnt] = plex[i*2  ]+Dx*Dthk;
		o_vy[o_vcnt] = pley[i*2  ]+Dy*Dthk;
		o_vz[o_vcnt] = plez[i*2  ]+Dz*Dthk;	o_fi[o_fcnt][3] = o_vcnt;	o_vcnt++;
		o_fcnt++;
	}
#endif
#if 0 // futa
	for( int i=0; i<plcnt; i++ ){
		o_fvcnt[o_fcnt] = MIN( plvcnt[i], 4 );
		for( int j=0; j<o_fvcnt[o_fcnt]; j++ ){
			o_vx[o_vcnt] = plx[i*4 + j] +Dx*Swth;	
			o_vy[o_vcnt] = ply[i*4 + j] +Dy*Swth;
			o_vz[o_vcnt] = plz[i*4 + j] +Dz*Swth;
			o_fi[o_fcnt][o_fvcnt[o_fcnt]-1-j] = o_vcnt;
			o_vcnt++;
		}
		o_fcnt++;
		o_fvcnt[o_fcnt] = MIN( plvcnt[i], 4 );
		for( int j=0; j<o_fvcnt[o_fcnt]; j++ ){
			o_vx[o_vcnt] = plx[i*4 + j] +Dx*Swth_Dthk;	
			o_vy[o_vcnt] = ply[i*4 + j] +Dy*Swth_Dthk;
			o_vz[o_vcnt] = plz[i*4 + j] +Dz*Swth_Dthk;
			o_fi[o_fcnt][j] = o_vcnt;
			o_vcnt++;
		}
		o_fcnt++;
	}
	for( int i=0; i<plecnt; i++ ){
		o_fvcnt[o_fcnt] = 4;
		o_vx[o_vcnt] = plex[i*2]+Dx*Swth_Dthk;
		o_vy[o_vcnt] = pley[i*2]+Dy*Swth_Dthk;
		o_vz[o_vcnt] = plez[i*2]+Dz*Swth_Dthk;	o_fi[o_fcnt][0] = o_vcnt;	o_vcnt++;
		o_vx[o_vcnt] = plex[i*2+1]+Dx*Swth_Dthk;
		o_vy[o_vcnt] = pley[i*2+1]+Dy*Swth_Dthk;
		o_vz[o_vcnt] = plez[i*2+1]+Dz*Swth_Dthk;	o_fi[o_fcnt][1] = o_vcnt;	o_vcnt++;
		o_vx[o_vcnt] = plex[i*2+1]+Dx*Swth;
		o_vy[o_vcnt] = pley[i*2+1]+Dy*Swth;
		o_vz[o_vcnt] = plez[i*2+1]+Dz*Swth;		o_fi[o_fcnt][2] = o_vcnt;	o_vcnt++;
		o_vx[o_vcnt] = plex[i*2  ]+Dx*Swth;
		o_vy[o_vcnt] = pley[i*2  ]+Dy*Swth;
		o_vz[o_vcnt] = plez[i*2  ]+Dz*Swth;		o_fi[o_fcnt][3] = o_vcnt;	o_vcnt++;
		o_fcnt++;
	}
#endif

end:
	return ret;
}
