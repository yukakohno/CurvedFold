#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include "crease.h"
#include "util.h"
#include "Bspline.h"

crease::crease()
{
	init();
}

crease::~crease()
{
}

void crease::init()
{
	rl=0;
	flg_org=0;
	flg_src_s1e1=0;

	org_idx = -1;
	org_cnt = 0;
	org_x = org_y = NULL;

	Pcnt = 0;
	memset( Px2d, 0, sizeof(double)*MAX_CPCNT );
	memset( Py2d, 0, sizeof(double)*MAX_CPCNT );
	memset( Px, 0, sizeof(double)*MAX_CPCNT );
	memset( Py, 0, sizeof(double)*MAX_CPCNT );
	memset( Pz, 0, sizeof(double)*MAX_CPCNT );
	memset( Pa, 0, sizeof(double)*MAX_CPCNT );

	unit_m44( m3 );
	unit_m33( m2 );

	initX();
	maxrlen = 400;

	Xsidx = Xeidx = -1;
	Xxs2d = Xys2d = Xxs = Xys = Xzs = Xxs2d0 = Xys2d0 = 0.0;
	Xxe2d = Xye2d = Xxe = Xye = Xze = Xxe2d0 = Xye2d0 = 0.0;
}

void crease::initX()
{
	Xcnt = 0;
	memset( Xx2d, 0, sizeof(double)*MAX_SPCNT );
	memset( Xy2d, 0, sizeof(double)*MAX_SPCNT );
	memset( Tx2d, 0, sizeof(double)*MAX_SPCNT );
	memset( Ty2d, 0, sizeof(double)*MAX_SPCNT );
	memset( Nx2d, 0, sizeof(double)*MAX_SPCNT );
	memset( Ny2d, 0, sizeof(double)*MAX_SPCNT );
	memset( d2d, 0, sizeof(double)*MAX_SPCNT );
	memset( k2d, 0, sizeof(double)*MAX_SPCNT );

	memset( Xx, 0, sizeof(double)*MAX_SPCNT );
	memset( Xy, 0, sizeof(double)*MAX_SPCNT );
	memset( Xz, 0, sizeof(double)*MAX_SPCNT );
	memset( Tx, 0, sizeof(double)*MAX_SPCNT );
	memset( Ty, 0, sizeof(double)*MAX_SPCNT );
	memset( Tz, 0, sizeof(double)*MAX_SPCNT );
	memset( Nx, 0, sizeof(double)*MAX_SPCNT );
	memset( Ny, 0, sizeof(double)*MAX_SPCNT );
	memset( Nz, 0, sizeof(double)*MAX_SPCNT );
	memset( Bx, 0, sizeof(double)*MAX_SPCNT );
	memset( By, 0, sizeof(double)*MAX_SPCNT );
	memset( Bz, 0, sizeof(double)*MAX_SPCNT );
	memset( dx, 0, sizeof(double)*MAX_SPCNT );
	memset( kv, 0, sizeof(double)*MAX_SPCNT );
	memset( tr, 0, sizeof(double)*MAX_SPCNT );

	memset( alpha, 0, sizeof(double)*MAX_SPCNT );
	memset( cosa, 0, sizeof(double)*MAX_SPCNT );
	memset( sina, 0, sizeof(double)*MAX_SPCNT );
	memset( tana, 0, sizeof(double)*MAX_SPCNT );
	memset( da, 0, sizeof(double)*MAX_SPCNT );
	memset( betal, 0, sizeof(double)*MAX_SPCNT );
	memset( betar, 0, sizeof(double)*MAX_SPCNT );
	memset( cotbl, 0, sizeof(double)*MAX_SPCNT );
	memset( cotbr, 0, sizeof(double)*MAX_SPCNT );
	memset( cosbl, 0, sizeof(double)*MAX_SPCNT );
	memset( cosbr, 0, sizeof(double)*MAX_SPCNT );
	memset( sinbl, 0, sizeof(double)*MAX_SPCNT );
	memset( sinbr, 0, sizeof(double)*MAX_SPCNT );
	memset( rlx, 0, sizeof(double)*MAX_SPCNT );
	memset( rly, 0, sizeof(double)*MAX_SPCNT );
	memset( rlz, 0, sizeof(double)*MAX_SPCNT );
	memset( rrx, 0, sizeof(double)*MAX_SPCNT );
	memset( rry, 0, sizeof(double)*MAX_SPCNT );
	memset( rrz, 0, sizeof(double)*MAX_SPCNT );
	memset( rlx_cp, 0, sizeof(double)*MAX_SPCNT );
	memset( rly_cp, 0, sizeof(double)*MAX_SPCNT );
	memset( rrx_cp, 0, sizeof(double)*MAX_SPCNT );
	memset( rry_cp, 0, sizeof(double)*MAX_SPCNT );
	memset( rllen, 0, sizeof(double)*MAX_SPCNT );
	memset( rrlen, 0, sizeof(double)*MAX_SPCNT );
	init_rlflg();
	init_rrflg();

	Xxs2d = Xys2d = Xxs = Xys = Xzs = 0.0;
	Xxe2d = Xye2d = Xxe = Xye = Xze = 0.0;
}

int crease::copy( crease *c )
{
	int ret=0;

	rl = c->rl;
	flg_org = c->flg_org;

	org_idx = c->org_idx;
	org_cnt = c->org_cnt;
	org_x = c->org_x;
	org_y = c->org_y;

	Pcnt = c->Pcnt;
	memcpy( Px2d, c->Px2d, sizeof(double)*MAX_CPCNT );
	memcpy( Py2d, c->Py2d, sizeof(double)*MAX_CPCNT );
	memcpy(  Px,  c->Px, sizeof(double)*MAX_CPCNT );
	memcpy(  Py,  c->Py, sizeof(double)*MAX_CPCNT );
	memcpy(  Pz,  c->Pz, sizeof(double)*MAX_CPCNT );
	memcpy(  Pa,  c->Pa, sizeof(double)*MAX_CPCNT );

	memcpy( m3, c->m3, sizeof(double)*16 );
	memcpy( m2, c->m2, sizeof(double)*9 );

	Xcnt = c->Xcnt;
	memcpy( Xx2d, c->Xx2d, sizeof(double)*MAX_SPCNT );
	memcpy( Xy2d, c->Xy2d, sizeof(double)*MAX_SPCNT );
	memcpy( Tx2d, c->Tx2d, sizeof(double)*MAX_SPCNT );
	memcpy( Ty2d, c->Ty2d, sizeof(double)*MAX_SPCNT );
	memcpy( Nx2d, c->Nx2d, sizeof(double)*MAX_SPCNT );
	memcpy( Ny2d, c->Ny2d, sizeof(double)*MAX_SPCNT );
	memcpy( d2d, c->d2d, sizeof(double)*MAX_SPCNT );
	memcpy( k2d, c->k2d, sizeof(double)*MAX_SPCNT );

	memcpy( Xx, c->Xx, sizeof(double)*MAX_SPCNT );
	memcpy( Xy, c->Xy, sizeof(double)*MAX_SPCNT );
	memcpy( Xz, c->Xz, sizeof(double)*MAX_SPCNT );
	memcpy( Tx, c->Tx, sizeof(double)*MAX_SPCNT );
	memcpy( Ty, c->Ty, sizeof(double)*MAX_SPCNT );
	memcpy( Tz, c->Tz, sizeof(double)*MAX_SPCNT );
	memcpy( Nx, c->Nx, sizeof(double)*MAX_SPCNT );
	memcpy( Ny, c->Ny, sizeof(double)*MAX_SPCNT );
	memcpy( Nz, c->Nz, sizeof(double)*MAX_SPCNT );
	memcpy( Bx, c->Bx, sizeof(double)*MAX_SPCNT );
	memcpy( By, c->By, sizeof(double)*MAX_SPCNT );
	memcpy( Bz, c->Bz, sizeof(double)*MAX_SPCNT );
	memcpy( dx, c->dx, sizeof(double)*MAX_SPCNT );
	memcpy( kv, c->kv, sizeof(double)*MAX_SPCNT );
	memcpy( tr, c->tr, sizeof(double)*MAX_SPCNT );

	memcpy( alpha, c->alpha, sizeof(double)*MAX_SPCNT );
	memcpy( cosa, c->cosa, sizeof(double)*MAX_SPCNT );
	memcpy( sina, c->sina, sizeof(double)*MAX_SPCNT );
	memcpy( tana, c->tana, sizeof(double)*MAX_SPCNT );
	memcpy( da, c->da, sizeof(double)*MAX_SPCNT );
	memcpy( betal, c->betal, sizeof(double)*MAX_SPCNT );
	memcpy( betar, c->betar, sizeof(double)*MAX_SPCNT );
	memcpy( cotbl, c->cotbl, sizeof(double)*MAX_SPCNT );
	memcpy( cotbr, c->cotbr, sizeof(double)*MAX_SPCNT );
	memcpy( cosbl, c->cosbl, sizeof(double)*MAX_SPCNT );
	memcpy( cosbr, c->cosbr, sizeof(double)*MAX_SPCNT );
	memcpy( sinbl, c->sinbl, sizeof(double)*MAX_SPCNT );
	memcpy( sinbr, c->sinbr, sizeof(double)*MAX_SPCNT );
	memcpy( rlx, c->rlx, sizeof(double)*MAX_SPCNT );
	memcpy( rly, c->rly, sizeof(double)*MAX_SPCNT );
	memcpy( rlz, c->rlz, sizeof(double)*MAX_SPCNT );
	memcpy( rrx, c->rrx, sizeof(double)*MAX_SPCNT );
	memcpy( rry, c->rry, sizeof(double)*MAX_SPCNT );
	memcpy( rrz, c->rrz, sizeof(double)*MAX_SPCNT );
	memcpy( rlx_cp, c->rlx_cp, sizeof(double)*MAX_SPCNT );
	memcpy( rly_cp, c->rly_cp, sizeof(double)*MAX_SPCNT );
	memcpy( rrx_cp, c->rrx_cp, sizeof(double)*MAX_SPCNT );
	memcpy( rry_cp, c->rry_cp, sizeof(double)*MAX_SPCNT );
	memcpy( rllen, c->rllen, sizeof(double)*MAX_SPCNT );
	memcpy( rrlen, c->rrlen, sizeof(double)*MAX_SPCNT );
	memcpy( rlflg, c->rlflg, sizeof(ruling_end_type)*MAX_SPCNT );
	memcpy( rrflg, c->rrflg, sizeof(ruling_end_type)*MAX_SPCNT );

	Xsidx = c->Xsidx;	Xeidx = c->Xeidx;
	Xxs2d = c->Xxs2d;	Xys2d = c->Xys2d;
	Xxe2d = c->Xxe2d;	Xye2d = c->Xye2d;
	Xxs = c->Xxs;		Xys = c->Xys;		Xzs = c->Xzs;
	Xxe = c->Xxe;		Xye = c->Xye;		Xze = c->Xze;
	Xstype = c->Xstype;	Xetype = c->Xetype;

	return ret;
}

void crease::init_rlflg()
{
	for( int i=0; i<MAX_SPCNT; i++ ){
		rlflg[i] = RETYPE_UNDEF;
	}
}

void crease::init_rrflg()
{
	for( int i=0; i<MAX_SPCNT; i++ ){
		rrflg[i] = RETYPE_UNDEF;
	}
}

