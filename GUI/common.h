#ifndef COMMON_H
#define COMMON_H

# define N_C_MODE 4
typedef enum _C_MODE { CMODE_A, CMODE_B, CMODE_C, CMODE_R } C_MODE;

# define N_P_IDX 10
typedef enum _P_IDX {
	P_CV2D, P_CV3D, P_TRSN, P_FLDA, P_RULL, P_RULR,
	P_TRSN1, P_FLDA1, P_RULL1, P_RULR1
} P_IDX;

//#define DEBUG_MOUSE_EVENT // マウスイベントの状態を出力

#endif
