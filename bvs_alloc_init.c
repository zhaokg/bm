#include <string.h>  //memset
#include <math.h>    //sqrt

#include "abc_000_warning.h"

#include "abc_vec.h"
#include "abc_ts_func.h"
#include "abc_mem.h"
#include "bvs_header.h"
#include "bvs_fun.h"

#define MODEL (*model)
//f32_fill_val_matrixdiag(MODEL.SIG2, opt->prior.sig2, q);
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#define max3(a,b,c)         max( max(a,b), c)
#define max4(a,b,c,d)       max( max3(a,b,c),d)
#define min(a,b)            (((a) < (b)) ? (a) : (b))
 
void AllocInit_Model(BVS_MODEL_PTR model, BVS_OPTIONS_PTR opt, MemPointers* MEM) {	

	// not needed bcz it is initialized to 0 when model is declared
	// memset(model, 0, sizeof(BVS_MODEL));
	I32 N      = opt->io.meta.Nobj;
	I32 Kmax   = opt->prior.Kmax;
	I32 q      = 1;

	// sig2 is intialized only once here. This default value is used only for the the first chain
	// of the first pixel. All other chains or pixels will be using the existing value left from 
	// previous runs. This may have some unintended consequences
 
	// BEASTV2
	MODEL.sig2       = &MODEL._sig2; // sig2 is needed only for generating beta,  
	MODEL.precVec    = &MODEL._precVal; // sig2 is needed only for generating beta,  
	MODEL.logPrecVec = &MODEL._logPrecVal; // sig2 is needed only for generating beta,  

	MODEL.sig2[0]       = opt->prior.sig2;
	MODEL.precVec[0]    = opt->prior.precValue;
	MODEL.logPrecVec[0] = logf(MODEL.precVec[0]);

	MODEL.Kterms = 0;

	MemNode nodes[100];
	int     nid = 0;

	nodes[nid++] = (MemNode){  &MODEL.terms,		  sizeof(BVS_TERMS) * Kmax,	.align = 64 }; 
	nodes[nid++] = (MemNode){  &MODEL.beta,			  sizeof(F32) * Kmax * q,	.align = 4 };  

	nodes[nid++] = (MemNode){  &MODEL.curr.XtX,		  sizeof(F32) * Kmax * Kmax,	.align = 4 };  
	nodes[nid++] = (MemNode){  &MODEL.curr.XtY ,	  sizeof(F32) * Kmax * q,	.align = 4 };  
	nodes[nid++] = (MemNode){  &MODEL.curr.cholXtX ,  sizeof(F32)* Kmax * Kmax,	.align = 4 };  
	nodes[nid++] = (MemNode){  &MODEL.curr.beta_mean, sizeof(I32)* Kmax * q,		.align = 4 };  
 
	nodes[nid++] = (MemNode){  &MODEL.prop.XtX,		   sizeof(F32) * Kmax * Kmax,	.align = 4 }; 
	nodes[nid++] = (MemNode){  &MODEL.prop.XtY ,	   sizeof(F32) * Kmax * q,	.align = 4 };  
	nodes[nid++] = (MemNode){  &MODEL.prop.cholXtX ,   sizeof(F32)* Kmax * Kmax,	.align = 4 };  
	nodes[nid++] = (MemNode){  &MODEL.prop.beta_mean,  sizeof(I32)* Kmax * q,		.align = 4 };  

	// Allocate 3 arrays altogether: Deviation (N*q) + avgDeivation (q) + extremePosVec (Npad16)
	nodes[nid++] = (MemNode){ &MODEL.deviation,		sizeof(F32) * N * q,		.align = 4 };    // changed for MRBEAST
	nodes[nid++] = (MemNode){ &MODEL.avgDeviation,	sizeof(F32) *q,				.align = 4 };  

	//	MODEL.extremePosVec needs to be 8-byte aligned; otherwise,there is an 
	//  run-time error in Debian Linux, R-devel, GCC ASAN/UBSAN.
	I32 Npad16 = (N + 15) / 16 * 16;
	nodes[nid++] = (MemNode){ &MODEL.extremePosVec,   sizeof(I08) * Npad16,	.align = 8 };  
	// MODEL.extremPosNum    = 0; // No need to initialize bcz it will be reset to yInfo.n in beast_core

	nodes[nid++] = (MemNode){ .addr = NULL, };	
	MEM->alloclist(MEM, nodes, AggregatedMemAlloc, NULL);


	Alloc_BasisState_All(&MODEL, opt, MEM);

	 
#undef MODEL
}

void Alloc_Xterms(F32PTR * Xmars, F32PTR*  Xnewterm,  BVS_OPTIONS_PTR opt, MemPointers * MEM) 
{	
	I32 N        =  opt->io.meta.Nobj;
	I32 Kmax     =  opt->prior.Kmax;
	I32 Knewterm =  10;  //TODO: just use 10 arbitraly
 	
	MemNode nodes[] = {
		{Xmars,		 sizeof(F32)*N*Kmax,	 .align=64},
		{Xnewterm,   sizeof(F32)*N*Knewterm + sizeof(I16)*opt->prior.Kmax, .align = 64}, //Extra bytes needed as a temp buf for Kposition in UpdateBasisAll
		{NULL, },
	};

	MEM->alloclist(MEM, nodes, AggregatedMemAlloc, NULL);
	
}

void Alloc_xInfo_yInfo(BVS_XINFO_PTR xInfo, BVS_YINFO_PTR yInfo, BVS_OPTIONS_PTR opt, I16PTR typedid, MemPointers* MEM) {
	I32  N        = opt->io.meta.Nobj;
	I32  Pobj     = opt->io.meta.Pobj;
	I32  Pall     = opt->prior.Pall;
	
	I32  Plin = opt->prior.Pbasis[LinearTerm];
	I32  Psqr = opt->prior.Pbasis[QuadraticTerm]; 

	I32 Nangle = hasTerm(Hinge2DRTerm) || hasTerm(Hinge1DRTerm) ? 40 : 0;

	I32 NhingeMEM = 0;  // size of mem allocated for xInfo.mem_GenOneSeg_hinege2d_2dr_1dr
	if (hasTerm(Hinge2DTerm))   NhingeMEM = N;
	if (hasTerm(Hinge1DRTerm))  NhingeMEM = N;
	if (hasTerm(Hinge2DRTerm))  NhingeMEM = 3*N;

	I32  q    = 1;  //q=1 for BEASTV2
	I32  qxq  = q*q;
	
	MemNode  nodes[100];
	VOIDPTR* badnodes[100];
	int      nid = 0, bid=0;	

	// Xorg contain the input X plus the derived X such as the sin/cos for BEAST and FFT
	nodes[nid++] = (MemNode){  &xInfo->Xorg,			 sizeof(F32) * N * Pall,	.align = 64 };
	nodes[nid++] = (MemNode){  &xInfo->Xones,			 sizeof(F32) * N ,			.align = 64 };
	nodes[nid++] = (MemNode){  &xInfo->X1norm ,			 sizeof(F32) * N * Plin,	.align = 64 };
	nodes[nid++] = (MemNode){  &xInfo->X2 ,				 sizeof(F32) * N * Psqr,	.align = 64 };
	nodes[nid++] = (MemNode){  &xInfo->Nunique ,		 sizeof(I32)* Pobj,			.align = 4 };
	nodes[nid++] = (MemNode){  &xInfo->SortedUnqiueIdx0,  sizeof(I32)*N* Pobj,		.align = 4 };
	nodes[nid++] = (MemNode){  &xInfo->RowsMissing ,	 sizeof(I32)*N,				.align = 4 };
	
	nodes[nid++] = (MemNode){  &xInfo->mean ,			 sizeof(F32) * Pall,		.align = 4 };
	nodes[nid++] = (MemNode){  &xInfo->sd ,				 sizeof(F32) * Pall,		.align = 4 };
	nodes[nid++] = (MemNode){  &xInfo->mean2 ,			 sizeof(F32) * Psqr,		.align = 4 };
	nodes[nid++] = (MemNode){  &xInfo->sd2 ,			 sizeof(F32) * Psqr,		.align = 4 };
	nodes[nid++] = (MemNode){  &xInfo->mem_GenOneSeg_hinege2d_2dr_1dr ,	max( sizeof(F32) * NhingeMEM, sizeof(I16)*opt->prior.Kmax),.align = 4};
	nodes[nid++] = (MemNode){  &xInfo->anglecoff ,	     sizeof(F32) * Nangle * 2,		.align = 4 };

	// Sizes: Y[Npad], sd[q], mean[q], YtY_plus_alpha2Q[q*q]
	nodes[nid++] = (MemNode){  &yInfo->Y ,					sizeof(F32) * N*q,  .align = 64 };
	nodes[nid++] = (MemNode){  &yInfo->YtY_plus_alpha2Q ,	sizeof(F32) * qxq,  .align = 4 };
	nodes[nid++] = (MemNode){  &yInfo->mean ,				sizeof(F32) * q,    .align = 4 };
	nodes[nid++] = (MemNode){  &yInfo->sd ,					sizeof(F32) * q,    .align = 4 };

	nodes[nid++] = (MemNode){ .addr = NULL, };
	
	if (!hasTerm(QuadraticTerm)) {
		badnodes[bid++] = &xInfo->X2;
		badnodes[bid++] = &xInfo->mean2;
		badnodes[bid++] = &xInfo->sd2;
	} 
	if (!hasTerm(StepTerm)) {
		badnodes[bid++] = &xInfo->Xones;
	}
	if (!hasTerm(StairTerm)) {}
	if (!hasTerm(PieceTerm)) {}
	if (!hasTerm(HingeTerm)) {	}
	if (!hasTerm(ChangeTerm)) {	}
	if (!hasTerm(Hinge2DTerm)) {}
	if (!hasTerm(Hinge2DRTerm)) {}

	if (! ( hasTerm(Hinge2DRTerm)|| hasTerm(Hinge1DRTerm))) {
		badnodes[bid++] = &xInfo->anglecoff;
	}

	if (! (hasTerm(Hinge2DRTerm) || hasTerm(Hinge2DTerm)|| hasTerm(Hinge1DRTerm))) {
		badnodes[bid++] = &xInfo->mem_GenOneSeg_hinege2d_2dr_1dr;
	}

	badnodes[bid++] = NULL;
	MEM->alloclist(MEM, nodes, AggregatedMemAlloc, badnodes);

   
	xInfo->Pobj    = Pobj;
	xInfo->Pall    = Pall;
	for (int i = 0; i < MAX_NUM_BASIS; i++) {
		xInfo->xi_Pbasis[i]    = opt->prior.Pbasis[i];
		xInfo->xi_Pindices0[i] = opt->prior.Pindices0[i];			
	}
	xInfo->xi_PchangeY0 = opt->prior.PchangeY0;
	xInfo->xi_PseasonY0 = opt->prior.PseasonY0;

 	xInfo->Nangle   = Nangle;
	if (Nangle > 0 ) {
		f32_seq(xInfo->anglecoff, 0, 2 * 3.1415926 / Nangle, Nangle);
		f32_copy(xInfo->anglecoff, xInfo->anglecoff + Nangle, Nangle);
		f32_sincos_vec_inplace(xInfo->anglecoff + Nangle, xInfo->anglecoff, Nangle);
	}
}

#include "abc_000_warning.h"