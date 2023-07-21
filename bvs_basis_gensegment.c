#include <math.h> //sqrtf
#include "abc_000_macro.h" 
#include "abc_000_warning.h"

#include "abc_vec.h"
#include "abc_blas_lapack_lib.h" //r_ippsSet_32 r_cblas_scopy
#include "abc_common.h"          // normalize_x_factor normalize
#include "bvs_header.h"

#include "bvs_fun.h"

/*
* 	r_ippsSet_32f(scale, X + (seg->R1) - 1, Nseg);
	f32_seq(X + seg->R1 - 1, a, b, Nseg);
	r_cblas_scopy(Nseg, TERMS, 1, X + seg->R1 - 1, 1);
	f32_normalize_x_factor_inplace(X + seg->R1 - 1, Nseg, scale);
	f32_seq(X + k * Npad + r1 - 1, 1, 1, segLength);
	f32_normalize_inplace(X + k * Npad + r1 - 1, segLength);
	f32_mul_val_inplace(scale, X + k * Npad + r1 - 1, segLength);	
	f32_mul_val_inplace(scale, X + k * Npad + r1 - 1, segLength);
	r_ippsSubC_32f_I(sum / Nseg, X + r1 - 1, Nseg); //centering the data vector
	r_cblas_sscal(Nseg, scalingFactor, X + seg->R1 - 1, 1L);
	r_ippsMulC_32f_I(scale, Xt + r1 - 1, segLength);
*/

static void GenOneTerm_const(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {
	r_ippsSet_32f(1.0, xcol, xinfo->N);	
}

static void GenOneTerm_linear(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {
	int vidx0  = term->var0;  // one-based
	//int pidx0 = xinfo->xi_Pindices0[LinearTerm][var0];  //  not used
	int N     = xinfo->N;
 	r_cblas_scopy(N, xinfo->X1norm + N * vidx0, 1L, xcol, 1L);
	term->mean = xinfo->mean[vidx0];
	term->sd   = xinfo->sd[vidx0];
	
}
static void GenOneTerm_quadratic(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {
	int    vidx0  = term->var0;  //  1-based 
	int    pidx0 = xinfo->xi_Pindices0[QuadraticTerm][vidx0];  //  1-based 
	int    N     = xinfo->N;	

	if (xinfo->X2) {
		r_cblas_scopy(N, xinfo->X2 + N * vidx0, 1L, xcol, 1L);
		term->mean = xinfo->mean2[vidx0];
		term->sd   = xinfo->sd2[vidx0];
	} else {
		r_cblas_scopy(N, xinfo->Xorg + N * pidx0, 1L, xcol, 1L);
		f32_mul_vec_inplace(xcol, xcol, N);
		f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
	} 	
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
static F32PTR GetKnot1D(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo,  F32PTR xknot)  { 	
	int    N       = xinfo->N;	
	int    pidx0   = xinfo->xi_Pindices0[term->type][term->var0];  	//// zero-based for basisState but 1-based for vars
	F32PTR X       = xinfo->Xorg + N * pidx0;
	int    knotidx = xinfo->SortedUnqiueIdx0[N * pidx0 + term->knot - 1L];	
	xknot[0]       = X[knotidx]; 
 	return X;
}
static void GetKnot2D(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xknot, F32PTR *xcols) {
	int    N       =  xinfo->N;	
	int    pidx0   = xinfo->xi_Pindices0[term->type][term->var0s[0]];	 //// zero-based for basisState but 1-based for vars
	F32PTR X       = xinfo->Xorg + N * pidx0;
	int    knotidx = xinfo->SortedUnqiueIdx0[N * pidx0 + term->knots[0] - 1L];
	xknot[0] = X[knotidx];
	xcols[0] = X;

	pidx0    = xinfo->xi_Pindices0[term->type][term->var0s[1]];	 //// zero-based for basisState but 1-based for vars
	X        = xinfo->Xorg + N * pidx0;
	knotidx  = xinfo->SortedUnqiueIdx0[N * pidx0 + term->knots[1] - 1L];
	xknot[1] = X[knotidx];
	xcols[1] = X;	
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

static void GenOneTerm_step(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) { 
	
	F32     knot ;
	F32PTR  X   = GetKnot1D(term, xinfo,&knot );
	I32     N   = xinfo->N;
	
	if (term->side > 0)
		f32_step_pos(X, xinfo->Xones, xcol, knot, N);
	else
		f32_step_neg(X, xinfo->Xones, xcol, knot, N);	

	f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
}

static void GenOneTerm_stair(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {

	int     N       = xinfo->N;
	int     pidx0   = xinfo->xi_Pindices0[term->type][term->var0];  // // zero-based for basisState but 1-based for vars
	int     Nunqiue = xinfo->Nunique[pidx0];
	F32PTR  X       = xinfo->Xorg + N * pidx0;

 	int k0 = term->knot;
	int kU = term->knotnext;

	F32    x0 = *(X + xinfo->SortedUnqiueIdx0[N * pidx0 + k0 - 1L]);
	F32    xup = *(X + xinfo->SortedUnqiueIdx0[N * pidx0 + kU - 1L]);

	
	F32 X0     = (k0 == 1) ? -1.e300      :  x0;
	F32 XUPPER = (kU == Nunqiue) ? 1.e300 : xup;
 
	for (int i = 0; i < N; i++) {
  	   xcol[i] = X[i] > X0 && X[i] <= XUPPER ? 1.f : 0.f;
	}
 
	F32 dotProduct = DOT(N, xcol, xcol);
	F32 scale = 1 / sqrtf(dotProduct / N);
	f32_mul_val_inplace(scale, xcol, N);
	term->mean = 0;
	term->sd = 1 / scale;

	//f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
}
static void GenOneTerm_piece(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {

	int     N       = xinfo->N;
	int     pidx0   = xinfo->xi_Pindices0[term->type][term->var0];  // // zero-based for basisState but 1-based for vars
	int     Nunqiue = xinfo->Nunique[pidx0];
	F32PTR  X       = xinfo->Xorg + N * pidx0;

	int kL = term->knotprev; 
	int k0 = term->knot;
	int kU = term->knotnext;
 
	F32    xlow  = *(X + xinfo->SortedUnqiueIdx0[N * pidx0 + kL - 1L]);
	F32    x0   = *(X + xinfo->SortedUnqiueIdx0[N * pidx0 + k0 - 1L]);
	F32    xup  = *(X + xinfo->SortedUnqiueIdx0[N * pidx0 + kU - 1L]);

	F32 slpL = 1. / (x0 - xlow);
	F32 slpU = 1. / (x0 - xup);
	F32 aL   = -xlow *slpL;
	F32 aU   = -xup * slpU;

	F32 XLOWER = (kL ==  1) ?      -1.e300 : xlow;
	F32 XUPPER = (kU == Nunqiue) ? 1.e300 : xup;
	if (k0 == 1) {		
		for (int i = 0; i < N; i++) {
			xcol[i] = X[i] <= XUPPER ? X[i] * slpU + aU : 0.f;
		}	
	} else {
		for (int i = 0; i < N; i++) {
			if (X[i] > XUPPER || X[i] < XLOWER) {
				xcol[i] = 0;
			} else {
				xcol[i] = X[i] < x0 ? X[i] * slpL + aL : X[i] * slpU + aU;
			}			
		}
	}

	F32 dotProduct = DOT(N, xcol, xcol);
	F32 scale      = 1 / sqrtf(dotProduct / N);
	f32_mul_val_inplace(scale, xcol, N);	
	term->mean = 0;
	term->sd   = 1 / scale;

 
	//f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
}

static void GenOneTerm_season(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {

	int     N       = xinfo->N;
	int     pidx0   = xinfo->xi_Pindices0[term->type][term->var0];  // // zero-based for basisState but 1-based for vars
	int     Nunqiue = xinfo->Nunique[pidx0];
	F32PTR  X       = xinfo->Xorg + N * pidx0;

	int k0 = term->knot;
	int kU = term->knotnext;

	F32    x0 = *(X + xinfo->SortedUnqiueIdx0[N * pidx0 + k0 - 1L]);
	F32    xup = *(X + xinfo->SortedUnqiueIdx0[N * pidx0 + kU - 1L]);

	F32 X0      = (k0 == 1) ?      -1.e300 : x0;
	F32 XUPPER  = (kU == Nunqiue) ? 1.e300 : xup;

	int    pidx0y = xinfo->xi_PseasonY0[term->var0];  //  1-based  
	F32PTR Y      = xinfo->Xorg + N * pidx0y;

	for (int i = 0; i < N; i++) {
		xcol[i] = X[i] > X0 && X[i] <= XUPPER ? Y[i] : 0.f;
	}

	F32 dotProduct = DOT(N, xcol, xcol);
	F32 scale = 1 / sqrtf(dotProduct / N);
	f32_mul_val_inplace(scale, xcol, N);
	term->mean = 0;
	term->sd = 1 / scale;
}

static void GenOneTerm_hinge(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {

	F32     knot ;
	F32PTR  X = GetKnot1D(term, xinfo, &knot);
	I32     N = xinfo->N;
	if ( term->side > 0) //term->sides[0]>0
		f32_hinge_pos(X, xcol,knot, N);
	else
		f32_hinge_neg(X, xcol, knot, N);

	//F32 dotProduct = DOT(N, xcol, xcol);
	//F32 scale      = 1 / sqrtf(dotProduct / N);
	//f32_mul_val_inplace(scale, xcol, N);
	//f32_mul_val_inplace(scale, xcol, N);
	//term->mean = 0;
	//term->sd   = 1 / scale;

	f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
}

static void GenOneTerm_hingepair(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {


	F32     knot;
	F32PTR  X = GetKnot1D(term, xinfo, &knot);
	I32     N = xinfo->N;
	if (term->side > 0) //term->sides[0]>0
		f32_hinge_pos(X, xcol, knot, N);
	else
		f32_hinge_neg(X, xcol, knot, N);

	f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
}

static void GenOneTerm_change(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {

	F32     knot;
	F32PTR  X = GetKnot1D(term, xinfo, &knot);

	int    vidx0  = term->var0;
	int    pidx0y = xinfo->xi_PchangeY0[vidx0];  //  1-based 
	int    N	 = xinfo->N;
	F32PTR Y	 = xinfo->Xorg + N * pidx0y;
 
	if (term->side > 0)
		f32_step_pos(X, Y, xcol, knot, N);
	else
		f32_step_neg(X, Y, xcol, knot, N);
	//f32_normalize_inplace(xcol, N);

	//f32_normalize_std_avg_inplace(xcol, xinfo->N, &term->mean, &term->sd);

	F32 dotProduct = DOT(N, xcol, xcol);
	F32 scale       = 1 / sqrtf(dotProduct / N);
	f32_mul_val_inplace(scale, xcol, N);
	term->mean = 0;
	term->sd  = 1 / scale; 
}

static void GenOneTerm_hinge2d_core(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {

	int      N        = xinfo->N;
	F32      knots[2];
	F32PTR   Xdata[2];
	GetKnot2D(term, xinfo, knots, Xdata);
	
 
	if (term->sides[0] > 0) //term->sides[0]>0
		f32_hinge_pos(Xdata[0], xcol, knots[0], N);
	else
		f32_hinge_neg(Xdata[0], xcol, knots[0], N);
	 
	if (term->sides[1] > 0) //term->sides[0]>0
		f32_hinge_pos(Xdata[1], xinfo->mem_GenOneSeg_hinege2d_2dr_1dr, knots[1], N);
	else
		f32_hinge_neg(Xdata[1], xinfo->mem_GenOneSeg_hinege2d_2dr_1dr, knots[1], N);

	f32_mul_vec_inplace(xinfo->mem_GenOneSeg_hinege2d_2dr_1dr, xcol, N);
	
}
static void GenOneTerm_hinge2d(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {
	GenOneTerm_hinge2d_core( term, xinfo,  xcol);
	f32_normalize_std_avg_inplace(xcol, xinfo->N, &term->mean, &term->sd);	
}


static void GenOneTerm_hinge2dr_core(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {
	int    N = xinfo->N;

	F32      knots[2];
	F32PTR   Xdata[2];
	GetKnot2D(term, xinfo, knots, Xdata);

	F32PTR t1 = xinfo->mem_GenOneSeg_hinege2d_2dr_1dr;
	F32PTR t2 = t1 + N;
	F32PTR t3 = t2 + N;

	F32   a0 = xinfo->anglecoff[term->sides[0]];
	F32   b0 = xinfo->anglecoff[term->sides[0]+xinfo->Nangle];

	F32   a1 = xinfo->anglecoff[term->sides[1]];
	F32   b1 = xinfo->anglecoff[term->sides[1] + xinfo->Nangle];

	f32_copy(Xdata[0], t1, N); 	f32_scale_inplace(a0, -a0 * knots[0], t1, N);
	f32_copy(Xdata[1], t2, N);	f32_scale_inplace(b0, -b0 * knots[1], t2, N);
	f32_add_vec_inplace(t2, t1, N);
	f32_hinge_pos(t1, t1,0, N);

 	
	f32_copy(Xdata[0], t2, N); 	f32_scale_inplace(a1, -a1 * knots[0], t2, N);
	f32_copy(Xdata[1], t3, N);	f32_scale_inplace(b1, -b1 * knots[1], t3, N);
	f32_add_vec_inplace(t3, t2, N);
	f32_hinge_pos(t2, xcol, 0, N);

	f32_mul_vec_inplace(t1, xcol, N);	
}
static void GenOneTerm_hinge2dr(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {
	GenOneTerm_hinge2dr_core( term, xinfo,  xcol);
	f32_normalize_std_avg_inplace(xcol, xinfo->N, &term->mean, &term->sd);	
}


static void GenOneTerm_hinge1dr_core(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {
	int    N        = xinfo->N;

	F32      knots[2];
	F32PTR   Xdata[2];
	GetKnot2D(term, xinfo, knots, Xdata);

	F32PTR t1 = xinfo->mem_GenOneSeg_hinege2d_2dr_1dr;
 
	F32   a0 = xinfo->anglecoff[term->sides[0]];
	F32   b0 = xinfo->anglecoff[term->sides[0]+xinfo->Nangle];

	f32_copy(Xdata[0], t1, N);   	f32_scale_inplace(a0, -a0 * knots[0], t1, N);
	f32_copy(Xdata[1], xcol, N);	f32_scale_inplace(b0, -b0 * knots[1], xcol, N);
	f32_add_vec_inplace(t1, xcol, N);
	f32_hinge_pos(xcol, xcol,0, N);	
}
static void GenOneTerm_hinge1dr(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {
	GenOneTerm_hinge1dr_core(term, xinfo, xcol);
	f32_normalize_std_avg_inplace(xcol, xinfo->N, &term->mean, &term->sd);
}


typedef void (*pfnGenOneTerm)(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol);

static pfnGenOneTerm __funcs[] = {
	GenOneTerm_const,
	GenOneTerm_linear,
	GenOneTerm_quadratic,
	GenOneTerm_step,
	GenOneTerm_stair,
	GenOneTerm_piece,
	GenOneTerm_season,
	GenOneTerm_hinge,
	GenOneTerm_hingepair,
	GenOneTerm_change,
	GenOneTerm_hinge2d,
	GenOneTerm_hinge2dr,
	GenOneTerm_hinge1dr,
};

static pfnGenOneTerm  *funcs = &__funcs[1];

void GenTerms(BVS_TERMS_PTR terms, BVS_XINFO_PTR xinfo, F32PTR xcols,I32 K) {
 
	int N = xinfo->N;
	for (int i = 0; i < K; ++i) {
		funcs[terms[i].type](terms + i, xinfo, xcols + i * N);
	}	 
}

int GenOneTerm_hinge2d_GetGoodDataNpts(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {
	GenOneTerm_hinge2d_core(term, xinfo, xcol);
	I32PTR indices = xinfo->mem_GenOneSeg_hinege2d_2dr_1dr;
	return f32_findindex(xcol, indices, 0, xinfo->N, CMP_GT);
}
int GenOneTerm_hinge2dr_GetGoodDataNpts(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {
	GenOneTerm_hinge2dr_core(term, xinfo, xcol);
	I32PTR indices = xinfo->mem_GenOneSeg_hinege2d_2dr_1dr;
	return f32_findindex(xcol, indices, 0, xinfo->N, CMP_GT);
}

int GenOneTerm_hinge1dr_GetGoodDataNpts(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol) {
	GenOneTerm_hinge1dr_core(term, xinfo, xcol);
	I32PTR indices = xinfo->mem_GenOneSeg_hinege2d_2dr_1dr;
	return f32_findindex(xcol, indices, 0, xinfo->N, CMP_GT);
}


#include "abc_000_warning.h"