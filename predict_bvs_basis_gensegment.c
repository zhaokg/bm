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

static void GenOneTerm_const(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol) {
	r_ippsSet_32f(1.0, xcol, newxinfo->N);	
}

static void GenOneTerm_linear(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol) {
	int vidx0 = term->var0;
	//int pidx0 = xinfo->xi_Pindices0[LinearTerm][term->vars[0] - 1L];  //  1-based 
	int Nnew  = newxinfo->N;
 	r_cblas_scopy(Nnew, newxinfo->Xnorm1 + Nnew * vidx0, 1L, xcol, 1L);
	//term->mean = xinfo->mean[var0];
   //term->sd   = xinfo->sd[var0];	
}
static void GenOneTerm_quadratic(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol) {
	int vidx0 = term->var0;
	//int pidx0 = xinfo->xi_Pindices0[QuadraticTerm][term->vars[0] - 1L];  //  1-based 
	int Nnew  = newxinfo->N;

	if (newxinfo->X2) {
		r_cblas_scopy(Nnew, newxinfo->X2 + Nnew * vidx0, 1L, xcol, 1L);
		//term->mean = xinfo->mean2[var0];
		//term->sd   = xinfo->sd2[var0];
	} else {
		//r_cblas_scopy(N, xinfo->Xorg + N * pidx0, 1L, xcol, 1L);
		//f32_mul_vec_inplace(xcol, xcol, N);
		//f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
	} 	
}


/////////////////////////////////////////////////////////////////////??????????????????????????????
static F32PTR GetKnot1DXnew(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR knot) {
	int    N       = xinfo->N;
	int    pidx0   = xinfo->xi_Pindices0[term->type][ term->var0 ];   
	F32PTR X       = xinfo->Xorg + N *pidx0;
	int    knotidx = xinfo->SortedUnqiueIdx0[N * pidx0 + term->knot - 1L];
	knot[0]        = X[knotidx];
	int    Nnew = newxinfo->N;
	F32PTR Xnew = newxinfo->Xorg + Nnew * pidx0;
    return Xnew;
}
static void GetKnot2DXnew(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR knots, F32PTR * Xnewdata) {
	int    N       = xinfo->N;
	int    pidx0   = xinfo->xi_Pindices0[term->type][ term->var0s[0]];  //  1-based 	
	F32PTR X       = xinfo->Xorg + N *pidx0;
	int    knotidx = xinfo->SortedUnqiueIdx0[N * pidx0 + term->knots[0] - 1L];
	knots[0]        = X[knotidx];

	int    Nnew = newxinfo->N;
	F32PTR Xnew = newxinfo->Xorg + Nnew * pidx0;
	Xnewdata[0] = Xnew;

	pidx0   = xinfo->xi_Pindices0[term->type][ term->var0s[1]];  //  1-based 	
	X       = xinfo->Xorg + N *pidx0;
	knotidx = xinfo->SortedUnqiueIdx0[N * pidx0 + term->knots[1] - 1L];
	knots[1]  = X[knotidx];

	Nnew = newxinfo->N;
	Xnew = newxinfo->Xorg + Nnew * pidx0;
	Xnewdata[1] = Xnew;
}
/////////////////////////////////////////////////////////////////////??????????????????????????????

static void GenOneTerm_step(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol) {

	F32    knot;
	F32PTR Xnew = GetKnot1DXnew(term, xinfo, newxinfo,&knot);
	 
	int    Nnew = newxinfo->N;
	if (term->side > 0)
		f32_step_pos(Xnew, newxinfo->Xones, xcol, knot, Nnew);
	else
		f32_step_neg(Xnew, newxinfo->Xones, xcol, knot, Nnew);

	f32_scale_inplace(1. / term->sd, -term->mean / term->sd, xcol, Nnew);
	//f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
}
static void GenOneTerm_stair(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol) {

	int     N       = xinfo->N;
	int     pidx0   = xinfo->xi_Pindices0[term->type][term->var0];  // // zero-based for basisState but 1-based for vars
	int     Nunqiue = xinfo->Nunique[pidx0];
	F32PTR  X       = xinfo->Xorg + N * pidx0;

 	int k0 = term->knot;
	int kU = term->knotnext;
	F32 x0  = *(X + xinfo->SortedUnqiueIdx0[N * pidx0 + k0 - 1L]);
	F32 xup = *(X + xinfo->SortedUnqiueIdx0[N * pidx0 + kU - 1L]);
	
	F32 X0     = (k0 == 1) ? -1.e300      :  x0;
	F32 XUPPER = (kU == Nunqiue) ? 1.e300 : xup;
 

	int    Nnew = newxinfo->N;
	F32PTR Xnew = newxinfo->Xorg + Nnew * pidx0;

	for (int i = 0; i < N; i++) {
  	   xcol[i] = X[i] > X0 && X[i] <= XUPPER ? 1.f : 0.f;
	}
 
	f32_scale_inplace(1 / term->sd, -term->mean / term->sd, xcol, Nnew);

	//f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
}
static void GenOneTerm_piece(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol) {

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

	int    Nnew = newxinfo->N;
	F32PTR Xnew = newxinfo->Xorg + Nnew * pidx0;

	if (k0 == 1) {		
		for (int i = 0; i < Nnew; i++) {
			xcol[i] = Xnew[i] <= XUPPER ? Xnew[i] * slpU + aU : 0.f;
		}	
	} else {
		for (int i = 0; i < Nnew; i++) {
			if (Xnew[i] > XUPPER || Xnew[i] < XLOWER) {
				xcol[i] = 0;
			} else {
				xcol[i] = Xnew[i] < x0 ? Xnew[i] * slpL + aL : Xnew[i] * slpU + aU;
			}			
		}
	}

	f32_scale_inplace(1 / term->sd, -term->mean / term->sd, xcol, Nnew);
	//f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
}


static void GenOneTerm_hinge(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol) {

	F32    knot;
	F32PTR Xnew = GetKnot1DXnew(term, xinfo, newxinfo, &knot);

	int    Nnew = newxinfo->N;
	if ( term->side>0) //term->sides[0]>0
		f32_hinge_pos(Xnew, xcol,knot, Nnew);
	else
		f32_hinge_neg(Xnew, xcol, knot, Nnew);
	//f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
	f32_scale_inplace(1 / term->sd, -term->mean / term->sd, xcol, Nnew);
}


static void GenOneTerm_change(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol) {

	F32    knot;
	F32PTR Xnew = GetKnot1DXnew(term, xinfo, newxinfo, &knot);

	int    vidx0  = term->var0;
	int    pidx0y = xinfo->xi_PchangeY0[vidx0];  //  1-based   
	int    Nnew   = newxinfo->N; 
	F32PTR Ynew   = newxinfo->Xnorm1 + Nnew * pidx0y;

	if (term->side > 0)
		f32_step_pos(Xnew, Ynew, xcol, knot, Nnew);
	else
		f32_step_neg(Xnew, Ynew, xcol, knot, Nnew);
	//f32_normalize_inplace(xcol, N);
	//term->mean = 0;
	//term->sd   = 1;
	f32_scale_inplace(1.f / term->sd, -term->mean / term->sd, xcol, Nnew);
}

static void GenOneTerm_hinge2d(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol) {
	
	F32    knots[2];
	F32PTR Xnewdata[2];
	GetKnot2DXnew(term, xinfo, newxinfo, knots, Xnewdata);
	
	int    Nnew = newxinfo->N;
	if (term->sides[0] > 0) //term->sides[0]>0
		f32_hinge_pos(Xnewdata[0], xcol, knots[0], Nnew);
	else
		f32_hinge_neg(Xnewdata[0], xcol, knots[0], Nnew);


	if (term->sides[1] > 0) //term->sides[0]>0
		f32_hinge_pos(Xnewdata[1], xinfo->mem_GenOneSeg_hinege2d_2dr_1dr, knots[1], Nnew);
	else
		f32_hinge_neg(Xnewdata[1], xinfo->mem_GenOneSeg_hinege2d_2dr_1dr, knots[1], Nnew);


	f32_mul_vec_inplace(xinfo->mem_GenOneSeg_hinege2d_2dr_1dr, xcol, Nnew);
	f32_scale_inplace(1 / term->sd, -term->mean / term->sd, xcol, Nnew);
	//f32_normalize_std_avg_inplace(xcol, N, &term->mean, &term->sd);
	 
}
static void GenOneTerm_hinge2dr(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol) {
	F32    knots[2];
	F32PTR Xnewdata[2];
	GetKnot2DXnew(term, xinfo, newxinfo, knots, Xnewdata);

	int    Nnew = newxinfo->N;
	F32PTR t1 = xinfo->mem_GenOneSeg_hinege2d_2dr_1dr;
	F32PTR t2 = t1 + Nnew;
	F32PTR t3 = t2 + Nnew;

	F32   a0 = xinfo->anglecoff[term->sides[0]];
	F32   b0 = xinfo->anglecoff[term->sides[0]+xinfo->Nangle];

	F32   a1 = xinfo->anglecoff[term->sides[1]];
	F32   b1 = xinfo->anglecoff[term->sides[1] + xinfo->Nangle];

	f32_copy(Xnewdata[0], t1, Nnew);
	f32_scale_inplace(a0, -a0 * knots[0], t1, Nnew);

	f32_copy(Xnewdata[1], t2, Nnew);
	f32_scale_inplace(b0, -b0 * knots[1], t2, Nnew);
	f32_add_vec_inplace(t2, t1, Nnew);
	f32_hinge_pos(t1, t1,0, Nnew);

 	
	f32_copy(Xnewdata[0], t2, Nnew);
	f32_scale_inplace(a1, -a1 * knots[0], t2, Nnew);

	f32_copy(Xnewdata[1], t3, Nnew);
	f32_scale_inplace(b1, -b1 * knots[1], t3, Nnew);
	f32_add_vec_inplace(t3, t2, Nnew);
	f32_hinge_pos(t2, xcol, 0, Nnew);

	f32_mul_vec_inplace(t1, xcol, Nnew);
	f32_scale_inplace(1 / term->sd, -term->mean / term->sd, xcol, Nnew);
 
}
static void GenOneTerm_hinge1dr(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol) {

	F32    knots[2];
	F32PTR Xnewdata[2];
	GetKnot2DXnew(term, xinfo, newxinfo, knots, Xnewdata);

	int    Nnew = newxinfo->N;

	F32PTR t1 = xinfo->mem_GenOneSeg_hinege2d_2dr_1dr;
 
	F32   a0 = xinfo->anglecoff[term->sides[0]];
	F32   b0 = xinfo->anglecoff[term->sides[0]+xinfo->Nangle];
 

	f32_copy(Xnewdata[0], t1, Nnew);
	f32_scale_inplace(a0, -a0 * knots[0], t1, Nnew);

	f32_copy(Xnewdata[1], xcol, Nnew);
	f32_scale_inplace(b0, -b0 * knots[1], xcol, Nnew);
	f32_add_vec_inplace(t1, xcol, Nnew);
	f32_hinge_pos(xcol, xcol,0, Nnew);

 	f32_scale_inplace(1 / term->sd, -term->mean / term->sd, xcol, Nnew);

}

typedef void (*pfnGenOneTerm_Predict)(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcol);

static pfnGenOneTerm_Predict __funcs[] = {
	GenOneTerm_const,
	GenOneTerm_linear,
	GenOneTerm_quadratic,
	GenOneTerm_step,
	GenOneTerm_stair,
	GenOneTerm_piece,
	GenOneTerm_hinge,
	GenOneTerm_hinge, // for HingePair
	GenOneTerm_change,
	GenOneTerm_hinge2d,
	GenOneTerm_hinge2dr,
	GenOneTerm_hinge1dr,
};

static pfnGenOneTerm_Predict*funcs = &__funcs[1];

void GenTermsPredict(BVS_TERMS_PTR terms, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcols,I32 K) {
 
	int N = xinfo->N;
	for (int i = 0; i < K; ++i) {
		funcs[terms[i].type](terms + i, xinfo, newxinfo, xcols + i * N);
	}	 
}
#include "abc_000_warning.h"