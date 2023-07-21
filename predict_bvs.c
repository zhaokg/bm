#include "abc_000_macro.h"
#include "abc_000_warning.h"

#if defined(MSVC_COMPILER)
#include "intrin.h"                //_rdstc
#endif

#include <stdio.h>	               //fprintf fopen FILE
#include <string.h>	               //memset memcpy
#include <time.h>
#include <math.h>

#include "abc_001_config.h"
#include "abc_mem.h"              // Independent of R/Matlab,  VC/GNU, or MY/MKL.
#include "abc_blas_lapack_lib.h"  // Slight dependence on the choice of VC/GNU. Dependence on MY/MKL. Independacne of R/Matlab.
#include "abc_ide_util.h"
#include "abc_ts_func.h"
#include "abc_timer.h"
#include "abc_mcmc.h"
#include "abc_mat.h"
#include "abc_rand.h"
#include "abc_vec.h"   // for f32_add_v_v2_vec_in_place, f32_diff_back,i32_increment_bycon_inplace i32_to_f32_scaelby_inplace, f32_sx_sxx_toavstd_inplace 
#include "abc_math.h"  // for fastexp, fastsqrt only

#include "globalvars.h"  

#include "bvs_header.h"
#include "bvs_io.h" 
#include "bvs_fun.h" 

#define LOCAL(...) do{ __VA_ARGS__ } while(0);
 
 
int bvs_predict_corev4(BVS_XINFO_PTR  xinfo, BVS_TERMS_ENCODER_PTR coder, 
	BVS_PRED_XINFO * newxinfo, BVS_PRIOR_PTR  prior, BVS_PRED_OPT *newopt){

	// A struct to track allocated pointers   
	// const MemPointers MEM;
	// do not use 'const' bcz Clang will asssume other fields as zeros (e.g., alloc, and alloc0).
	MemPointers MEM = (MemPointers){.init = mem_init,};
	MEM.init(&MEM);
	MEM.checkHeader = 1;
	//mem = &MEM;

	// Get the Option parameters from the global pointer GLOBAL_OPTIONS, which was already 
	// set prior to this point (e.g., time-series length N, and number of pixels).	
 
	
	// Pre-allocate memory to save samples for calculating credibile intervals	
	CI_PARAM     ciParam = {0,};
	CI_RESULT    *ci     ;
 

	// Allocate MEMORY FOR BASIS VARIABLE: Initialzie two pointers to BASIS

	
	// Allocate mem for current covariate/design matrix (Xt_mars), proposed new terms (Xnewterm),
	// and subset matrix corresponding to rows of missing values.		
	int  N      = newxinfo->N;
	int  Ktrue  = coder->Kmax_actual;
	 F32PTR Xcur= MyALLOC(MEM, N*Ktrue, U32, 64);
	 F32PTR Xpre= MyALLOC(MEM, N*Ktrue, U32, 64);	
	
	xinfo->mem_GenOneSeg_hinege2d_2dr_1dr = MyALLOC(MEM, N*3,  U32, 64); // A temporary mem used in GenOneSeg_hinege2d
	coder->PreList                = MyALLOC(MEM, Ktrue, BVS_TERMS, 64); // Used in ExtractTerms
	coder->mem                    = MyALLOC(MEM, Ktrue+Ktrue+1L, I16, 64); // Used in EnocdeTremList
	// Allocate the output memory for a single chain (resultChain) and the averaged
	  
	int        P = xinfo->Pall;
	BVS_PRED_RESULT result;
	result.P    = MyALLOC0(MEM,  P, F32, 64); // A temporary mem used in GenOneSeg_hinege2d
	result.Y    =  MyALLOC0(MEM, N, F32, 64); // A temporary mem used in GenOneSeg_hinege2d
	result.SD   =  MyALLOC0(MEM, N, F32, 64); // A temporary mem used in GenOneSeg_hinege2d
	//result.Ychains = MyALLOC0(MEM, N*coder->nModels, F32, 64); 
	// Do not allocate this huge mem; instead use the acllocated from IDE
	F32PTR  Ychains= newopt->bComputeYchains? newxinfo->mat->Ychains: NULL;
	
	result.Yparts     = MyALLOC0(MEM, N*P, F32, 64); // A temporary mem used in GenOneSeg_hinege2d
	result.BetaLinear = MyALLOC0(MEM, P, F32, 64); // A temporary mem used in GenOneSeg_hinege2d

#define  hasTerm(x) (prior->type2id[x] >= 0)

	I32 maxCptNumHinge = prior->maxKnotNum_Max[HingeTerm];
	if ( hasTerm(HingeTerm)  ) {
		int  nelem =  (1+ maxCptNumHinge)*xinfo->xi_Pbasis[HingeTerm];
		result.Khinge = MyALLOC0(MEM, nelem, I32, 64); // A temporary mem used in GenOneSeg_hinege2d
	}
  
	VOIDPTR buf = MyALLOC0(MEM, P*N, F32, 64); // A temporary mem used in GenOneSeg_hinege2d
	/****************************************************************************/
	//		THE OUTERMOST LOOP: Loop through all the time series one by one
	/****************************************************************************/
	// Get conversion factors from counts to seceonds
 
	BVS_TERMS *out = MyALLOC(MEM, Ktrue, BVS_TERMS, 64);

	F64  beta_const = 0;

	coder->curBetaID  = 0;
	coder->curModelID = 0;
	coder->curTermID  = 0;
	coder->Kpre       = 0;
	for (int idModel = 0; idModel < coder->nModels; idModel++) {

		int					 Kcode     = coder->KcodeVec[coder->curModelID];
		BVS_TERMS_ANY_PTR    termscode = &coder->TERMS[coder->curTermID];
		F32PTR               beta      = &coder->BETA[coder->curBetaID];
		int Kout                       = ExtractNextTermFromList(coder, out);

		beta_const += beta[0];
		beta++;         //move forward to skip the consts


		BVS_TERMS* terms = out;
		terms--;       // move backward one step to have a ghost const term at index=0
		{
			I08PTR alreadyIncluded = buf;
			memset(alreadyIncluded, 0, P);

			for (I32 j = 1; j < (Kout + 1); j++) {
				I32 vidx0 = terms[j].var0s[0];
				I32 pidx0 = xinfo->xi_Pindices0[terms[j].type][vidx0];
				if (!alreadyIncluded[pidx0]) {
					result.P[pidx0]++;
					alreadyIncluded[pidx0] = 1;
				}
				if (terms[j].type == Hinge2DTerm || terms[j].type == Hinge2DRTerm || terms[j].type == Hinge1DRTerm) {
					I32 vidx1 = terms[j].var0s[1];
					I32 pidx1 = xinfo->xi_Pindices0[terms[j].type][vidx1];
					if (!alreadyIncluded[pidx1]) {
						result.P[pidx1]++;
						alreadyIncluded[pidx1] = 1;
					}
				}
			}


			if (hasTerm(LinearTerm)) {
				I32PTR Pindices0 = xinfo->xi_Pindices0[LinearTerm];
				for (I32 j = 1; j < (Kout + 1); j++) {
					if (terms[j].type == LinearTerm) {
						I32 vidx0 = terms[j].var0;
						I32 pidx0 = Pindices0[vidx0];
						result.BetaLinear[pidx0] += beta[j - 1];
					}					
				}
			}

		}
		///////////////////////////////////////////////////////////////////////////////////////////////////
		int Kterms = 0;
		for (int i = 0; i < Kcode; i++) {
			if (termscode[i].type >= 0) {
				GenTermsPredict(&termscode[i].tm, xinfo, newxinfo, Xcur + Kterms * N, 1);
				//GenTerms(&termscode[i].tm, xinfo, Xcur + Kterms * N,  1);
				Kterms++;
			}
			else if (termscode[i].type == CMPLIST_RANGE) {
				for (int j = 0; j < termscode[i].num; j++) {
					int s = termscode[i].range[j].s;
					int e = termscode[i].range[j].e;
					int kse = e - s + 1;
					memcpy(Xcur + Kterms * N, Xpre + (s - 1) * N, sizeof(F32) * kse * N);
					Kterms += kse;
				}
			}
			else {
				assert(0);
			}
		}
		assert(Kterms == Kout);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (Kterms > 0) {
			r_cblas_sgemv(CblasColMajor, CblasNoTrans, N, Kterms, 1.f, Xcur, N, beta, 1L, 0.f, buf, 1L);
			f32_add_v_v2_vec_inplace(buf, result.Y, result.SD, N);		
		}

		if (Ychains) {
			f32_copy(buf, Ychains + N * idModel, N); //TODO: buf must be pre-filled with zeros if Keterms==0
		}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (hasTerm(HingeTerm)) {		
			I32PTR count = buf;
			int    Phinge = xinfo->xi_Pbasis[HingeTerm];
			memset(count, 0, Phinge * sizeof(I32));

			for (I32 j = 1; j < (Kout + 1); j++) {
				if (terms[j].type != HingeTerm) continue;
				I32 vidx0 = terms[j].var0s[0];
				count[vidx0]++;				
			}

			for (int i = 0; i < Phinge; i++) {
				result.Khinge[(1 + maxCptNumHinge) * i + count[i]]++;
			}		

		}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		for (I32 i = 1; i <= Kterms; i++) {
			int type = terms[i].type;
			if (type == Hinge2DTerm || type == Hinge2DRTerm || type == Hinge1DRTerm) {
				continue;
			}
			I32 vidx0  = terms[i].var0s[0];
			I32 pidx0  = xinfo->xi_Pindices0[type][vidx0];

			f32_axpy_inplace(beta[i-1L], Xcur + N * (i-1), result.Yparts + pidx0 * N, N);
		}
		/*
		for (I32 i = 0; i < Kterms; i++) {
			if (MODEL.terms[i].type == Hinge2DTerm) {
				continue;
			}
			I32 pidx0 = MODEL.terms[i].vars[0] - 1L;
			f32_axpy_inplace(BETA[i], Xt_mars + N * i, resultChain.Ys + pidx0 * N, N);
		}
		*/
		VOIDPTR tmp;
		tmp  = Xcur;
		Xcur = Xpre;
		Xpre = tmp;
	}

	/**/
 	///////////////////////////////////////////////////////////////////////////////////
	i32_to_f32_scaleby_inplace(result.P, P, 1.0/coder->nModels);
	if (hasTerm(HingeTerm)) {
		i32_to_f32_scaleby_inplace(result.Khinge,(1+maxCptNumHinge)* xinfo->xi_Pbasis[HingeTerm],	1.0 / coder->nModels);
	}
	 
	F32 Ystd = prior->alpha2;
	F32 Ymean = prior->alpha1;
	beta_const = beta_const / coder->nModels,
	f32_sx_sxx_to_avgstd_inplace(result.Y, result.SD, coder->nModels, Ystd, Ymean, N); 
	f32_add_val_inplace(Ystd* beta_const, result.Y, N);
 
	f32_mul_val_inplace(1.0 / coder->nModels, result.BetaLinear, P);
	if (prior->type2id[LinearTerm] >= 0) {
		for (I32 j = 0; j < xinfo->Pall; j++) {
 			result.BetaLinear[j] = result.BetaLinear[j] /xinfo->sd[j]* Ystd;
		}
	}

	f32_mul_val_inplace(Ystd/ coder->nModels, result.Yparts, N*P);

	if (Ychains) {
		f32_scale_inplace(Ystd, Ymean + (Ystd * beta_const), Ychains, N*coder->nModels);
		if (newopt->odtype == DATA_DOUBLE) {
			f32_to_f64_inplace(Ychains, N* coder->nModels); // Ychains points to mat.Ychains
		}		
	} 
	
	BVS_PRED_RESULT_PTR mat=newxinfo->mat;
	int len = P;
	f32_to_strided_mem(result.P, mat[0].P, len, 1, 0, DATA_DOUBLE);
	if (hasTerm(HingeTerm)) {
		len = (1 + maxCptNumHinge) * xinfo->xi_Pbasis[HingeTerm] ;
		f32_to_strided_mem(result.Khinge, mat[0].Khinge, len, 1, 0, DATA_DOUBLE);
	}

	len = N;
	f32_to_strided_mem(result.Y, mat[0].Y, len, 1, 0, DATA_DOUBLE);
	len = N;
	f32_to_strided_mem(result.SD, mat[0].SD, len, 1, 0, DATA_DOUBLE);


	len = P;
	f32_to_strided_mem(result.BetaLinear, mat[0].BetaLinear, len,1, 0, DATA_DOUBLE);

	len = P*N;
	f32_to_strided_mem(result.Yparts, mat[0].Yparts, len, 1, 0, DATA_DOUBLE);
 
 /*
 
		{
			F32PTR MEMBUF1 = Xnewterm;

			//Counting probability of being breakpoints	




			if (MODEL.type2id[QuadraticTerm] >= 0) {
				I32    K = MODEL.basisState[MODEL.type2id[QuadraticTerm]].K;
				I16PTR Kposition = MODEL.basisState[MODEL.type2id[QuadraticTerm]].Kposition;
				for (I32 i = 0; i < K; i++) {
					I32 pidx0 = MODEL.terms[Kposition[i] - 1L].vars[0] - 1L;
					resultChain.inProbSqr[pidx0]++;
				}

			}

			if (MODEL.type2id[StepTerm] >= 0) {
				I32    K = MODEL.basisState[MODEL.type2id[StepTerm]].K;
				I16PTR Kposition = MODEL.basisState[MODEL.type2id[StepTerm]].Kposition;
				for (I32 i = 0; i < K; i++) {
					I32 pidx0 = MODEL.terms[Kposition[i] - 1L].vars[0] - 1L;
					I32 knot0 = MODEL.terms[Kposition[i] - 1L].knots[0] - 1L;
					resultChain.probStep[pidx0 * N + knot0]++;
				}

			}
			if (MODEL.type2id[HingeTerm] >= 0) {
				I32    K = MODEL.basisState[MODEL.type2id[HingeTerm]].K;
				I16PTR Kposition = MODEL.basisState[MODEL.type2id[HingeTerm]].Kposition;
				for (I32 i = 0; i < K; i++) {
					I32 pidx0 = MODEL.terms[Kposition[i] - 1L].vars[0] - 1L;
					I32 knot0 = MODEL.terms[Kposition[i] - 1L].knots[0] - 1L;
					resultChain.probHinge[pidx0 * N + knot0]++;
				}
			}

			if (MODEL.type2id[ChangeTerm] >= 0) {
				I32    K = MODEL.basisState[MODEL.type2id[ChangeTerm]].K;
				I16PTR Kposition = MODEL.basisState[MODEL.type2id[ChangeTerm]].Kposition;
				for (I32 i = 0; i < K; i++) {
					I32 pidx0 = MODEL.terms[Kposition[i] - 1L].vars[0] - 1L;
					I32 knot0 = MODEL.terms[Kposition[i] - 1L].knots[0] - 1L;
					resultChain.probChange[pidx0 * N + knot0]++;
				}
			}

			I32 K = MODEL.Kterms;

			resultChain.KProb[K - 1L]++;

			if (opt->prior.isBeast) {
				f32_axpy_inplace(BETA[0], Xt_mars, resultChain.Ystair, N);
				for (int i = 1; i < K; i++) {
					F32PTR Ydst;
					if (MODEL.terms[i].type == LinearTerm) {
						if (MODEL.terms[i].vars[0] == 1) { Ydst = resultChain.Ylin; }
						else { Ydst = resultChain.Yseason; }
					}
					if (MODEL.terms[i].type == StepTerm) { Ydst = resultChain.Ystair; }
					if (MODEL.terms[i].type == HingeTerm) { Ydst = resultChain.Ylin; }
					if (MODEL.terms[i].type == ChangeTerm) { Ydst = resultChain.Yseason; }
					f32_axpy_inplace(BETA[i], Xt_mars + N * i, Ydst, N);
				}

			}

 

		 

			//Compute the averaged  signals
			//r_cblas_sgemv(CblasColMajor, CblasNoTrans, Npad, K, 1.f, X, Npad,beta, 1L, 0.f,	Y, 1L);
			r_cblas_sgemv(CblasColMajor, CblasNoTrans, Npad, K, 1.f, Xt_mars, Npad, MODEL.curr.beta_mean, 1L, 0.f, MEMBUF1, 1L);
			//basis->ComputeY(Xt_mars, BETA, MEMBUF1, basis, Npad);
			f32_add_v_v2_vec_inplace(MEMBUF1, resultChain.Y, resultChain.SD, N);
			MEMBUF1 += Npad;
		}

		*/
 

  
	MEM.free_all(&MEM);
 
	return 1;
} /* End of beastST() */


#include "abc_000_warning.h"