#include <string.h>  //memset
#include <math.h>    //sqrt

#include "abc_000_warning.h"

#include "abc_vec.h"
#include "abc_ts_func.h"
#include "abc_mem.h"
#include "abc_vec.h"
#include "bvs_header.h"
#include "abc_ide_util.h"  //printf
#include "abc_common.h"  //strcicmp
#include "abc_blas_lapack_lib.h"  //strcicmp

static int _IsTwoTermsEqual(BVS_TERMS_PTR A, BVS_TERMS_PTR B) {
	#define a (*A)
	#define b (*B)
	
	// var0s[0] is aliased to var0 through a union

	if (a.type != b.type || a.var0s[0]!= b.var0s[0]) return 0;

	//type==type. &&  v0=v1
	I32 type = a.type;
	if (type == LinearTerm || type == QuadraticTerm) return 1;

	// Not LinearTerm or QuadraticTerm
	if (type == StepTerm || type == HingeTerm || type == ChangeTerm) {
		//assert(a.sides[0] == 99);
		return a.knot == b.knot && a.side == b.side;
	}
	if (type == PieceTerm || type == StairTerm || type == SeasonTerm) {
		return a.knotnext == b.knotnext && a.knot == b.knot;
	}
	if (type == Hinge2DTerm || type == Hinge1DRTerm|| type == Hinge2DRTerm) {
		return a.var0s[1] == b.var0s[1] && *((I64PTR)a.knots) == *((I64PTR)b.knots)
			&& *((I16PTR)a.sides) == *((I16PTR)b.sides);
	}

	assert(0);
	#undef a
	#undef b
}


int EncodeTermsList(BVS_TERMS_PTR pre, BVS_TERMS_PTR cur, BVS_TERMS_ANY_PTR code, int Npre, int Ncur, I16PTR mem) {
	// a+b=>c
	if (Ncur == 0) {
		// If the list is empty		
		int Ncode = 0;
		code[Ncode].type = CMPLIST_RANGE;
		code[Ncode].num = 0;
		Ncode++;
		return Ncode;
	}

	if (Npre == 0) {
		// If the prev list is empty, and then just copy the current list to c
		int Ncode = Ncur;
		memcpy(code, cur, sizeof(BVS_TERMS_ANY) * Ncode);// BVS_TERMS_CTRL must be smaller or equal to BVS_TERMS in sizee;
		return Ncode;
	}


	/****************************/
	// Match a (pre) and b (cur)
	/****************************/
	memset(mem, 0, sizeof(I16) * (Npre + Ncur));
	I16PTR statPre = mem;
	I16PTR statCur = mem + Npre;
	statCur[Ncur] = -1110L; // Add an extra one to indicate the end: To be used in the encoding process

	int first_NoMatch_in_Pre = 0;

	for (int i = 0; i < Ncur; ++i) {
		BVS_TERMS_PTR curterm = &cur[i];
		I32           bingo = _False_;
		for (int j = first_NoMatch_in_Pre; j < Npre; ++j) {
			if (statPre[j] == 0 && _IsTwoTermsEqual(&pre[j], curterm)) {
				statCur[i] = j + 1L;
				statPre[j] = 1L;
				bingo = _True_;
				break;
			}
		}

		// Update first_nomatch_in_a
		if (!bingo) {
			statCur[i] = -1L; // New terms and No match
			first_NoMatch_in_Pre = first_NoMatch_in_Pre; // the same as before
		}
		else {
			int new_first_NoMatch = Npre;
			for (int j = first_NoMatch_in_Pre; j < Npre; ++j) {
				if (statPre[j] == 0) {
					new_first_NoMatch = j;
					break;
				}
			}
			first_NoMatch_in_Pre = new_first_NoMatch;
		}
	}

	/****************************/
	// Start to Encode b;
	/****************************/
	int Ncode = 0;

	int start_idx = statCur[0];
	int end_idx = statCur[0];
	BVS_TRANGE rangeSE[5];
	int        nrange = 0;
	for (int i = 0; i < Ncur; ++i) {

		int idx = statCur[i];
		if (idx <0) {	

			// First to check if the prev range should be terminated
			if (end_idx > 0) {
			 // We MUST terminat the range group				
				rangeSE[nrange] = (BVS_TRANGE){start_idx, end_idx};
				nrange++;

				code[Ncode].type = CMPLIST_RANGE;
				code[Ncode].num  = nrange;
				memcpy(code[Ncode].range, rangeSE, sizeof(I16) * 2 * nrange);

				Ncode++;
			}
			end_idx =start_idx = statCur[i + 1];
			// if (i+1)=Nb;  start_idx=end_idx=-2L;  We will check if it is equal to -2 to see if the loop quit from this branch
			nrange    = 0;

			code[Ncode++] = *(BVS_TERMS_ANY_PTR)(cur + i);
			continue;
		}	
		else {
		//(idx > 0)
			assert(start_idx > 0 && end_idx >= start_idx);			

			if (idx - end_idx == 1 || idx - end_idx == 0) {
				// idx-end_idx ==  0  (1)  for the first ite if statB[0] > 0 
				// or (2) from the prev ite is from the "if ((idx <0)" part
				end_idx = idx;
			} else {		
				rangeSE[nrange] = (BVS_TRANGE){ start_idx, end_idx };
				nrange++;

				end_idx =start_idx = idx;	   

				if (nrange == 5) {			
					code[Ncode].type = CMPLIST_RANGE;
					code[Ncode].num = nrange;
					memcpy(code[Ncode].range, rangeSE, sizeof(I16) * 2 * nrange);

					Ncode++;
					nrange = 0;
				}
			}

		}////(idx > 0)


	}//for (int i = 0; i < Nb; ++i)

	if (end_idx != -1110L) {
		rangeSE[nrange] = (BVS_TRANGE){ start_idx, end_idx }; 
		nrange++;

		code[Ncode].type = CMPLIST_RANGE;
		code[Ncode].num = nrange;
		memcpy(code[Ncode].range, rangeSE, sizeof(I16) * 2 * nrange);

		Ncode++;
	}

	return Ncode;
}

int DecodeTermsList(BVS_TERMS_PTR pre,  BVS_TERMS_ANY_PTR code, BVS_TERMS_PTR out, int Npre, int Ncode) {
	// a+b=>c
	int Nout = 0;

	if (code[0].type == CMPLIST_RANGE && code[0].num==0) {
		assert(Ncode == 1);		
		return Nout = 0;
	}

	if (Npre == 0) {
	//  All the elments of c are not encoded; the case without any terms is handled above
		memcpy(out, code, sizeof(BVS_TERMS) * Ncode);		
		return Nout=Ncode;;
	}

	for (int i = 0; i < Ncode; i++) {
		if (code[i].type >= 0) {
			out[Nout++] = code[i].tm;
		} 
		else if (code[i].type == CMPLIST_RANGE) {
			for (int j = 0; j < code[i].num; j++) {				
				int s = code[i].range[j].s;
				int e = code[i].range[j].e;
				memcpy(out + Nout, pre + s - 1, sizeof(BVS_TERMS) * (e-s + 1)); 
				Nout += (e - s) + 1;				
			}			
		}
	}

	return Nout;
}


void AllocInit_TermEncoder(BVS_TERMS_ENCODER* coder, MemPointers * MEM, BVS_OPTIONS_PTR opt) {

	I64  nSamples = opt->mcmc.chainNumber * opt->mcmc.samples;
	I64  Kmax     = opt->prior.Kmax;
	coder->TERMS    = MyALLOC(*MEM, nSamples*(Kmax-1), BVS_TERMS, 8); // a better
	coder->BETA     = MyALLOC(*MEM, nSamples * Kmax, F32, 4); // a better
	coder->KcodeVec = MyALLOC(*MEM, nSamples , I16,2); // a better
	coder->KorgVec  = MyALLOC(*MEM, nSamples, I16, 2); // a better
	coder->PreList  = MyALLOC(*MEM,  Kmax, BVS_TERMS, 8); // a better
	coder->mem      = MyALLOC(*MEM, Kmax+Kmax+1, I16, 2); // a better

	coder->nTerms       = 0;
	coder->nTermsMax   =  nSamples * (Kmax - 1);
	coder->Kmax        = Kmax;
	coder->Kmax_actual = 0;
	coder->Kpre   =0;
	coder->nModels = 0;

	coder->nBeta = 0;

}
void InsertNewTermList(BVS_TERMS_PTR NewTermList, I32 Knewterm, BVS_TERMS_ENCODER* coder, F32PTR beta) {

	int Kcode= EncodeTermsList(coder->PreList, NewTermList, coder->TERMS + coder->nTerms, coder->Kpre, Knewterm, coder->mem);

	memcpy(coder->PreList, NewTermList, sizeof(BVS_TERMS) * Knewterm);
	coder->Kpre    = Knewterm;
	coder->nTerms += Kcode;

	coder->KorgVec[coder->nModels]   = Knewterm;
	coder->KcodeVec[coder->nModels]  = Kcode;
	coder->nModels++; 

	memcpy(coder->BETA + coder->nBeta, beta, sizeof(F32) * (Knewterm + 1));
	coder->nBeta += Knewterm + 1;

	coder->Kmax_actual = max(coder->Kmax_actual, Knewterm);
}


int ExtractNextTermFromList(BVS_TERMS_ENCODER* coder, BVS_TERMS_PTR out) {

	// before the first run, initilize coder following:
	// curBetaID=0;
	// curModelID=0;
	// curTermID=0
	// Kpre = 0;

	if (coder->curModelID >= coder->nModels) {
		return -1L;
	}

	int Kcode  = coder->KcodeVec[coder->curModelID];
	int Kout   = DecodeTermsList(coder->PreList, coder->TERMS + coder->curTermID, out, coder->Kpre, Kcode);

	assert(Kout == coder->KorgVec[coder->curModelID]);

	coder->Kpre      = Kout;
	memcpy(coder->PreList, out, sizeof(BVS_TERMS) * Kout);

	coder->curModelID+=1L;
	coder->curBetaID += Kout + 1;
	coder->curTermID += Kout;

	return Kout;
}

#define RoundTo4(t)    ((t+3)/4*4)

void SaveOutput(BVS_XINFO* xinfo, BVS_TERMS_ENCODER* coder, BVS_OPTIONS_PTR opt) {

#define x (*xinfo)
	I64  n =
		RoundTo4(sizeof(BVS_XINFO))/4 +
		x.N * x.Pobj /*xorg*/
		+ x.Pall         /*mean*/
		+ x.Pall         /*sd*/
		+ x.xi_Pbasis[QuadraticTerm]	/*mean2*/
		+ x.xi_Pbasis[QuadraticTerm]	/*sd2*/
		+ x.xi_Pbasis[LinearTerm]	/*Plinear*/
		+ x.xi_Pbasis[QuadraticTerm]	/*Psqrind0*/
		+ x.xi_Pbasis[StepTerm]	/*Xstepind0*/
		+ x.xi_Pbasis[HingeTerm]	/*Xhingeind0*/
		+ x.xi_Pbasis[ChangeTerm]	/*Xhingeind0*/
		+ x.xi_Pbasis[ChangeTerm]	/*PchangeY0*/
		+ x.xi_Pbasis[Hinge2DTerm]		/*Xhingeind0*/
		+ x.xi_Pbasis[Hinge2DRTerm]	/*Xhingeind0*/
		+ x.xi_Pbasis[Hinge1DRTerm]	/*Xhingeind0*/
		+ x.Nangle*2	/*anglecoff*/
		+ x.Pobj	    /*Nunique*/
		+ x.Pobj * x.N /* SortedUnqiueIdx0*/;

	int     len = n;
	I32PTR  out;
	VOIDPTR odata = PROTECT(  CreateNumVar(DATA_INT32,  &len,1, &out) );

	I64 csum = 0;
	len = RoundTo4(sizeof(BVS_XINFO)) / 4; 	f32_copy(xinfo, out + csum, len); csum += len;
	len = x.N * x.Pobj; 					f32_copy(x.Xorg, out + csum, len); csum += len;
	len = x.Pall;          					f32_copy(x.mean, out + csum, len); csum += len;
	len = x.Pall;          					f32_copy(x.sd , out + csum, len); csum += len;
	len = x.xi_Pbasis[QuadraticTerm];       f32_copy(x.mean2, out + csum, len); csum += len;
	len = x.xi_Pbasis[QuadraticTerm];       f32_copy(x.sd2, out + csum, len); csum += len;
	len = x.xi_Pbasis[LinearTerm];       	f32_copy(x.xi_Pindices0[LinearTerm], out + csum, len); csum += len;
	len = x.xi_Pbasis[QuadraticTerm];    	f32_copy(x.xi_Pindices0[QuadraticTerm], out + csum, len); csum += len;
	len = x.xi_Pbasis[StepTerm];			f32_copy(x.xi_Pindices0[StepTerm], out + csum, len); csum += len;
	len = x.xi_Pbasis[HingeTerm];			f32_copy(x.xi_Pindices0[HingeTerm], out + csum, len); csum += len;
	len = x.xi_Pbasis[ChangeTerm];			f32_copy(x.xi_Pindices0[ChangeTerm], out + csum, len); csum += len;
	len = x.xi_Pbasis[ChangeTerm];			f32_copy(x.xi_PchangeY0, out + csum, len); csum += len;
	len = x.xi_Pbasis[Hinge2DTerm];			f32_copy(x.xi_Pindices0[Hinge2DTerm], out + csum, len); csum += len;
	len = x.xi_Pbasis[Hinge2DRTerm];		f32_copy(x.xi_Pindices0[Hinge2DRTerm], out + csum, len); csum += len;
	len = x.xi_Pbasis[Hinge1DRTerm];		f32_copy(x.xi_Pindices0[Hinge1DRTerm], out + csum, len); csum += len;
	len = x.Nangle * 2;						f32_copy(x.anglecoff, out + csum, len); csum += len;
	len = x.Pobj;          					f32_copy(x.Nunique, out + csum, len); csum += len;
	len = x.Pobj * x.N;					 	f32_copy(x.SortedUnqiueIdx0, out + csum, len); csum += len;


	ReplaceStructField(opt->ans, "xinfo", odata);

	/*****************************************************************************************/
	/*
	I16PTR            KorgVec;  // size: nSamples*nChains
	I16PTR            KcodeVec; // size: nSamples*nChains

	BVS_TERMS_ANY_PTR TERMS;     // size: nSamples*nChains*(Kmax-1)
	F32PTR            BETA;      // size: nSamples*nChains*(Kmax)
	*/

	n = RoundTo4(sizeof(BVS_TERMS_ENCODER))  +
		coder->nModels * sizeof(I16) + coder->nModels * sizeof(I16)
		+ coder->nTerms * sizeof(BVS_TERMS_ANY)
		+ coder->nBeta * sizeof(F32);

	I08PTR  out8;
	len = n;
	odata = PROTECT( CreateNumVar(DATA_INT32, &len, 1, &out8) );

	csum = 0;
	len = RoundTo4(sizeof(BVS_TERMS_ENCODER)); 		memcpy(out8 + csum, coder, len); csum += len;
	len = coder->nModels * sizeof(I16); 			memcpy(out8 + csum, coder->KorgVec, len); csum += len;
	len = coder->nModels * sizeof(I16); 			memcpy(out8 + csum, coder->KcodeVec, len); csum += len;
	len = coder->nTerms * sizeof(BVS_TERMS_ANY); 	memcpy(out8 + csum, coder->TERMS, len); csum += len;
	len = coder->nBeta * sizeof(F32); 			    memcpy(out8 + csum, coder->BETA, len); csum += len;

	ReplaceStructField(opt->ans, "terms", odata);

	/*****************************************************************************************/
	/*****************************************************************************************/

	len = sizeof(BVS_OPTIONS);
	len = (len + 3) / 4 * 4;
	len = len / 4;
	odata = PROTECT(CreateNumVar(DATA_INT32, &len, 1, &out8));

	csum = 0;
	len = sizeof(BVS_OPTIONS); 	memcpy(out8 + csum, &opt->prior, len); csum += len;
	
	ReplaceStructField(opt->ans, "opt", odata);


	UNPROTECT(3);

#undef x
}

#define hasTerm(x)  (prior->type2id[x] >=0)
int  GetArg_Predict(VOIDPTR prhs[], int nrhs, BVS_XINFO_PTR xinfo, BVS_TERMS_ENCODER_PTR coder,
	                 BVS_PRED_XINFO *newxinfo, BVS_PRIOR_PTR* prior_ptr, BVS_PRED_OPT *newopt)

{
	// Check the number of inputs:  the first arg is the algorithm name.
	if (nrhs < 2L) {
		r_error("ERROR: At least one input argument is needed!\n");
		return 0;
	}

	VOIDPTR S;

	////////////////////////////////////////////////////////
	newopt->bComputeYchains = 0;
	if (nrhs >= 4) {		
		S = prhs[3];
		if (IsNumeric(S)) {
			newopt->bComputeYchains = GetScalar(S);
		} else if (IsStruct(S)) {

		} else {
			newopt->bComputeYchains = 0;
		}
	}	
#if   R_INTERFACE==1
	newopt->odtype = DATA_DOUBLE;
#elif M_INTERFACE==1
	newopt->odtype = DATA_FLOAT;
#endif
 
	////////////////////////////////////////////////////////////////
	
	S = prhs[1];
	if (!IsStruct(S)) {
		r_error("ERROR: the input must be a struct variable.\n");
		return 0;
	}
	VOIDPTR tmp = GetField(S,"xinfo");
	I08PTR  out = GetData(tmp);

    

	#define x (*xinfo)

	I64 csum = 0;
	I64 len;

	len = RoundTo4(sizeof(BVS_XINFO));  *xinfo=*(BVS_XINFO_PTR)out; csum += len;
	len = sizeof(I32) *x.N * x.Pobj; 				   x.Xorg  = len==0? NULL: out + csum; csum += len;
	len = sizeof(I32) * x.Pall;          				x.mean  = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.Pall;          				x.sd    = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.xi_Pbasis[QuadraticTerm];     x.mean2 = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.xi_Pbasis[QuadraticTerm];     x.sd2   = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.xi_Pbasis[LinearTerm];        x.xi_Pindices0[LinearTerm]   = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.xi_Pbasis[QuadraticTerm];     x.xi_Pindices0[QuadraticTerm] = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.xi_Pbasis[StepTerm];			x.xi_Pindices0[StepTerm]     = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.xi_Pbasis[HingeTerm];			x.xi_Pindices0[HingeTerm]    = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.xi_Pbasis[ChangeTerm];		x.xi_Pindices0[ChangeTerm]   = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.xi_Pbasis[ChangeTerm];		x.xi_PchangeY0               = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.xi_Pbasis[Hinge2DTerm];		x.xi_Pindices0[Hinge2DTerm]  = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.xi_Pbasis[Hinge2DRTerm];		x.xi_Pindices0[Hinge2DRTerm] = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.xi_Pbasis[Hinge1DRTerm];		x.xi_Pindices0[Hinge1DRTerm] = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.Nangle*2;     x.anglecoff = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.Pobj;         x.Nunique = len == 0 ? NULL : out + csum; csum += len;
	len = sizeof(I32) * x.Pobj * x.N;   x.SortedUnqiueIdx0 = len == 0 ? NULL : out + csum; csum += len;

	/////////////////////////////////////////////////////////////////
	BVS_PRIOR_PTR  prior;
	tmp = GetField(S, "opt");	
	prior_ptr[0] = prior = GetData(tmp);

	/////////////////////////////////////////////////////////////////

	tmp = GetField(S, "terms");
	out = GetData(tmp);
	
	csum = 0;
	len = RoundTo4(sizeof(BVS_TERMS_ENCODER)); 		*coder=*(BVS_TERMS_ENCODER_PTR)out; csum += len;
	len = coder->nModels * sizeof(I16); 			coder->KorgVec = out + csum; csum += len;
	len = coder->nModels * sizeof(I16); 			coder->KcodeVec = out + csum; csum += len; 
	len = coder->nTerms  * sizeof(BVS_TERMS_ANY); 	coder->TERMS  = out + csum; csum += len;
	len = coder->nBeta * sizeof(F32); 			    coder->BETA    = out + csum; csum += len;  

	
    #undef x

	/*
	 !(IsDouble(DATA) && numel > 2) && 	!(IsSingle(DATA) && numel > 2) &&
			!(IsInt32(DATA)  && numel > 2) &&	!(IsInt16(DATA) && numel > 2) &&
			!(IsStruct(DATA) && numel >= 1 ) &&
			!(IsCell(DATA) )
	*/
 

	/*************************************************/
	/* Get the dimension of the input              */
	/*************************************************/
	VOID_PTR  XobjInput = prhs[2];
	I32       ndimx  = GetNumOfDim(XobjInput);

	if (ndimx == 0) {
		// the input is a vector: this branch is possible only for R vectors (not  Matlab)
		ndimx = 1L;
		
		newxinfo->Pobj = GetNumberOfElements(XobjInput);
		newxinfo->N     = 1L;

		if (xinfo->Pobj == 1) {
			newxinfo->N    = newxinfo->Pobj;
			newxinfo->Pobj = 1L;
		}
	} //ndims is impossible to be 1L
	else if (ndimx == 2) {
		// Matlab: a vector or matrix;
		// R:      a matrix or a vector with a dim attribute.
		newxinfo->N    = GetDim1(XobjInput);
		newxinfo->Pobj = GetDim2(XobjInput);

	}

	int N    = newxinfo->N;
	int Pobj = newxinfo->Pobj;
	int Pall   = xinfo->Pall;

	newxinfo->Xorg   = malloc(sizeof(F32) * N * Pall);	
	CopyNumericObjToF32Arr(newxinfo->Xorg, XobjInput, N * Pobj);
 

	/******************************************************/
	// For FFT, generate sin cons terms
	/******************************************************/
	if (prior->isFFT) {
	
	 	F32    freq_factor    = 2. * 3.141592653589793f / prior->fft_period;
		I32    maxSeasonOrder = prior->beast_maxSeasonOrder;
		F32PTR ptr            = newxinfo->Xorg + N*(prior->fft_Pbase-1);
		F32PTR xcol           = newxinfo->Xorg + N*(prior->fft_Ptime- 1);

		for (int order = 1; order <= maxSeasonOrder; order++) {
			f32_copy(xcol, ptr, N);
			f32_mul_val_inplace(freq_factor * (F32)order, ptr, N);
			f32_copy(ptr, ptr + N, N);

			f32_sincos_vec_inplace(ptr + N, ptr, N); //(sin,cos): we want the cos term first and then the sin term.
		   //because for the highest order period/2, sins are all zeros.
			ptr += 2 * N;
		}
 	}

	/******************************************************/
	// For time series decomposition, generate sin cons terms
	/******************************************************/
	if (prior->isBeast &&  prior->beast_Pextra > 0) {

		F32    freq_factor = 2.0f * 3.141592653589793f / prior->beast_period;
		F32PTR ptr         = newxinfo->Xorg + N * (prior->beast_Pbase - 1L);
		F32PTR xcol        = newxinfo->Xorg + N * (prior->beast_Ptime - 1L);
		for (int order = 1; order <= prior->beast_maxSeasonOrder; order++) {
			f32_copy(xcol, ptr, N);
			f32_mul_val_inplace(freq_factor * (F32)order, ptr, N);
			f32_copy(ptr, ptr + N, N);

			f32_sincos_vec_inplace(ptr + N, ptr, N); //(sin,cos): we want the cos term first and then the sin term.
		   //because for the highest order period/2, sins are all zeros.
			ptr += 2 * N;
		}
	}

	//CopyStrideMEMToF32Arr(Y + i * N/*dst*/, io->pdata[i] /*src*/, N, stride, offset, io->dtype[i]);
	//f32_set_nan_by_value(Y + i * N, N, io->meta.missingValue);

	//	int  nMissing = f32_find_nans(Y, N, badRowIdx);
	//f32_mat_multirows_extract_set_by_scalar(Y, N, 1L, Ycopy, badRowIdx, nMissing, 0);

	//f32_mat_multirows_extract_set_by_scalar(Xtmp, N, K + 1, Xcopy, badRowIdx, nMissing, 0);
	//	linear_regression(Y, X, N, N, K, B, Yfit, NULL, XtX);
	//	f32_mat_multirows_set_by_submat(Xtmp, N, K + 1, Xcopy, badRowIdx, nMissing);
	 
	int     Plin      = xinfo->xi_Pbasis[LinearTerm];
	I32PTR  Pindices0 = xinfo->xi_Pindices0[LinearTerm];
	newxinfo->Xnorm1  = malloc(sizeof(F32) * N * Plin); 		
	for (int p = 0; p < Plin; p++) {
		f32_copy( newxinfo->Xorg + N* Pindices0[p], newxinfo->Xnorm1 + p * N, N);
		f32_scale_inplace(1./xinfo->sd[p], -xinfo->mean[p]/xinfo->sd[p], newxinfo->Xnorm1+p*N, N);
	}
	

	if (hasTerm(StepTerm)) {
		int   P = xinfo->xi_Pbasis[StepTerm];
		assert(P > 0);
		newxinfo->Xones = malloc(sizeof(F32) * N);
		f32_fill_val(1.0f, newxinfo->Xones, N);
	}

	if ( hasTerm(QuadraticTerm) ) {
		int     P         = xinfo->xi_Pbasis[QuadraticTerm];
		I32PTR  Pindices0 = xinfo->xi_Pindices0[QuadraticTerm];
		assert(P > 0);
		newxinfo->X2      = malloc(sizeof(F32) * N * P);
		for (int i = 0; i < P; i++) {
			f32_copy(newxinfo->Xorg + Pindices0[i] * N, newxinfo->X2 + i * N,  N);
			f32_mul_vec_inplace(newxinfo->X2 + i * N, newxinfo->X2 + i * N, N);
			f32_scale_inplace(1. / xinfo->sd2[i], -xinfo->mean2[i] / xinfo->sd2[i], newxinfo->X2 + i * N, N);
		}
	}
 
	return 1;
}

VOIDPTR Predict_Alloc_Output(BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, BVS_PRED_RESULT_PTR mat, BVS_PRIOR_PTR prior, int nModels, BVS_PRED_OPT *newopt) {

	int N     = newxinfo->N;
	int P     = xinfo->Pall;
	int dtype = newopt->odtype;
 
	FIELD_ITEM  fieldList[] = {
			{"Y",            dtype,		2,   {N,1,},       &mat->Y}, // Needed to be changed to reflect the output dim
			{"SD",           dtype,		2,   {N,1,},       &mat->SD}, // Needed to be changed to reflect the output dim
			{"Ychains",      dtype,		2,   {N,nModels,}, &mat->Ychains}, // Needed to be changed to reflect the output dim
			{"Khinge",       dtype,		2,   {(1+prior->maxKnotNum_Max[HingeTerm]), xinfo->xi_Pbasis[HingeTerm],}, &mat->Khinge}, // Needed to be changed to reflect the output dim
			{"Yparts",       dtype,		2,   {N,P,},       &mat->Yparts}, // Needed to be changed to reflect the output dim
			{"P",            dtype,		2,   {P,1,},       &mat->P}, // Needed to be changed to reflect the output dim
			{"BetaLinear",   dtype,		2,   {P,1,},       &mat->BetaLinear}, // Needed to be changed to reflect the output dim			
			{0,},
	};
	I32    nfields = 99999;
	if (!hasTerm(LinearTerm)) {
		RemoveField(fieldList, nfields, "BetaLinear");
	}
	if (!newopt->bComputeYchains) RemoveField(fieldList, nfields, "Ychains");
	if (!hasTerm(HingeTerm)) RemoveField(fieldList, nfields, "Khinge");

 
	return CreateStructVar(fieldList, nfields); // no need to protect it bvz it is returned immediately
}
#include "abc_000_warning.h"