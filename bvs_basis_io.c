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
 
 
int ReadPreprocessInputData(BVS_OPTIONS_PTR opt, BVS_XINFO_PTR xinfo, BVS_YINFO_PTR yinfo, F32PTR MEMBUF) {

	int Nobj = xinfo->Nobj = opt->io.meta.Nobj;
	int Pobj = xinfo->Pobj = opt->io.meta.Pobj;
	CopyNumericObjToF32Arr(xinfo->Xorg, opt->io.meta.Xobj, Nobj * Pobj);
	CopyNumericObjToF32Arr(yinfo->Y, opt->io.meta.Yobj, Nobj);

	/******************************************************/
	// Find missing indices
	/******************************************************/
	int    nMissing = 0;
	{  // Detemine bad Rows
		I32PTR isBadRow = MEMBUF;
		memset(isBadRow, 0, sizeof(I32) * Nobj);
		F32PTR  Xorg = xinfo->Xorg;
		for (int p = 0; p < Pobj; p++) {
			for (int i = 0; i < Nobj; i++) {
				isBadRow[i] |= IsNaN(Xorg[i]);
			}
			Xorg += Nobj;
		}
		F32PTR  Yraw = yinfo->Y;
		for (int i = 0; i < Nobj; i++) {
			isBadRow[i] |= IsNaN(Yraw[i]);
		}

		for (int i = 0; i < Nobj; i++) {
			if (isBadRow[i]) {
				xinfo->RowsMissing[nMissing++] = i;
			}
		}

	}
	xinfo->Nmissing = nMissing;
	xinfo->N = Nobj - nMissing;
	yinfo->N = Nobj - nMissing;


	/******************************************************/
	// Remove NANs
	/******************************************************/
	int N = xinfo->N;
	if (nMissing > 0) {

		I32PTR  RowsMissing = xinfo->RowsMissing;
		F32PTR  Xraw = xinfo->Xorg;
		F32PTR  Xorg = xinfo->Xorg;
		for (int p = 0; p < Pobj; p++) {

			int j = 0;
			int nOmit = 0;
			for (int i = 0; i < Nobj; i++) {
				if (nOmit < nMissing && i == RowsMissing[nOmit]) {
					nOmit++;
				}
				else {
					Xorg[j++] = Xraw[i];
				}
			}
			Xraw += Nobj;
			Xorg += N;

		}

		F32PTR  Yraw = yinfo->Y;
		F32PTR  Y = yinfo->Y;
		int j = 0;
		int nOmit = 0;
		for (int i = 0; i < Nobj; i++) {
			if (nOmit < nMissing && i == RowsMissing[nOmit]) {
				nOmit++;
			}
			else {
				Y[j++] = Yraw[i];
			}
		}

	}// Remove NANs


	//******************************************************/
	// Sort X
	/******************************************************/
	F32PTR Xorg = xinfo->Xorg;
	F32PTR Xcopy = MEMBUF;
	I32PTR Indices = Xcopy + N;
	for (int p = 0; p < Pobj; p++) {

		i32_seq(Indices, 0, 1, N);
		f32_copy(Xorg + N * p, Xcopy, N);
		QuickSortA(Xcopy, Indices, 0, N - 1);

		int Nunique = 1;
		F32 Xcur = Xcopy[0];
		for (int i = 1; i < N; ++i) {
			if (Xcopy[i] > Xcur) {
				Xcur = Xcopy[i];
				Indices[Nunique] = Indices[i];
				Nunique++;
			}
		}

		xinfo->Nunique[p] = Nunique;
		memcpy(xinfo->SortedUnqiueIdx0 + N * p, Indices, sizeof(I32) * Nunique);
	}


	/******************************************************/
	// For FFT, generate sin cons terms
	/******************************************************/
	if (opt->prior.isFFT) {

		int    ptime0 = opt->prior.fft_Ptime - 1;
		F32PTR Xtime  = Xorg + N * ptime0;
		I32    Nt     = xinfo->Nunique[ptime0];

		opt->prior.fft_period = (Xtime[xinfo->SortedUnqiueIdx0[N * ptime0 + Nt - 1]] - Xtime[xinfo->SortedUnqiueIdx0[N * ptime0 + 0]]) / (Nt - 1.0) * Nt;
 
		F32    freq_factor = 2. * 3.141592653589793f / opt->prior.fft_period;		
		// If maxSeaonOder is not set or invalid, choose the largest possible value N/2
		I32    maxFFTOrder = opt->prior.fft_maxOrder;
		F32PTR ptr         = Xorg + N * (opt->prior.fft_Pbase - 1);
		for (int order = 1; order <= maxFFTOrder; order++) {
			f32_copy(Xtime, ptr, N);
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
	if (opt->prior.isBeast && opt->prior.beast_Pextra > 0) {

		assert(!IsNaN(opt->prior.beast_period));

		int    ptime0 = opt->prior.beast_Ptime - 1;
		F32PTR Xtime  = Xorg + N * ptime0;
		F32    freq_factor = 2.0f * 3.141592653589793f / opt->prior.beast_period;
		F32PTR ptr    = Xorg + N * (opt->prior.beast_Pbase - 1);

		for (int order = 1; order <= opt->prior.beast_maxSeasonOrder; order++) {
			f32_copy(Xtime, ptr, N);
			f32_mul_val_inplace(freq_factor * (F32)order, ptr, N);
			f32_copy(ptr, ptr + N, N);

			f32_sincos_vec_inplace(ptr + N, ptr, N); //(sin,cos): we want the cos term first and then the sin term.
		   //because for the highest order period/2, sins are all zeros.
			ptr += 2 * N;
		}
	}


	/******************************************************/
	// Normalization
	/******************************************************/

	f32_normalize_std_avg_inplace(yinfo->Y, N, yinfo->mean, yinfo->sd);

 
	int    P    = opt->prior.Pall;
	int    Plin = opt->prior.Pbasis[LinearTerm];

	for (int i = 0; i < Plin; i++) {
		int    pidx0   = xinfo->xi_Pindices0[LinearTerm][i];
		F32PTR xsrc    = xinfo->Xorg +N*pidx0;
		F32PTR xdst   = xinfo->X1norm+N* i;
		memcpy(xdst, xsrc, N * sizeof(F32));

		if ( (pidx0+1) <= xinfo->Pobj) {
			f32_normalize_std_avg_inplace(xdst, N, xinfo->mean + i, xinfo->sd + i);
		}	else {
			// For BEAST and FFT the cos and sin terms are not centerned and only 
			// scaled to have a sum-proudct of N
		 
			F32 dotProduct = DOT(N, xdst, xdst);
			F32 scale      = 1 / sqrtf(dotProduct / N);
			f32_mul_val_inplace(scale, xdst, N);
			xinfo->mean[i] = 0;
			xinfo->sd[i]   = 1 / scale;
		}
	}
	

	//******************************************************/
	// Pre-compute YtY
	/******************************************************/
	yinfo->YtY_plus_alpha2Q[0] = DOT(N, yinfo->Y, yinfo->Y) + opt->prior.alpha2;


	/******************************************************/
	// Handle for the quadratic termsAlloc_Xterms
	/******************************************************/	
	if (hasTerm(QuadraticTerm)) { //xinfo->X2 && xinfo->mean2 && xinfo->sd2
		F32PTR X2   = xinfo->X2;
		I32    Psqr = xinfo->xi_Pbasis[QuadraticTerm];
		for (int p = 0; p < Psqr; p++) {
			F32PTR xcol = X2 + p * N;
			memcpy(xcol, xinfo->Xorg+xinfo->xi_Pindices0[QuadraticTerm][p]*N, N * sizeof(F32));
			f32_mul_vec_inplace(xcol, xcol, N);
			f32_normalize_std_avg_inplace(xcol, N, xinfo->mean2 + p, xinfo->sd2 + p);
		}
	}

	//nMissing = f32_normalize_multicols_zeroout_nans(xinfo->X, rowsMissing, N, N, P,  xinfo->mean, xinfo->sd);
    //nMissing = f32_normalize_multicols_zeroout_nans(yinfo->Y, rowsMissing, N, N, 1L, yinfo->mean, yinfo->sd);

	return 1;
}
#include "abc_000_warning.h"