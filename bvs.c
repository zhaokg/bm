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

//#include <crtdbg.h>

#define LOCAL(...) do{ __VA_ARGS__ } while(0);
//extern MemPointers* mem;
//time_t start, end; 

void  GenarateRandomBasis(BVS_MODEL_PTR  model, BVS_RNDSTREAM * PRAND, BVS_XINFO_PTR xinfo, BVS_OPTIONS_PTR opt) {

	// Initit_BasisState_all must be called before this func

	I32 Kfixed = 0;
	// return a cosnt term only
	model->terms[0].type  = ConstTerm;	
	model->terms[0].var0  = model->terms[0].knot = -1L; // For CompareTerms
	model->Kterms = 1;
	Kfixed++;

	NewPropTerm    NEW;     // moved here bvz its xcols has two fixed lements, N and Npad
	NEWCOLINFOv2   NewCol; // moved here bvz its xcols has two fixed lements, N and Npad
	NEW.xcols     = &NewCol;
	NewCol.N      = xinfo->N;
	NewCol.Nlda   = xinfo->N;
	if (opt->prior.FixedLinearTerms0 && model->type2id[LinearTerm] >=0) {
	    
		BSTATE_LINEAR *blin = &model->basisState[model->type2id[LinearTerm]];	    
	       
		for (int i = 0; i < opt->prior.nFixedLinearTerms; i++) {
			int p0 = opt->prior.FixedLinearTerms0[i];
			
			if ( (p0+1) >= 1 && (p0+1) <= xinfo->Pall) {
				
				NEW.terms[0].type  = LinearTerm;
				NEW.terms[0].var0   = blin->pidx0_to_vidx0[p0];
				NEW.terms[0].knot  = -1L; // For CompareTerms
				NEW.jumpType   = BIRTH;
				NEW.nTerms     = 1;

				NewCol.K       = model->Kterms;
				NewCol.nbands  = 1L;
				NewCol.ks_x[0]     = model->Kterms + 1L;
				NewCol.kterms_x[0] = 0;
				NewCol.ks_xnewterm[0]     = 1L;
				NewCol.kterms_xnewterm[0] = 1L;		 
				get_parts_for_newinfo(&NewCol);

				Update_BasisState_All(&NEW, model, xinfo); //xInfo.Nunique is nseed inside

				Kfixed++;
			}	
		}
		
	}

	model->Kfixed = Kfixed;
 
}

void IncreasePrecValues(BVS_MODEL_PTR model) {
	// 1:UniformPrec
	model->precVec[0]    = 2*model->precVec[0];
	model->logPrecVec[0] = logf(model->precVec[0]);
}  

void ComputeMargLik(BVS_MODELDATA_PTR data, BVS_MODEL_PTR model,	BVS_YINFO_PTR yInfo,  BVS_HyperPar_PTR hyperPar)
{
	 I32 K = data->K;
	 solve_U_as_LU_invdiag_sqrmat(data->cholXtX, data->XtY, data->beta_mean, K);
	//Compute alpha2_star: Two ways to compute alpha2_star: YtY-m*V*m
	/* ---THE FIRST WAY---
			//GlobalMEMBuf_1st = Xnewterm + K_newTerm*Npad;
			r_cblas_scopy(KNEW, beta_mean_prop, 1, GlobalMEMBuf_2nd, 1);
			//cblas_dtrmv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag, const MKL_INT n, const double *a, const MKL_INT lda, double *x, const MKL_INT incx);
			r_cblas_strmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, KNEW, cholXtX_prop, KNEW, GlobalMEMBuf_2nd, 1);
			//double cblas_ddot (const MKL_INT n, const double *x, const MKL_INT incx, const double *y, const MKL_INT incy);
			basis_prop->alpha2_star = yInfo.YtY - DOT(KNEW, GlobalMEMBuf_2nd, GlobalMEMBuf_2nd);
	*/
	/*---THE SECOND WAY-- */
	F32 alpha2_star			=  (yInfo->YtY_plus_alpha2Q[0] - DOT(K, data->XtY,data->beta_mean))*0.5;

	//half_log_det_post  = sum(log(1. / diag(U)))	 
	//half_log_det_post  = -sum_log_diagv2(MODEL.prop.cholXtX, KNEW);
	F32 half_log_det_post = sum_log_diagv2(data->cholXtX, K);

	////Copy the diagonal of P_U to buf
	// vmsLn(KNEW, GlobalMEMBuf_1st, GlobalMEMBuf_1st, VML_HA);
	// ippsSum_32f(GlobalMEMBuf_1st, KNEW, &FLOAT_SHARE.half_log_det_post, ippAlgHintFast);
	// r_ippsSumLn_32f(GlobalMEMBuf_2nd, KNEW, &half_log_det_post);	

	//half_log_det_prior = -0.5 * (k_SN * log(prec(1)) + length(k_const_terms) * log(prec(2)) + length(linear_terms) * log(prec(3)));
	// F32 half_log_det_prior = -0.5f * KNEW * *modelPar.LOG_PREC[0]

	F32 half_log_det_prior = -0.5f*model->logPrecVec[0] * (K-1);  // K-1 because the first term is the constant term and its prec is set to 0.

	//log_ML = half_log_det_post - half_log_det_prior - (n / 2 + b) * log(a + a_star / 2);
	F32 marg_lik   = half_log_det_post  -half_log_det_prior	- yInfo->alpha1_star * fastlog(alpha2_star);

	data->alpha2_star = alpha2_star;
	data->marg_lik    = marg_lik;
}
 
void BVS_EvaluateModel( BVS_MODELDATA *curmodel, F32PTR Xt_mars,	
	BVS_TERMS_PTR terms, I32 Kterm,
	BVS_XINFO_PTR xInfo,BVS_YINFO_PTR  yInfo, 
	BVS_HyperPar *hyperPar, F32 precVal, VOID_PTR stream )
{
	//TO RUN THIS FUNCTION, MUST FIRST RUN CONVERTBAIS so to GET NUMBERS OF BASIS TERMS	
	I32 N    = xInfo->N;
	I32 Npad = N;
	I32 K    = Kterm;	 
	curmodel->K = K;
	GenTerms(terms, xInfo, Xt_mars, K);

	//Calcuate X'*X and Calcuate X'*Y
	//DGEMM('T', 'N', &m, &n, &K, &alpha, X_mars, &K, X_mars, &K, &beta, XtX, &m);	
	//cblas_dgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const double alpha, const double *a, const MKL_INT lda, const double *b, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc);

	F32PTR XtX = curmodel->XtX;
	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K, K, N, 1.f, Xt_mars, Npad, Xt_mars, Npad, 0.f, XtX, K);
	//cblas_dgemmt (const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT n, const MKL_INT k, const double alpha, const double *a, const MKL_INT lda, const double *b, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc);
	//cblas_dgemmt(CblasColMajor, CblasLower, CblasTrans, CblasNoTrans, k, N, 1, X_mars, N, X_mars, N, 0, XtY, k);				

	F32PTR XtY = curmodel->XtY;
	//cblas_dgemv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n, const double alpha, const double *a, const MKL_INT lda, const double *x, const MKL_INT incx, const double beta, double *y, const MKL_INT incy);
	r_cblas_sgemv(CblasColMajor, CblasTrans, Npad, K, 1, Xt_mars, Npad, yInfo->Y, 1, 0, XtY, 1);

	//X_mars has been constructed. Now use it to calcuate its margin likelihood
	//Add precison values to the diagonal of XtX: post_P=XtX +diag(prec)	
	F32PTR cholXtX = curmodel->cholXtX;
	//f32_copy( XtX,  cholXtX, K * K);
	//f32_add_val_matrixdiag(cholXtX, precVal, K);
	chol_addCol_skipleadingzeros_prec_nostartprec_invdiag(XtX, cholXtX, &precVal, K, 1, K);

	//Solve inv(Post_P)*XtY using  Post_P*b=XtY to get beta_mean
	F32PTR beta_mean = curmodel->beta_mean;	
	solve_U_as_LU_invdiag_sqrmat(cholXtX, XtY, beta_mean, K);
	/*	
		//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, double * a, lapack_int lda);
		r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', K, cholXtX, K); // Choleskey decomposition; only the upper triagnle elements are used
		//LAPACKE_spotrs (int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
		r_cblas_scopy(K, XtY, 1, beta_mean, 1);
		r_LAPACKE_spotrs(LAPACK_COL_MAJOR, 'U', K, 1, cholXtX, K, beta_mean, K);	
	*/
	
	//Compute beta = beta_mean + Rsig2 * randn(p, 1);
	//Usig2 = (1 / sqrt(sig2)) * U; 		beta = beta_mean + linsolve(Usig2, randn(p, 1), opts);
	//status = vdRngGaussian( method, stream, n, r, a, sigma );

	/**********************************************************/
	// Sample beta from beta_mean and cholXtX
	/********************************************************/
	/*
	F32PTR beta = model->beta;
	r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, (*(VSLStreamStatePtr *)stream), K, beta, 0, 1);

	// LAPACKE_strtrs (int matrix_layout , char uplo , char trans , char diag , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
	{
		//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, 1, post_P_U, K, beta, K);
		solve_U_as_U_invdiag(cholXtX, beta, K, K);
	}
	r_ippsMulC_32f_I(sqrtf(model->sig2), beta, K);
	r_ippsAdd_32f_I(beta_mean, beta, K);
	*/
	/**********************************************************/
	// Sample beta from beta_mean and cholXtX
	/********************************************************/

	//Compute alpha2_star	 
	/*
	r_cblas_scopy(K, beta_mean, 1, GlobalMEMBuf_1st, 1);
	//cblas_dtrmv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag, const MKL_INT n, const double *a, const MKL_INT lda, double *x, const MKL_INT incx);
	r_cblas_strmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, K, cholXtX, K, GlobalMEMBuf_1st, 1);
	F32 alpha2_star = pyInfo->YtY - DOT(K, GlobalMEMBuf_1st, GlobalMEMBuf_1st);
	 */
	F32 alpha2_star      = (yInfo->YtY_plus_alpha2Q[0] - DOT(K,  XtY,  beta_mean)) * 0.5;

	//half_log_det_post; = sum(log(1. / diag(U)))
	F32 half_log_det_post = sum_log_diagv2(cholXtX, K);
	
	//half_log_det_prior = -0.5 * (k_SN * log(prec(1)) + length(k_const_terms) * log(prec(2)) + length(linear_terms) * log(prec(3)));		
	F32 half_log_det_prior	= -0.5f * (K-1)*logf(precVal);  // K-1 because the first term is constant and its prec is set to 0

	//log_ML    = half_log_det_post - half_log_det_prior - (n / 2 + b) * log(a + a_star / 2);
	F32 marg_lik = half_log_det_post - half_log_det_prior - yInfo->alpha1_star * fastlog(alpha2_star);

	curmodel->alpha2_star = alpha2_star;
	curmodel->marg_lik    = marg_lik;

	return;
}
 
 
// F32PTR Ys = Xnewterm + N * 2; // the first two cols of Xnewtemrs were used to store new terms
void ComputeMarginalY_ProposedModel( BVS_MODEL_PTR model, NewPropTerm  *new, 
                    NEWCOLINFOv2  *newCol, 	I32PTR * Pindices0,  I32 pidx0, F32PTR Y, int N) 
{
	 
		f32_fill_val(0, Y, N);
		int Kidx_newterm = 1;
		for (int i = 0; i < (newCol->nbands * 2 + 1L); i++) {
		
			F32PTR X    = newCol->parts[i].X;
			int    Ksrc = newCol->parts[i].ks_src;
			int    Kdst = newCol->parts[i].ks_dst;
			BVS_TERMS* terms = (X == newCol->X)?  model->terms: new->terms;
		
			for (int j = 0; j < newCol->parts[i].kterms; j++) {			
				BVS_TERMS* term = &terms[Ksrc + j - 1];
				int type = term->type;
				if (type == LinearTerm || type == HingeTerm || type == StepTerm) {			 
					if (pidx0 == Pindices0[type][term->var0] ) {						
						F32  beta = model->prop.beta_mean[Kdst + j - 1];
						f32_axpy_inplace(beta, X + (Ksrc + j - 1) * N, Y, N);
					}
				}
			}		
		}
 
}

int CheckCurveProperties_Increasing(F32PTR Y, int N, int p0, BVS_XINFO_PTR xinfo) {

	int     Nlen     = xinfo->Nunique[p0];
	I32PTR  indices  = xinfo->SortedUnqiueIdx0 + N * p0;


	int   bingo = 1;
	F32   Yprev = Y[indices[0]];
	for (int i = 1; i < Nlen; ++i) {
		F32 Ycur = Y[indices[i]];
		if (Ycur < Yprev) {
			bingo = 0;
			break;
		}
		Yprev = Ycur;
	}

	return bingo;
}
	
int CheckCurveProperties_Decreasing(F32PTR Y, int N, int p0, BVS_XINFO_PTR xinfo) {

	int     Nlen = xinfo->Nunique[p0];
	I32PTR  indices = xinfo->SortedUnqiueIdx0 + N * p0;


	int   bingo = 1;
	F32   Yprev = Y[indices[0]];
	for (int i = 1; i < Nlen; ++i) {
		F32 Ycur = Y[indices[i]];
		if (Ycur > Yprev) {
			bingo = 0;
			break;
		}
		Yprev = Ycur;
	}

	return bingo;
}

int CheckCurveProperties_n_shape(F32PTR Y, int N, int p0, BVS_XINFO_PTR xinfo) {
	int     Nlen = xinfo->Nunique[p0];
	I32PTR  indices = xinfo->SortedUnqiueIdx0 + N * p0;

	int bingo = 1;	
	for (int i = 1; i <= (Nlen - 2); i++) {
		int ipre = indices[i - 1];
		int icur = indices[i];
		int inxt = indices[i + 1];
		if (Y[icur] < Y[ipre] && Y[icur] < Y[inxt]) {
			bingo = 0;
			break;
		}
	}

	return bingo;
}

int CheckCurveProperties_u_shape(F32PTR Y, int N, int p0, BVS_XINFO_PTR xinfo) {

	int     Nlen = xinfo->Nunique[p0];
	I32PTR  indices = xinfo->SortedUnqiueIdx0 + N * p0;

	int bingo = 1;
	for (int i = 1; i <= (Nlen - 2); i++) {
		int ipre = indices[i - 1];
		int icur = indices[i];
		int inxt = indices[i + 1];
		if (Y[icur] > Y[ipre] && Y[icur] > Y[inxt]) {
			bingo = 0;
			break;
		}
	}

	return bingo;
}


#define VerifyMemHeader() MEM.verify_header(&MEM);
//#define VerifyMemHeader() ;


int beast2_main_corev4()  {

	// A struct to track allocated pointers   
	// const MemPointers MEM;
	// do not use 'const' bcz Clang will asssume other fields as zeros (e.g., alloc, and alloc0).
	MemPointers MEM = (MemPointers){.init = mem_init,};
	MEM.init(&MEM);
	MEM.checkHeader = 1;
	//mem = &MEM;
	GLOBAL_TMP = &MEM;

	// Get the Option parameters from the global pointer GLOBAL_OPTIONS, which was already 
	// set prior to this point (e.g., time-series length N, and number of pixels).	
	const BVS_OPTIONS_PTR	opt   = GLOBAL_OPTIONS;
	BVS_EXTRA         extra = opt->extra;
	const int  q = 1L;
	
	// Pre-allocate memory to save samples for calculating credibile intervals	
	CI_PARAM     ciParam = {0,};
	CI_RESULT    *ci     ;
	if (extra.computeCredible) {
		ci = MyALLOC(MEM, opt->prior.numBasis+1L, CI_RESULT, 0);
		ConstructCIStruct(opt->mcmc.credIntervalAlphaLevel, opt->mcmc.samples, opt->io.meta.Nobj *q,  //for MRBEAST
							opt->prior.numBasis + 1L,&MEM, &extra.fastCIComputation, &ciParam, ci );
	}

	// Allocate MEMORY FOR BASIS VARIABLE: Initialzie two pointers to BASIS
	BVS_MODEL  MODEL = {0,};	
	AllocInit_Model(&MODEL, opt, &MEM);

	VSLStreamStatePtr stream;
	//Initializing the random number generaotr	
	LOCAL( 	
		U64 seed = (opt->mcmc.seed == 0) ? TimerGetTickCount() : (opt->mcmc.seed+0x4f352a3dc);
		r_vslNewStream(&stream, VSL_BRNG_MT19937, seed);   	
	)

	const U32PTR  RND32        = MyALLOC(MEM, MAX_RAND_NUM,     U32, 64);
	const U16PTR  RND16        = MyALLOC(MEM, MAX_RAND_NUM * 2, U16, 64);
	const U08PTR  RND08        = MyALLOC(MEM, MAX_RAND_NUM * 4, U08, 64);	
	const F32PTR  RNDGAMMA     = MyALLOC(MEM, MAX_RAND_NUM,     F32, 64);	
	const U32PTR  RND32_END    = RND32		+ MAX_RAND_NUM - 50;
	const U16PTR  RND16_END    = RND16		+ MAX_RAND_NUM * 2 - 50;
	const U08PTR  RND08_END    = RND08		+ MAX_RAND_NUM * 4 - 50 -3;     //-3 bcz GenRandomBasis will also consume the stream bits
	const F32PTR  RNDGAMMA_END = RNDGAMMA	+ MAX_RAND_NUM - MODEL.nPrec-1L;
		
	// Allocate mem for current covariate/design matrix (Xt_mars), proposed new terms (Xnewterm),
	// and subset matrix corresponding to rows of missing values.		
	const F32PTR Xt_mars;
	const F32PTR Xnewterm;      //will be re-used as a temp mem for multiple purposes		
	Alloc_Xterms(&Xt_mars, &Xnewterm,opt,&MEM);

	

	// yInfo used to save the current time series to be processed
	BVS_XINFO     xInfo;
	BVS_YINFO     yInfo;
	Alloc_xInfo_yInfo(&xInfo, &yInfo, opt, MODEL.type2id, &MEM);
 

	// Allocate the output memory for a single chain (resultChain) and the averaged
	// result of all chains ( result)
	BVS_RESULT resultChain = { NULL,}, result={ NULL, };
	Alloc_Result(&resultChain, opt, &MEM); 	
	Alloc_Result(&result,      opt, &MEM);
	
	if (extra.computeCredible) {
		I32  N          = opt->io.meta.Nobj;//Correct for the inconsitency of X and Y in gemm and gemv		
		// the last componet is the sum of all the bases
		for (I32 i = 0; i < MODEL.NUMBASIS+1L; i++) {	 
				ci[i].result     = resultChain.CI+ 2*N*i;  
			    ci[i].newDataRow = resultChain.Y +  N * i;
		}
	} //NUMVAR_FOR_CI=3
 
 
	const BVS_HyperPar  hyperPar = {.alpha_1=opt->prior.alpha1, .alpha_2 = opt->prior.alpha2,.del_1 = opt->prior.delta1,  .del_2 = opt->prior.delta2};

	/****************************************************************************/
	//		THE OUTERMOST LOOP: Loop through all the time series one by one
	/****************************************************************************/
	// Get conversion factors from counts to seceonds
	InitTimerFunc();
	StartTimer();
	SetBreakPointForStartedTimer();
	

	//#define __DEBUG__
	#undef  __DEBUG__ 

	#ifdef __DEBUG__
		// Allocate a mem block and memset it to zero
		I32    N          = opt->io.meta.Nobj;
		I32    Npad       = (N + 7) / 8 * 8; Npad =  N;//Correct for the inconsitency of X and Y in gemm and gemv
		F32PTR flagSat    = MyALLOC0(MEM, N, I32, 64);
		F32PTR Xdebug     = MyALLOC0(MEM, Npad*(opt->prior.Kmax+ opt->prior.Kmax), I32, 64); // one Kmax for Xt_mars, and another for Xtbackup
	#endif


	//#define XtX_ByGroup XtX_ByGroup_FULL
	//#define MatxMat     MatxMat_FULL
	//#define MatxVec     MatxVec_FULL

	// Calculate the total number of time sereies to be processed
 
	const U32  MCMC_SAMPLES  = opt->mcmc.samples;
	const U32  MCMC_THINNING = opt->mcmc.thinningFactor;
	const U32  MCMC_BURNIN   = opt->mcmc.burnin;
	const U32  MCMC_CHAINNUM = opt->mcmc.chainNumber;

 
	//for (U32 pixelIndex = 1; pixelIndex <= NUM_PIXELS; pixelIndex++)
	{  

		// Fecth a new time-series: set up Y, nMissing,  n, rowsMissing		
		F32PTR MEMBUF           = Xnewterm; // Xnewterm is a temp mem buf.
		
		r_printf("Importing data from R and pre-processing the data ... \n");
		if (!ReadPreprocessInputData(opt, &xInfo, &yInfo, MEMBUF)) {
			r_vslDeleteStream(&stream);
			MEM.free_all(&MEM);
			r_printf("Error: Cannot read the input data ... \n");
			return 0;
		}		
		r_printf("Data imported and preprocessed... \nStart running the MCMC algorithm ... \n");

		// Print a blank line to be backspaced by the follow
		if (extra.printProgressBar) {
			F32 frac = 0.0; I32 firstTimeRun = 1;
			printProgress(frac, extra.consoleWidth, Xnewterm, firstTimeRun);
		}

		int skipCurrentPixel = 0;

		if (q == 1) {
				// alpha2_star  = alpha_2 + 0.5(YtY-beta*X'Y) = 0.5* (  [YtY+2*alpha_2] - beta*X'Y  )
				// YtY+2*alpha_2  is pre-cacluated here. THe "2" before alpha_2 is to undo the division
				// later in the calcution of alpha2_star
				yInfo.YtY_plus_alpha2Q[0] = yInfo.YtY_plus_alpha2Q[0] + 2 * hyperPar.alpha_2;
				//Pre-compute alpha1_star, which depends only on n.
				yInfo.alpha1_star         = yInfo.N * 0.5 + hyperPar.alpha_1;
		}	 
		/****************************************************************************/
		// GENERATE A STREAM OF RANDOM NUMBERS FOR FUTURE USE
		/****************************************************************************/
		BVS_RNDSTREAM  RND;
		{
			RND.rnd32 = RND32, RND.rnd16 = RND16, RND.rnd08 = RND08, RND.rndgamma = RNDGAMMA;
			//vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE,  stream, MAX_RAND_NUM, rnd, 0, 1.0);		
			r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, MAX_RAND_NUM, (U32PTR)RND32);
			r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, MAX_RAND_NUM, (U32PTR)RND16);
			r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, MAX_RAND_NUM, (U32PTR)RND08);
			// Depends on n and alpha_1 to generate rndgamma
			// Needed to resmaple Sig2	
			// Cann't be used to sample prec bcz the degree of freedom there depends on the number of terms Kterms
			r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, MAX_RAND_NUM, RNDGAMMA, ( hyperPar.alpha_1 + yInfo.N * 0.5f), 0, 1);
		}
		VerifyMemHeader();
		 
		// Clear up and zero out A(RESULT) for initialization	 
		Init_Result(&result, opt, 0);		

		VerifyMemHeader();
		BVS_TERMS_ENCODER  coder;
		AllocInit_TermEncoder(&coder, &MEM, opt);
		VerifyMemHeader();
		/****************************************************************************/
		//Iterate all the chains. The individual chain result will be saved into resultChain
		/****************************************************************************/
		for ( U32 chainNumber =0;  chainNumber < MCMC_CHAINNUM; chainNumber++)
		{
			const I32  N      = xInfo.N; 
			//const I32  Npad   = (N + 7) / 8 * 8; 
			const I32  Npad   = N;//Correct for the inconsitency of X and Y in gemm and gemv
			const I32  Npad16 = (N + 15) /16 * 16;	
			/****************************************************************************/
			//                 GENERATE AN INITIAL MODEL
			/****************************************************************************/		
			VerifyMemHeader();
			{   
				Init_BasisState_All(&MODEL, opt, &xInfo);
				// Generate random knots (numknot,KNOTs and ORDERS)								
				GenarateRandomBasis(&MODEL,&RND,&xInfo,opt);//CHANGE: nunKnot, ORDER, KNOT, K, KBase, Ks, Ke, or TERM_TYPE										
	
			 	BVS_EvaluateModel(&MODEL.curr, Xt_mars, MODEL.terms, MODEL.Kterms, &xInfo,&yInfo, &hyperPar, opt->prior.precValue, &stream);
			
				//MODEL.basisState[MODEL.type2id[LinearTerm]]->K  =  1;
				//MODEL.basisState[MODEL.type2id[LinearTerm]]->Kposition[0]=1L;

			}		
			VerifyMemHeader();
			{
				// Clear up and zero out resultChain for initialization
				Init_Result(&resultChain, opt, 0);

				// Reset all nots to ones bcz the real extrem positiosn won'tY be updated till samples > 1
				//memset(MODEL.extremePosVec, 1, N); 
			}
			VerifyMemHeader();
			/**********************************************************************************************/
			// PREPARE FOR THE START OF THE MAIN LOOP
			// The proposed basis needs XtX and XtY of the current basis to compute XtX_prop and XtY_prop.
			// Other varialbes (e.g., beta and cholXtX) are not cross-used by the basis and basis_prop
			/**********************************************************************************************/
			U32 ite            = 0;
			U32 sample         = 0;
			U32 subSampleIndex = 0;

			PROP_DATA PROPINFO = { .samples=&sample, .mem = Xnewterm, .model = &MODEL,.pRND =&RND, .nSample_ExtremVecNeedUpdate =1L,       
								   .sigFactor = opt->prior.sigFactor,  .Kmax  = opt->prior.Kmax, .xinfo=&xInfo, .prior=&opt->prior};

			I32            numBadIterations = 0;
			NewPropTerm    NEW;     // moved here bvz its xcols has two fixed lements, N and Npad
			NEWCOLINFOv2   NewCol; // moved here bvz its xcols has two fixed lements, N and Npad
			NewCol.N    = N;
			NewCol.Nlda = Npad;
			NewCol.X    = Xt_mars;
			NewCol.Xnewterm = Xnewterm;
			NEW.xcols   = &NewCol;
 
			while (sample < MCMC_SAMPLES)
			{
				ite++;
				/**********************************************************************/
				/*     Re-generate a pool of random numbers if almost used up          */
				/***********************************************************************/
				//vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, stream, MAX_RAND_NUM, rnd32, 0.f, 1.0f);
				if (RND.rnd32    >= RND32_END)    {r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, (RND.rnd32 - RND32), (U32PTR)RND32);	                               RND.rnd32    = RND32;   }
				if (RND.rnd16    >= RND16_END)    {r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, ((char*)RND.rnd16 - (char*)RND16 + 3) / sizeof(U32), (U32PTR)RND16); RND.rnd16    = RND16;   }
				if (RND.rnd08    >= RND08_END)    {r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, ((char*)RND.rnd08 - (char*)RND08 + 3) / sizeof(U32), (U32PTR)RND08); RND.rnd08    = RND08;   }
				if (RND.rndgamma >= RNDGAMMA_END) {r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE,      stream, MAX_RAND_NUM, RNDGAMMA, (hyperPar.alpha_1+yInfo.N*0.5f), 0.f, 1.f);    RND.rndgamma = RNDGAMMA;}
			
				// IMPLEMENT THE NEW PROPOSED BASIS		       
				VerifyMemHeader();
				NewCol.K = MODEL.curr.K;                      // Needed inside ProposeMove
				ProposeMove(&MODEL, &NEW, &NewCol, &PROPINFO);	
				get_parts_for_newinfo(&NewCol);
				VerifyMemHeader();
 
				//r_printf("%2d| %2d %2d %3d  %2d\n", NEW.term.type,NEW.jumpType,NEW.term.vars[0], NEW.term.knots[0], MODEL.curr.K);
					
				/**********************************************************************/
				// Generate the new terms for the propsed step: To add or move a bk, two segments 
				// are re-generated; to remove a bk or merge two bks, one segment are re-generated.	
				/**********************************************************************/								
				I32 Knewterm = 0;				
				for (I32 i = 0; i < NEW.nTerms; i++) { 
					// NEW.numSeg is 0 if removing terms for ChORDER(newOrder<oldeTerm)					
					GenTerms(NEW.terms+i, &xInfo, Xnewterm+i*N, 1L);
					Knewterm    += 1L;
				} // Iterate through all the new segments
				assert(Knewterm == NewCol.Knewterm); // Knewterm is the number of terms created not the extra col number
				VerifyMemHeader();

				/*************************************************************************/
				// Get XtX_prop: Copy parts of XtX to XtT_prop and fill new components 
				/*************************************************************************/		
		 
		 		update_XtX_from_Xnewterm_v2(MODEL.curr.XtX, MODEL.prop.XtX, &NewCol); 
				update_XtY_from_Xnewterm_v2(MODEL.curr.XtY, MODEL.prop.XtY, yInfo.Y, &NewCol, q);

	 			/*************************************************************/
				// XtX_prop has been constructed. Now use it to get the marg lik
				/***************************************************************/			
				if (1L) {		
					//Add precison values to the diagonal of XtX: post_P=XtX +diag(prec) 
					/*
					//Solve inv(Post_P)*XtY using  Post_P*b=XtY to get beta_mean
					//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, double * a, lapack_int lda);
					r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', KNEW, MODEL.prop.cholXtX, KNEW); // Choleskey decomposition; only the upper triagnle elements are used
					*/
					int KOLD = NewCol.K;
					int KNEW = NewCol.Knew;
					for (I32 i=1; i<NewCol.Kchol; ++i) 
						SCPY(i, MODEL.curr.cholXtX+(i-1)*KOLD, MODEL.prop.cholXtX+(i-1)*KNEW);
					VerifyMemHeader();
					//precFunc.UpdateXtXPrec_nTermsPerGrp(&MODEL, basis, &NEW, &NewCol); //&NEW is used only for OrderWise
					chol_addCol_skipleadingzeros_prec_nostartprec_invdiag(
										   MODEL.prop.XtX  + (NewCol.Kchol-1)*KNEW,
							               MODEL.prop.cholXtX,
							               MODEL.precVec,  KNEW, NewCol.Kchol, KNEW);
					VerifyMemHeader();
					//chol_full_v2(MODEL.prop.XtX, MODEL.prop.cholXtX, KNEW, KNEW);
			       /*
					for (rI32 i = 1; i <= (NEW.k1_new - 1L); i++) 	r_cblas_scopy(i, MODEL.curr.cholXtX + (i - 1L) * KOLD, 1L, MODEL.prop.cholXtX + (i - 1L) * KNEW, 1L);
					chol_addCol(MODEL.prop.cholXtX+ (NEW.k1_new - 1L)*KNEW, MODEL.prop.cholXtX, KNEW, NEW.k1_new, KNEW);
					//chol_addCol(MODEL.prop.cholXtX + (1 - 1L) * KNEW, MODEL.prop.cholXtX, KNEW, 1, KNEW);
					*/

					/*{
					//LAPACKE_dpotrs (int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );			
					SCPY(KNEW, MODEL.prop.XtY, MODEL.prop.beta_mean);
					r_LAPACKE_spotrs(LAPACK_COL_MAJOR, 'U', KNEW, 1, MODEL.prop.cholXtX, KNEW, MODEL.prop.beta_mean, KNEW);
					}*/

					MODEL.prop.K = KNEW;  // K is needed for computing marg_Lik
					// In MODEL, basis's K is still the old ones. They will be updated only if the proposal is accetped
					// via basis->CalcBasisKsKeK_TermType(basis). The Ks of the proposed bases are stored in MODEL.prop.ntermP				   
				   VerifyMemHeader();
				   ComputeMargLik(&MODEL.prop, &MODEL, &yInfo, &hyperPar);
				   //if (MODEL.prop.marg_lik != MODEL.prop.marg_lik || fabs(MODEL.prop.marg_lik )>FLOAT_MAX || MODEL.prop.alpha2_star <0.f) {					   				  
				   VerifyMemHeader();
				   if ( IsNaN(MODEL.prop.marg_lik) || IsInf(MODEL.prop.marg_lik ) ) {
					   if (++numBadIterations < 20) {					    	   
						   IncreasePrecValues(&MODEL);
						   chol_addCol_skipleadingzeros_prec_nostartprec_invdiag(MODEL.curr.XtX, MODEL.curr.cholXtX, MODEL.precVec, MODEL.curr.K, 1L, MODEL.curr.K);
						   ComputeMargLik(&MODEL.curr, &MODEL, &yInfo, &hyperPar); 
						   continue;
					   }  else {
						   skipCurrentPixel = 2;
						   break;
					   }					   
				   }  else {
					   numBadIterations = 0;
				   } //if (marg_lik_prop != marg_lik_prop || alpha2_star_prop <0.f) 

				   if (q == 1) {// added for MRBEAST
					   MODEL.prop.alpha2_star = max(MODEL.prop.alpha2_star,  MIN_ALPHA2_VALUE);
				   }
				}
			   /****************************************************************************************/
			   /*    DETERMINE WHETHER OR NOT TO ACCCEPT THE PROPSOED STEP                              */
			   /****************************************************************************************/
			
	 

				// First, calcuate a factor adjusting the likelihood change		 

				F32 delta_lik = MODEL.prop.marg_lik - MODEL.curr.marg_lik ;
				
				//acceptTheProposal = *(RND.rnd16)++ < fastexp(delta_lik) * 65535.0f;				
				I08     acceptTheProposal;
				if      (delta_lik >   0)   acceptTheProposal = 1;
				else if (delta_lik < -23) 	acceptTheProposal = 0;				
				else {				 
					F32    expValue = fastexp(delta_lik);
					if     (delta_lik > -0.5) 	acceptTheProposal = *(RND.rnd08)++ < expValue * 255.0f;
					else if(delta_lik > -5  )   acceptTheProposal = *(RND.rnd16)++ < expValue * 65535.0f;
					else						acceptTheProposal = *(RND.rnd32)++ < expValue * 4.294967296e+09;					 
				}

				if(acceptTheProposal){

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					int bingo = 1;
					F32PTR Y = Xnewterm + N * 2; // The firest two colums still contain the newly-generated terms
					for (int i = 0; i < opt->prior.num_con1d_Decrease; i++) {
						int pidx0 = opt->prior.con1d_Decrease[i]-1L;
						ComputeMarginalY_ProposedModel(&MODEL, &NEW, &NewCol, xInfo.xi_Pindices0, pidx0, Y, N);
						bingo =CheckCurveProperties_Decreasing(Y, N, pidx0, &xInfo);
						if (bingo == 0) break;
					}
					if (bingo == 0) continue;

					for (int i = 0; i < opt->prior.num_con1d_Increase; i++) {
						int pidx0 = opt->prior.con1d_Increase[i] - 1L;
						ComputeMarginalY_ProposedModel(&MODEL, &NEW, &NewCol, xInfo.xi_Pindices0, pidx0, Y, N);
						bingo = CheckCurveProperties_Increasing(Y, N, pidx0, &xInfo);
						if (bingo == 0) break;
					}
					if (bingo == 0) continue;

					for (int i = 0; i < opt->prior.num_con1d_nShape; i++) {
						int pidx0 = opt->prior.con1d_nShape[i] - 1L;
						ComputeMarginalY_ProposedModel(&MODEL, &NEW, &NewCol, xInfo.xi_Pindices0, pidx0, Y, N);
						bingo = CheckCurveProperties_n_shape(Y, N, pidx0, &xInfo);
						if (bingo == 0) break;
					}
					if (bingo == 0) continue;

					for (int i = 0; i < opt->prior.num_con1d_uShape; i++) {
						int pidx0 = opt->prior.con1d_uShape[i] - 1L;
						ComputeMarginalY_ProposedModel(&MODEL, &NEW, &NewCol, xInfo.xi_Pindices0,  pidx0, Y, N);
						bingo = CheckCurveProperties_u_shape(Y, N, pidx0, &xInfo);
						if (bingo == 0) break;
					}
					if (bingo == 0) continue;

					for (int i = 0; i < opt->prior.num_con1d_NoWiggle; i++) {
						int pidx0 = opt->prior.con1d_NoWiggle[i] - 1L;
						ComputeMarginalY_ProposedModel(&MODEL, &NEW, &NewCol, xInfo.xi_Pindices0, pidx0, Y, N);
						bingo = CheckCurveProperties_u_shape(Y, N, pidx0, &xInfo) ||
							    CheckCurveProperties_n_shape(Y, N, pidx0, &xInfo);
						if (bingo == 0) break;
					}
					if (bingo == 0) continue;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					VerifyMemHeader();
					// Inserting XnewTerms into Xt_mars
					swap_cols_bands_within_matrx(&NewCol);

					/****************************************************/
					//Find the good positions of the proposed MOVE
					//Then update the knotLists and order
					/****************************************************/		
					VerifyMemHeader();
					NEW.KpositionBuf = Xnewterm + NewCol.Knewterm * N;
					Update_BasisState_All(&NEW, &MODEL,&xInfo); //xInfo.Nunique is nseed inside: MODEL.ktermks is updated here
					VerifyMemHeader();

					//basis->CalcBasisKsKeK_TermType(basis);
					//UpdateBasisKbase(MODEL.b, MODEL.NUMBASIS, basis-MODEL.b);//basisIdx=basis-b Re-compute the K indices of the bases after the basisID 
	
					//Switching between Basis and Basis_prop
					{
						//http: //stackoverflow.com/questions/3647331/how-to-swap-two-numbers-without-using-temp-variables-or-arithmetic-operations
						//basis  = ((I64)basis ^ (I64)basis_prop);basis_prop = ((I64)basis ^ (I64)basis_prop); basis      =  ((I64)basis ^ (I64)basis_prop);						
					
						#define Exchange(x,y)   {void * _restrict tmp; tmp=MODEL.x;  MODEL.x=MODEL.y; MODEL.y=tmp;}
						Exchange(curr.XtX,       prop.XtX);
						Exchange(curr.XtY,       prop.XtY);
						Exchange(curr.beta_mean, prop.beta_mean);
						//Exchange(curr.beta,      prop.beta); //no need to exchange
						Exchange(curr.cholXtX,          prop.cholXtX);
						//Exchange(curr.precXtXDiag,      prop.precXtXDiag);
						//Exchange(curr.nTermsPerPrecGrp, prop.nTermsPerPrecGrp);   //needed for compontwise and orderwise

						if (q == 1) MODEL.curr.alpha2_star = MODEL.prop.alpha2_star; //BEASTV4
						else 	    Exchange(curr.alphaQ_star, prop.alphaQ_star);    //MRBEAST 
						
						MODEL.curr.marg_lik    = MODEL.prop.marg_lik;
						MODEL.curr.K		   = MODEL.prop.K;  //GetNumOfXmarCols(&MODEL): this function should also give KNEW; if not, there must be something wrong!
						#undef Exchange											
					}
					#ifdef __DEBUG__	
					if (q == 1) {
						BVS_EvaluateModel(&MODEL.prop, Xdebug, MODEL.terms, MODEL.Kterms, &xInfo,
							&yInfo, &hyperPar, MODEL.precVec[0], &stream);	

						for (int i = 1; i < MODEL.Kterms; i++) {
							F32 sum=0;
							for (int j = 0; j < N; j++) {
								sum += fabsf(Xt_mars[i * N + j] - Xdebug[i * N + j]);
							}
							if (sum > 0.000001) {
								int a = 1;
							}
						}
						if (MODEL.curr.marg_lik != MODEL.prop.marg_lik) {
							int a = 1;
						}
						r_printf("ite:%d K:%d|%f|%f|diff:%f\n", ite, MODEL.curr.K, MODEL.curr.marg_lik, MODEL.prop.marg_lik, MODEL.prop.marg_lik - MODEL.curr.marg_lik);
					    r_printf(" %f[%f]-%f %f\n", (MODEL.curr.alpha2_star), (MODEL.prop.alpha2_star),	yInfo.alpha1_star* (log(MODEL.curr.alpha2_star) - log(MODEL.prop.alpha2_star)), yInfo.alpha1_star);
 					}
					#endif


				} //(*rnd32++ < exp(marg_lik_prop - basis->marg_lik))

				VerifyMemHeader();

				/****************************************************************************************/
				//
				// For buin-in iterations, no need to  make posterior inference.
				//
				/****************************************************************************************/

				U08 bResampleParameter  = (ite % 20 == 0) && ite > 200;
				U08 bStoreCurrentSample = (ite % MCMC_THINNING == 0);

				/****************************************************************************************/
				//    First, Re-SAMPLING SIG2
				/****************************************************************************************/
				if (bResampleParameter || bStoreCurrentSample) {
					//vdRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, &sig2, (alpha_1+n/2), 0, 1.0/(alpha_1+basis->alpha2_star *0.5));
					F32  sig2_inv = (*RND.rndgamma++) * 1.f / MODEL.curr.alpha2_star;
					MODEL.sig2[0] = 1.0f / sig2_inv;
					MODEL.sig2[0] = max(MODEL.sig2[0], MIN_SIG2_VALUE);
					//r_printf("ite-%d SIG %f %f\n", ite, MODEL.sig2,  MODEL.sig2*yInfo.sd*yInfo.sd);			 
				}
		
				/****************************************************************************************/
				//  Re-sample beta to be used for either re-sampling prec (ite%20=0) or predicting Y (ite%thiningFactor=0)
				/****************************************************************************************/
				if (bResampleParameter || (bStoreCurrentSample && extra.useMeanOrRndBeta)) {

					//Compute beta = beta_mean + Rsig2 * randn(p, 1);
					//Usig2 = (1 / sqrt(sig2)) * U; 		beta = beta_mean + linsolve(Usig2, randn(p, 1), opts);
					//status = vdRngGaussian( method, stream, n, r, a, sigma );
					I32 K = MODEL.curr.K;
					r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, K, MODEL.beta, 0, 1);
					//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, 1, MODEL.curr.cholXtX, K, MODEL.curr.beta, K); // LAPACKE_strtrs (int matrix_layout , char uplo , char trans , char diag , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
					solve_U_as_U_invdiag(MODEL.curr.cholXtX, MODEL.beta, K, K);
					r_ippsMulC_32f_I(fastsqrt(MODEL.sig2[0]), MODEL.beta, K);
					r_ippsAdd_32f_I(MODEL.curr.beta_mean, MODEL.beta, K);

				}

				/************************************************************/
				// Re-sample the precison parameters and re-calcuate marg_lik and beta
				/************************************************************/
				if (bResampleParameter && q==1 && MODEL.curr.K >1) 	{
 
					I32 K    = MODEL.curr.K;
					F32 sumq = DOT(K-1, MODEL.beta+1L, MODEL.beta+1L); // skip the first beta correspoding to the const term
					F32 newPrecVal;
					r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1L, &newPrecVal, (hyperPar.del_1+(K-1)*0.5f), 0.f, 1.f);
					newPrecVal          = newPrecVal / (hyperPar.del_2 + 0.5f * sumq / MODEL.sig2[0]);
					MODEL.precVec[0]    = newPrecVal > MIN_PREC_VALUE ? newPrecVal : MODEL.precVec[0];
					MODEL.logPrecVec[0] = logf(MODEL.precVec[0]);
			
					I32 ntries = 0;
					do {
						if (ntries++ > 0){IncreasePrecValues(&MODEL);	}						
						chol_addCol_skipleadingzeros_prec_nostartprec_invdiag(MODEL.curr.XtX, MODEL.curr.cholXtX, MODEL.precVec, MODEL.curr.K, 1L, MODEL.curr.K);
						ComputeMargLik( &MODEL.curr, &MODEL, &yInfo, &hyperPar);
					} while (  IsNaN(MODEL.curr.marg_lik) && ntries < 20 );

					if ( IsNaN(MODEL.curr.marg_lik) ) { 
						skipCurrentPixel = 3;
						break;
					} 

					/* No need to re-sample beta because it is not really used
					r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, K, beta, 0, 1);
					// LAPACKE_strtrs (int matrix_layout , char uplo , char trans , char diag , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
					r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, 1, cholXtX, K, beta, K);
					r_ippsMulC_32f_I(sqrtf(modelPar.sig2), beta, K);
					r_ippsAdd_32f_I(beta_mean, beta, K);
					*/
					/* /FInally, re-sample sig2 based on the lastest alpha2_star
					//vdRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, &sig2, (alpha_2+n/2), 0, 1.0/(alpha_1+basis->alpha2_star *0.5));
					modelPar.sig2 = (*rndgamma++)*1.0f / (modelPar.alpha_1 + basis->alpha2_star *0.5f);
					modelPar.sig2 = 1.f / modelPar.sig2;
					*/
				}
				VerifyMemHeader();
				if (extra.printProgressBar && ite % 1000 == 0) {
					I64 IterationPerChain = (MCMC_SAMPLES * MCMC_THINNING + MCMC_BURNIN) ;
					F32 frac = (F32)(chainNumber * IterationPerChain + ite) / (IterationPerChain * MCMC_CHAINNUM);
					printProgress(frac, extra.consoleWidth, Xnewterm, 0);
				}
				VerifyMemHeader();

				if (ite <= MCMC_BURNIN)   continue;
				if (!bStoreCurrentSample) continue;

				sample++;



				/**********************************************/
				//
				//      Start to compute final results
				//
				/**********************************************/

				*resultChain.marg_lik += MODEL.curr.marg_lik;
			 	*resultChain.sig2     += MODEL.sig2[0];
				*((I32PTR)resultChain.Kmean)  += MODEL.curr.K;
			

				{

					////////////////////////////////////////////////////////////////////////////////////
					F32PTR MEMBUF1   = Xnewterm;					  						
					//Counting probability of being breakpoints	
					I08PTR included = MEMBUF1;
					memset(included, 0, xInfo.Pall);
					for (I32 i = 1; i < MODEL.Kterms; i++) {
						I32 vidx0 = MODEL.terms[i].var0s[0];
						I32 type  = MODEL.terms[i].type;
						I32 pidx0 = xInfo.xi_Pindices0[type][vidx0];
						included[pidx0] = 1;
						if (type >= Hinge2DTerm) {
							vidx0 = MODEL.terms[i].var0s[1];
							pidx0 = xInfo.xi_Pindices0[type][vidx0];
							included[pidx0] = 1;
						}
					}
					for (int i = 0; i < xInfo.Pall; i++) {
						resultChain.probAll[i] += included[i];
					}
					////////////////////////////////////////////////////////////////////////////////////
  

					F32PTR BETA = extra.useMeanOrRndBeta == 0 ? MODEL.curr.beta_mean : MODEL.beta;
					int type;
					VerifyMemHeader();
					////////////////////////////////////////////////////////////////////////////////////
					type = LinearTerm;
					if (MODEL.type2id[type] >= 0) {
						I32    K         = MODEL.basisState[MODEL.type2id[type]].K;
						I16PTR Kposition = MODEL.basisState[MODEL.type2id[type]].Kposition;
						for (I32 i = 0; i < K; i++) {
							I32 vidx0 = MODEL.terms[Kposition[i] - 1L].var0;							
							resultChain.probLin[vidx0]++;
						}
		 						/*************************************************************************/
						F32PTR Ydst = resultChain.Ylin;
						for (int i = 0; i < K; i++) {
							I32 k0    =  Kposition[i] - 1L;
							I32 vidx0 = MODEL.terms[k0].var0;
							*(Ydst + vidx0) += BETA[k0];							
						}
 
					}
					VerifyMemHeader();
					type = StepTerm;
					if (MODEL.type2id[type] >= 0) {
						I32    K         = MODEL.basisState[MODEL.type2id[type]].K;
						I16PTR Kposition = MODEL.basisState[MODEL.type2id[type]].Kposition;
						for (I32 i = 0; i < K; i++) {
							I32 vidx0 = MODEL.terms[Kposition[i] - 1L].var0;
							I32 knot0 = MODEL.terms[Kposition[i] - 1L].knot- 1L;
							if (knot0 == 0) continue;
							resultChain.probStep[vidx0 *N+knot0]++;
						}
						/*************************************************************************/
						F32PTR Ydst = resultChain.Ystep;
						for (int i = 0; i < K; i++) {
							I32 k0 = Kposition[i] - 1L;
							I32 vidx0 = MODEL.terms[k0].var0;
							f32_axpy_inplace(BETA[k0], Xt_mars + N * k0, Ydst + vidx0 * N, N);
						}

					}
					type = StairTerm;
					if (MODEL.type2id[type] >= 0) {
						I32    K         = MODEL.basisState[MODEL.type2id[type]].K;
						I16PTR Kposition = MODEL.basisState[MODEL.type2id[type]].Kposition;
						for (I32 i = 0; i < K; i++) {
							I32 vidx0 = MODEL.terms[Kposition[i] - 1L].var0;
							I32 knot0 = MODEL.terms[Kposition[i] - 1L].knot- 1L;
							if (knot0 == 0) continue;
							resultChain.probStair[vidx0 *N+knot0]++;
						}
						/*************************************************************************/
						F32PTR Ydst = resultChain.Ystair;
						for (int i = 0; i < K; i++) {
							I32 k0 = Kposition[i] - 1L;
							I32 vidx0 = MODEL.terms[k0].var0;
							f32_axpy_inplace(BETA[k0], Xt_mars + N * k0, Ydst + vidx0*N, N);
						}

					}
					type = PieceTerm;
					if (MODEL.type2id[type] >= 0) {
						I32    K         = MODEL.basisState[MODEL.type2id[type]].K;
						I16PTR Kposition = MODEL.basisState[MODEL.type2id[type]].Kposition;
						for (I32 i = 0; i < K; i++) {
							I32 vidx0 = MODEL.terms[Kposition[i] - 1L].var0;
							I32 knot0 = MODEL.terms[Kposition[i] - 1L].knot- 1L;
							if (knot0 == 0) continue;
							resultChain.probPiece[vidx0 *N+knot0]++;
						}
						/*************************************************************************/
						F32PTR Ydst = resultChain.Ypiece;
						for (int i = 0; i < K; i++) {
							I32 k0 = Kposition[i] - 1L;
							I32 vidx0 = MODEL.terms[k0].var0;
							f32_axpy_inplace(BETA[k0], Xt_mars + N * k0, Ydst + vidx0 * N, N);
						}

					}
					type = SeasonTerm;
					if (MODEL.type2id[type] >= 0) {
						I32    K         = MODEL.basisState[MODEL.type2id[type]].K;
						I16PTR Kposition = MODEL.basisState[MODEL.type2id[type]].Kposition;
						for (I32 i = 0; i < K; i++) {
							I32 vidx0 = MODEL.terms[Kposition[i] - 1L].var0;
							I32 knot0 = MODEL.terms[Kposition[i] - 1L].knot- 1L;
							if (knot0 == 0) continue;
							resultChain.probSeason[vidx0 *N+knot0]++;
						}
						/*************************************************************************/
						F32PTR Ydst = resultChain.Yseason;
						for (int i = 0; i < K; i++) {
							I32 k0 = Kposition[i] - 1L;
							I32 vidx0 = MODEL.terms[k0].var0;
							f32_axpy_inplace(BETA[k0], Xt_mars + N * k0, Ydst + vidx0 * N, N);
						}

					}
					VerifyMemHeader();
					type = HingeTerm;
					if (MODEL.type2id[type] >= 0  ) {
						I32    K         = MODEL.basisState[MODEL.type2id[type]].K;
						I16PTR Kposition = MODEL.basisState[MODEL.type2id[type]].Kposition;
						for (I32 i = 0; i < K; i++) {
							I32 vidx0 = MODEL.terms[Kposition[i] - 1L].var0;
							I32 knot0 = MODEL.terms[Kposition[i] - 1L].knot- 1L;
							resultChain.probHinge[vidx0 *N+knot0]++;
						}
						/*************************************************************************/
						F32PTR Ydst = resultChain.Yhinge;
						for (int i = 0; i < K; i++) {
							I32 k0 = Kposition[i] - 1L;
							I32 vidx0 = MODEL.terms[k0].var0;
							f32_axpy_inplace(BETA[k0], Xt_mars + N * k0, Ydst + vidx0 * N, N);
						}

					}
					VerifyMemHeader();
					type = HingePairTerm;
					if (MODEL.type2id[type] >= 0) {
						I32    K = MODEL.basisState[MODEL.type2id[type]].K;
						I16PTR Kposition = MODEL.basisState[MODEL.type2id[type]].Kposition;
						for (I32 i = 0; i < K; i++) {
							I32 vidx0 = MODEL.terms[Kposition[i] - 1L].var0;
							I32 knot0 = MODEL.terms[Kposition[i] - 1L].knot - 1L;
							resultChain.probHingePair[vidx0 * N + knot0]++;
						}
						/*************************************************************************/
						F32PTR Ydst = resultChain.YhingePair;
						for (int i = 0; i < K; i++) {
							I32 k0 = Kposition[i] - 1L;
							I32 vidx0 = MODEL.terms[k0].var0;
							f32_axpy_inplace(BETA[k0], Xt_mars + N * k0, Ydst + vidx0 * N, N);
						}

					}
					VerifyMemHeader();
					type = ChangeTerm;
					if (MODEL.type2id[type] >= 0) {
						I32    K = MODEL.basisState[MODEL.type2id[type]].K;
						I16PTR Kposition = MODEL.basisState[MODEL.type2id[type]].Kposition;
						for (I32 i = 0; i < K; i++) {
							I32 vidx0 = MODEL.terms[Kposition[i] - 1L].var0;
							I32 knot0 = MODEL.terms[Kposition[i] - 1L].knot - 1L;
							resultChain.probChange[vidx0 * N + knot0]++;
						}
						/*************************************************************************/
						F32PTR Ydst = resultChain.Ychange;
						for (int i = 0; i < K; i++) {
							I32 k0 = Kposition[i] - 1L;
							I32 vidx0 = MODEL.terms[k0].var0;
							f32_axpy_inplace(BETA[k0], Xt_mars + N * k0, Ydst + vidx0 * N, N);
						}

					}

					I32 K = MODEL.Kterms;			 
					resultChain.KProb[K - 1L]++;
	 
					//Compute the averaged  signals
					//r_cblas_sgemv(CblasColMajor, CblasNoTrans, Npad, K, 1.f, X, Npad,beta, 1L, 0.f,	Y, 1L);
					r_cblas_sgemv(CblasColMajor, CblasNoTrans, Npad, K, 1.f, Xt_mars, Npad, MODEL.curr.beta_mean, 1L, 0.f, MEMBUF1, 1L);
					//basis->ComputeY(Xt_mars, BETA, MEMBUF1, basis, Npad);
					f32_add_v_v2_vec_inplace(MEMBUF1, resultChain.Y, resultChain.SD, N);
					MEMBUF1 += Npad;			 
				}
				
				VerifyMemHeader();
				/******************************************************************/
				//Inser the terms into the term List (the first constant term is skipped via MODEL.terms+1L)
				/******************************************************************/
				InsertNewTermList(MODEL.terms+1L, MODEL.Kterms-1, &coder, MODEL.curr.beta_mean);
		
				VerifyMemHeader();
	
			    /*************************************************/
				// Compute ci intervals: new row of data have already calculated
				// and saved in Xnewterm, Xnewterm+Npad, and Xnweterm+2*Npad
				/*************************************************/
				if (extra.computeCredible && 0x000000000) {

					// when  *RND.rnd16++ <= ciParam.subsampleFraction_x_INT16MAX, samples are included;					
					if (extra.fastCIComputation && !(*RND.rnd16++ < ciParam.subsampleFraction_x_INT16MAX)) {
						//if (*rnd32++ < subsampleFraction*4.294967296000000e+09)
						// The current sample not included. No need to insert it into the ci strips.
						// So, just skip to the next iteration	
					} else {
						// New row of data for slope, seasonal, and trend components: MEMBUF1=slope over time	 					    
						for (int i = 0; i < MODEL.NUMBASIS; i++)
							InsertNewRowToUpdateCI(&ciParam, &ci[i]);
					}

				} // if (extra.computeCredible)

			}//WHILE(sample<SAMPLE)

			/*****************************************************************************/
			//
			// One chain is done; then start post-processing the result of the chain     
			//
			/*****************************************************************************/

			if (!skipCurrentPixel)	{

				F32   inv_sample = 1.f / sample;
				int   sum;
				#define GetSum(arr) (r_ippsSum_32s_Sfs(arr, N, &sum, 0), sum) // parenthesis operator
				int   P = xInfo.Pall;
				int   N = xInfo.N;

				*resultChain.marg_lik  *= inv_sample;
				// FOR MRBEAST
				for (int col = 0; col < q; col++)	{ 
					for (int i = 0; i < q; i++) {
						resultChain.sig2[col*q+i] = resultChain.sig2[col*q+i] * inv_sample * yInfo.sd[col] * yInfo.sd[i];
					}
				}
				VerifyMemHeader();
				i32_to_f32_scaleby_inplace(resultChain.Kmean,  1, inv_sample);
				i32_to_f32_scaleby_inplace(resultChain.KProb,  opt->prior.Kmax , inv_sample);
				
		 
				I32  Psqr = opt->prior.Pbasis[QuadraticTerm]; 
				I32  Plin = opt->prior.Pbasis[LinearTerm];
				I32  Pstep = opt->prior.Pbasis[StepTerm];
				I32  Pstair = opt->prior.Pbasis[StairTerm];
				I32  Ppiece = opt->prior.Pbasis[PieceTerm];
				I32  Pseason = opt->prior.Pbasis[SeasonTerm];
				I32  Phinge = opt->prior.Pbasis[HingeTerm];
				I32  Phingepair = opt->prior.Pbasis[HingePairTerm];
				I32  Pchange = opt->prior.Pbasis[ChangeTerm];


				i32_to_f32_scaleby_inplace(resultChain.probAll, P, inv_sample);

				if (MODEL.type2id[LinearTerm] >= 0) {
					i32_to_f32_scaleby_inplace(resultChain.probLin, 1 * Plin, inv_sample);
					f32_mul_val_inplace(yInfo.sd[0] * inv_sample, resultChain.Ylin, Plin);
				}
				if (MODEL.type2id[StepTerm] >= 0) {
					i32_to_f32_scaleby_inplace(resultChain.probStep, N*Pstep, inv_sample);
					f32_scale_inplace(yInfo.sd[0] * inv_sample, yInfo.mean[0], resultChain.Ystep, N* Pstep);
				}
				if (MODEL.type2id[StairTerm] >= 0) {
					i32_to_f32_scaleby_inplace(resultChain.probStair, N * Pstair, inv_sample);
					f32_mul_val_inplace(yInfo.sd[0] * inv_sample, resultChain.Ystair, N* Pstair);
				}
				if (MODEL.type2id[PieceTerm] >= 0) {
					i32_to_f32_scaleby_inplace(resultChain.probPiece, N * Ppiece, inv_sample);
					f32_mul_val_inplace(yInfo.sd[0] * inv_sample, resultChain.Ypiece, N* Ppiece);
				}
				if (MODEL.type2id[SeasonTerm] >= 0) {
					i32_to_f32_scaleby_inplace(resultChain.probSeason, N * Pseason, inv_sample);
					f32_mul_val_inplace(yInfo.sd[0] * inv_sample, resultChain.Yseason, N * Pseason);
				}
				if (MODEL.type2id[HingeTerm] >= 0) {
					i32_to_f32_scaleby_inplace(resultChain.probHinge, N * Phinge, inv_sample);
					f32_mul_val_inplace(yInfo.sd[0] * inv_sample, resultChain.Yhinge, N* Phinge);
				}
				if (MODEL.type2id[HingePairTerm] >= 0) {
					i32_to_f32_scaleby_inplace(resultChain.probHingePair, N * Phingepair, inv_sample);
					f32_mul_val_inplace(yInfo.sd[0] * inv_sample, resultChain.YhingePair, N* Phingepair);
				}
				if (MODEL.type2id[ChangeTerm] >= 0) {
					i32_to_f32_scaleby_inplace(resultChain.probChange, N * Pchange, inv_sample);
					f32_mul_val_inplace(yInfo.sd[0] * inv_sample, resultChain.Ychange, N* Pchange);
				}

		
				//FOR MRBEAST
				for (int i = 0; i < q; i++) {
					F32 offset = 0.0f;
					f32_sx_sxx_to_avgstd_inplace(resultChain.Y + i * N, resultChain.SD + i * N, sample, yInfo.sd[i], yInfo.mean[i], N);
				}

				if (extra.computeCredible &0x000000000000) {
					//FOR MRBEAST
					for (int i = 0; i < q; i++) {
						f32_scale_inplace(yInfo.sd[i], yInfo.mean[i], resultChain.CI + N * i, N);
						f32_scale_inplace(yInfo.sd[i], yInfo.mean[i], resultChain.CI + N * q + N * i, N);
						//r_ippsMulC_32f_I(,    resultChain.tCI+(2*N)*i, N + N),
						//r_ippsSubC_32f_I(-yInfo.mean[i], ); //ippsAddC_32f_I(yInfo.mean, result.tCI, N + N);
					}
				}
 

			}// Finish computing the result of the single chiain

			VerifyMemHeader();
			/**************************************************/
			//   Add up the individual chain to the Result
			/**************************************************/
			if (!skipCurrentPixel) 
			{
				//https://stackoverflow.com/questions/13216423/error-pasting-and-red-does-not-give-a-valid-preprocessing-token
				//#define _1(x)      *(result.##x) += *(resultChain.##x) // Working with MSVC but not GCC

				#define _1(x)      *(result.x) += *(resultChain.x)
				#define _N(x)      r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N)
				#define _P(x)      r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, xInfo.Pall)
				#define _Nq(x)     r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N*q)
				#define _q(x)      r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, q)
				#define _q2(x)     r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, q*q)
				#define _2N(x)     r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N+N)
				#define _2Nq(x)    r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N*q+N*q)
				#define _kmax(x)   r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, opt->prior.Kmax)
				#define _NP(x)     r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N*xInfo.Pall)

				_1(marg_lik);
				_q2(sig2);  //Fpr MRBEAST
				_1(Kmean);
								
				_N(Y); 
				_N(SD);
				_kmax(KProb);
				
				I32 P         = opt->prior.Pall;
		 
				I32  Psqr = opt->prior.Pbasis[QuadraticTerm]; 
				I32  Plin = opt->prior.Pbasis[LinearTerm];
				I32  Pstep = opt->prior.Pbasis[StepTerm];
				I32  Pstair = opt->prior.Pbasis[StairTerm];
				I32  Ppiece = opt->prior.Pbasis[PieceTerm];
				I32  Pseason = opt->prior.Pbasis[SeasonTerm];
				I32  Phinge = opt->prior.Pbasis[HingeTerm];
				I32  Phingepair = opt->prior.Pbasis[HingePairTerm];
				I32  Pchange = opt->prior.Pbasis[ChangeTerm];

				_P(probAll);
  
				if (MODEL.type2id[LinearTerm] >= 0) {
					r_ippsAdd_32f_I((F32PTR)resultChain.probLin, (F32PTR)result.probLin, 1 * Plin);
					r_ippsAdd_32f_I((F32PTR)resultChain.Ylin, (F32PTR)result.Ylin, Plin);
 
				}

				if (MODEL.type2id[StepTerm] >= 0) {
					r_ippsAdd_32f_I((F32PTR)resultChain.probStep, (F32PTR)result.probStep, N*Pstep);		
					r_ippsAdd_32f_I((F32PTR)resultChain.Ystep, (F32PTR)result.Ystep, N* Pstep);
				}
				if (MODEL.type2id[StairTerm] >= 0) {
					r_ippsAdd_32f_I((F32PTR)resultChain.probStair, (F32PTR)result.probStair, N * Pstair);
					r_ippsAdd_32f_I((F32PTR)resultChain.Ystair, (F32PTR)result.Ystair, N* Pstair);
				}
				if (MODEL.type2id[PieceTerm] >= 0) {
					r_ippsAdd_32f_I((F32PTR)resultChain.probPiece, (F32PTR)result.probPiece, N * Ppiece);
					r_ippsAdd_32f_I((F32PTR)resultChain.Ypiece, (F32PTR)result.Ypiece, N* Ppiece);
				}
				if (MODEL.type2id[SeasonTerm] >= 0) {
					r_ippsAdd_32f_I((F32PTR)resultChain.probSeason, (F32PTR)result.probSeason, N * Pseason);
					r_ippsAdd_32f_I((F32PTR)resultChain.Yseason, (F32PTR)result.Yseason, N * Pseason);
				}
				if (MODEL.type2id[HingeTerm] >= 0) {
					r_ippsAdd_32f_I((F32PTR)resultChain.probHinge, (F32PTR)result.probHinge, N* Phinge);	
					r_ippsAdd_32f_I((F32PTR)resultChain.Yhinge, (F32PTR)result.Yhinge, N* Phinge);
				}
				if (MODEL.type2id[HingePairTerm] >= 0) {
					r_ippsAdd_32f_I((F32PTR)resultChain.probHingePair, (F32PTR)result.probHingePair, N * Phingepair);
					r_ippsAdd_32f_I((F32PTR)resultChain.YhingePair, (F32PTR)result.YhingePair, N* Phingepair);
				}
				if (MODEL.type2id[ChangeTerm] >= 0) {
					r_ippsAdd_32f_I((F32PTR)resultChain.probChange, (F32PTR)result.probChange, N* Pchange);					
					r_ippsAdd_32f_I((F32PTR)resultChain.Ychange, (F32PTR)result.Ychange, N* Pchange);
				}
			 


				#undef _1
				#undef _N
				#undef _Nq 
				#undef _q 
				#undef _q2
				#undef _2N
				#undef _2Nq
 
			}
			VerifyMemHeader();
			// Jump out of the chainumber loop
			if (skipCurrentPixel) {
				r_warning("\nWARNING(#%d):The max number of bad iterations exceeded. Can't decompose the current time series\n", skipCurrentPixel);
				break;
			}
			VerifyMemHeader();
		}
		/*********************************/
		// WHILE(chainNumber<chainNumber)
		/*********************************/
 

	/******************************************************/
	//
	// Finish all the chains and now acverage all of them
	//
	/******************************************************/

		// Average the results from multiple chains
		if (MCMC_CHAINNUM >= 2 && !skipCurrentPixel) 	{

			I32  N              = xInfo.N;						
			F32  invChainNumber = 1.f /(F32)MCMC_CHAINNUM;

			#define _1(x)      *((F32PTR)result.x)*=invChainNumber
			#define _N(x)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, N)
			#define _Nq(x)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x,N*q)
			#define _P(x)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, xInfo.Pall)
			#define _q(x)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x,q)
			#define _q2(x)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x,q*q)
			#define _2N(x)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, N+N)
			#define _2Nq(x)    r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, N*q+N*q)
			#define _kmax(x)   r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, opt->prior.Kmax)
			#define _NP(x)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, N*xInfo.Pall)

			F32 maxncpProb;	 
			_1(marg_lik);			
			_q2(sig2); 
			_1(Kmean);
			_N(Y);
			_N(SD);
			_kmax(KProb);

			_P(probAll);	

			if (MODEL.type2id[LinearTerm] >= 0)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.probLin, xInfo.xi_Pbasis[LinearTerm] * 1);
			if (MODEL.type2id[StepTerm] >= 0)       r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.probStep,  xInfo.xi_Pbasis[StepTerm] *N);
			if (MODEL.type2id[StairTerm] >= 0)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.probStair, xInfo.xi_Pbasis[StairTerm] * N);
			if (MODEL.type2id[PieceTerm] >= 0)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.probPiece, xInfo.xi_Pbasis[PieceTerm] * N);
			if (MODEL.type2id[SeasonTerm] >= 0)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.probSeason, xInfo.xi_Pbasis[SeasonTerm] * N);
			if (MODEL.type2id[HingeTerm] >= 0)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.probHinge, xInfo.xi_Pbasis[HingeTerm] * N);
			if (MODEL.type2id[HingePairTerm] >= 0)  r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.probHingePair, xInfo.xi_Pbasis[HingePairTerm] * N);
			if (MODEL.type2id[ChangeTerm] >= 0)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.probChange, xInfo.xi_Pbasis[ChangeTerm] * N);
		 
			if (MODEL.type2id[LinearTerm] >= 0)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.Ylin, xInfo.xi_Pbasis[LinearTerm] * 1);
			if (MODEL.type2id[StepTerm] >= 0)       r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.Ystep,  xInfo.xi_Pbasis[StepTerm] *N);
			if (MODEL.type2id[StairTerm] >= 0)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.Ystair, xInfo.xi_Pbasis[StairTerm] * N);
			if (MODEL.type2id[PieceTerm] >= 0)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.Ypiece, xInfo.xi_Pbasis[PieceTerm] * N);
			if (MODEL.type2id[SeasonTerm] >= 0)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.Yseason, xInfo.xi_Pbasis[SeasonTerm] * N);
			if (MODEL.type2id[HingeTerm] >= 0)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.Yhinge, xInfo.xi_Pbasis[HingeTerm] * N);
			if (MODEL.type2id[HingePairTerm] >= 0)  r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.YhingePair, xInfo.xi_Pbasis[HingePairTerm] * N);
			if (MODEL.type2id[ChangeTerm] >= 0)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.Ychange, xInfo.xi_Pbasis[ChangeTerm] * N);
 
 
			#undef _1
			#undef _N
			#undef _2N
			#undef _skn_1
			#undef _tkn_1
			#undef _okn_1
		}
 
 
		// Compute R2 and RMSE
		// At this point, yInfo.Y is not used any longer, so is re-used
		// to stote the sume of fitted S,T and/or O.
		// Xnewterm still contains the orignal data and shiuldn't not be touched
  		if (!skipCurrentPixel) { 
			
			/*
			I32  N  = opt->io.N ;
			I32  Nq = N * q;  // For MRBEAST

			I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
			I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
			I08 hasTrendCmpnt   = 1;
			//Xnewterm is used at this point and should not be touched!
			F32PTR BUF      = yInfo.Y;
			f32_fill_val(0., BUF, Nq);
	
			if (hasTrendCmpnt)   f32_add_vec_inplace(result.tY, BUF, Nq);
			if (hasSeasonCmpnt)  f32_add_vec_inplace(result.sY, BUF, Nq);
			if (hasOutlierCmpnt) f32_add_vec_inplace(result.oY, BUF, Nq);
		
			for (int j = 0; j < q; ++j) {
				//For MRBEAST
				F32 r = f32_corr_rmse_nan(BUF+N*j, Xnewterm + N*j, N, &result.RMSE[j]);
				result.R2[j] = r * r;
			}
			*/
		}
 
	    if (skipCurrentPixel) Init_Result(&result, opt, getNaN());
		/*********************************************/
	 	 
		BVS_WriteOutput(opt, &result, &xInfo);
		 
		/*
		BVS_TERMS out[200];
		coder.curBetaID = 0;
		coder.curModelID = 0;
		coder.curTermID = 0;
		coder.Kpre      = 0;
		for (int i = 0; i < coder.nModels; i++) {
			ExtractNextTermFromList(&coder, out);	 
		}
	   */
		opt->prior.alpha1 = yInfo.mean[0];
		opt->prior.alpha2 = yInfo.sd[0];

 		SaveOutput(&xInfo, &coder, opt);

	} //for (U32 pixelIndex = 1; pixelIndex <= TOTALNUMPIXELS; pixelIndex++)


 

	/***********************************************************/
	// This is the ending bracekt of the iteration through pixels
	/***********************************************************/

	r_vslDeleteStream(&stream);
	MEM.free_all(&MEM);
 
	return 1;
} /* End of beastST() */


#include "abc_000_warning.h"