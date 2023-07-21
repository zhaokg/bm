#include "abc_000_warning.h"

#include "abc_001_config.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
 
#ifndef ARM64_OS
	#include <immintrin.h> //https://stackoverflow.com/questions/56049110/including-the-correct-intrinsic-header
#endif

#include "abc_datatype.h"
#include "abc_blas_lapack_lib.h"
#include "abc_common.h"
#include "abc_ide_util.h"
#include "abc_pthread.h"
#include "abc_timer.h"
#include "abc_vec.h"
#include "abc_cpu.h"
#include "abc_timer.h"
#include "abc_rand.h"

//#include "abc_dir.h"
 
#include "globalvars.h"

#if defined(WIN64_OS) 
/*---------WINDOW-------------------------*/
//extern void DllExport WinMainDemoST(BEAST_OPTIONS_PTR  option);
//extern void DllExport WinMainDemoTrend(BEAST_OPTIONS_PTR  option);
/*---------WINDOW-------------------------*/
#endif


//#include "..\TinyTiff\tinytiffreader.h"

#ifndef ARM64_OS
	#include "abc_math_avx.h"
#endif


#include "bvs_io.h"
#include "bvs_fun.h"

#define IS_STRING_EQUAL(a, b)  (strcicmp(a, b) == 0)



void * mainFunction(void *prhs[], int nrhs) {

	//r_printf("\033[1m\033[32m" "\033[4m" "%cERROR: Essential input paramaters are missing!\n",149);
	if (nrhs == 0) {
		r_error("ERROR: Essential input paramaters are missing!\n");
		return IDE_NULL;
	}
	if (!IsChar(prhs[0])) {
		r_error("ERROR: The very first parameter must be a string specifying the algorithm name!\n");
		return IDE_NULL;
	}
	 
	if (nrhs >= 7) {
		int avxOption = GetScalar(prhs[nrhs - 1]);
		SetupRoutines_UserChoice(avxOption);
	}
	else {
		int printlevel = 0;
		SetupRoutines_ByCPU(printlevel); // prnitlevel=1: only print hhe libray choice result
	}

	//SetupVectorFunction_Generic();
	//SetupPCG_GENERIC(); 	 	
	//print_funcs();

	/*
	int i, rc;
	long t1 = 1, t2 = 2, t3 = 3;
	pthread_t      threads[3];
	pthread_attr_t attr;


	pthread_mutex_init(&count_mutex, NULL); //Initialize mutex and condition variable objects
	pthread_cond_init(&count_threshold_cv, NULL);

	// For portability, explicitly create threads in a joinable state 
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	pthread_create(&threads[0], &attr, watch_count, (void *)t1);
	pthread_create(&threads[1], &attr, inc_count, (void *)t2);
	pthread_create(&threads[2], &attr, inc_count, (void *)t3);

	//Wait for all threads to complete
	for (i = 0; i<NUM_THREADS; i++) {
	   pthread_join(threads[i], NULL);
	}
	printf("Main(): Waited on %d  threads. Done.\n", NUM_THREADS);


	//Clean up and exit
	pthread_join(threads[i], NULL);
	pthread_attr_destroy(&attr);
	pthread_mutex_destroy(&count_mutex);
	pthread_cond_destroy(&count_threshold_cv);
	pthread_exit(NULL);
	return;
	*/
 
 
	#define STRING_LEN 20

	char  algorithm[STRING_LEN +1];
	GetCharArray(prhs[0], algorithm, STRING_LEN);

	// Now 'algorithm' should be eitehr the default algorith or the value from Opt.

	void * ANS  = NULL;
	int    nptr = 0;
	if (IS_STRING_EQUAL(algorithm, "bvstrain"))  {
		/*
		// Initialize mutex and condition variable objects
		pthread_mutex_init(&MUTEX_WRITE, NULL);
		pthread_cond_init(&CONDITION_WRITE, NULL);
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		DATA_AVAILABLE_WRITE = false;
		pthread_create(&THREADID_WRITE, &attr, WRITE, (VOID_PTR )&threadParWrite);
		pthread_attr_destroy(&attr); 	
		*/

	    //stackoverflow.com/questions/10828294/c-and-c-partial-initialization-of-automatic-structure
	    //stackoverflow.com/questions/11152160/initializing-a-struct-to-0
				
		//Warning from MacOS: suggest braces around initialization of subobject [-Wmissing-braces]
		BVS_OPTIONS      option = { {{0,},}, }; 
	
		if (BVS_GetArgs(prhs, nrhs, &option) == 0) {	
			BVS_DeallocatePriorPts(&option.prior);
			return IDE_NULL;
		}		
		
		option.io.out.result		 = NULL;		 	// result to be allocated in OUput_allocMEM
		int q = 1;
		if (q == 1) {
			ANS = PROTECT(Alloc_OutPut(&option)); nptr++;	
			option.ans = ANS; // used in SaveOutput
		}  
		/**********************************/
	
		GLOBAL_OPTIONS = (BVS_OPTIONS_PTR)&option;
		beast2_main_corev4();		
		r_printf("\n");
	 
		BVS_DeallocatePriorPts(&option.prior);
 
	} else if (IS_STRING_EQUAL(algorithm, "bvspredict")) {
		/*
		// Initialize mutex and condition variable objects
		pthread_mutex_init(&MUTEX_WRITE, NULL);
		pthread_cond_init(&CONDITION_WRITE, NULL);
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		DATA_AVAILABLE_WRITE = false;
		pthread_create(&THREADID_WRITE, &attr, WRITE, (VOID_PTR )&threadParWrite);
		pthread_attr_destroy(&attr);
		*/

		/**********************************/
		//stackoverflow.com/questions/10828294/c-and-c-partial-initialization-of-automatic-structure
		//stackoverflow.com/questions/11152160/initializing-a-struct-to-0
		//BVS_OPTIONS      option = {0,};		
		BVS_PRED_XINFO    newxinfo = { 0, };		
		BVS_XINFO         xinfo    = { 0, }; // may contain thread-specific pointers
		BVS_TERMS_ENCODER coder    = { 0, };
		BVS_PRIOR_PTR     prior;
		BVS_PRED_OPT      newopt;

		GetArg_Predict(prhs, nrhs, &xinfo, &coder, &newxinfo, &prior , &newopt);
 
		BVS_PRED_RESULT mat;
		newxinfo.mat = &mat;
		int q = 1;
		if (q == 1) {
			//ANS = PROTECT(Alloc_OutPut(&option)); nptr++;
			//option.ans = ANS; // used in SaveOutput	 
  
			ANS = PROTECT( Predict_Alloc_Output(&xinfo, &newxinfo, &mat, prior, coder.nModels, &newopt) );
			nptr++;
			//AddIntegerAttribute(out, "hasOutlier", hasOutlierCmpnt);	  
		}
		/**********************************/

		//GLOBAL_OPTIONS = (BVS_OPTIONS_PTR)&option;		
		 bvs_predict_corev4(&xinfo, &coder,  &newxinfo, prior, &newopt);
		 
		free0(newxinfo.Xnorm1);
		free0(newxinfo.X2);
		free0(newxinfo.Xones);
		free0(newxinfo.Xorg);
	}

 
 
	
	/*
		// Initialize mutex and condition variable objects
		pthread_mutex_init(&MUTEX_WRITE, NULL);
		pthread_cond_init(&CONDITION_WRITE, NULL);
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		DATA_AVAILABLE_WRITE = false;
		pthread_create(&THREADID_WRITE, &attr, WRITE, (VOID_PTR )&threadParWrite);
		pthread_attr_destroy(&attr); 	
	*/
	/*
	else if (IS_STRING_EQUAL(algorithm, "beastMT"))
	{
	

		BEAST_OPTIONS         option;
		BEAST_IO              io;
		BEAST_RESULT          result;
		memset(&io,     0, sizeof(BEAST_IO));
		memset(&result, 0, sizeof(BEAST_RESULT));
		option.io = &io;
		option.io->out.result  = &result;
		option.io->isRegularOrdered = 1;


		if (!BEAST_Get1stArg_Data(prhs, nrhs, &option) ||
			!BEAST_Get2ndArg_MetaData(prhs, nrhs, &option) ||
			!BEAST_Get3rdArg_Prior(prhs, nrhs, &option) ||
			!BEAST_Get4thArg_MCMC(prhs, nrhs, &option) ||
			!BEAST_Get5thArg_FLAGS(prhs, nrhs, &option)) {
			return NULL;
		}

		BEAST_print_options(&option);
		ANS = PROTECT(BEAST_AllocateOutput(&option)); nptr++;
		GLOBAL_OPTIONS = (BEAST_OPTIONS_PTR)&option;
	 
		
 
		pthread_mutex_init(&mutex, NULL); //Initialize mutex and condition variable objects
		pthread_cond_init(&condVar, NULL);

		// For portability, explicitly create threads in a joinable state 
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
 
		cpu_set_t cpus;

		
		I32 NUM_THREADS;

		I32 nCores   = GetNumCores();
		nCores       = max(nCores, 1);

		NUM_THREADS  = option.extra.numCPUCoresToUse;
		if (NUM_THREADS == 0) {
			NUM_THREADS = nCores - 1;
		} else if (NUM_THREADS < 0)	{
			NUM_THREADS = nCores + NUM_THREADS;
		}
		NUM_THREADS = min(nCores- 1, NUM_THREADS);
		NUM_THREADS = max(NUM_THREADS, 1);
		NUM_THREADS = min(NUM_THREADS, option.io->numOfPixels);	

		thread_id    = malloc(sizeof(pthread_t) * NUM_THREADS);

		NEXT_PIXEL_INDEX = 1;
		for (I32 i = 0; i < NUM_THREADS; i++) 	{	
			 CPU_ZERO(&cpus);
			 CPU_SET(i, &cpus);
			 pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);

			 extern 	int beast_main_mthrd(void *);
			 pthread_create(&thread_id[i], &attr, beast_main_mthrd, (void *)NULL);
		}
		pthread_attr_destroy(&attr);

		if (option.extra.printProgressBar) {
			// Print a blank line to be backspaced by the follow
			void* BUF;
			BUF = malloc(option.extra.consoleWidth * 3);
			printProgress2(0, 0, option.extra.consoleWidth, BUF, 1);

			PERCENT_COMPLETED = 0;
			REMAINING_TIME    = 10000;
			// https://www.mathworks.com/matlabcentral/answers/101658-is-it-possible-to-start-new-threads-from-a-c-mex-file
			// https://stackoverflow.com/questions/28227313/multithreaded-pthreads-matlab-mex-function-causes-matlab-to-crash-after-exitin
			// https://stackoverflow.com/questions/54010898/how-do-you-print-to-console-in-a-multi-threaded-mex-function
			while (PERCENT_COMPLETED < 1.f && NEXT_PIXEL_INDEX < option.io->numOfPixels) {			
				printProgress2(PERCENT_COMPLETED, REMAINING_TIME, option.extra.consoleWidth, BUF, 0);
				Sleep_ms(2 * 1000);
			}
			printProgress2(1.0, 0, option.extra.consoleWidth, BUF, 0);
			free(BUF);
		} else	{
			r_printf("\nRbeast: Waiting on %d threads...\n", NUM_THREADS);
		}


		// Wait for all threads to complete
		for (I32 i = 0; i<NUM_THREADS; i++) {
			pthread_join(thread_id[i], NULL);
		}
		r_printf("\nRbeast: Waited on %d threads. Done.\n", NUM_THREADS);

		//Clean up and exit		
		pthread_mutex_destroy(&mutex);
		pthread_cond_destroy(&condVar);
		//pthread_exit(NULL);	
		free(thread_id);
	}
	 */

 
	UNPROTECT(nptr);

	return ANS == NULL ? IDE_NULL : ANS;	
}
 

#if R_INTERFACE==1
#include <R_ext/libextern.h>
#include "Rembedded.h"

//	R_FlushConsole(): a R functin to flush the print
#if defined(MSVC_COMPILER)
SEXP DllExport rexFunction1(SEXP rList, SEXP dummy)
#else
SEXP DllExport rexFunction(SEXP rList, SEXP dummy)
#endif
{
  	
 /*
	FILELIST_PTR flist = GetFlist("a:\\", "tif");
	PrintFlist(flist);
	FreeFlist(flist);
	return R_NilValue;
	
	SEXP IMG = R_NilValue;
 
	struct TinyTIFFReaderFile* tiffr = TIFF_open("a:\\landsat_postfire_2012.tif");
	//tiffr = TIFF_open("G:\\land\\ndvi\\LT05_018033_19840327.2001-07-14.tif");
	
	if (!tiffr) {
		printf("XXX\n");
		return IMG;
	}
	
	TIFF_readIFD(tiffr, 0);
	TIFF_getCurrentFrameInfo(tiffr);

	uint32_t width  = TIFF_getWidth(tiffr);
	uint32_t height = TIFF_getHeight(tiffr);
	double* image   = (uint16_t*)malloc(width*height*sizeof(double));

	int sample = GetScalar(VECTOR_ELT(rList, 0));

	TIFF_getSampleData(tiffr, image, sample);
	TIFF_close(tiffr);
	
	///////////////////////////////////////////////////////////////////
	// HERE WE CAN DO SOMETHING WITH THE IMAGE IN image (ROW-MAJOR!)
	///////////////////////////////////////////////////////////////////
	IMG = PROTECT(Rf_allocMatrix(REALSXP, width, height));
	//IMG = PROTECT(Rf_alloc3DArray(REALSXP, width, height, 6));
	double *tmp=REAL(IMG);
	for (I32 i = 0; i < height*width; i++)
		tmp[i] = *((double*)image + i);
	free(image);
	UNPROTECT(1);
	return IMG; 
 */

	if (!isNewList(rList)) 	return R_NilValue;
 
	// stat.ethz.ch/pipermail/r-help/2001-September/015081.html
	//SEXP   prhs = VECTOR_DATA(inList);
	SEXP   prhs[10];
	int    nrhs = length(rList);
	nrhs = min(10L, nrhs);
	for (int i = 0; i < nrhs; i++) 	
		prhs[i] = VECTOR_ELT(rList, i);

	SEXP ans;
	PROTECT(ans=mainFunction(prhs, nrhs));
	UNPROTECT(1);
	return ans != NULL ? ans : R_NilValue;
}


#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}


#if (defined(WIN64_OS) || defined(WIN32_OS)) 
	SEXP TetrisSetTimer(SEXP action, SEXP seconds, SEXP envior);
	static const R_CallMethodDef CallEntries[] = {
		#if defined(MSVC_COMPILER)
		CALLDEF(rexFunction1,    2),
		#else
		CALLDEF(rexFunction,    2),
		#endif
		//CALLDEF(TetrisSetTimer, 3),
		{ NULL, NULL, 0 }
	};
#else
static const R_CallMethodDef CallEntries[] = {
				CALLDEF(rexFunction, 2),
				{ NULL, NULL, 0 }
			};
#endif



void  R_init_Rbeast(DllInfo *dll) {
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	//R_forceSymbols(dll, TRUE);
}

#elif M_INTERFACE ==1

/***************************************************************
* function: mexFunction() - Matlab interface function.         *
* INPUTS:                                                      *
*   nlhs - The number of output variables to be assigned by    *
*     the mex function.                                        *
*   plhs[] - Empty array of MATLAB matricies, of size nlhs.    *
*     This is to be filled in by this application.             *
*   nrhs - Number of input arguments to this mex funciton.     *
*   prhs[] - Array of MATLAB matricies from which input data   *
*     is taken.                                                *
***************************************************************/

#include "abc_date.h"
//void DllExport mexFunction(int nlhs, mxArray * _restrict plhs[],   int nrhs, const mxArray * _restrict prhs[]) {
/* 
   Restirct is really pf no use here. More important, the decleration of mexFunction in mex.h has no restrict keyword;
   MVSC is OK with the difference; gcc, however, failed with a complaint of conflicting types for ‘mexFunction’. 
   What flag can circuvume this? The best thing I could find is
   https://stackoverflow.com/questions/66951324/conflicting-types-compiling-a-ld-preload-wrapper
*/
void DllExport mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
		
	//https://stackoverflow.com/questions/19813718/mex-files-how-to-return-an-already-allocated-matlab-array
	//https://www.mathworks.com/matlabcentral/answers/77048-return-large-unchange-mxarray-from-mex
	//https://stackoverflow.com/questions/18847833/is-it-possible-return-cell-array-that-contains-one-instance-in-several-cells/18849127#18849127
	//https://undocumentedmatlab.com/articles/matlabs-internal-memory-representation/
	// NOT ALWAYS WORK
	// plhs[0] = prhs[0];
	// mxCreateSharedDataCopy
	
	mxArray * ans = mainFunction(prhs, nrhs);
	plhs[0] = ans;
	return;
}
#endif

 

#include "abc_000_warning.h"