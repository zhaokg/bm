#include "abc_000_warning.h"
#include "abc_001_config.h"



#include "abc_ide_util.h"
#include "abc_vec.h"
#include "abc_mem.h"
#include "bvs_io.h"

//#define  _(x)						RemoveField(fieldList, nfields, #x), mat->s##x = NULL

void* Alloc_OutPut(A(OPTIONS_PTR)  opt) {	
	
	const BVS_IO_PTR  io = &opt->io;
	
	// Moved from glue_code to accomodate multivriate time series
	// Added to accomodate the beast_thread() for Win GUI where Output_AllocMem is called
	// each time a thread is created.
	free0(io->out.result);
	
	I32 q = 1;
	const BVS_RESULT_PTR  mat =io->out.result = malloc0(sizeof(BVS_RESULT) * q);
 
#if   R_INTERFACE==1
	io->out.dtype = DATA_DOUBLE;
#elif M_INTERFACE==1 || P_INTERFACE==1
	io->out.dtype = DATA_FLOAT;
#endif

	DATA_TYPE    dtype  = io->out.dtype; // DATA_FLOAT or DATA_DOUBLE                            
 
	// stackoverflow.com/questions/2124339/c-preprocessor-va-args-number-of-arguments
	#define NUMARGS(...)                (sizeof((int[]){__VA_ARGS__})/sizeof(int))
	#define NARGS(...)                  (sizeof((int[]){0, ##__VA_ARGS__})/sizeof(int)-1)
	#define _(name, ...)                {#name, dtype, NUMARGS(__VA_ARGS__),{__VA_ARGS__},(void ** )&mat->tY##name }
	#define _1(name, ...)               _(name, __VA_ARGS__)  

	 
	int  Kmax         = opt->prior.Kmax;
	I32  P            = opt->prior.Pall;
	I32  Plin         = opt->prior.Pbasis[LinearTerm];
	I32  Pstep        = opt->prior.Pbasis[StepTerm];
	I32  Pstair       = opt->prior.Pbasis[StairTerm];
	I32  Ppiece       = opt->prior.Pbasis[PieceTerm];
	I32  Pseason      = opt->prior.Pbasis[SeasonTerm];
	I32  Phinge       = opt->prior.Pbasis[HingeTerm];
	I32  Phingepair   = opt->prior.Pbasis[HingePairTerm];
	I32  Pchange      = opt->prior.Pbasis[ChangeTerm];
	

	//#define __DEBUG__  
	
	 int   N = io->meta.Nobj;
	int    isMultiVariate = (q == 1) ? 0L : 1L;
	int ROW = 1;
	int COL = 1;
	I32PTR nrows, ncols;
	FIELD_ITEM  fieldList[ ] ={
	 	{"marg_lik",  dtype,		2,   {ROW,COL,},  &mat->marg_lik},
		#ifdef __DEBUG__
        {"R2",      dtype,		2,   {(N + 7) / 8 * 8,300,},  &mat->R2},
		{"RMSE",    dtype,		2,   {(N + 7) / 8 * 8,300,},  &mat->RMSE},
	    #else
		{"R2",        dtype,		2,   {ROW,COL,},  &mat->R2,  .extra = isMultiVariate},
		{"RMSE",      dtype,		2,   {ROW,COL,},  &mat->RMSE,.extra = isMultiVariate},
		#endif
	    {"sig2",      dtype,		2,   {q,q,},       &mat->sig2}, //Needed to be changed to reflect the output dim
		// {"sig2",      dtype,		4,   {ROW,COL,q,q},  &mat->sig2},
		//{"trend",     DATA_STRUCT,	0,	 {0,},       (void**)trend},
		//{"season",    DATA_STRUCT,	0,	 {0,},       (void**)season},
		{"Y",            dtype,		2,   {N,1,},           &mat->Y}, //Needed to be changed to reflect the output dim
		{"SD",           dtype,		2,   {N,1,},           &mat->SD}, //Needed to be changed to reflect the output dim
		{"KProb",        dtype,		2,   {Kmax,1,},        &mat->KProb}, //Needed to be changed to reflect the output dim
		{"Kmean",        dtype,		2,   {1,1,},           &mat->Kmean}, //Needed to be changed to reflect the output dim
		{"Ylin",         dtype,		2,   {Plin,1,},        &mat->Ylin}, //Needed to be changed to reflect the output dim
		{"Ystep",      dtype,		2,   {N,Pstep,},      &mat->Ystep}, //Needed to be changed to reflect the output dim
		{"Ystair",       dtype,		2,   {N,Pstair,},      &mat->Ystair}, //Needed to be changed to reflect the output dim
		{"Ypiece",       dtype,		2,   {N,Ppiece,},      &mat->Ypiece}, //Needed to be changed to reflect the output dim
		{"Yseason",      dtype,		2,   {N,Pseason,},      &mat->Yseason}, //Needed to be changed to reflect the output dim
		{"Yhinge",      dtype,		2,   {N,Phinge,},& mat->Yhinge }, //Needed to be changed to reflect the output dim
		{"YhingePair",      dtype,		2,   {N,Phingepair,},&mat->YhingePair }, //Needed to be changed to reflect the output dim
		{"Ychange",      dtype,		2,   {N,Pchange,},&mat->Ychange }, //Needed to be changed to reflect the output dim
		{"probAll",       dtype,	2,   {P,1,},	   &mat->probAll}, //Needed to be changed to reflect the output dim
		{"probLin",       dtype,	2,   {Plin,1,},	   &mat->probLin}, //Needed to be changed to reflect the output dim
		{"probStep",     dtype,	    2,   {N,Pstep,  },     &mat->probStep}, //Needed to be changed to reflect the output dim
		{"probStair",    dtype,	    2,   {N,Pstair,  },    &mat->probStair}, //Needed to be changed to reflect the output dim
		{"probPiece",    dtype,	    2,   {N,Ppiece,  },    &mat->probPiece}, //Needed to be changed to reflect the output dim
 	    {"probSeason",    dtype,	 2,   {N,Pseason,  },   &mat->probSeason}, //Needed to be changed to reflect the output dim
		{"probHinge",    dtype,	    2,   {N,Phinge, },     &mat->probHinge}, //Needed to be changed to reflect the output dim
		{"probHingePair",  dtype,	2,   {N,Phingepair, }, &mat->probHingePair}, //Needed to be changed to reflect the output dim
		{"probChange", dtype,	    2,   {N,Pchange,},     &mat->probChange}, //Needed to be changed to reflect the output dim
		{"xinfo",        DATA_INT32,	 2,   {1,1,},    NULL}, //Needed to be changed to reflect the output dim
		{"terms",        DATA_INT32,	 2,   {1,1,},    NULL}, //Needed to be changed to reflect the output dim
		{"opt",          DATA_INT32,	 2,   {1,1,},    NULL}, //Needed to be changed to reflect the output dim
		{"nrows",        DATA_INT32,  2,  {1,1,},   &nrows },
		{"ncols",        DATA_INT32,  2,  {1,1,},   &ncols },
		{0,},
	};

	I32    nfields = 99999;
		  
	if (opt->extra.removeSingletonDims) {
		RemoveSingltonDims(fieldList, nfields);
	}


	if (!hasTerm(LinearTerm)) {
		RemoveField(fieldList, nfields, "Ylin");
		RemoveField(fieldList, nfields, "probLin");
	}
	if (!hasTerm(StepTerm)) {
		RemoveField(fieldList, nfields, "Ystep");
		RemoveField(fieldList, nfields, "probStep");
	}
	if (!hasTerm(StairTerm)) {
		RemoveField(fieldList, nfields, "Ystair");
		RemoveField(fieldList, nfields, "probStair");
	}
	if (!hasTerm(PieceTerm)) {
		RemoveField(fieldList, nfields, "Ypiece");
		RemoveField(fieldList, nfields, "probPiece");
	}
	if (!hasTerm(HingeTerm)) {
		RemoveField(fieldList, nfields, "Yhinge");
		RemoveField(fieldList, nfields, "probHinge");
	}
	if (!hasTerm(HingePairTerm)) {
		RemoveField(fieldList, nfields, "YhingePair");
		RemoveField(fieldList, nfields, "probHingePair");
	}
	if (!hasTerm(ChangeTerm)) {
		RemoveField(fieldList, nfields, "Ychange");
		RemoveField(fieldList, nfields, "probChange");
	}

 	I32       nprt = 0;
	VOID_PTR  out  = PROTECT(CreateStructVar(fieldList, nfields));                
	nprt++;

	nrows[0] = ROW;
	ncols[0] = COL; 
	AddStringAttribute(out,  "class",    "bvs");	 
	AddIntegerAttribute(out, "intattr",  9999);	  
	 
	UNPROTECT(nprt);
	return out;

	#undef NUMARGS
    #undef NARGS
    #undef _
	#undef _2
	#undef _3
	#undef _4   
}

void  Alloc_Result(A(RESULT_PTR)  result, A(OPTIONS_PTR)  opt, MemPointers* _restrict MEM) {		 

	const I32 N    = opt->io.meta.Nobj;
	const I32 q    = 1;
	const I32 Nq   = N * q;

	int  Kmax         = opt->prior.Kmax;
	I32  P            = opt->prior.Pall;
	I32  Plin         = opt->prior.Pbasis[LinearTerm];
	I32  Pstep        = opt->prior.Pbasis[StepTerm];
	I32  Pstair       = opt->prior.Pbasis[StairTerm];
	I32  Ppiece       = opt->prior.Pbasis[PieceTerm];
	I32  Pseason      = opt->prior.Pbasis[SeasonTerm];
	I32  Phinge       = opt->prior.Pbasis[HingeTerm];
	I32  Phingepair   = opt->prior.Pbasis[HingePairTerm];
	I32  Pchange      = opt->prior.Pbasis[ChangeTerm];
	
	
	*result = (BVS_RESULT){0, }; //memset result to zero


	result->marg_lik = MEM->alloc(MEM, sizeof(F32) * 1,   0);
	result->sig2	 = MEM->alloc(MEM, sizeof(F32) * q*q, 0);
	result->R2		 = MEM->alloc(MEM, sizeof(F32) * q,   0);
	result->RMSE	 = MEM->alloc(MEM, sizeof(F32) * q,   0);

 
	result->Y       = MEM->alloc(MEM, sizeof(F32) * N, 0);
	result->SD      = MEM->alloc(MEM, sizeof(F32) * N, 0);
	result->KProb   = MEM->alloc(MEM, sizeof(F32) * Kmax, 0);
	result->Kmean   = MEM->alloc(MEM, sizeof(F32) * 1, 0);
 

	result->probAll = MEM->alloc(MEM, sizeof(F32) * P, 0);
	if (hasTerm(LinearTerm)) 	     result->probLin   = MEM->alloc(MEM, sizeof(F32) * 1 * Plin, 0);
	if (hasTerm(StepTerm)   ) 	 result->probStep      = MEM->alloc(MEM, sizeof(F32) * N*Pstep, 0);   
	if (hasTerm(StairTerm)) 	 result->probStair     = MEM->alloc(MEM, sizeof(F32) * N * Pstair, 0);
	if (hasTerm(PieceTerm)) 	 result->probPiece     = MEM->alloc(MEM, sizeof(F32) * N * Ppiece, 0);
	if (hasTerm(SeasonTerm)) 	 result->probSeason    = MEM->alloc(MEM, sizeof(F32) * N * Pseason, 0);
	if (hasTerm(HingeTerm) ) 	 result->probHinge     = MEM->alloc(MEM, sizeof(F32) * N * Phinge, 0); 
	if (hasTerm(HingePairTerm))  result->probHingePair = MEM->alloc(MEM, sizeof(F32) * N * Phingepair, 0);
	if (hasTerm(ChangeTerm)) 	 result->probChange    = MEM->alloc(MEM, sizeof(F32) * N * Pchange, 0);

	if (hasTerm(LinearTerm)) 	  result->Ylin      = MEM->alloc(MEM, sizeof(F32) * 1 * Plin, 0);
	if (hasTerm(StepTerm)   ) 	 result->Ystep      = MEM->alloc(MEM, sizeof(F32) * N* Pstep, 0);   
	if (hasTerm(StairTerm)) 	 result->Ystair     = MEM->alloc(MEM, sizeof(F32) * N * Pstair, 0);
	if (hasTerm(PieceTerm)) 	 result->Ypiece     = MEM->alloc(MEM, sizeof(F32) * N * Ppiece, 0);
	if (hasTerm(SeasonTerm)) 	 result->Yseason    = MEM->alloc(MEM, sizeof(F32) * N * Pseason, 0);
	if (hasTerm(HingeTerm) ) 	 result->Yhinge     = MEM->alloc(MEM, sizeof(F32) * N * Phinge, 0); 
	if (hasTerm(HingePairTerm))  result->YhingePair = MEM->alloc(MEM, sizeof(F32) * N * Phingepair, 0);
	if (hasTerm(ChangeTerm)) 	 result->Ychange    = MEM->alloc(MEM, sizeof(F32) * N * Pchange, 0);


}

void  Init_Result(A(RESULT_PTR)  result, A(OPTIONS_PTR)  opt, const F32 nan)
{

	const I32 N = opt->io.meta.Nobj;
	const I32 q = 1;
	const I32 Nq = N * q;
	const I32 Kmax = opt->prior.Kmax;
 
	I32  P            = opt->prior.Pall;
	I32  Plin         = opt->prior.Pbasis[LinearTerm];
	I32  Pstep        = opt->prior.Pbasis[StepTerm];
	I32  Pstair       = opt->prior.Pbasis[StairTerm];
	I32  Ppiece       = opt->prior.Pbasis[PieceTerm];
	I32  Pseason      = opt->prior.Pbasis[SeasonTerm];
	I32  Phinge       = opt->prior.Pbasis[HingeTerm];
	I32  Phingepair   = opt->prior.Pbasis[HingePairTerm];
	I32  Pchange      = opt->prior.Pbasis[ChangeTerm];
	

	*result->marg_lik = nan;
	f32_fill_val(nan, result->sig2, q * q);
	f32_fill_val(nan, result->R2, q);
	f32_fill_val(nan, result->RMSE, q);

	f32_fill_val(nan, result->Y, N);
	f32_fill_val(nan, result->SD, N);
	f32_fill_val(nan, result->KProb, Kmax);

	*result->Kmean = nan;

	f32_fill_val(nan, result->probAll, P);

	if (hasTerm(LinearTerm)) 	  f32_fill_val(nan, result->probLin,  1 * Plin);
	if (hasTerm(StepTerm)  ) 	  f32_fill_val(nan, result->probStep, N*Pstep);
	if (hasTerm(StairTerm)) 	  f32_fill_val(nan, result->probStair, N * Pstair);
	if (hasTerm(PieceTerm)) 	  f32_fill_val(nan, result->probPiece, N * Ppiece);
	if (hasTerm(SeasonTerm)) 	  f32_fill_val(nan, result->probSeason, N * Pseason);
	if (hasTerm(HingeTerm) ) 	  f32_fill_val(nan, result->probHinge, N * Phinge);
	if (hasTerm(HingePairTerm))   f32_fill_val(nan, result->probHingePair, N * Phingepair);
	if (hasTerm(ChangeTerm)) 	  f32_fill_val(nan, result->probChange, N * Pchange);

	if (hasTerm(LinearTerm)) 	  f32_fill_val(nan, result->Ylin,  1 * Plin);
	if (hasTerm(StepTerm)  ) 	  f32_fill_val(nan, result->Ystep, N*Pstep);
	if (hasTerm(StairTerm)) 	  f32_fill_val(nan, result->Ystair, N * Pstair);
	if (hasTerm(PieceTerm)) 	  f32_fill_val(nan, result->Ypiece, N * Ppiece);
	if (hasTerm(SeasonTerm)) 	  f32_fill_val(nan, result->Yseason, N * Pseason);
	if (hasTerm(HingeTerm) ) 	  f32_fill_val(nan, result->Yhinge, N * Phinge);
	if (hasTerm(HingePairTerm))   f32_fill_val(nan, result->YhingePair, N * Phingepair);
	if (hasTerm(ChangeTerm)) 	  f32_fill_val(nan, result->Ychange, N * Pchange);
 
}


 

#include "abc_000_warning.h"