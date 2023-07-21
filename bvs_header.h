#pragma once

#include "abc_000_macro.h"
#include "abc_datatype.h"   //#include <inttypes.h>  #include <stdint.h>
#include "abc_ide_util.h"
#include "abc_mat.h"

#ifndef SOLARIS_COMPILER

#if R_INTERFACE==1 
	#ifdef beta
		#undef beta  //beta is refined to Rf_beta in Rmath.h, which causes errors in expanding model->beta
	#endif
#endif

#define _IN_ 
#define _OUT_

#define MAX_RAND_NUM	  5000L
#define MIN_PREC_VALUE    0.00001
#define MIN_SIG2_VALUE    0.00001
#define MIN_ALPHA2_VALUE  0.000001

#define PROB_SAMPLE_EXTREME_VECTOR 0.5

#define SEASONID   0
#define TRENDID    1
#define OUTLIERID  2  
#define DUMMYID    3
#define SVDID      4
#define A(xxx)     BVS_##xxx

typedef U32 TKNOT,  *_restrict TKNOT_PTR;
typedef U08 TORDER, *_restrict TORDER_PTR;
#define rTKNOT_PTR   register  TKNOT_PTR
#define rTORDER_PTR  register  TORDER_PTR
#define rTKNOT       register  TKNOT
#define rTORDER      register  TORDER

/*************************************************/
// Declarations of data types for input options
/*************************************************/

typedef enum { 
	ConstTerm=-1,
	LinearTerm,QuadraticTerm, StepTerm, StairTerm, PieceTerm, SeasonTerm, HingeTerm, HingePairTerm,ChangeTerm,
	Hinge2DTerm,Hinge2DRTerm, Hinge1DRTerm,
	LAST_BASIS_VALUE
} TEMR_TYPE;

#define MAX_NUM_BASIS  (Hinge1DRTerm+1)
#define MAX_NUM_BASIS  LAST_BASIS_VALUE

// 0: Use a constant precison value
// 1: Use a uniform precision value for all terms
typedef enum { ConstPrec, UniformPrec, ComponentWise, OrderWise} PRECPRIOR_TYPE;

typedef struct BVS_METADATA {

	VOID_PTR  Xobj, Yobj;	
	I32       Nobj, Pobj;
	DATA_TYPE XobjDtype, YobjDtype;
	I08		  nrhs;
	I08       ndimx;

	////////////////////////////////
	I08	   isMetaStruct; 
	F32    missingValue; 
	F32    maxMissingRate;

} BVS_METADATA, * _restrict BVS_METADATA_PTR;

struct BVS_RESULT;
typedef struct BVS_RESULT BVS_RESULT, * _restrict BVS_RESULT_PTR;

typedef struct BVS_IO {
  	BVS_METADATA	meta;  
	I08				ndim;
	I32				dims[3];
	// q is the number of time series;q=1 is for univaraite TS in BEASTv4.
	// q is added to ensure a consistent API inteface btw BEAST and MRBEAST
	// For irregular time seires, Nraw is obtained from ndims[whichDimIsTime=1]
	struct {
		BVS_RESULT* result;
		DATA_TYPE   dtype;		
	} out;
} BVS_IO, * _restrict BVS_IO_PTR;

typedef struct BVS_PRIOR {
	I16	  basisType[MAX_NUM_BASIS];			
	I16   type2id[MAX_NUM_BASIS];
	I16	  numBasis;

	I32    Pobj;
	I32    Pall;
	I32    Pbasis[MAX_NUM_BASIS];
	I32PTR Pindices0[MAX_NUM_BASIS];
	I32PTR PchangeY0;
	I32PTR PseasonY0;
	I32PTR FixedLinearTerms0;
	I32    nFixedLinearTerms;

	I32    nGrpHinge2d, nGrpHinge2dr, nGrpHinge1dr;
	I32PTR sGrpHinge2d, sGrpHinge2dr, sGrpHinge1dr;
	I32PTR eGrpHinge2d, eGrpHinge2dr, eGrpHinge1dr;

	I08   isBeast;
	I08   isFFT;
	I32   beast_Pbase;
	I32   beast_Pextra;
	I32   beast_Ptime;
	F32   beast_period;
	F32   beast_maxSeasonOrder;

	I32   fft_Pbase;
	I32   fft_Pnew;
	I32   fft_Ptime;
	F32   fft_period;
	F32   fft_maxOrder;

	I32      maxKnotNum_Max[MAX_NUM_BASIS];
 
	I16PTR   maxKnotNum_Vec[MAX_NUM_BASIS];
	I16PTR   minSepDist_Vec[MAX_NUM_BASIS];

  
	int  num_con1d_Decrease, num_con1d_Increase, num_con1d_nShape;
	int  num_con1d_uShape, num_con1d_NoWiggle;
	I32PTR con1d_Decrease;
	I32PTR con1d_Increase;
	I32PTR con1d_nShape;
	I32PTR con1d_uShape;
	I32PTR con1d_NoWiggle;

	U08	  modelPriorType;
	U08   precPriorType;

	U16   Kmax; 

	F32   sigFactor;
	F32   sig2;
	F32	  precValue;
	F32   alpha1, alpha2, delta1, delta2;
	

} BVS_PRIOR, * _restrict BVS_PRIOR_PTR;

typedef struct BVS_MCMC {
	U64   seed;	                  // Unsigned long long seed;
	F32   credIntervalAlphaLevel;
	F32   ridgeFactor;
	F32   trendResamplingOrderProb;
	F32   seasonResamplingOrderProb;
	U16   maxMoveStepSize;
	U32   burnin, samples, chainNumber;
	U16   thinningFactor;
} BVS_MCMC, * _restrict BVS_MCMC_PTR;
 
typedef struct BVS_EXTRA {
	I08   useMeanOrRndBeta;
	I08   dumpInputData;
	U08   numThreadsPerCPU;
	U16   numParThreads;
	U16   numCPUCoresToUse;
	U16   consoleWidth;
	I08   whichOutputDimIsTime;
	I08   removeSingletonDims;

	I08   ncpStatMethod;
	Bool  computeCredible;
	Bool  fastCIComputation;
	Bool  computeSeasonOrder;
	Bool  computeTrendOrder;
 

	Bool  printOptions;
	Bool  printProgressBar;

} BVS_EXTRA, * _restrict BVS_EXTRA_PTR;

typedef struct BVS_OPTIONS {
	BVS_IO		    io;
	BVS_MCMC		mcmc;
	BVS_EXTRA		extra;
	BVS_PRIOR		prior;
	VOIDPTR         ans; // used only in SaveOutput
} BVS_OPTIONS, * _restrict BVS_OPTIONS_PTR;

typedef struct BVS_XINFO {
	I32PTR RowsMissing;
	int    Nmissing;
	int    Nobj, Pobj;
	int    N,    Pall ;
	F32PTR Xorg, Xones;
	F32PTR X1norm, mean, sd, X2, mean2, sd2;
	

	I32    xi_Pbasis[MAX_NUM_BASIS];
	I32PTR xi_Pindices0[MAX_NUM_BASIS];
	I32PTR xi_PchangeY0;
	I32PTR xi_PseasonY0;

	I32PTR Nunique;
	I32PTR SortedUnqiueIdx0;

	F32PTR mem_GenOneSeg_hinege2d_2dr_1dr;// temp mem used in GenOneTerm_hinge2d 	
	I32    Nangle;   // for hinge2dr
	F32PTR anglecoff;// for hinge2dr

} BVS_XINFO, * _restrict BVS_XINFO_PTR;

typedef struct BVS_YINFO {
	F32PTR     Y;
	F32PTR     mean, sd;
	int        N;
	F32PTR     YtY_plus_alpha2Q;
	F32        alpha1_star; 
	
} BVS_YINFO, * _restrict BVS_YINFO_PTR;

typedef struct BVS_RESULT {		
	F32PTR  marg_lik, sig2, R2, RMSE;
	F32PTR  Y,  SD, CI; 
	F32PTR  Ylin, Ystep, Ystair, Ypiece, Yseason, Yhinge, YhingePair, Ychange;
	I32PTR  KProb;
	F32PTR  Kmean;
	I32PTR  probAll;
	I32PTR  probLin, probStep, probStair, probPiece, probSeason, probHinge, probHingePair,probChange; 
} BVS_RESULT, * _restrict BVS_RESULT_PTR;

 
typedef struct BVS_RNDSTREAM {
	F32PTR  rndgamma;
	U32PTR  rnd32;
	U16PTR  rnd16;
	U08PTR  rnd08;	
} BVS_RNDSTREAM, * _restrict BVS_RANDSEEDPTR;


struct BVS_BASIS;
struct BVS_MODEL;
typedef struct BVS_BASIS  * _restrict BVS_BASIS_PTR;
typedef struct BVS_MODEL BVS_MODEL, * _restrict BVS_MODEL_PTR;
typedef struct PROPOSE_STRUCT {	
	//CORESULT_PTR       keyresult;
	I32PTR              samples;	
	F32PTR				mem;
	BVS_MODEL_PTR	    model;
	BVS_XINFO_PTR       xinfo;
	BVS_RANDSEEDPTR		pRND;
	BVS_PRIOR_PTR       prior;
	I32                 Kmax;
	I32                 nSample_ExtremVecNeedUpdate;
	F32                 sigFactor;	
	F32                 outlierSigFactor;
} PROP_DATA, *_restrict PROP_DATA_PTR;


/************************************************/
typedef struct  __BVS_TERMS {
	TEMR_TYPE type;
	union {
		struct {
			I16 var0s[2];  // vars[0] is aliased to var, as used in IsTermsEqual()
			I32 knots[2];
			I08 sides[2];
		};
		struct {
			I16 var0;
			I32 knotprev;
			I32 knot;
			I32 knotnext;
			I08 side;
		};

	};
	
	F32 mean, sd;
} BVS_TERMS, * _restrict BVS_TERMS_PTR; ;

 
#define  CMPLIST_RANGE -2
 
typedef struct __BVS_TRANGE{I16 s,e;}BVS_TRANGE;
typedef struct  __BVS_TERMS_CTRL {
	TEMR_TYPE      type;
	BVS_TRANGE     range[5];
	I08           num;
} BVS_TERMS_CTRL, * _restrict BVS_TERMS_CTRL_PTR;

typedef union  __BVS_TERMS_ANY {
	struct {
		TEMR_TYPE    type;
		BVS_TRANGE   range[5];
		I08          num;
	};
	BVS_TERMS      tm;
	BVS_TERMS_CTRL ctrl;
} BVS_TERMS_ANY, * _restrict BVS_TERMS_ANY_PTR;

typedef struct __TERMS_ENCODER {	
	I64     nTermsMax, nTerms;	
	I16     Kmax, Kmax_actual;
	I32     nModels;
	I64     nBeta;
	I16PTR  mem; // has at least Kmax+(Kmax+1) elements

	I16               Kpre;
	BVS_TERMS_PTR     PreList;  // size: Kmax

	I16PTR            KorgVec;  // size: nSamples*nChains
	I16PTR            KcodeVec; // size: nSamples*nChains

	BVS_TERMS_ANY_PTR TERMS;     // size: nSamples*nChains*(Kmax-1)
	F32PTR            BETA;      // size: nSamples*nChains*(Kmax)

	//////////////////////////////////
	// FOR EXTRACTION ONLY
	I64 curTermID;
	I64 curBetaID;
	I32 curModelID;
	
}BVS_TERMS_ENCODER, * _restrict BVS_TERMS_ENCODER_PTR;
/************************************************/

typedef struct ___BVS_PRED_RESULT{
	F32PTR Y,SD;
	F32PTR Ychains;
	I32PTR Khinge;
	F32PTR Yparts, SDparts;
	I32PTR P;	
	I32PTR Plinear;
	I32PTR Psquare;
	I32PTR Phinge;
	I32PTR Pchange;
	F32PTR BetaLinear;
} BVS_PRED_RESULT, * _restrict BVS_PRED_RESULT_PTR;

/*
	F32PTR Xorg;
	F32PTR mean, sd, X2, mean2, sd2;
	int    N, P, Praw;
	I32PTR Xsqrind0, Xstepind0, Xhingeind0, Xchangeind0, Xhinge2dind0;
	int    Psqr, Pstep, Phinge, Pchange, Pchange_timevar0, Phinge2d;
	I32PTR Nunique;
	I32PTR SortedUnqiueIdx0;
*/

typedef struct ___BVS_NEWXINO {
	F32PTR          Xorg, Xnorm1, X2, Xones;
	I32             N, Pobj;
	BVS_PRED_RESULT* mat;
} BVS_PRED_XINFO, * _restrict BVS_PRED_XINFO_PTR;


typedef struct __BVS_PRED_OPT {
	int bComputeYchains;
	int odtype;
}BVS_PRED_OPT ;
//////////////////////////////////////////////////

typedef struct _NEWTERM {
	//I32           R1, R2; // start and end rows of the non-zero rowss: not used here
	I16PTR        KpositionBuf; // needed in the UpdateBasisAll as a buf of temporay Kpositon
	I32           idx0_in_cpList; // FOR MERGE ONLY 
	I32           knot_new; // for birth and move
	I32           knot_old; // for death and move
	BVS_TERMS     terms[3]; // piecewise may need three terms; hingpair needs two; others need only one
	NEWCOLINFOv2* xcols;	
	enum MOVETYPE jumpType;
	I08           nTerms;
} NewPropTerm, * _restrict NewPropTerm_PTR;


//period is used as "basis->bConst.dummy.period" in the following functons: computeXY,InitPrecPriorMEM,DD_CalcBasisKsKeK_prec0123
//TERMS needed only if meta->io.deseasonlaized=TRUE
 
 //https://gcc.gnu.org/onlinedocs/gcc/Unnamed-Fields.html

typedef struct PROP_PROB_STRUCT {U08 birth, death, merge, move;} PROP_PROB_STRUCT;
/*
typedef struct BVS_BASIS {
	struct {
  	    void        (*Propose)(BVS_BASIS_PTR, NewPropTerm *, PROP_DATA*);
		int         (*CalcBasisKsKeK_TermType)(BVS_BASIS_PTR  basis);
		void        (*UpdateGoodVec)(BVS_BASIS_PTR basis, NewPropTerm * new, I32 Npad16_not_used);
		void        (*ComputeY)(F32PTR X, F32PTR beta, F32PTR Y, BVS_BASIS_PTR basis, I32 Npad);
		F32         (*ModelPrior)(BVS_BASIS_PTR basis, NewPropTerm * new, I32 N);
	};	

	struct {
		I16    maxKnotNum;
		TORDER minOrder, maxOrder;
	} prior;

	PROP_PROB_STRUCT	propprob;
	I16					mcmc_Kstopping;
 
	I16      nPrec;
	I16      offsetPrec;		

} BVS_BASIS, * _restrict BVS_BASIS_PTR;
 */
typedef struct {	
	I16PTR    Kposition; // pmax+1: the first term fixed at the constant term
	I16       K;
	TEMR_TYPE type;
	I32       goodVarNum;
	I08PTR    goodVarBinVec;   // ((pmax+15)/16*16;	
	I32       P, Ppad16;
	I32       maxMovStep;
	I32PTR    pidx0_to_vidx0;
	/*--------------------------*/
} BSTATE_LINEAR;
typedef struct {	
	I16PTR  Kposition; // pmax+1: the first term fixed at the constant term
	I16     K;
	TEMR_TYPE type;
	I32     goodVarNum;
	I08PTR  goodVarBinVec;   // ((pmax+15)/16*16;	
	I32     P, Ppad16;
	I32     maxMovStep;
	I32PTR  pidx0_to_vidx0; //unused
	/*--------------------------*/
} BSTATE_QUADRATIC;

typedef struct {	
	I16PTR  Kposition;  //	
	I16     K;
	TEMR_TYPE type;
	I32     goodVarNum;
	I08PTR  goodVarBinVec;
	I32     P, Ppad16;
	I32     maxMovStep;
	I32PTR pidx0_to_vidx0; //unused
	/*--------------------------*/

	I08PTR  goodKnotBinMat;     //(N+15)/16*16 * pmax
	I32PTR  goodKnotNumVec;  // pmax	
	I32     Npad16; //(N+15)/16*16

	I32    MAXCPTNUM;
	I32PTR cptNumVec;
	I32PTR cptMat;
	I32    sbad, ebad;
	I16PTR maxCptNumVec;
	I16PTR minSepDistVec;

} BSTATE_STEP;

typedef struct {
	I16PTR  Kposition;  //	
	I16     K;
	TEMR_TYPE type;
	I32     goodVarNum;
	I08PTR  goodVarBinVec;
	I32     P, Ppad16;
	I32     maxMovStep;
	I32PTR  pidx0_to_vidx0;
	/*--------------------------*/

	I08PTR  goodKnotBinMat;   //(N+15)/16*16 * pmax
	I32PTR  goodKnotNumVec;   // pmax	
	I32     Npad16;           //(N+15)/16*16

	I32    MAXCPTNUM;
	I32PTR cptNumVec;
	I32PTR cptMat;
	I32    sbad, ebad;
	I16PTR maxCptNumVec;
	I16PTR minSepDistVec;
	
} BSTATE_STAIR, BSTATE_PIECE;

typedef struct {
	I16PTR  Kposition;  //	
	I16     K;
	TEMR_TYPE type;
	I32     goodVarNum;
	I08PTR  goodVarBinVec;
	I32     P, Ppad16;
	I32     maxMovStep;
	I32PTR  pidx0_to_vidx0;
	/*--------------------------*/

	I08PTR  goodKnotBinMat;     //(N+15)/16*16 * pmax
	I32PTR  goodKnotNumVec;  // pmax	
	I32     Npad16; //(N+15)/16*16

	I32    MAXCPTNUM;
	I32PTR cptNumVec;
	I32PTR cptMat;
	I32    sbad, ebad;
	I16PTR maxCptNumVec;
	I16PTR minSepDistVec;
	
} BSTATE_HINGE;

typedef struct {
	I16PTR  Kposition;  //	
	I16     K;
	TEMR_TYPE type;
	I32     goodVarNum;
	I08PTR  goodVarBinVec;
	I32     P, Ppad16;
	I32     maxMovStep;
	I32PTR pidx0_to_vidx0; //unused
	/*--------------------------*/

	I08PTR  goodKnotBinMat;     //(N+15)/16*16 * pmax
	I32PTR  goodKnotNumVec;  // pmax	
	I32     Npad16; //(N+15)/16*16

	I32    MAXCPTNUM;
	I32PTR cptNumVec;
	I32PTR cptMat;
	I32    sbad, ebad;
	I16PTR maxCptNumVec;
	I16PTR minSepDistVec;
 
} BSTATE_HINGEPAIR;

typedef struct {
	I16PTR  Kposition;  //	
	I16     K;
	TEMR_TYPE type;
	I32     goodVarNum;
	I08PTR  goodVarBinVec;
	I32     P, Ppad16;
	I32     maxMovStep;
	I32PTR  pidx0_to_vidx0; 
	/*--------------------------*/

	I08PTR  goodKnotBinMat;     //(N+15)/16*16 * pmax
	I32PTR  goodKnotNumVec;  // pmax	
	I32     Npad16; //(N+15)/16*16

	I32    MAXCPTNUM;
	I32PTR cptNumVec;
	I32PTR cptMat;
	I32    sbad, ebad;
	I16PTR maxCptNumVec;
	I16PTR minSepDistVec; 
} BSTATE_CHANGE;
typedef struct {
	I16PTR  Kposition;  //	
	I16     K;
	TEMR_TYPE type;
	I32     goodVarNum_unused;
	I08PTR  goodVarBinVec_unused;
	I32     P, Ppad16_unused;
	I32     maxMovStep_unused;
	I32PTR pidx0_to_vidx0; //unused
	/*--------------------------*/


	I08PTR  goodKnotBinMat_unused;     //(N+15)/16*16 * pmax
	I32PTR  goodKnotNumVec_unused;  // pmax	
	I32     Npad16_unused; //(N+15)/16*16

	I32    maxCptNum_unused;
	I32PTR cptNumVec_unused;
	I32PTR cptMat_unused;
	I32 sbad, ebad;// minSeptDist_unused;
	I32 minpts;

} BSTATE_HINGE2D;

typedef struct {
	I16PTR  Kposition;  //	
	I16     K;
	TEMR_TYPE type;
	I32     goodVarNum_unused;
	I08PTR  goodVarBinVec_unused;
	I32     P, Ppad16_unused;
	I32     maxMovStep_unused;
	I32PTR pidx0_to_vidx0; //unused
	/*--------------------------*/


	I08PTR  goodKnotBinMat_unused;     //(N+15)/16*16 * pmax
	I32PTR  goodKnotNumVec_unused;  // pmax	
	I32     Npad16_unused; //(N+15)/16*16

	I32    maxCptNum_unused;
	I32PTR cptNumVec_unused;
	I32PTR cptMat_unused;
	I32 sbad, ebad;// minSeptDist_unused;
	I32 minpts;

} BSTATE_HINGE2DR;

typedef union {
	struct {
		I16PTR  Kposition;  //	
		I16     K;
		TEMR_TYPE type;
		I32     goodVarNum;
		I08PTR  goodVarBinVec;
		I32     P, Ppad16;
		I32     maxMovStep;
		I32PTR pidx0_to_vidx0; //unused
	};

	BSTATE_LINEAR    lin;
	BSTATE_QUADRATIC sqr;
	BSTATE_STEP      step;
	BSTATE_PIECE     piece;
	BSTATE_HINGE     hinge;
	BSTATE_CHANGE    change;
} BSTATE_ANY;

typedef struct {
	F32PTR XtX, XtY, cholXtX, beta_mean;
	//F32PTR precXtXDiag;
	//I08PTR nTermsPerPrecGrp;
	union {
		F32      alpha2_star;
		F32PTR   alphaQ_star; // added for MRBEAST
	};	
	F32    marg_lik;
	I32    K;
} BVS_MODELDATA, *_restrict  BVS_MODELDATA_PTR;


typedef struct BVS_MODEL {
	F32 _sig2, _precVal, _logPrecVal;

	BVS_MODELDATA curr, prop;
	F32PTR  beta; 		

	F32PTR  sig2;	
	F32PTR  precVec; 
	F32PTR  logPrecVec;
	I16     nPrec;


	I08PTR08 extremePosVec;
	F32PTR   deviation;
	F32PTR   avgDeviation;  // Changed to a pointer for a consistent API with MRBEAST
	I32      extremPosNum;

	BVS_TERMS_PTR terms;
	I16           Kterms, Kfixed;

	///////XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	I16       NUMBASIS;
	I16       NumBasisPad16;
	I32       goodBasisBirthableNum;
	I08PTR    basisBirthable;
	I16PTR    basisType; //set to point to opt->basiTypes
	I16PTR    type2id;    //set to point to opt->type2id
	BSTATE_ANY* basisState;


} BVS_MODEL, * _restrict BVS_MODEL_PTR;

/*
typedef struct {
	I32           minSepDist;
	BVS_YINFO_PTR yInfo;
} KNOT2BINVEC, * _restrict KNOT2BINVEC_PTR;
*/

typedef struct BVS_HyperPar {
	F32 alpha_1, alpha_2, del_1, del_2;
} BVS_HyperPar, * _restrict BVS_HyperPar_PTR;


#define RoundTo64(N)  ((N+63)/64*64)
#define RoundTo16(N)  ((N+15)/16*16)
#define hasTerm(x)     (opt->prior.type2id[x] >=0 )
#define AggregatedMemAlloc 0
 

#include "assert.h"

#elif defined(SOLARIS_COMPILER)
#include "beastv2_header_solaris.h"
#endif

