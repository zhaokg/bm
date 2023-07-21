#include <string.h>  //memset
#include <math.h>    //sqrt

#include "abc_000_warning.h"

#include "abc_vec.h"
#include "abc_ts_func.h"
#include "abc_mem.h"
#include "abc_vec.h"
#include "bvs_header.h"

//f32_fill_val_matrixdiag(MODEL.SIG2, opt->prior.sig2, q);

#define MODEL               (*model)
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#define max3(a,b,c)         max( max(a,b), c)
#define max4(a,b,c,d)       max( max3(a,b,c),d)
#define min(a,b)            (((a) < (b)) ? (a) : (b))


static void Alloc_BasisState_LinearQuadratic(BSTATE_LINEAR *b, BVS_OPTIONS_PTR opt, MemPointers* MEM) { 
	int P              = opt->prior.Pbasis[b->type];
    int Ppad16         = RoundTo16(P);
	b->Kposition       = MyALLOC(*MEM,  P+1L,     I16, 2); // a better size is min(Kmax, P+1)
	b->goodVarBinVec   = MyALLOC(*MEM,  Ppad16,   I08, 8); // a better size is min(Kmax, P+1)
	if (b->type == LinearTerm) {
		b->pidx0_to_vidx0 = MyALLOC(*MEM, opt->prior.Pall, I32, 4);
	}
}

// Called after xinfo is filled with the true information
static int  Init_BasisState_LinearQuadratic(BSTATE_LINEAR* b, BVS_XINFO* xinfo) {
		
	int P        = xinfo->xi_Pbasis[b->type];
	int Ppad16   = RoundTo16(P);

	memset(b->goodVarBinVec, 0, Ppad16);
	memset(b->goodVarBinVec, 1, P);
	b->goodVarNum = P;
	b->K          = 0;
	b->P          = P;
	b->Ppad16     = Ppad16;

	if (b->type == LinearTerm) {
		// Fill the array to get the lndices0 for a linear term: to be used for checking DEATH
		
		for (int i = 0; i < xinfo->Pall; i++) {
			b->pidx0_to_vidx0[i] = -9999999L;
		}

		I32PTR Pindices0 = xinfo->xi_Pindices0[LinearTerm];
		for (int i = 0; i < P; i++) {
			// Whe pidx0_to_vidx0 was allocated, it has been initialized to zeros
			// So, zeros mean bad values
			b->pidx0_to_vidx0[Pindices0[i]] = i;
		}
	}

	int basisBrithaable = b->goodVarNum > 0;
	return basisBrithaable;
	
}

static void Alloc_BasisState_StepStairPieceSeasonHingeChange_HingePair(BSTATE_HINGE* b,  BVS_OPTIONS_PTR opt, MemPointers* MEM) {

	int N      = opt->io.meta.Nobj;
	int Npad16 = RoundTo16(N);
	int Kmax   = opt->prior.Kmax;
	
	int type = b->type;

	int P            = opt->prior.Pbasis[type];
	int Ppad16       = RoundTo16(P);
	b->Kposition     = MyALLOC(*MEM, Kmax,  I16, 2);
	b->goodVarBinVec = MyALLOC(*MEM, Ppad16, I08, 8); // a better size is min(Kmax, P+1)

	b->goodKnotBinMat = MyALLOC(*MEM, Npad16 * P, I08, 8); // a better size is min(Kmax, P+1)
	b->goodKnotNumVec = MyALLOC(*MEM, P, I32, 4); // a better size is min(Kmax, P+1)

	
	// Used to get the indices0 for a given linear term
    // Needed to eheck the linkage between linearterm and hinge/change term
	b->pidx0_to_vidx0 = (type == HingeTerm || type == ChangeTerm)
						? MyALLOC0(*MEM, opt->prior.Pall, I32, 4) 
		                : NULL;

	b->maxCptNumVec   = opt->prior.maxKnotNum_Vec[type];
	b->minSepDistVec  = opt->prior.minSepDist_Vec[type];
	b->MAXCPTNUM      = opt->prior.maxKnotNum_Max[type];
  
	b->cptNumVec = MyALLOC(*MEM, P, I32, 4); // a better size is min(Kmax, P+1)
	if (type == PieceTerm || type == StairTerm || type == SeasonTerm) {
		// the 1st and last elemants (1 and n) are added as fixed chpt
		b->MAXCPTNUM += 2;
		b->cptMat    = MyALLOC(*MEM, P * b->MAXCPTNUM, I32, 4);  
		b->cptMat++; // move ahead bye one so cptMat[-1]=1;
	} else {
		b->cptMat = MyALLOC(*MEM, P * b->MAXCPTNUM, I32, 4);  
	}
 
}

#define MIN_GOOD_POINTS 4
static int  Init_BasisState_StepStairPieceSeasonHingeChange_HingePair(BSTATE_HINGE* b, BVS_XINFO* xinfo) {

	int    type = b->type;
	int    P    = xinfo->xi_Pbasis[type];

	if (type == HingeTerm || type == ChangeTerm) {
	// Fill the array to get the lndices0 for a linear term: to be used for checking DEATH
		
		for (int i = 0; i < xinfo->Pall; i++) {
			b->pidx0_to_vidx0[i] = -9999999L;
		}
		I32PTR Pindices0 = xinfo->xi_Pindices0[b->type];
		for (int i = 0; i < P; i++) {
			// Whe pidx0_to_vidx0 was allocated, it has been initialized to zeros
			// So, zeros mean bad values
			b->pidx0_to_vidx0[Pindices0[i]] = i;
		}
	}


	if (type == StepTerm)      { b->sbad = 1;	b->ebad = 0;  }
	if (type == StairTerm)     { b->sbad = 1;	b->ebad = 1; }
	if (type == PieceTerm)     { b->sbad = 2;	b->ebad = 2;  }
	if (type == SeasonTerm)   { b->sbad = 2;	b->ebad = 2; }
	if (type == HingeTerm)     { b->sbad = 2;	b->ebad = 2;  }
	if (type == ChangeTerm)    { b->sbad = 2;	b->ebad = 2;  }
	if (type == HingePairTerm) { b->sbad = 2;	b->ebad = 2;  }
	
	int N      = xinfo->N;  // Not Nobj	
	int Ppad16 = RoundTo16(P);
	int Npad16 = RoundTo16(N);

	b->Npad16 = Npad16;
	b->P      = P;
	b->Ppad16 = Ppad16;
	b->K      = 0;	

	// Fill xinfo->Xones
	if (xinfo->Xones) {
		f32_fill_val(1.0, xinfo->Xones, N);
	}

	memset(b->goodVarBinVec,    0, Ppad16);
	memset(b->goodKnotBinMat,   0, Npad16 * P);
	memset(b->goodKnotNumVec,   0, sizeof(I32) * P);


	if (type == PieceTerm || type == StairTerm || type == SeasonTerm) {
		//cptMat has one extra eleme at the start and end of each col;
		//cptMat has also been move ahead by one elem so that cptMat[-1]=1;		
		I32PTR  cptMat = b->cptMat;
		cptMat--;
		memset(cptMat, 0, sizeof(I32) * P * b->MAXCPTNUM);
		for (int i = 0; i < P; i++) {
			cptMat[b->MAXCPTNUM*i] = 1L;                       // pre-fill the element at the negative 1 with 1
			if (type == PieceTerm) b->cptNumVec[i]  = -999999L; // a negative value indicating the variable not chosen at all
			if (type == StairTerm) b->cptNumVec[i]  = 0; 
			if (type == SeasonTerm) b->cptNumVec[i] = -999999L;
		}

	} else {
		memset(b->cptMat, 0, sizeof(I32) * P * b->MAXCPTNUM);
		memset(b->cptNumVec, 0, sizeof(I32) * P);
	}
	
	
	int s = b->sbad +1;
	int e = N - b->ebad;
 
	I08PTR goodKnotBinMat = b->goodKnotBinMat;
	I32PTR Pindices0      = xinfo->xi_Pindices0[type];
	int goodVarNum = 0;
	for (int i = 0; i < P; i++) {
		I32 pidx0  = Pindices0[i];
		I32 Nlen   = xinfo->Nunique[pidx0];
		I32 e      = Nlen - b->ebad;
		I32 segLen = e-s+1L;
		if (segLen > 0) {			
			memset(goodKnotBinMat + s - 1L, 1L, segLen);  // -1 converted to uint32 will be a super large number
			b->goodKnotNumVec[i] = segLen;
		} else {
			b->goodKnotNumVec[i] = 0;
		}
				
		goodKnotBinMat      += Npad16;
		b->goodVarBinVec[i]  = b->goodKnotNumVec[i] > 0 && b->cptNumVec[i] < b->maxCptNumVec[i];
		goodVarNum          += b->goodVarBinVec[i];
	}
	b->goodVarNum = goodVarNum;
	
	int basisBrithaable = b->goodVarNum > 0;
	return basisBrithaable;
}

 static void Alloc_BasisState_Hinge2d(BSTATE_HINGE2D* b,  BVS_OPTIONS_PTR opt, MemPointers* MEM) {

	int N      = opt->io.meta.Nobj;
	int Npad16 = RoundTo16(N);
	int Kmax   = opt->prior.Kmax;

	int P            = opt->prior.Pbasis[b->type];	 //b->type assigned in  Alloc_BasisState_All
	int Ppad16       = RoundTo16(P);
	b->Kposition     = MyALLOC(*MEM, Kmax,  I16, 2);	
	
}
 static int  Init_BasisState_Hinge2d(BSTATE_HINGE2D* b, BVS_XINFO* xinfo) {

	 int    P = xinfo->xi_Pbasis[b->type]; //b->type assigned in  Alloc_BasisState_All
	 b->sbad = 2;
	 b->ebad = 2;

 	 int N = xinfo->N;  // Not Nobj	
	 int Ppad16 = RoundTo16(P);
	 int Npad16 = RoundTo16(N);

	 b->K = 0;
	 b->P = P;
 
	 //int s = b->sbad + 1;
	// int e = N - b->ebad;

	 b->minpts = MIN_GOOD_POINTS;

	 int basisBrithaable = 1L; //always brithable
	 return basisBrithaable;
 
}
  
 static void Alloc_BasisState_Hinge2dr(BSTATE_HINGE2D* b,  BVS_OPTIONS_PTR opt, MemPointers* MEM) {

	int N      = opt->io.meta.Nobj;
	int Npad16 = RoundTo16(N);
	int Kmax   = opt->prior.Kmax;

	int P            = opt->prior.Pbasis[b->type];	 //b->type assigned in  Alloc_BasisState_All
	int Ppad16       = RoundTo16(P);
	b->Kposition     = MyALLOC(*MEM, Kmax,  I16, 2);	
 
}
 static int  Init_BasisState_Hinge2dr(BSTATE_HINGE2DR* b, BVS_XINFO* xinfo) {

	 int    P = xinfo->xi_Pbasis[b->type]; //b->type assigned in  Alloc_BasisState_All
	 b->sbad = 2;
	 b->ebad = 2;
 
	 int N       = xinfo->N;  // Not Nobj	
	 int Ppad16 = RoundTo16(P);
	 int Npad16 = RoundTo16(N);

	 b->K = 0; // number of terms
	 b->P = P; // number of allowed hinge variabls
 
	 b->minpts = MIN_GOOD_POINTS;
	 //int s = b->sbad + 1;
	// int e = N - b->ebad;

	 int basisBrithaable = 1L; //always brithable
	 return basisBrithaable;
}
 
  static void Alloc_BasisState_Hinge1dr(BSTATE_HINGE2D* b,  BVS_OPTIONS_PTR opt, MemPointers* MEM) {

	int N      = opt->io.meta.Nobj;
	int Npad16 = RoundTo16(N);
	int Kmax   = opt->prior.Kmax;

	int P            = opt->prior.Pbasis[b->type];	 //b->type assigned in  Alloc_BasisState_All
	int Ppad16       = RoundTo16(P);
	b->Kposition     = MyALLOC(*MEM, Kmax,  I16, 2);	
 
}
 static int  Init_BasisState_Hinge1dr(BSTATE_HINGE2DR* b, BVS_XINFO* xinfo) {

	 int    P = xinfo->xi_Pbasis[b->type]; //b->type assigned in  Alloc_BasisState_All
	 b->sbad = 2;
	 b->ebad = 2;
 
	 int N       = xinfo->N;  // Not Nobj	
	 int Ppad16 = RoundTo16(P);
	 int Npad16 = RoundTo16(N);

	 b->K = 0; // number of terms
	 b->P = P; // number of allowed hinge variabls
 
	 b->minpts = MIN_GOOD_POINTS;
	 //int s = b->sbad + 1;
	// int e = N - b->ebad;


	 int basisBrithaable = 1L; //always brithable
	 return basisBrithaable;
}



typedef void(*pfunAllocBasis)(
	BSTATE_HINGE* b, BVS_OPTIONS_PTR opt, MemPointers* MEM
	);
static pfunAllocBasis AllocBasisFuncs[MAX_NUM_BASIS] = {NULL,};

typedef int (*pfunInitBasis)(
	BSTATE_HINGE* b, BVS_XINFO* xinfo
	);
static pfunInitBasis InitBasisFuncs[MAX_NUM_BASIS] = {NULL, };

static void SetUp_BasisState_Functions(void) {
	if (AllocBasisFuncs[LinearTerm]) {
		return;
	}
	AllocBasisFuncs[LinearTerm]    = Alloc_BasisState_LinearQuadratic;
	AllocBasisFuncs[QuadraticTerm] = Alloc_BasisState_LinearQuadratic;
	AllocBasisFuncs[StepTerm]      = Alloc_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	AllocBasisFuncs[StairTerm]     = Alloc_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	AllocBasisFuncs[PieceTerm]     = Alloc_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	AllocBasisFuncs[SeasonTerm]    = Alloc_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	AllocBasisFuncs[HingeTerm]     = Alloc_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	AllocBasisFuncs[ChangeTerm]    = Alloc_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	AllocBasisFuncs[HingePairTerm] = Alloc_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	AllocBasisFuncs[Hinge2DTerm]   = Alloc_BasisState_Hinge2d;
	AllocBasisFuncs[Hinge2DRTerm]  = Alloc_BasisState_Hinge2dr;
	AllocBasisFuncs[Hinge1DRTerm]  = Alloc_BasisState_Hinge1dr;

	InitBasisFuncs[LinearTerm]    = Init_BasisState_LinearQuadratic;
	InitBasisFuncs[QuadraticTerm] = Init_BasisState_LinearQuadratic;
	InitBasisFuncs[StepTerm]      = Init_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	InitBasisFuncs[StairTerm]     = Init_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	InitBasisFuncs[PieceTerm]     = Init_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	InitBasisFuncs[SeasonTerm]    = Init_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	InitBasisFuncs[HingeTerm]     = Init_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	InitBasisFuncs[ChangeTerm]    = Init_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	InitBasisFuncs[HingePairTerm] = Init_BasisState_StepStairPieceSeasonHingeChange_HingePair;
	InitBasisFuncs[Hinge2DTerm]   = Init_BasisState_Hinge2d;
	InitBasisFuncs[Hinge2DRTerm]  = Init_BasisState_Hinge2dr;
	InitBasisFuncs[Hinge1DRTerm]  = Init_BasisState_Hinge1dr;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
void Alloc_BasisState_All(BVS_MODEL_PTR model, BVS_OPTIONS_PTR opt, MemPointers* MEM) {

	SetUp_BasisState_Functions();


	MODEL.basisType  = opt->prior.basisType;
	MODEL.type2id    = opt->prior.type2id;

	I32   NumBasis   = MODEL.NUMBASIS      = opt->prior.numBasis;
	I32   NumBasis16 = MODEL.NumBasisPad16 = RoundTo16(NumBasis);

	MODEL.basisState     = MyALLOC0(*MEM, NumBasis,  BSTATE_ANY, 64);
	MODEL.basisBirthable = MyALLOC0(*MEM, NumBasis16, I08,      8);// must be 8byte-aligned
	 	 
 
	for (int i = 0; i < model->NUMBASIS; ++i) {					
		I32 type = MODEL.basisState[i].type= model->basisType[i];
		AllocBasisFuncs[type](model->basisState+i,  opt, MEM);
	}
 
}
void Init_BasisState_All(BVS_MODEL_PTR model, BVS_OPTIONS_PTR opt, BVS_XINFO *xinfo) {
	
	memset(model->basisBirthable, 0, sizeof(I08)*model->NumBasisPad16);

	int goodBasisBirthableNum = 0;
	for (int i = 0; i < model->NUMBASIS; ++i) { 
		model->basisBirthable[i] = InitBasisFuncs[model->basisType[i]](model->basisState+i, xinfo);
		goodBasisBirthableNum += model->basisBirthable[i];
		 
	}
	model->goodBasisBirthableNum = goodBasisBirthableNum;
}

#include "abc_000_warning.h"