#include <string.h>
#include <math.h>
#include "abc_000_warning.h"
#include "abc_mcmc.h"  //RANDINT RANDINT_SKIPINE
#include "bvs_header.h"
#include "abc_sort.h"  //insertionSort
#include "abc_vec.h"   // for i08_sum_binvec only
#include "bvs_fun.h"   // GenTerms_hinge2d_getgoodpts

#define type2bs(type) &model->basisState[model->type2id[type]]
#define hasBasis(x)   (model->type2id[x] >= 0)

static int __i32_find_nth_nonzero_elem(I32PTR good, I32 N, I32 nth) {
	I32  idx0       = -999999;
	I32  csumPosLoc = 0;
	for (int i = 0; i < N; i++) {
		csumPosLoc += (good[i] > 0);
		if (csumPosLoc == nth) {
			idx0 = i;
			break;
		}
	}
	I32 idx1 = idx0 + 1L;
	return idx1;
}

static int  __pick_rand_varidx0(BSTATE_LINEAR* b, BVS_RANDSEEDPTR	PRND) {  
	// Randomly choose a good loc: newKnot is the newly chosen breakpoint	
    // NewVarIdx retunrd is 1-based 
	assert(b->goodVarNum > 0);
	I32  randLoc    = RANDINT(1, (I32)b->goodVarNum, *(PRND->rnd32)++);
	I32  newVarIdx  = i08_find_nth_onebyte_binvec(b->goodVarBinVec, (I32)b->Ppad16, randLoc);
	return newVarIdx-1L;
}
static int  __pick_rand_knot(BSTATE_STEP* b, int vidxsel0, BVS_RANDSEEDPTR	PRND) {
	I32  goodKnotNum = b->goodKnotNumVec[vidxsel0];
	I32  randLoc     = RANDINT(1, (I32)goodKnotNum, *(PRND->rnd32)++);
	I32  newKnotIdx1 = i08_find_nth_onebyte_binvec(b->goodKnotBinMat + vidxsel0 * b->Npad16, (I32)b->Npad16, randLoc);
	assert(goodKnotNum > 0);
	return newKnotIdx1;
}

void  Birth_LinearQuadratic(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, PROP_DATA_PTR info) {
	BVS_TERMS       *term = NEW->terms;
	BVS_RANDSEEDPTR	 PRND = info->pRND;
	BSTATE_LINEAR   * b   = type2bs(type);
	term->var0  = __pick_rand_varidx0(b, PRND);
	term->type   = type;	

	newcol->nbands = 1L;
	newcol->ks_x[0]        = newcol->K +1L;
	newcol->kterms_x[0]    = 0;
	newcol->ks_xnewterm[0]     = 1L;
	newcol->kterms_xnewterm[0] = 1L;
 
	NEW->nTerms= 1L;  // Number of terms: two terms for HingeParTErm
}

void   Birth_Step(NewPropTerm* NEW,  BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, PROP_DATA_PTR info) {
	BVS_TERMS       *term = NEW->terms;
	BVS_RANDSEEDPTR	PRND = info->pRND;
	BSTATE_HINGE*   b    = type2bs(type); 	
	term->var0   = __pick_rand_varidx0(b, PRND);		 // NewVarIdx retunrd above is 1-based but terms needs a 0-based index  
	term->knot   = __pick_rand_knot(b, term->var0, PRND);
	term->side   = 1;
	term->type   = type;

	newcol->nbands          = 1L;
	newcol->ks_x[0]            = newcol->K +1L;
	newcol->kterms_x[0]        = 0;
	newcol->ks_xnewterm[0]     = 1L;
	newcol->kterms_xnewterm[0] = 1L;
 
	NEW->knot_new  = term->knot;
	NEW->nTerms    = 1L;  // Number of terms: two terms for HingeParTErm	
}

static int __FindModelTermByKnot(BVS_MODEL_PTR model, BSTATE_HINGE* b, I32 vidx0, I32 knot ) {

	int kidx = -9999999L; 
	BVS_TERMS_PTR terms = model->terms; 
	for (int i = 0; i < b->K; i++) {
		int Ktmp = b->Kposition[i];
		if (terms[Ktmp - 1].var0 == vidx0 && terms[Ktmp - 1].knot == knot) {
			kidx = Ktmp;
			break;
		}
	}
 
	return kidx;
}

void  Birth_Stair(NewPropTerm* NEW,  BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, PROP_DATA_PTR info) {

	// If HInge and CHange are included, the linear term must be also included
	BVS_RANDSEEDPTR	PRND      = info->pRND;
	BVS_TERMS       *terms    = NEW->terms;
	BSTATE_HINGE    *b        = type2bs(type);	
 	 
	 I32    vidx0    = __pick_rand_varidx0(b, PRND);
	 I32    pidx0    = info->xinfo->xi_Pindices0[type][vidx0];
	 I32    Nlen     = info->xinfo->Nunique[pidx0];

	 I32     cptNum   = b->cptNumVec[vidx0];
	 I32PTR  cptList  = b->cptMat + vidx0* b->MAXCPTNUM;

	 terms[0].type = type;


	 /***************************************/
	 // b->cptNumVec[vidx0] >= 0
	 // Add a new knot and first pick up a random knot
	 /***************************************/
	 cptList[cptNum] = Nlen;  // make sure there is one extra element of N

	 I32 knotnew  = __pick_rand_knot(b, vidx0, PRND);


	 if (cptNum == 0) {
		 // if cptNum >0
		 terms[0].type = type;
		 terms[0].var0 = vidx0;
		 terms[0].knot     = 1;
		 terms[0].knotnext = knotnew;

		 terms[1].type = type;
		 terms[1].var0 = vidx0;
		 terms[1].knot = knotnew;
		 terms[1].knotnext = Nlen;

		 newcol->nbands    = 2L;
		 newcol->ks_x[0]     = newcol->K + 1L;;
		 newcol->ks_x[1]     = newcol->K + 1L;
		 newcol->kterms_x[0] = 0;
		 newcol->kterms_x[1] = 0;

		 newcol->ks_xnewterm[0] = 1L;
		 newcol->ks_xnewterm[1] = 2L;
		 newcol->kterms_xnewterm[0] = 1L;
		 newcol->kterms_xnewterm[1] = 1L;

		 NEW->knot_new = knotnew;
		 NEW->nTerms = 2L;  // Number of terms: two terms for HingeParTErm
		 return;

		 return;
	 }

	 // if cptNum >0
	 {
		 I32 knewidx0 = 0;
		 for (; knewidx0 < cptNum && knotnew > cptList[knewidx0]; knewidx0++);

		 // Add one term and change one term
		 int K1 = __FindModelTermByKnot(model, b, vidx0, cptList[knewidx0 - 1]);

		 terms[0] = model->terms[K1 - 1];
		 terms[0].knotnext = knotnew;

		 terms[1].type = type;
		 terms[1].var0 = vidx0;
		 terms[1].knot = knotnew;
		 terms[1].knotnext = cptList[knewidx0];

		 newcol->nbands = 2L;
		 newcol->ks_x[0] = K1;
		 newcol->ks_x[1] = newcol->K + 1L;
		 newcol->kterms_x[0] = 1;
		 newcol->kterms_x[1] = 0;

		 newcol->ks_xnewterm[0] = 1L;
		 newcol->ks_xnewterm[1] = 2L;
		 newcol->kterms_xnewterm[0] = 1L;
		 newcol->kterms_xnewterm[1] = 1L;

		 NEW->knot_new = knotnew;
		 NEW->nTerms = 2L;  // Number of terms: two terms for HingeParTErm
		 return;

	 }

}

void  Birth_PieceSeason(NewPropTerm* NEW,  BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, PROP_DATA_PTR info) {

	// If HInge and CHange are included, the linear term must be also included
	BVS_RANDSEEDPTR	PRND      = info->pRND;
	BVS_TERMS       *terms    = NEW->terms;
	BSTATE_HINGE    *b        = type2bs(type);	
 
 
	 I32    vidx0    = __pick_rand_varidx0(b, PRND);
	 I32    pidx0 = info->xinfo->xi_Pindices0[type][vidx0];
	 I32    Nlen  = info->xinfo->Nunique[pidx0];
 

	 I32     cptNum  = b->cptNumVec[vidx0];
	 I32PTR  cptList  = b->cptMat + vidx0* b->MAXCPTNUM;

	 terms[0].type = type;
	 
	 /***************************************/
	 // The variable not included at all. Then add a linear term
	 /***************************************/
	 if (cptNum < 0) {
		 terms[0].var0      = vidx0;
		 terms[0].knotprev  = 1;
		 terms[0].knot      = 1;
		 terms[0].knotnext = Nlen;
		 
		 newcol->nbands      = 1L;
		 newcol->ks_x[0]     = newcol->K + 1L;
		 newcol->kterms_x[0] = 0;
		 newcol->ks_xnewterm[0]     = 1L;
		 newcol->kterms_xnewterm[0] = 1L;

		 NEW->knot_new  = 1;
		 NEW->nTerms    = 1L;  // Number of terms: two terms for HingeParTErm
		 return;
	 }


	 /***************************************/
	 // b->cptNumVec[vidx0] >= 0
	 // Add a new knot and first pick up a random knot
	 /***************************************/
	 cptList[cptNum] = Nlen;  // make sure there is one extra element of N

	 I32 knotnew  = __pick_rand_knot(b, vidx0, PRND);
	 I32 knewidx0 = 0;	 
	 for (; knewidx0 < cptNum && knotnew > cptList[knewidx0]; knewidx0++);

	 // Add one term and change one term
	 if (knewidx0 == cptNum || type == SeasonTerm) {
		 int K1 = __FindModelTermByKnot(model, b, vidx0, cptList[knewidx0 - 1]);

		 terms[0]          = model->terms[K1-1];
		 terms[0].knotnext = knotnew;

		 terms[1].type     = type;
		 terms[1].var0     = vidx0;
		 terms[1].knotprev = terms[0].knot;
		 terms[1].knot     = knotnew;
		 terms[1].knotnext = cptList[knewidx0]; // is Nlen for PieceTerm
	 
		 newcol->nbands  = 2L;
		 newcol->ks_x[0] = K1;
		 newcol->ks_x[1] = newcol->K + 1L;
		 newcol->kterms_x[0] = 1;
		 newcol->kterms_x[1] = 0;

		 newcol->ks_xnewterm[0] = 1L;
		 newcol->ks_xnewterm[1] = 2L;
		 newcol->kterms_xnewterm[0] = 1L;
		 newcol->kterms_xnewterm[1] = 1L;

		 NEW->knot_new  = knotnew;
		 NEW->nTerms    = 2L;  // Number of terms: two terms for HingeParTErm
		 return;
	 }


	 // Add one term and change two terms
	 {
		 int K1 = __FindModelTermByKnot(model, b, vidx0, cptList[knewidx0 - 1]);
		 int K2 = __FindModelTermByKnot(model, b, vidx0, cptList[knewidx0]);

		 terms[0]          = model->terms[K1 - 1];
		 terms[0].knotnext = knotnew;

		 terms[1]          = model->terms[K2 - 1];
		 terms[1].knotprev = knotnew;


		 terms[2].type     = type;
		 terms[2].var0     = vidx0;
		 terms[2].knotprev = terms[0].knot;
		 terms[2].knot     = knotnew;
		 terms[2].knotnext = terms[1].knot;

		 if (K1 > K2) {
			 int tmp = K1;
			 K1 = K2;
			 K2 = tmp;
		 }

		 newcol->nbands = 3L;
		 newcol->ks_x[0] = K1;
		 newcol->ks_x[1] = K2;
		 newcol->ks_x[2] = newcol->K + 1L;
		 newcol->kterms_x[0] = 1;
		 newcol->kterms_x[1] = 1;
		 newcol->kterms_x[2] = 0;

		 newcol->ks_xnewterm[0] = 1L;
		 newcol->ks_xnewterm[1] = 2L;
		 newcol->ks_xnewterm[2] = 3L;
		 newcol->kterms_xnewterm[0] = 1L;
		 newcol->kterms_xnewterm[1] = 1L;
		 newcol->kterms_xnewterm[2] = 1L;

		 NEW->knot_new = knotnew;
		 NEW->nTerms   = 3L;  // Number of terms: two terms for HingeParTErm
		 return;	 
	 }

}

 

static int ___DetermineSide_HingeChange(BVS_TERMS* term, BVS_XINFO_PTR xinfo, BVS_RANDSEEDPTR	PRND) {
	
	return   *PRND->rnd08++ > 128;

	// IT IS A BAD STRAGETY TO USE THE RELATIVE POSITION OF KNOT TO DETERMINE THE SIDE
	// THIS WILL CREATE AN ARITIFICAL BREAKPOINT in THE MODDEL

	int    pidx0 = xinfo->xi_Pindices0[term->type][term->var0];  // // zero-based for basisState but 1-based for vars
	int    N     = xinfo->N;
	F32PTR X     = xinfo->Xorg + N * pidx0;

	int    knotidx = xinfo->SortedUnqiueIdx0[N * pidx0 + term->knot - 1L];

	F32 knot = X[knotidx];
	F32 mean = xinfo->mean[pidx0]; 
	F32 sd   = xinfo->sd[pidx0];
	F32 diff = (knot - mean);     
	if (diff > sd/2 || diff < sd / 2) {
		return diff > 0 ? 99:-99;
	}
	return   *PRND->rnd08++ > 128;

}

void  Birth_HingeChange(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, PROP_DATA_PTR info) {

	// If HInge and CHange are included, the linear term must be also included
	BVS_TERMS       *term = NEW->terms;
	BVS_RANDSEEDPTR	PRND = info->pRND;
	BSTATE_HINGE    *b    = type2bs(type);
	BSTATE_HINGE    *blin = type2bs(LinearTerm);  
	
	I32PTR          Pindices0 = info->xinfo->xi_Pindices0[type];

	#define IsLinearTermIncluded(p0)  (blin->goodVarBinVec[blin->pidx0_to_vidx0[p0]] == 0)

	I32 vidx0;
	// (1) FIrst try our luck by picking up some random var index from b->goodVarBinVec	
	I32  foundIncludedVar = 0;
	for (int i = 0; i < 3; ++i) {			
		vidx0 = __pick_rand_varidx0(b, PRND);
		if ( IsLinearTermIncluded(Pindices0[vidx0]) ) {
			// has been inlcuded in the model
			foundIncludedVar = 1;
			break;
		}
	}


	// (2) If not finding one, do a more intensive search by considering only the included linear terms
	if (!foundIncludedVar) {
	
		// Check which var has the linar term included
		I32PTR linVarInUse   = info->mem;
		I32    goodLinVarNum = 0;
		for (int i = 0; i < b->P; ++i) {			
			if (IsLinearTermIncluded(Pindices0[i]) && b->goodVarBinVec[i]) {
				linVarInUse[goodLinVarNum++] = i;
			}
		}

		if (goodLinVarNum > 0) {			 
			vidx0            = linVarInUse[ RANDINT(0, goodLinVarNum-1L, *(PRND->rnd32)++) ];
			foundIncludedVar = 1;
		}
	}

	if (foundIncludedVar) {			
			term->var0 = vidx0;
			term->knot = __pick_rand_knot(b, vidx0, PRND);
			term->type = type;
			term->side = ___DetermineSide_HingeChange(term, info->xinfo, PRND); //  *(PRND->rnd08)++ > 128; //TODO
			
			newcol->nbands = 1L;
			newcol->ks_x[0] = newcol->K + 1L;
			newcol->kterms_x[0] = 0;
			newcol->ks_xnewterm[0] = 1L;
			newcol->kterms_xnewterm[0] = 1L;
			NEW->knot_new = term->knot;
			NEW->nTerms   = 1L;  // Number of terms: two terms for HingeParTErm
			return;
	} 

	// (3) OUt of luck ,finding no good varbiable at all
	// That is, no linear terms included or the inluclued linar terms have no free knots
	Birth_LinearQuadratic(NEW, model, LinearTerm, newcol, info);

}

void  Birth_HingePair(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, PROP_DATA_PTR info) {

	// If HInge and CHange are included, the linear term must be also included
	BVS_TERMS       *term = NEW->terms;
	BVS_RANDSEEDPTR	PRND = info->pRND;
	BSTATE_HINGE*   b    = type2bs(type);
 
	I32  vidx0= __pick_rand_varidx0(b, PRND); ;

	term[0].var0  = vidx0;
	term[0].knot = __pick_rand_knot(b, vidx0, PRND);
	term[0].type = type;
	term[0].side = 0; //  *(PRND->rnd08)++ > 128; //TODO

	// TODO: terms[1] points to the NEW.terms2 field
	term[1].var0   = term[0].var0;
	term[1].knot  = term[0].knot;
	term[1].type  = type;
	term[1].side  = 1; //  *(PRND->rnd08)++ > 128; //TODO

	newcol->nbands       = 1L;
	newcol->ks_x[0]      = newcol->K + 1L;
	newcol->kterms_x[0]        = 0;
	newcol->ks_xnewterm[0]     = 1L;
	newcol->kterms_xnewterm[0] = 2L;
	NEW->knot_new = term->knot;
	NEW->nTerms      = 2L;  // Number of terms: two terms for HingeParTErm	
	
}

void  Birth_Hinge2d(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, PROP_DATA_PTR info) {
	// Randomly choose a good loc: newKnot is the newly chosen breakpoint	

	// NewVarIdx retunrd above is 1-based 
	BVS_TERMS       *term = NEW->terms;
	BVS_RANDSEEDPTR	 PRND = info->pRND;
	BSTATE_HINGE2D*  b    = type2bs(type);

	I32  P  = b->P;
	I32  v0 = 1;
	I32  v1 = 2;

	if (P > 2) {		
		BVS_PRIOR_PTR prior = info->prior;
		int nGrp   = RANDINT(1, prior->nGrpHinge2d, *PRND->rnd16++);  
		int start = prior->sGrpHinge2d[nGrp - 1];
		int end   = prior->eGrpHinge2d[nGrp - 1];		
		if ( end-start+1 == 2) {
			v0 = start;
			v1 = end;
		}	else {
			v0 = RANDINT(start, end, *PRND->rnd16++);
			RANDINT_SKIPONE(v1, start, v0, end, *PRND->rnd16++);
		}     	
	}

	int Vidx0_1 = (min(v0, v1)) - 1L;
	int Vidx0_2 = (max(v0, v1)) - 1L;

	term->var0s[0] = Vidx0_1;
	term->var0s[1] = Vidx0_2;

	BVS_XINFO_PTR xinfo = info->xinfo;
	int P0 = xinfo->xi_Pindices0[type][Vidx0_1];
	int P1 = xinfo->xi_Pindices0[type][Vidx0_2];

	int sSkip = b->sbad;
	int eSkip = b->ebad;

	I32        bingo       = 0;
	I32        bestGoodPts = -1L;
	BVS_TERMS  bestTerms;
	for (int i = 0; i < 10; i++) {

		term->knots[0] = RANDINT(sSkip+1, xinfo->Nunique[P0]-eSkip, *PRND->rnd32++); //TODO: buggy if no knot are available
		term->knots[1] = RANDINT(sSkip+1, xinfo->Nunique[P1]-eSkip, *PRND->rnd32++);

		term->sides[0] = *PRND->rnd08++ > 128;
		term->sides[1] = *PRND->rnd08++ > 128;

		term->type = type;// GenOneTerm_goodpoints needs to call GetKnot which requires type
		int goodPts = GenOneTerm_hinge2d_GetGoodDataNpts(term, info->xinfo, info->mem);
		if (goodPts > b->minpts) {
			bingo = 1;
			break;
		}	else {
			if (goodPts > bestGoodPts) {
				bestGoodPts = goodPts;
				bestTerms.knots[0] = term->knots[0];
				bestTerms.knots[1] = term->knots[1];
				bestTerms.sides[0] = term->sides[0];
				bestTerms.sides[1] = term->sides[1];
			}
		}
	}
	
	if (!bingo) {
		term->knots[0] = bestTerms.knots[0];
		term->knots[1] = bestTerms.knots[1];
		term->sides[0]= bestTerms.sides[0]  ;
		term->sides[1] = bestTerms.sides[1];
	}

	term->type    = type;

	newcol->nbands       = 1L;
	newcol->ks_x[0]      = newcol->K + 1L;
	newcol->kterms_x[0]        = 0;
	newcol->ks_xnewterm[0]     = 1L;
	newcol->kterms_xnewterm[0] = 1L;

	NEW->nTerms = 1L;  // Number of terms: two terms for HingeParTErm	
}
void  Birth_Hinge2dr(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, PROP_DATA_PTR info) {
	// Randomly choose a good loc: newKnot is the newly chosen breakpoint	

	// NewVarIdx retunrd above is 1-based 
	BVS_TERMS       *term = NEW->terms;
	BVS_RANDSEEDPTR	 PRND = info->pRND;
	BSTATE_HINGE2D*  b    = type2bs(type);

	I32  P  = b->P;
	I32  v0 = 1, v1 = 2;

	if (P > 2) {		
		BVS_PRIOR_PTR prior = info->prior;
		int nGrp   = RANDINT(1, prior->nGrpHinge2dr, *PRND->rnd16++);  
		int start = prior->sGrpHinge2dr[nGrp - 1];
		int end   = prior->eGrpHinge2dr[nGrp - 1];		
		if ( end-start+1 == 2) {
			v0 = start;
			v1 = end;
		}	else {
			v0 = RANDINT(start, end, *PRND->rnd16++);
			RANDINT_SKIPONE(v1, start, v0, end, *PRND->rnd16++);
		}     	
	}

	int Vidx0_1 = (min(v0, v1)) - 1L;
	int Vidx0_2 = (max(v0, v1)) - 1L;

	term->var0s[0] = Vidx0_1;
	term->var0s[1] = Vidx0_2;

	BVS_XINFO_PTR xinfo = info->xinfo;
	int P0 = xinfo->xi_Pindices0[type][Vidx0_1];
	int P1 = xinfo->xi_Pindices0[type][Vidx0_2];

	int sSkip = b->sbad;
	int eSkip = b->ebad;
	 
	I32        bingo = 0;
	I32        bestGoodPts = -1L;
	BVS_TERMS  bestTerms;
	for (int i = 0; i < 10; i++) {
		term->knots[0] = RANDINT(sSkip+1, xinfo->Nunique[P0]-eSkip, *PRND->rnd32++);
		term->knots[1] = RANDINT(sSkip+1, xinfo->Nunique[P1]-eSkip, *PRND->rnd32++);

		term->sides[0] = RANDINT(0, xinfo->Nangle-1, *PRND->rnd16++);
		term->sides[1] = RANDINT(0, xinfo->Nangle-1, *PRND->rnd16++);
 
		term->type = type;// GenOneTerm_goodpoints needs to call GetKnot which requires type
		int goodPts = GenOneTerm_hinge2dr_GetGoodDataNpts(term, info->xinfo, info->mem);
		if (goodPts > b->minpts) {
			bingo = 1;
			break;
		} else {
			if (goodPts > bestGoodPts) {
				bestGoodPts = goodPts;
				bestTerms.knots[0] = term->knots[0];
				bestTerms.knots[1] = term->knots[1];
				bestTerms.sides[0] = term->sides[0];
				bestTerms.sides[1] = term->sides[1];
			}
		}
	}

	if (!bingo) {
		term->knots[0] = bestTerms.knots[0];
		term->knots[1] = bestTerms.knots[1];
		term->sides[0] = bestTerms.sides[0];
		term->sides[1] = bestTerms.sides[1];
	}

	term->type    = type;

	newcol->nbands = 1L;
	newcol->ks_x[0] = newcol->K + 1L;
	newcol->kterms_x[0] = 0;
	newcol->ks_xnewterm[0] = 1L;
	newcol->kterms_xnewterm[0] = 1L;

	NEW->nTerms = 1L;  // Number of terms: two terms for HingeParTErm
}

void  Birth_Hinge1dr(NewPropTerm* NEW,  BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, PROP_DATA_PTR info) {
	// Randomly choose a good loc: newKnot is the newly chosen breakpoint	

	// NewVarIdx retunrd above is 1-based 
	BVS_TERMS       *term = NEW->terms;
	BVS_RANDSEEDPTR	 PRND = info->pRND;
	BSTATE_HINGE2D*  b    = type2bs(type);

	I32  P  = b->P;
	I32  v0 = 1, v1 = 2;

	if (P > 2) {		
		BVS_PRIOR_PTR prior = info->prior;
		int nGrp   = RANDINT(1, prior->nGrpHinge1dr, *PRND->rnd16++);  
		int start = prior->sGrpHinge1dr[nGrp - 1];
		int end   = prior->eGrpHinge1dr[nGrp - 1];		
		if ( end-start+1 == 2) {
			v0 = start;
			v1 = end;
		}	else {
			v0 = RANDINT(start, end, *PRND->rnd16++);
			RANDINT_SKIPONE(v1, start, v0, end, *PRND->rnd16++);
		}     	
	}

	int Vidx0_1 = (min(v0, v1)) - 1L;
	int Vidx0_2 = (max(v0, v1)) - 1L;

	term->var0s[0] = Vidx0_1;
	term->var0s[1] = Vidx0_2;

	BVS_XINFO_PTR xinfo = info->xinfo;
	int P0 = xinfo->xi_Pindices0[type][Vidx0_1];
	int P1 = xinfo->xi_Pindices0[type][Vidx0_2];

	int sSkip = b->sbad;
	int eSkip = b->ebad;

	I32        bingo       = 0;
	I32        bestGoodPts = -1L;
	BVS_TERMS  bestTerms;
	for (int i = 0; i < 10; i++) {

		term->knots[0] = RANDINT(sSkip+1, xinfo->Nunique[P0]-eSkip, *PRND->rnd32++); //TODO: buggy if no knot are available
		term->knots[1] = RANDINT(sSkip+1, xinfo->Nunique[P1]-eSkip, *PRND->rnd32++);

		term->sides[0] = RANDINT(0, xinfo->Nangle - 1, *PRND->rnd16++);
		term->sides[1] = 0;// filled to make sure the terms comparison operators work

		term->type = type;// GenOneTerm_goodpoints needs to call GetKnot which requires type
		int goodPts = GenOneTerm_hinge1dr_GetGoodDataNpts(term, info->xinfo, info->mem);
		if (goodPts > b->minpts) {
			bingo = 1;
			break;
		}	else {
			if (goodPts > bestGoodPts) {
				bestGoodPts = goodPts;
				bestTerms.knots[0] = term->knots[0];
				bestTerms.knots[1] = term->knots[1];
				bestTerms.sides[0] = term->sides[0];
				bestTerms.sides[1] = term->sides[1];
			}
		}
	}

 	if (!bingo) {
		term->knots[0] = bestTerms.knots[0];
		term->knots[1] = bestTerms.knots[1];
		term->sides[0] = bestTerms.sides[0];
		term->sides[1] = bestTerms.sides[1];
	}

	term->type    = type;

	newcol->nbands       = 1L;
	newcol->ks_x[0]      = newcol->K + 1L;
	newcol->kterms_x[0]        = 0;
	newcol->ks_xnewterm[0]     = 1L;
	newcol->kterms_xnewterm[0] = 1L;

	NEW->nTerms = 1L;  // Number of terms: two terms for HingeParTErm
}


/*******************************************************************/
static int IsLinearTermKillable(BVS_MODEL_PTR model, BVS_TERMS_PTR term, PROP_DATA_PTR info) {
	// Check if a liinear term is linked to a Hinge or Change term

	I32   vidx0  = term->var0;
	I32   pidx0  = info->xinfo->xi_Pindices0[LinearTerm][vidx0];
	
	if (hasBasis(HingeTerm)) {
	// Need to make sure the linear term has no associated HingTerm and ChangeTerm
		BSTATE_HINGE* b     = type2bs(HingeTerm);   
		I32           vidx0 = b->pidx0_to_vidx0[pidx0];
		if (vidx0 >=0 && b->cptNumVec[vidx0] > 0) {
			// there is a hinge term assoacited with the var being in use.
			// vidx0 = -1 when the linear term is not used for Hingterm
			return 0;
		}
	}

	if (hasBasis(ChangeTerm)) {
		// Need to make sure the linear term has no associated HingTerm and ChangeTerm
		BSTATE_HINGE* b      = type2bs(ChangeTerm);   
		I32           vidx0  = b->pidx0_to_vidx0[pidx0] ;
		if (vidx0 >= 0 && b->cptNumVec[vidx0] > 0) {
			// there is a hinge term assoacited with the var being in use.
			// vidx0 = -1 when the linear term is not used for Hingterm
			return 0;
		}
	}

	return 1L;
}

I32 CheckMovalbe_Linear(BVS_MODEL_PTR model, BVS_TERMS_PTR term, PROP_DATA_PTR info) {

	I32   bid = model->type2id[LinearTerm];

	if (!model->basisBirthable[bid]) {
	// No more terms available to add/switch, so it is not movable, either.
		return 0;
	}

	if (!IsLinearTermKillable(model, term, info)) {
	// If the current term is associated with some in-use Hinge/Change term, do not kill or move it out
		return 0;
	}
	
	// No hingerterm and changeterm basis are present
	// basisBirthable[bid]must be TRUE up to this point
	return model->basisBirthable[bid];	
	// if basisBirthable[bid]=TRUE, no need to check b->goodNum
}

I32 CheckMovalbe_Quadratic(BVS_MODEL_PTR model, BVS_TERMS_PTR term, PROP_DATA_PTR info_unused) {
	I32   bid = model->type2id[QuadraticTerm];
	return model->basisBirthable[bid];
	// if basisBirthable[bid]=TRUE, no need to check b->goodNum
}

I32 CheckMovalbe_StepStairHingeChange(BVS_MODEL_PTR model, BVS_TERMS_PTR term, PROP_DATA_PTR info_unused) {

	I32           bid = model->type2id[term->type];
	BSTATE_HINGE* bs = &model->basisState[bid];
	//BUGGY: a step/hinge/change basis is not birthable due to ether no good knot or the max number of knots reached
	// Even if it is not birhtable, it may be still movabale, as long as there are more good knots available

	// if (!model->basisBirthable[bid]) return 0;		

	int           vidx0 = term->var0;	
	return bs->goodKnotNumVec[vidx0] > 0; // No need to check b->knotNumVec[pidx]  > maxKnotNum
	 
}

I32 CheckMovalbe_PieceSeason(BVS_MODEL_PTR model, BVS_TERMS_PTR term, PROP_DATA_PTR info_unused) {
	I32           bid = model->type2id[term->type];
	BSTATE_HINGE* bs  = &model->basisState[bid];
	//BUGGY: a step/hinge/change basis is not birthable due to ether no good knot or the max number of knots reached
	// Even if it is not birhtable, it may be still movabale, as long as there are more good knots available
	
	// if (!model->basisBirthable[bid]) return 0;		

	int   vidx0 = term->var0;
	if (bs->cptNumVec[vidx0] <= 0) {
		assert(bs->cptNumVec[vidx0] == 0); // if there is a term, cptNum must be >-0
		return 0;  // if it is zero, it is just a global linear term
	} else {
		// this is not accurate. The right way is to check the bin vec within a window around the current knot
		return bs->goodKnotNumVec[vidx0] > 0; // No need to check b->knotNumVec[pidx]  > maxKnotNum
	}
	

}
static int __IsKnotTwoSided_HingePair(BVS_TERMS_PTR term , BSTATE_HINGEPAIR* bs, BVS_MODEL_PTR model) {
	int vidx0 = term->var0;
	int knot  = term->knot;
	int sides = 1-term->side;

	BVS_TERMS* terms = model->terms;
	for (int i = 0; i < bs->K; i++) {
		int k = bs->Kposition[i];
		if (terms[k - 1].var0 == vidx0 && terms[k - 1].knot == knot && terms[k - 1].side == sides) {
			return i+1;
		}			
	}
	return 0;
}

I32 CheckMovalbe_HingePair(BVS_MODEL_PTR model, BVS_TERMS_PTR term, PROP_DATA_PTR info_unused) {

	I32  bid = model->type2id[term->type];
	BSTATE_HINGEPAIR* bs = &model->basisState[bid];
	//BUGGY: a step/hinge/change basis is not birthable due to ether no good knot or the max number of knots reached
	// Even if it is not birhtable, it may be still movabale, as long as there are more good knots available
		
	// DO NOT MOVE a two-sided knot
	if (__IsKnotTwoSided_HingePair(term, bs, model)) {
		return 0;
	}

	int    vidx0 = term->var0;	
	return bs->goodKnotNumVec[vidx0] > 0; // No need to check b->knotNumVec[pidx]  > maxKnotNum

}

I32 CheckMovalbe_Hinge2d_Hinge2dr_Hinge1dr(BVS_MODEL_PTR model, BVS_TERMS_PTR term, PROP_DATA_PTR info_unused) {
	return 1L;
}

// Move a term
#define mwin  5

int  MoveBirth_LinearQuadraticTerm(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, int kmove, PROP_DATA_PTR info) {
	// Randomly choose a good loc: newKnot is the newly chosen breakpoint	

	newcol->nbands = 1;
	newcol->ks_x[0]            = kmove;
	newcol->kterms_x[0]        = 1L;
	newcol->ks_xnewterm[0]     = 1L;
	newcol->kterms_xnewterm[0] = 1L;

	BVS_TERMS_PTR    term  = NEW->terms;
	BVS_RANDSEEDPTR	 PRND = info->pRND;
	BSTATE_LINEAR   *b    = type2bs(type);
	 
	term->var0  = __pick_rand_varidx0(b, PRND);;
	term->type  = type;

	NEW->nTerms = 1;
	return 1;

	
	I32  varidx1  = term->var0 + 1L;
	
	// First, try to find a neighoring var
	I32 s = varidx1 - mwin;
	I32 e = varidx1 + mwin;
	s = max(s, 1);
	e = min(e, b->P);
	if (s == e) {
		// there is only one var, b->P==1
		return 0; 
	}

	I32 goodVarInSeg = 0;
	for (int i = s; i <= e; ++i) {
		goodVarInSeg += b->goodVarBinVec[i-1];
	}

	I32  varidx0_new;
	if (goodVarInSeg > 0) {
	// If there is any var chosable within the moving segment
		I32 rndLoc = RANDINT(1, goodVarInSeg, *PRND->rnd16++);
		I32 cumsum = 0;
		for (int i = s; i <= e; ++i) {
			cumsum += b->goodVarBinVec[i - 1L];
			if (cumsum == rndLoc) {
				varidx0_new = i-1L;
				break;
			}
		}
	}	
	else {
	// if not found, then choose any of them.
		varidx0_new = __pick_rand_varidx0(b, PRND);
	}

	
	term->var0  = varidx0_new;
	term->type  = type;
	NEW->nTerms = 1;
	return 1;
}

int  MoveBirth_StepHingeTermChange(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, int kmove, PROP_DATA_PTR info) {
	 
	newcol->nbands             = 1;
	newcol->ks_x[0]            = kmove;
	newcol->kterms_x[0]        = 1L;
	newcol->ks_xnewterm[0]     = 1L;
	newcol->kterms_xnewterm[0] = 1L;

	BVS_TERMS_PTR    term  = NEW->terms;
	term[0] = model->terms[kmove - 1];

	//(111) First, if the term is flippable, it is ranonly flopped to the other sides 
    #define IsSideFippable(x)  ((x==0) || (x==1))
	if ( IsSideFippable(term->side)) {
	// if the knot is around the mean, flip it randomly
		if (*info->pRND->rnd08++  < (255*0.3) ) {
			//TODO: randomly flip the sides
			term->side    = 1 - term->side;
			NEW->knot_old  = term->knot;
			NEW->knot_new  = term->knot;
			NEW->nTerms    = 1;
			return 1L;
		}
	}
	 

	//(222) If not flipped, randomly choose a good loc: newKnot is the newly chosen breakpoint	
 
	BSTATE_HINGE  *b = type2bs(type);

	I32PTR Pindices0 =  info->xinfo->xi_Pindices0[b->type];
	I32    vidx0     =  term->var0;
	I32    pidx0     = Pindices0[vidx0];
	I32    N         = info->xinfo->Nunique[pidx0];
	
	//I08PTR goodKnotBin = b->goodKnotBinMat + varidx0 * b->Npad16;
	I32PTR  knotList = b->cptMat + vidx0 * b->MAXCPTNUM;
	I32     KnotNum  = b->cptNumVec[vidx0];

	I32   knotold = term->knot;
	I32   idx0;
	for (idx0 = 0; idx0 < KnotNum && knotList[idx0] != knotold; idx0++);
	// idxo is the index of knot wrt KnotList

	I32 msep      = b->minSepDistVec[vidx0];
	I32 prevStart = idx0 == 0         ? b->sbad + 1 : (knotList[idx0 - 1] + msep+1L);
	I32 nextEnd   = idx0 == KnotNum-1 ? N- b->ebad  : (knotList[idx0 + 1] - msep-1L);

	I32 s = max(knotold - mwin, prevStart);
	I32 e = min(knotold + mwin, nextEnd);

	if (s == e) {
		return 0;
	}

	assert(e>s);
	I32 knotnew;
	RANDINT_SKIPONE(knotnew, s, knotold, e, *(info->pRND->rnd32)++);

	term->type  = type;
	term->var0   = vidx0;
	term->knot  = knotnew;
	//term->sides[0] = ___DetermineSide_HingeChange(term, info->xinfo, info->pRND); //*(info->pRND->rnd08)++ > 128; //TODO
	term->side    = term->side; // keep the sides the same when moving the knot
	NEW->knot_new = knotnew;
	NEW->knot_old = knotold;
	NEW->nTerms   = 1;

	return 1;
}

int  MoveBirth_Stair(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, int kmove, PROP_DATA_PTR info) {

	// We assume there must be at least one cpt.
	BSTATE_HINGE* b = type2bs(type);
	BVS_TERMS* term = &model->terms[kmove - 1];
	I32PTR Pindices0 = info->xinfo->xi_Pindices0[type];
	I32    vidx0    = term->var0;
	I32    pidx0    = Pindices0[vidx0];
	I32    N       = info->xinfo->Nunique[pidx0];

	I32     cptNum  = b->cptNumVec[vidx0];
	I32PTR  cptList = b->cptMat + vidx0 * b->MAXCPTNUM;
	cptList[cptNum] = N; // make sure there is an extra N at the end

	I32 knotold = term->knot;
	int idx0;
	if (knotold == 1L) {
		// Illegal to move the first segment, we have to randomly choose another knot
		idx0    = RANDINT(0L, (cptNum - 1), *(info->pRND->rnd32)++);
		knotold = cptList[idx0];
		// kmove must be updated
		kmove = __FindModelTermByKnot(model, b, vidx0, knotold);
		assert(kmove > 1);
	}
	else {
		for (idx0 = 0; idx0 < cptNum && knotold != cptList[idx0]; idx0++);
		// kidx0 is the index of knot wrt KnotList
		assert(idx0 < cptNum);
	}


	I32 msep = b->minSepDistVec[vidx0];
	I32 prevStart = idx0 == 0 ? b->sbad + 1 : (cptList[idx0 - 1] + msep + 1L);
	I32 nextEnd = idx0 == cptNum - 1 ? N - b->ebad : (cptList[idx0 + 1] - msep - 1L);
	I32 s = max(knotold - mwin, prevStart);
	I32 e = min(knotold + mwin, nextEnd);

	if (s == e) {
		return 0;
	}
	assert(e > s);

	I32 knotnew;
	RANDINT_SKIPONE(knotnew, s, knotold, e, *(info->pRND->rnd32)++);

	I32 K1;
	K1 = __FindModelTermByKnot(model, b, vidx0, cptList[idx0 - 1]);

	BVS_TERMS* terms = NEW->terms;
	terms[0]          = model->terms[K1 - 1];
	terms[0].knotnext = knotnew;

	terms[1].type     = type;
	terms[1].var0     = vidx0;
	terms[1].knot     = knotnew;
	terms[1].knotnext = cptList[idx0 + 1];

	// if the sleected knot is the last one in the cptList: Only two terms need to changed
	int k1 = min(K1, kmove);
	int k2 = max(K1, kmove);
	newcol->nbands = 2;
	newcol->ks_x[0] = k1;
	newcol->ks_x[1] = k2;
	newcol->kterms_x[0] = 1L;
	newcol->kterms_x[1] = 1L;
	newcol->ks_xnewterm[0] = 1L;
	newcol->ks_xnewterm[1] = 2L;
	newcol->kterms_xnewterm[0] = 1L;
	newcol->kterms_xnewterm[1] = 1L;

	NEW->knot_new = knotnew;
	NEW->knot_old = knotold;
	NEW->nTerms   = 2;
 
	return 1;

}

int  MoveBirth_PieceSeason(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, int kmove, PROP_DATA_PTR info) {
	
	// We assume there must be at least one cpt.
	BSTATE_HINGE*   b    = type2bs(type);
	BVS_TERMS   *   term = &model->terms[kmove - 1];
	I32PTR Pindices0 = info->xinfo->xi_Pindices0[type];
	I32    vidx0     = term->var0;
	I32    pidx0     = Pindices0[vidx0];
	I32    N         = info->xinfo->Nunique[pidx0];

	I32     cptNum  = b->cptNumVec[vidx0];
	I32PTR  cptList = b->cptMat + vidx0 * b->MAXCPTNUM;
	cptList[cptNum] = N; // make sure there is an extra N at the end

	I32 knotold  = term->knot;
	int idx0;
	if (knotold == 1L) {
		// Illegal to move the first segment, we have to randomly choose another knot
		idx0    = RANDINT(0L, (cptNum-1), *(info->pRND->rnd32)++);
		knotold = cptList[idx0];
		// kmove must be updated
		kmove = __FindModelTermByKnot(model, b, vidx0, knotold);
		assert(kmove > 1);
	} else {
		for (idx0=0; idx0 < cptNum && knotold != cptList[idx0]; idx0++);
		// kidx0 is the index of knot wrt KnotList
		assert(idx0 < cptNum);
	}


	I32 msep       = b->minSepDistVec[vidx0];
	I32 prevStart  = idx0 == 0          ? b->sbad + 1 : (cptList[idx0 - 1] + msep + 1L);
	I32 nextEnd    = idx0 == cptNum - 1 ? N - b->ebad : (cptList[idx0 + 1] - msep - 1L);
	I32 s = max(knotold - mwin, prevStart);
	I32 e = min(knotold + mwin, nextEnd);

	if (s == e) {
		return 0;
	}
	assert(e > s);

	I32 knotnew;
	RANDINT_SKIPONE(knotnew, s, knotold, e, *(info->pRND->rnd32)++);

	I32 K1;
	K1  = __FindModelTermByKnot(model, b, vidx0 , cptList[idx0 - 1]);

	BVS_TERMS* terms = NEW->terms;
	terms[0]          = model->terms[K1 - 1];
	terms[0].knotnext = knotnew;

	terms[1].type     = type;
	terms[1].var0     = vidx0;
	terms[1].knotprev = cptList[idx0 - 1];
	terms[1].knot     = knotnew;
	terms[1].knotnext = cptList[idx0 + 1];


	if (idx0 == cptNum - 1 || type == SeasonTerm) {
	// if the sleected knot is the last one in the cptList: Only two terms need to changed
	// For seasonterm, knotprev is ignored.
		int k1 = min(K1, kmove);
		int k2 = max(K1, kmove);
		newcol->nbands  =  2;
		newcol->ks_x[0] = k1;
		newcol->ks_x[1] = k2;
		newcol->kterms_x[0] = 1L;
		newcol->kterms_x[1] = 1L;

		newcol->ks_xnewterm[0]  = 1L;
		newcol->ks_xnewterm[1]  = 2L;
		newcol->kterms_xnewterm[0] = 1L;
		newcol->kterms_xnewterm[1] = 1L;

		NEW->knot_new = knotnew;
		NEW->knot_old = knotold;
		NEW->nTerms =2;
	} else {
	  // knot is not the last one and then three terms will be changed
      // This is the rightmost term

		int K2 = __FindModelTermByKnot(model, b, vidx0, cptList[idx0 + 1]);

		terms[2]          = model->terms[K2 - 1];
		terms[2].knotprev = knotnew;

		int m1 = min(min(K1, kmove), K2);
		int m3 = max(max(K1, kmove), K2);
		int m2 = (K1 > m1 && K1 < m3)   ? K1 :
			     (K2 > m1 && K2 < m3)   ? K2 :
			     kmove;
		newcol->nbands  = 3;
		newcol->ks_x[0] = m1;
		newcol->ks_x[1] = m2;
		newcol->ks_x[2] = m3;
		newcol->kterms_x[0] = 1L;
		newcol->kterms_x[1] = 1L;
		newcol->kterms_x[2] = 1L;

		newcol->ks_xnewterm[0]  = 1L;
		newcol->ks_xnewterm[1]  = 2L;
		newcol->ks_xnewterm[2]  = 3L;
		newcol->kterms_xnewterm[0] = 1L;
		newcol->kterms_xnewterm[1] = 1L;
		newcol->kterms_xnewterm[2] = 1L;

		NEW->knot_new = knotnew;
		NEW->knot_old = knotold;
		NEW->nTerms = 3;
	}

	return 1;
	
}

int  MoveBirth_HingePair(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, int kmove, PROP_DATA_PTR info) {

	newcol->nbands = 1;
	newcol->ks_x[0]            = kmove;
	newcol->kterms_x[0]        = 1L;
	newcol->ks_xnewterm[0]     = 1L;
	newcol->kterms_xnewterm[0] = 1L;


	BVS_TERMS_PTR    term = NEW->terms;
	term[0] = model->terms[kmove - 1];

	// ONLY ONE-SIDED KNOT CAN BE MOVED
	// At this point, __IsKnotTwoSided_HingePar has been called, and the knot is movable
	
	 //(111) First, if the term is flippable, it is ranonly flopped to the other sides
	
	// if the knot is around the mean, flip it randomly
	if (*info->pRND->rnd08++ < (255 * 0.3)) {
		//TODO: randomly flip the sides
		term->side = 1 - term->side;
		NEW->knot_new = term->knot;
		NEW->knot_old = term->knot;
		NEW->nTerms = 1;
		return 1L;	
	}
	 

	//(222) If not flipped, randomly choose a good loc: newKnot is the newly chosen breakpoint	
 
	BSTATE_HINGE  *b = type2bs(type);

	I32  vidx0  = term->var0;
	I32  pidx0    = info->xinfo->xi_Pindices0[type][vidx0];
	I32  N        = info->xinfo->Nunique[pidx0];
	
	//I08PTR goodKnotBin = b->goodKnotBinMat + varidx0 * b->Npad16;
	I32PTR  knotList = b->cptMat + vidx0 * b->MAXCPTNUM;
	I32     KnotNum  = b->cptNumVec[vidx0];

	I32   knotold = term->knot;
	I32   idx0;
	for (idx0 = 0; idx0 < KnotNum && knotList[idx0] != knotold; idx0++);
	// idxo is the index of knot wrt KnotList

	I32 msep      = b->minSepDistVec[vidx0];
	I32 prevStart = idx0 == 0         ? b->sbad + 1 : (knotList[idx0 - 1] + msep+1L);
	I32 nextEnd   = idx0 == KnotNum-1 ? N- b->ebad  : (knotList[idx0 + 1] - msep-1L);

	I32 s = max(knotold - mwin, prevStart);
	I32 e = min(knotold + mwin, nextEnd);

	if (s == e) {
		return 0;
	}

	assert(e>s);
	I32 knotnew;
	RANDINT_SKIPONE(knotnew, s, knotold, e, *(info->pRND->rnd32)++);

	term->type  =  type;
	term->var0  = vidx0;
	term->knot  = knotnew;
	term->side  = term->side; //Keep the same sides
	NEW->knot_new = knotnew;
	NEW->knot_old = knotold;
	NEW->nTerms   = 1;
	return 1;
}

int  MoveBirth_Hinge2d(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, int kmove, PROP_DATA_PTR info) {

	Birth_Hinge2d(NEW, model, type, newcol, info);
 
	// Brith_hinge2d filled newcol, assusming an insertion  at the end of X
	newcol->ks_x[0]     = kmove;
	newcol->kterms_x[0] = 1L;
	return 1L;
}
	// Randomly choose a good loc: newKnot is the newly chosen breakpoint	

int  MoveBirth_Hinge2dr(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, int kmove, PROP_DATA_PTR info) {

	Birth_Hinge2dr(NEW, model, type, newcol, info);

	// Brith_hinge2d filled newcol, assusming an insertion  at the end of X
	newcol->ks_x[0]     = kmove;
	newcol->kterms_x[0] = 1L;	
	return 1L;
}

int  MoveBirth_Hinge1dr(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, int kmove, PROP_DATA_PTR info) {

	Birth_Hinge1dr(NEW, model, type, newcol,info);

	// Brith_hinge2d filled newcol, assusming an insertion  at the end of X
	newcol->ks_x[0]     = kmove;
	newcol->kterms_x[0] = 1L;	

	return 1L;
}

typedef void (*pfuncBrith)(
	NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, PROP_DATA_PTR info
	);
static pfuncBrith BirthFuncs[] = {
	Birth_LinearQuadratic,
	Birth_LinearQuadratic,
	Birth_Step,
	Birth_Stair,
	Birth_PieceSeason, //piece
	Birth_PieceSeason, //season
	Birth_HingeChange,
	Birth_HingePair,
	Birth_HingeChange,
	Birth_Hinge2d,
	Birth_Hinge2dr,
	Birth_Hinge1dr,
};

// info is needed only for IsLinearTermKillable in CheckMovalbe_Linear.
typedef int(*pfuncCheckMove)(
	BVS_MODEL_PTR model, BVS_TERMS_PTR term, PROP_DATA_PTR info
	);
static pfuncCheckMove CheckMoveFuncs[] = {
	CheckMovalbe_Linear,
	CheckMovalbe_Quadratic,
	CheckMovalbe_StepStairHingeChange,
	CheckMovalbe_StepStairHingeChange,
	CheckMovalbe_PieceSeason, //piece
	CheckMovalbe_PieceSeason, //season
	CheckMovalbe_StepStairHingeChange,
	CheckMovalbe_HingePair,
	CheckMovalbe_StepStairHingeChange,
	CheckMovalbe_Hinge2d_Hinge2dr_Hinge1dr,
	CheckMovalbe_Hinge2d_Hinge2dr_Hinge1dr,
	CheckMovalbe_Hinge2d_Hinge2dr_Hinge1dr,
};

typedef int(*pfuncMoveBirth)(
	NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, int kmove,PROP_DATA_PTR info
	);
static pfuncMoveBirth MoveFuncs[] = {
	MoveBirth_LinearQuadraticTerm,
	MoveBirth_LinearQuadraticTerm,
	MoveBirth_StepHingeTermChange,
	MoveBirth_Stair,
	MoveBirth_PieceSeason, //piece
	MoveBirth_PieceSeason, //season
	MoveBirth_StepHingeTermChange,
	MoveBirth_HingePair,
	MoveBirth_StepHingeTermChange,
	MoveBirth_Hinge2d,
	MoveBirth_Hinge2dr,
	MoveBirth_Hinge1dr,
};

static int __WhichMove__(BVS_MODEL_PTR model, PROP_DATA_PTR info){
	
	I32 Kcur = model->curr.K;	
	if (Kcur == model->Kfixed) {
		return BIRTH;
	}

	BVS_RANDSEEDPTR	PRND = info->pRND;
	U08             rnd = *(PRND->rnd08)++;
	

	if (model->goodBasisBirthableNum == 0) {
		// return DEATH; BUGGY: it can get stuck to the death forever, although the true knot can be still reached via MOVE
		return (rnd > 128) ? DEATH : MOVE;
	}


	I32 Kmax = info->Kmax;
	if ( (Kcur+2) > Kmax) { // if ( Kcur == Kmax)  HingeParTerm will add two terms in the birht Proposal
		return (rnd>128)? DEATH:MOVE;
	}

	if (rnd < 255 * 0.333) {
		return BIRTH;
	} else if (rnd < 255 * 0.66) {
		return MOVE;
	}   else  {
		return DEATH;
	}
}

static int IsTermMergible(BVS_MODEL_PTR model, BVS_TERMS_PTR term, BVS_RANDSEEDPTR	PRND, I32 * idx0_in_cpList, I32* pKdel2_knotnext) {

	int type = term->type;
	if (!(type == StepTerm || type == HingeTerm || type == ChangeTerm || type==PieceTerm ||
		 type == StairTerm)) {
	// Only these three terms have the merge propose
		return 0; //return 0 for HINGE2D and LinearTerm
	}

	BSTATE_STEP  *b = type2bs(type);
	int vidx0  = term->var0;
	int cptNum = b->cptNumVec[vidx0];
	if (cptNum < 2)
	// there are at least two changepoints present to do the merge 
		return 0;
	
	I32    knot    = term->knot;
	I32PTR cptList = b->cptMat + vidx0 * b->MAXCPTNUM;	
	if (knot == cptList[cptNum-1]  || /*PieceTerm/StairTerm: idx0=-1L*/ knot==1L)
	// If the selected knot is the last changepoint
		return 0;
 
	U08 rnd = *(PRND->rnd08)++;
	if (rnd > 128)
		return 0;

	// Find the idx of knot
	I32 idx0 = 0; for (; cptList[idx0] != knot; idx0++);

	I32 knotnext = cptList[idx0 + 1L];
	if (knot == knotnext) // no knots bewteen the two changepoints
		return 0;
	I32 newknot = RANDINT(knot + 1L, knotnext - 1L, *PRND->rnd32++);

	I32 Kdel2_knotnext = __FindModelTermByKnot(model, b, vidx0 , knotnext);

    idx0_in_cpList[0]   = idx0;
	pKdel2_knotnext[0]  = Kdel2_knotnext;
	return newknot;
}

int  Death_Stair(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, int kdel, PROP_DATA_PTR info) {

	BVS_TERMS*    terms = NEW->terms;
	BSTATE_PIECE* b     = type2bs(type);

	I32    vidx0 = model->terms[kdel - 1].var0;
	I32    pidx0 = info->xinfo->xi_Pindices0[type][vidx0];
	I32    Nlen = info->xinfo->Nunique[pidx0];

	I32    cptNum = b->cptNumVec[vidx0];
	I32PTR cptList = b->cptMat + b->MAXCPTNUM * vidx0;

	//There are at least one changepoint
   //Selet a ranomd knot
	int idx1 = RANDINT(1L, cptNum, *(info->pRND->rnd16)++);
	int idx0 = idx1 - 1L;

	int K1 = __FindModelTermByKnot(model, b, vidx0, cptList[idx0 - 1]);
	int K2 = __FindModelTermByKnot(model, b, vidx0, cptList[idx0]);
	

	cptList[cptNum] = Nlen;

	if (cptNum == 1) {
		//It is a linearr temr and simply remove it
		terms[0] = model->terms[K1 - 1];
	 
		newcol->nbands      = 2;
		newcol->ks_x[0]     = min(K1, K2);
		newcol->ks_x[1]     = max(K1, K2);
		newcol->kterms_x[0] = 1L;
		newcol->kterms_x[1] = 1L;

		newcol->ks_xnewterm[0]     = 1L;
		newcol->ks_xnewterm[1]     = 1L;
		newcol->kterms_xnewterm[0] = 0L;
		newcol->kterms_xnewterm[1] = 0L;

		NEW->knot_old = cptList[idx0];
		NEW->nTerms   = 0;                  //Used to indicate no new term is generated
		NEW->jumpType = DEATH;
		return;
	}


	// REmove one term and change one term
  
	if (cptNum >= 1) {

		terms[0]          = model->terms[K1 - 1];
		terms[0].knotnext = cptList[idx0+1L];

		int k1 = min(K1, K2);
		int k2 = max(K1, K2);

		newcol->nbands      = 2;
		newcol->ks_x[0]     = k1;
		newcol->ks_x[1]     = k2;
		newcol->kterms_x[0] = 1L;
		newcol->kterms_x[1] = 1L;

		newcol->ks_xnewterm[0] = 1;
		newcol->ks_xnewterm[1] = 2;
		newcol->kterms_xnewterm[0] = 1L;
		newcol->kterms_xnewterm[1] = 0L;

		NEW->knot_old = cptList[idx0];
		NEW->nTerms    = 1L;
		NEW->jumpType = DEATH;
		return;
	}
 
	assert(0);
}

int  Merge_Stair(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, 
	             int kdel, int kdel2_knotnext, int idx0_in_cplist, I32 knotnew, PROP_DATA_PTR info) 
{

	BVS_TERMS    * terms = NEW->terms;
	BSTATE_PIECE * b   = type2bs(type);

	I32    vidx0   = model->terms[kdel - 1].var0;
	I32    pidx0 = info->xinfo->xi_Pindices0[type][vidx0];
	I32    Nlen = info->xinfo->Nunique[pidx0];

	I32    cptNum = b->cptNumVec[vidx0];
	I32PTR cptList = b->cptMat + b->MAXCPTNUM * vidx0;
	
	//There are at least twp changepoints 
	int idx0 = idx0_in_cplist;

	cptList[cptNum] = Nlen;
	// Remove two terms and generate two terms
		// 1, p0, p1, ...,[p_cn-2], p_cn_1, N

	int K1 = __FindModelTermByKnot(model, b, vidx0, cptList[idx0 - 1]);

	terms[0]          = model->terms[K1 - 1];
	terms[0].knotnext = knotnew;

	terms[1].type     = terms[0].type;
	terms[1].var0     = vidx0;
	terms[1].knot     = knotnew;
	terms[1].knotnext = cptList[idx0+2L];

	// K1, Kdel, Kdel2

	int m1 = min(min(K1, kdel), kdel2_knotnext);
	int m3 = max(max(K1, kdel), kdel2_knotnext);
	int m2 = K1 > m1 && K1 < m3 ? K1 :
		kdel > m1 && kdel < m3 ? kdel :
		kdel2_knotnext;
	newcol->nbands = 3;
	newcol->ks_x[0] = m1;
	newcol->ks_x[1] = m2;
	newcol->ks_x[2] = m3;
	newcol->kterms_x[0] = 1L;
	newcol->kterms_x[1] = 1L;
	newcol->kterms_x[2] = 1L;

	newcol->ks_xnewterm[0] = 1;
	newcol->ks_xnewterm[1] = 2;
	newcol->ks_xnewterm[2] = 3;
	newcol->kterms_xnewterm[0] = 1L;
	newcol->kterms_xnewterm[1] = 1L;
	newcol->kterms_xnewterm[2] = 0L;

	NEW->knot_new = knotnew;
	NEW->idx0_in_cpList = idx0;
	NEW->nTerms = 2L;
	NEW->jumpType = MERGE;
	return; 
}

int  Death_PieceSeason(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, int kdel, PROP_DATA_PTR info) {

	BVS_TERMS* terms = NEW->terms;
	BSTATE_PIECE* b = type2bs(type);

	I32    vidx0  = model->terms[kdel - 1].var0;
	I32    pidx0 = info->xinfo->xi_Pindices0[type][vidx0];
	I32    Nlen   = info->xinfo->Nunique[pidx0];

	I32    cptNum = b->cptNumVec[vidx0];
	I32PTR cptList = b->cptMat + b->MAXCPTNUM * vidx0;

	if (cptNum == 0) {
		//It is a linearr temr and simply remove it
		terms[0] = model->terms[kdel - 1];
		assert(terms[0].knot == 1);

		newcol->nbands       = 1;
		newcol->ks_x[0]      = kdel;
		newcol->kterms_x[0]  = 1L;
		newcol->ks_xnewterm[0]     = 1L;
		newcol->kterms_xnewterm[0] = 0L;
		NEW->knot_old = 1L;
		NEW->nTerms   = 0;                  //Used to indicate no new term is generated
		NEW->jumpType = DEATH;
		return;
	}

	//There are at least one changepoint
	//Selet a ranomd knot
	int idx0 = RANDINT(0L, cptNum-1L, *(info->pRND->rnd16)++);

	cptList[cptNum] = Nlen;
	// REmove one term and change one term
	if (idx0 == cptNum - 1 || type == SeasonTerm) {

		int K1 = __FindModelTermByKnot(model, b, vidx0, cptList[idx0 - 1]);
		int K2 = __FindModelTermByKnot(model, b, vidx0, cptList[idx0]);

		terms[0]          = model->terms[K1 - 1];
		terms[0].knotnext = cptList[idx0+1]; // is Nlen for PieceTerm

		int k1 = min(K1, K2);
		int k2 = max(K1, K2);

		newcol->nbands      = 2;
		newcol->ks_x[0]     = k1;
		newcol->ks_x[1]     = k2;
		newcol->kterms_x[0] = 1L;
		newcol->kterms_x[1] = 1L;

		newcol->ks_xnewterm[0] = 1;
		newcol->ks_xnewterm[1] = 2;
		newcol->kterms_xnewterm[0] = 1L;
		newcol->kterms_xnewterm[1] = 0L;

		NEW->knot_old = cptList[idx0];
		NEW->nTerms    = 1L;
		NEW->jumpType = DEATH;
		return;
	}

	// REmove one term and change two term
	if (idx0 != cptNum - 1) {

		int K1 = __FindModelTermByKnot(model, b, vidx0 , cptList[idx0 - 1]);
		int K2 = __FindModelTermByKnot(model, b, vidx0 , cptList[idx0]);
		int K3 = __FindModelTermByKnot(model, b, vidx0 , cptList[idx0 + 1L]);

		terms[0] = model->terms[K1 - 1];
		terms[0].knotnext = model->terms[K3 - 1].knot;

		terms[1] = model->terms[K3 - 1];
		terms[1].knotprev = model->terms[K1 - 1].knot;

		int m1 = min(min(K1, K2), K3);
		int m3 = max(max(K1, K2), K3);
		int m2 = K1 > m1 && K1 < m3 ? K1 :
			K2>m1 && K2 < m3 ? K2 :
			K3;

		newcol->nbands = 3;
		newcol->ks_x[0] = m1;
		newcol->ks_x[1] = m2;
		newcol->ks_x[2] = m3;
		newcol->kterms_x[0] = 1L;
		newcol->kterms_x[1] = 1L;
		newcol->kterms_x[2] = 1L;

		newcol->ks_xnewterm[0] = 1;
		newcol->ks_xnewterm[1] = 2;
		newcol->ks_xnewterm[2] = 3;
		newcol->kterms_xnewterm[0] = 1L;
		newcol->kterms_xnewterm[1] = 1L;
		newcol->kterms_xnewterm[2] = 0L;

		NEW->knot_old = cptList[idx0];
		NEW->nTerms   = 2L;
		NEW->jumpType = DEATH;
		return;
	}
}

int  Merge_PieceSeason(NewPropTerm* NEW, BVS_MODEL_PTR model, int type, NEWCOLINFOv2* newcol, 
	             int kdel, int kdel2_knotnext, int idx0_in_cplist, I32 knotnew, PROP_DATA_PTR info) 
{

	BVS_TERMS    * terms = NEW->terms;
	BSTATE_PIECE * b   = type2bs(type);

	I32    vidx0   = model->terms[kdel - 1].var0;
	I32    pidx0   = info->xinfo->xi_Pindices0[type][vidx0];
	I32    Nlen    = info->xinfo->Nunique[pidx0];

	I32    cptNum = b->cptNumVec[vidx0];
	I32PTR cptList = b->cptMat + b->MAXCPTNUM * vidx0;
	
	//There are at least twp changepoints 
	int idx0 = idx0_in_cplist;

	cptList[cptNum] = Nlen;
	// Remove two terms and generate two terms
	if (idx0 == (cptNum - 1)-1L || type==SeasonTerm) {
		// 1, p0, p1, ...,[p_cn-2], p_cn_1, N

		int K1 = __FindModelTermByKnot(model, b, vidx0 , cptList[idx0 - 1]);
	
		terms[0]          = model->terms[K1 - 1];
		terms[0].knotnext = knotnew;

		terms[1].type     = terms[0].type;
		terms[1].var0     = vidx0;
		terms[1].knotprev = cptList[idx0 - 1];
		terms[1].knot     = knotnew;
		terms[1].knotnext = cptList[idx0+2L]; // is Nlen for PieceTerm

		// K1, Kdel, Kdel2
		
		int m1 = min(min(K1, kdel), kdel2_knotnext);
		int m3 = max(max(K1, kdel), kdel2_knotnext);
		int m2 = K1 > m1 && K1 < m3     ? K1  :
				 kdel > m1 && kdel < m3 ? kdel:
			                              kdel2_knotnext;

		newcol->nbands      = 3;
		newcol->ks_x[0]     = m1;
		newcol->ks_x[1]     = m2;
		newcol->ks_x[2]     = m3;
		newcol->kterms_x[0] = 1L;
		newcol->kterms_x[1] = 1L;
		newcol->kterms_x[2] = 1L;

		newcol->ks_xnewterm[0] = 1;
		newcol->ks_xnewterm[1] = 2;
		newcol->ks_xnewterm[2] = 3;
		newcol->kterms_xnewterm[0] = 1L;
		newcol->kterms_xnewterm[1] = 1L;
		newcol->kterms_xnewterm[2] = 0L;

		NEW->knot_new       = knotnew;
		NEW->idx0_in_cpList = idx0;
		NEW->nTerms        = 2L;
		NEW->jumpType      = MERGE;
		return;
	}

	// REmove 2 term and change three term
	if (idx0 != cptNum - 2) {
		assert(idx0 < cptNum - 2);

		// 1, p0, p1, ...,[p],p_cn-2, p_cn_1, N

		int K1 = __FindModelTermByKnot(model, b, vidx0, cptList[idx0 - 1]);		
		int K4 = __FindModelTermByKnot(model, b, vidx0, cptList[idx0 + 2L]);

	

		terms[0]          = model->terms[K1 - 1];
		terms[0].knotnext = knotnew;

		terms[1].type     = terms[0].type;
		terms[1].var0      = vidx0;
		terms[1].knotprev = cptList[idx0 - 1];
		terms[1].knot     = knotnew;
		terms[1].knotnext = cptList[idx0 + 2L];;


		terms[2]          = model->terms[K4 - 1];
		terms[2].knotprev = knotnew;
 

		I32 Karr[4] = { K1, kdel, kdel2_knotnext, K4 };
		I32 index[4];
		InsertionSort_I32(Karr, index, 4);

		newcol->nbands    = 4;
		newcol->ks_x[0]     = Karr[0];
		newcol->ks_x[1]     = Karr[1];
		newcol->ks_x[2]     = Karr[2];
		newcol->ks_x[3]     = Karr[3];
		newcol->kterms_x[0] = 1L;
		newcol->kterms_x[1] = 1L;
		newcol->kterms_x[2] = 1L;
		newcol->kterms_x[3] = 1L;

		newcol->ks_xnewterm[0] = 1;
		newcol->ks_xnewterm[1] = 2;
		newcol->ks_xnewterm[2] = 3;
		newcol->ks_xnewterm[3] = 4;
		newcol->kterms_xnewterm[0] = 1L;
		newcol->kterms_xnewterm[1] = 1L;
		newcol->kterms_xnewterm[2] = 1L;
		newcol->kterms_xnewterm[3] = 0L;

		NEW->knot_new       = knotnew;
		NEW->idx0_in_cpList = idx0;
		NEW->nTerms        = 3L;
		NEW->jumpType      = MERGE;
		return;
	}
}

void ProposeMove( BVS_MODEL_PTR model, NewPropTerm* NEW, NEWCOLINFOv2* newcol,PROP_DATA_PTR info) {

	BVS_RANDSEEDPTR	PRND = info->pRND;

	I32 jumpType  = __WhichMove__(model, info);

	if (jumpType == BIRTH) {
 	// Geneate a new term: goodBasisBirthableNum >0 is checked in WhichMove
		I32  goodBasisNum = model->goodBasisBirthableNum;
		I32  randLoc      = RANDINT(1, (I32)goodBasisNum, *(PRND->rnd16)++);
		I32  newBasisIdx  = i08_find_nth_onebyte_binvec(model->basisBirthable, (I32)model->NumBasisPad16, randLoc);
		int  newtype      = model->basisType[newBasisIdx - 1L];

		BirthFuncs[newtype](NEW, model, newtype, newcol, info);					   
		NEW->jumpType = jumpType;
		return;
	}

	int Kcur = model->curr.K;

	if (jumpType == MOVE) {
		int trytimes  = 0;
		int isMovalbe = 0;
		I32 ksel;
		// Find a term whose basis can generate a new term
		while (trytimes++ < 10L && isMovalbe==0 ) {
			ksel      = RANDINT(model->Kfixed + 1L, Kcur, *(PRND->rnd16)++); // term 1 is fixed at CONST 
			isMovalbe = CheckMoveFuncs[model->terms[ksel-1].type](model, &model->terms[ksel - 1] , info);
		}

		if (isMovalbe) {				
			int oldtype = model->terms[ksel - 1].type;
			// the old term value is used in MoveBirth	
			int status=MoveFuncs[oldtype](NEW, model, oldtype, newcol, ksel, info);
			if (status > 0) {				
				NEW->jumpType  = MOVE;
				return;
			} 			
		} 

		// If ismoveable==0 or (isMoveable=1 && result==0)			  
		jumpType = DEATH;	 
	}


	if (jumpType == DEATH ) {
		int kdel;
		int trytimes = 0;		
		int bingo    = 0;
		//(1)  First a nonlinear term or a linear term taht is kilable
		while ( trytimes++ < 10 && bingo==0 ) {
			kdel = RANDINT(model->Kfixed+1L, Kcur, *(PRND->rnd16)++); // term 1 is fixed at CONST 
			if ( model->terms[kdel-1].type != LinearTerm  || IsLinearTermKillable(model, &model->terms[kdel -1], info)) {
				bingo = 1;
			}					
		}

		if ( bingo == 0 ) {
		//(2) If not found a term above, there must be some Hinge and Change terms in the terms
	    // (2a) If having the HingeTerm
			assert(hasBasis(HingeTerm) || hasBasis(ChangeTerm));

			if ( hasBasis(HingeTerm) ) {
				BSTATE_HINGE* b = type2bs(HingeTerm); 
				if (b->K > 0) {					
					int idx = RANDINT(1L, b->K, *(PRND->rnd16)++);
					kdel    = b->Kposition[idx - 1];
					bingo   = 1;
				}
			}
		
	     // (2b)  If having the ChangeTerm	and not finding a valid term above
			if (hasBasis(ChangeTerm) && bingo == 0 ) {
				BSTATE_HINGE* b = type2bs(ChangeTerm);
				if (b->K > 0) {
					int idx = RANDINT(1L, b->K, *(PRND->rnd16)++);
					kdel    = b->Kposition[idx - 1];
					bingo   = 1;
				}
			}

		} 
		assert(bingo == 1);
 
		//Used in the UpdatebasisState to determine the basis type
		NEW->terms[0] = model->terms[kdel- 1];  

		// A term is chosen for delation. TWo choice are possible
				
		I32 idx0_in_cpList;
		I32 Kdel2_nextknot;
		I32 newMergeknot = IsTermMergible(model, NEW->terms, PRND, &idx0_in_cpList, &Kdel2_nextknot);
		if (newMergeknot >0 ) {
		// (1) The term can be merged with its nehigoring knots
		//   Make a MERGE move if the chosen term is mergible
		//   TODO: currenlty, no merging for HingePariTerm
		//   the term is mergiable
			int type = model->terms[kdel - 1].type;
			if (type == PieceTerm || type == SeasonTerm) {
				Merge_PieceSeason(NEW, model, type, newcol, kdel, Kdel2_nextknot, idx0_in_cpList, newMergeknot, info);
				return;
			} 
			else if (type == StairTerm) {
				Merge_Stair(NEW, model, type, newcol, kdel, Kdel2_nextknot, idx0_in_cpList, newMergeknot, info);
				return;
			}
			else {
				NEW->terms[0].knot = newMergeknot;
				NEW->terms[0].side = ___DetermineSide_HingeChange(NEW->terms, info->xinfo, PRND);// *PRND->rnd08++ > 128; //TODO

				int K1 = min(kdel, Kdel2_nextknot);
				int K2 = max(kdel, Kdel2_nextknot);

				// Insert the new term at the lower K position
				newcol->nbands = 2;
				newcol->ks_x[0] = K1;   newcol->ks_x[1] = K2;
				newcol->kterms_x[0] = 1L;	newcol->kterms_x[1] = 1L;

				newcol->ks_xnewterm[0] = 1;    newcol->ks_xnewterm[1] = 2;
				newcol->kterms_xnewterm[0] = 1L;   newcol->kterms_xnewterm[1] = 0L;

				NEW->idx0_in_cpList = idx0_in_cpList;
				NEW->knot_new = newMergeknot;
				NEW->nTerms = 1;                  //Used to indicate no new term is generated
				NEW->jumpType = MERGE;
				return;
			} 

		} 
		else {			
		// (2) Simply delete the term	    
			int type = model->terms[kdel - 1].type;
			if (type == PieceTerm || type == SeasonTerm) {
				// Piecewiset Term
				Death_PieceSeason(NEW, model, type, newcol, kdel, info);
				return;
			} 
			else if (type == StairTerm) {
				Death_Stair(NEW, model, type, newcol, kdel, info);
			}
			else {
				newcol->nbands = 1;
				newcol->ks_x[0] = kdel;
				newcol->kterms_x[0] = 1L;
				newcol->ks_xnewterm[0] = 1L;
				newcol->kterms_xnewterm[0] = 0L;
				NEW->nTerms = 0;                  //Used to indicate no new term is generated
				NEW->knot_old = NEW->terms[0].knot;
				NEW->jumpType = DEATH;
				return;

			}

		} //if (newknot >0) {


	}

}
 

#include "abc_000_warning.h"
