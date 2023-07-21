#include "abc_000_warning.h"  //put before abc_blas_lapack_lib bcz teh header contains static code
#include "abc_mcmc.h"
#include "abc_vec.h"               //for i08_sum_binvec only
#include "abc_mem.h"               //
#include "globalvars.h"    
#include "abc_blas_lapack_lib.h"   //r_ippsSet_8u
#include "bvs_header.h"


// Indices are 1-based
#define InsertNewElem(dst,n,newIdx, newValue, T)                for(int i=(n); i>=(newIdx);i--) {*((T*)dst+i)=*((T*)dst+i-1);} \
                                                                 *((T*)dst+newIdx-1)=newValue;
#define RepeatElem(dst,n,newIdx, T)                             for(int i=(n); i>=(newIdx);i--) *((T*)dst+i)=*((T*)dst+i-1); 
#define ReplaceElem(dst,n,newIdx, newValue, T)                  *((T*)dst+newIdx-1)=newValue;
#define DeleteElem(dst,n,newIdx, T)                             memcpy((T*)dst+newIdx-1,(T*)dst+newIdx+1-1,sizeof(T)*(n-(newIdx)));
#define MergeTwoElemWithNewValue(dst,n,newIdx,newValue, T)      *((T*)dst+newIdx-1)=newValue;\
			                                                     memcpy((T*)dst+(newIdx+1) - 1, (T*)dst + (newIdx +2) - 1, sizeof(T)* (n - (newIdx+2L)+1L));
#define DoNothing(dst,n, T) ;

static void Update_LinearQudratic(BVS_MODEL * model, int bid, int flag, NewPropTerm *  new , BVS_XINFO_PTR xinfo) {

	BSTATE_LINEAR* bs            = &model->basisState[bid];
	I32            vidx0         = new->terms[0].var0; //vars[0] should never to be zero (the constant linear term is never touched)
	I08PTR         goodVarBinVec = bs->goodVarBinVec;
 

	if (flag == BIRTH) {
		goodVarBinVec[vidx0] = 0;
		bs->goodVarNum--;

		// Update basisBirthable
		if (bs->goodVarNum == 0) {
			model->basisBirthable[bid] = 0;
			model->goodBasisBirthableNum--;
		}
	}
	else if (flag == DEATH) {
		// new->term must be pre-load with the deleted term
		goodVarBinVec[vidx0] = 1; 
		bs->goodVarNum++;

		// Update basisBirthable
		if (bs->goodVarNum == 1) {
			model->basisBirthable[bid] = 1;
			model->goodBasisBirthableNum++;
		}
	}
	else if (flag == MOVE) {//move
		// MOVE MUST BE MADE WITHIN THE SAME BAIS TYPE< NOT ACROSS TYPES
		int        Kdel     = new->xcols->ks_x[0];		
		goodVarBinVec[model->terms[Kdel-1].var0] = 1;//delete
		goodVarBinVec[vidx0]                    = 0;//birth
	}

 
}

static void __cvt_knot_to_binvec(I32PTR knotList, I32 knotNum, I08PTR binvec,
	I32 Ntrue, I32 Npad, BSTATE_HINGE* bs, int vidx0) {
	
	memset(binvec,       1, Npad);
	memset(binvec+Ntrue, 0, Npad-Ntrue);
	for (int i = 0; i < knotNum; i++) {
		int knot = knotList[i];
		int s = knot - bs->minSepDistVec[vidx0];
		int e = knot + bs->minSepDistVec[vidx0];
		s = max(s, 1 + bs->sbad);
		e = min(e, Ntrue - bs->ebad);
		memset(binvec + s - 1, 0, e - s + 1);
	}
	memset(binvec, 0, bs->sbad);
	memset(binvec+Ntrue-(bs->ebad)-1, 0, bs->ebad);
	
}


#define i08_zerofill(x, s, e)    memset(x+s-1, 0, e-s+1);
#define i08_onefill(x, s, e)    memset(x+s-1,1L, e-s+1);
#define i08_sumbin(x, s, e)      i08_sum_binvec(x+s-1, e-s+1);

static void Update_StepStairPieceSeasonHingeChange(BVS_MODEL * model, int bid, int flag, NewPropTerm *  new, BVS_XINFO_PTR xinfo ) {
	
 	BSTATE_HINGE*  bs = &model->basisState[bid];
	int type     = bs->type;
	int Npad     = bs->Npad16;
	int vidx0    = new->terms[0].var0; //vars[0] should never to be zero (the constant linear term is never touched)
	I32PTR  Pindices0      = xinfo->xi_Pindices0[bs->type];
	I32     N              = xinfo->Nunique[Pindices0[vidx0]];
	I08PTR  goodKnotBinVec = bs->goodKnotBinMat + vidx0*Npad;
	I32PTR  goodKnotNumVec = bs->goodKnotNumVec;
	I32     minsep         = bs->minSepDistVec[vidx0];
	I32PTR  cptList = bs->cptMat + vidx0 * bs->MAXCPTNUM;
	I32     cptNum  = bs->cptNumVec[vidx0];

	if (flag == BIRTH) {
 
	 	if (cptNum >= 0) {
			int knot = new->knot_new;
			int s = knot - minsep;  s = max(s, 1 + bs->sbad);
			int e = knot + minsep;  e = min(e, N - bs->ebad);				
			int nGoodKnotNumInSeg = i08_sumbin(goodKnotBinVec, s,e);
			i08_zerofill(goodKnotBinVec,  s, e);
			goodKnotNumVec[vidx0] -= nGoodKnotNumInSeg;
			/***********************************/
			// Add the knot to knotList;
			/***********************************/		
			int newIdx;
			for (newIdx = 0; newIdx < cptNum && cptList[newIdx] < knot; ++newIdx);
			InsertNewElem(cptList, cptNum, newIdx + 1L, knot, I32);
			bs->cptNumVec[vidx0]++;
		}	else {
			// this branch is for PieceTerm
			bs->cptNumVec[vidx0]=0;
		}

		/* Update goodVarNum*/
		if (goodKnotNumVec[vidx0] == 0 || bs->cptNumVec[vidx0] == bs->maxCptNumVec[vidx0]) {
			bs->goodVarBinVec[vidx0] = 0;
			bs->goodVarNum--;
		}

		/* Update basisBirthable*/
		if (bs->goodVarNum == 0) {
			model->basisBirthable[bid] = 0;
			model->goodBasisBirthableNum--;
		}

	}
	else if (flag == DEATH) {//death	

 		if (cptNum > 0) 	{	
			int  knot = new->knot_old;
			int  oldKnotIdx;
			for (oldKnotIdx = 0; oldKnotIdx < cptNum && cptList[oldKnotIdx] != knot; oldKnotIdx++);
			assert(oldKnotIdx < cptNum);
			int LeftEnd  = oldKnotIdx == 0          ? 1 + bs->sbad : cptList[oldKnotIdx - 1] + minsep + 1L;
			int RightEnd = oldKnotIdx == cptNum - 1 ? N - bs->ebad : cptList[oldKnotIdx + 1] - minsep - 1L;			
			int s = max(knot - minsep, LeftEnd);
			int e = min(knot + minsep, RightEnd);

			int segLen = e - s + 1;
			int nBadKnotsNumInSeg = segLen - i08_sumbin(goodKnotBinVec, s, e);
			i08_onefill(goodKnotBinVec, s, e);
			goodKnotNumVec[vidx0] = goodKnotNumVec[vidx0] + nBadKnotsNumInSeg;

			/***********************/
			//Remove the knot from knotList;		
			/***********************/
			DeleteElem(cptList, cptNum, oldKnotIdx + 1L, I32);
			bs->cptNumVec[vidx0]--;
		} else {
			// this branch is for PieceTerm
			bs->cptNumVec[vidx0] = -1;
		}
		
		/* Update goodVarNum*/
		if (bs->goodVarBinVec[vidx0] == 0 && goodKnotNumVec[vidx0] > 0 && bs->cptNumVec[vidx0] < bs->maxCptNumVec[vidx0]) {
			bs->goodVarBinVec[vidx0] = 1L;
			bs->goodVarNum++;
		}
 
		/* Update basisBirthable*/
		if (model->basisBirthable[bid] == 0 && bs->goodVarNum == 1) {  //goodVarNum may be always 1 if there is only 1 covariate
			model->basisBirthable[bid] = 1;
			model->goodBasisBirthableNum++;
		}

	}
	else if (flag == MERGE) {//death	

		int oldknot1 = cptList[new->idx0_in_cpList];
		int oldknot2 = cptList[new->idx0_in_cpList + 1L];

		int idx0_k1 = new->idx0_in_cpList, idx0_k2 = new->idx0_in_cpList + 1L;
		int LeftEnd  = idx0_k1 == 0          ? 1 + bs->sbad : cptList[idx0_k1 - 1] + minsep + 1L;
		int RightEnd = idx0_k2 == cptNum - 1 ? N - bs->ebad : cptList[idx0_k2 + 1] - minsep - 1L;
 
		int s = max(oldknot1 - minsep, LeftEnd);
		int e = min(oldknot2 + minsep, RightEnd);
		int nGoodKnotNumInSeg = i08_sumbin(goodKnotBinVec,s,e); 
		i08_onefill(goodKnotBinVec, s, e);

		int knotnew = new->knot_new;
		int sNew = max(knotnew - minsep, s);
		int eNew = min(knotnew + minsep, e);
		i08_zerofill(goodKnotBinVec, sNew, eNew);
 
		I32 nGoodInSeg_new = (e - s) - (eNew - sNew);
		goodKnotNumVec[vidx0] += (nGoodInSeg_new - nGoodKnotNumInSeg);

		assert(nGoodInSeg_new - nGoodKnotNumInSeg >= 0);

		/***********************/
		//Remove the knot from knotList;		
		/***********************/
		cptList[idx0_k1] = knotnew;
		DeleteElem(cptList, cptNum, idx0_k2 + 1L, I32);
		bs->cptNumVec[vidx0]--;

		/* Update goodVarNum*/
		if (bs->goodVarBinVec[vidx0] == 0 && goodKnotNumVec[vidx0] > 0 && bs->cptNumVec[vidx0] < bs->maxCptNumVec[vidx0]) {
			bs->goodVarBinVec[vidx0] = 1L;
			bs->goodVarNum++;
		}

		/* Update basisBirthable*/
		if (model->basisBirthable[bid] == 0 && bs->goodVarNum == 1) {  //goodVarNum may be always 1 if there is only 1 covariate
			model->basisBirthable[bid] = 1;
			model->goodBasisBirthableNum++;
		}

	}
	else if (flag == MOVE) {//move
		// MOVE MUST BE MADE WITHIN THE SAME BAIS TYPE< NOT ACROSS TYPES
		
		int Knotold = new->knot_old;
		int Knotnew = new->knot_new;

		// (1) If only the sides is flpped, do nothing
		if (Knotold == Knotnew) {
			return;
		}

		// (2) Otherwise, replace the old knot

		int sOld = Knotold - minsep;
		int eOld = Knotold + minsep;

		int sNew = Knotnew - minsep;
		int eNew = Knotnew + minsep;
		
		// Not touchinb boundaries
		/* This doesnot work
		if (sOld >= 1 + bs->sbad && eOld <= N - bs->ebad && sNew >= 1 + bs->sbad && eNew <= N - bs->ebad) {
			// no chanbge to goodKnotNum			
			memset(goodKnotBinVec+sOld-1, 1L, eOld - sOld + 1);
			memset(goodKnotBinVec+sNew-1, 0L, eNew - sNew + 1);
			goodKnotNumVec[vidx0] = goodKnotNumVec[vidx0]; 
		}
	*/

		int oldIdx;	for (oldIdx = 0; oldIdx < cptNum && cptList[oldIdx] != Knotold; oldIdx++);
		int LeftEnd = (oldIdx == 0)           ? 1 + bs->sbad : cptList[oldIdx - 1] + minsep + 1L;
		int RightEnd = (oldIdx == cptNum - 1)  ? N - bs->ebad : cptList[oldIdx + 1] - minsep - 1L;

		sOld = max(sOld, LeftEnd);
		eOld = min(eOld, RightEnd);

		sNew = max(sNew, LeftEnd);
		eNew = min(eNew, RightEnd);

		int s = min(sOld, sNew);
		int e = max(eOld, eNew);

		int oldGoodNum = i08_sumbin(goodKnotBinVec, s, e);						
		i08_onefill(goodKnotBinVec,  sOld, eOld);
		i08_zerofill(goodKnotBinVec, sNew, eNew);		
		int newGoodNum = i08_sumbin(goodKnotBinVec,s,e);
		goodKnotNumVec[vidx0] += newGoodNum - oldGoodNum;
	 
		// Replace the knot from knotList;	
		ReplaceElem(cptList, knotNum, oldIdx + 1L, Knotnew, I32);

		bs->cptNumVec[vidx0]= bs->cptNumVec[vidx0];

		/* Update goodVarNum*/
		if (bs->goodVarBinVec[vidx0] == 0 && goodKnotNumVec[vidx0] > 0 && bs->cptNumVec[vidx0] < bs->maxCptNumVec[vidx0]) {
			bs->goodVarBinVec[vidx0] = 1L;
			bs->goodVarNum++;
		}
		else if (bs->goodVarBinVec[vidx0] == 1 && (goodKnotNumVec[vidx0] == 0 || bs->cptNumVec[vidx0] == bs->maxCptNumVec[vidx0])) {
			bs->goodVarBinVec[vidx0] = 0L;
			bs->goodVarNum--;
		}

		/* Update basisBirthable*/
		if (model->basisBirthable[bid] == 0 && bs->goodVarNum == 1) {
			//goodVarNum may be always 1 if there is only 1 covariate
			model->basisBirthable[bid] = 1;
			model->goodBasisBirthableNum++;
		}
		else if (model->basisBirthable[bid] == 1L && bs->goodVarNum == 0) {
			model->basisBirthable[bid] = 0;
			model->goodBasisBirthableNum--;
		}
	 
	}


	/*
	int bb=i08_sum_binvec(goodKnotBinVec, Npad);
	if (bb != goodKnotNumVec[vidx0]) {
		int ttt = 1;
	}


	for (int i = 0; i < bs->cptNumVec[vidx0]; i++) {
		if (cptList[i] < 1) {
			int aa = 111;
		}
	}
	
	if (bs->type == 2) {	
		I08 binvec[10000];
		__cvt_knot_to_binvec(knotList, bs->knotNumVec[vidx0], binvec, N, Npad, bs , vidx0 );
		if (memcmp(binvec, goodKnotBinVec, Npad)!=0) {
			char * errormsg="Error Occured!"
		}
	}
	*/
}


 
static int __IsKnotTwoSided_HingePar(BVS_TERMS_PTR term, BSTATE_HINGEPAIR* bs, BVS_MODEL_PTR model) {
	int vidx0  = term->var0;
	int knot   = term->knot;
	int sides = 1 - term->side;

	BVS_TERMS* terms = model->terms;
	for (int i = 0; i < bs->K; i++) {
		int k = bs->Kposition[i];
		if (terms[k - 1].var0 == vidx0 && terms[k - 1].knot == knot && terms[k - 1].side == sides) {
			return i + 1;
		}
	}
	return 0;
}

static void Update_HingePair(BVS_MODEL * model, int bid, int flag, NewPropTerm *  new, BVS_XINFO_PTR xinfo ) {
	
 	BSTATE_HINGE*  bs = &model->basisState[bid];
	
	int Npad     = bs->Npad16;
	int vidx0    = new->terms[0].var0; //vars[0] should never to be zero (the constant linear term is never touched)
	int knot     = new->terms[0].knot;

	I32PTR Pindices0 =  xinfo->xi_Pindices0[bs->type];
	int    N         =  xinfo->Nunique[Pindices0[vidx0]];

	I08PTR  goodKnotBinVec = bs->goodKnotBinMat + vidx0*Npad;
	I32PTR  goodKnotNumVec = bs->goodKnotNumVec;
	I32PTR  cptList        = bs->cptMat + vidx0 * bs->MAXCPTNUM;

	int minsep = bs->minSepDistVec[vidx0];
	if (flag == BIRTH) {
	//THIS BRANCH IS THE SAME AS Update_StepHingeChange
		int s = knot - minsep; s = max(s, 1 + bs->sbad);
		int e = knot + minsep; e = min(e, N - bs->ebad);

		int nGoodKnotNumInSeg = i08_sum_binvec(goodKnotBinVec + s - 1, e - s + 1L);
		memset(goodKnotBinVec + s - 1, 0L, e - s + 1);
		goodKnotNumVec[vidx0] = goodKnotNumVec[vidx0] - nGoodKnotNumInSeg;


		/***********************************/
		//Add the knot to knotList;
		/***********************************/
		int cptNum = bs->cptNumVec[vidx0];
		int newIdx;
		for (newIdx = 0; newIdx < cptNum && knot>cptList[newIdx]; newIdx++);
		InsertNewElem(cptList, cptNum, newIdx + 1L, knot, I32);
		bs->cptNumVec[vidx0]++;

		/* Update goodVarNum*/
		if (goodKnotNumVec[vidx0] == 0 || bs->cptNumVec[vidx0] == bs->maxCptNumVec[vidx0]) {
			bs->goodVarBinVec[vidx0] = 0;
			bs->goodVarNum--;
		}

		/* Update basisBirthable*/
		if (bs->goodVarNum == 0) {
			model->basisBirthable[bid] = 0;
			model->goodBasisBirthableNum--;
		}
	}
	else if (flag == DEATH) {//death	

		if (__IsKnotTwoSided_HingePar(new->terms, bs, model)) {
			//Do NOthing
			return;
		}

		int s, e, cptNum, oldKnotIdx;
		{
			
			cptNum = bs->cptNumVec[vidx0];
			for (oldKnotIdx = 0; oldKnotIdx < cptNum && cptList[oldKnotIdx] != knot; oldKnotIdx++);
			assert(oldKnotIdx < cptNum);
			int LeftEnd  = oldKnotIdx == 0          ? 1 + bs->sbad : cptList[oldKnotIdx - 1] + minsep + 1L;
			int RightEnd = oldKnotIdx == cptNum - 1 ? N - bs->ebad : cptList[oldKnotIdx + 1] - minsep - 1L;

			s = knot - minsep;
			e = knot + minsep;
			s = max(s, LeftEnd);
			e = min(e, RightEnd);		
		}
		
		int segLen = e - s + 1;
		int nBadKnotsNumInSeg = segLen - i08_sum_binvec(goodKnotBinVec + s - 1, segLen);
		memset(goodKnotBinVec + s - 1, 1L, segLen);
		goodKnotNumVec[vidx0] = goodKnotNumVec[vidx0] + nBadKnotsNumInSeg;
 
		/***********************/
		//Remove the knot from knotList;		
		/***********************/
		DeleteElem(cptList, cptNum, oldKnotIdx + 1L, I32);
		bs->cptNumVec[vidx0]--;

		/* Update goodVarNum*/
		if (bs->goodVarBinVec[vidx0] == 0 && goodKnotNumVec[vidx0] > 0 && bs->cptNumVec[vidx0] < bs->maxCptNumVec[vidx0]) {
			bs->goodVarBinVec[vidx0] = 1L;
			bs->goodVarNum++;
		}

		/* Update basisBirthable*/
		if (model->basisBirthable[bid] == 0 && bs->goodVarNum == 1) {  
			//goodVarNum may be always 1 if there is only 1 covariate
			model->basisBirthable[bid] = 1;
			model->goodBasisBirthableNum++;
		}
	}
	else if (flag == MERGE) {//death	
    // NO MREGE FOR HingePairTerm
 
	}
	else if (flag == MOVE) {//move
	  // MOVE MUST BE MADE WITHIN THE SAME BAIS TYPE< NOT ACROSS TYPES
	  // termsdel->vars[0] EQQUAL pidx
		int Kdel            = new->xcols->ks_x[0];
		BVS_TERMS* termsdel = &model->terms[Kdel - 1];

		I32 oldKnot = termsdel->knot;

		// (1) If only the sides is flpped, do nothing
		if (knot == oldKnot) {
			return;
		}

		// (2) Otherwise, replace the old knot

		int sOld = oldKnot - minsep;
		int eOld = oldKnot + minsep;

		int sNew = knot - minsep;
		int eNew = knot + minsep;
		
		// Not touchinb boundaries
		/* This doesnot work
		if (sOld >= 1 + bs->sbad && eOld <= N - bs->ebad && sNew >= 1 + bs->sbad && eNew <= N - bs->ebad) {
			// no chanbge to goodKnotNum			
			memset(goodKnotBinVec+sOld-1, 1L, eOld - sOld + 1);
			memset(goodKnotBinVec+sNew-1, 0L, eNew - sNew + 1);
			goodKnotNumVec[vidx0] = goodKnotNumVec[vidx0];
 
		}
		*/

		int cptNum = bs->cptNumVec[vidx0];
		int oldIdx;	for (oldIdx = 0; oldIdx < cptNum && cptList[oldIdx] != oldKnot; oldIdx++);
		int LeftEnd = (oldIdx == 0)           ? 1 + bs->sbad  : cptList[oldIdx - 1] + minsep + 1L;
		int RightEnd = (oldIdx == cptNum - 1)  ? N - bs->ebad : cptList[oldIdx + 1] - minsep - 1L;

		sOld = max(sOld, LeftEnd);
		eOld = min(eOld, RightEnd);

		sNew = max(sNew, LeftEnd);
		eNew = min(eNew, RightEnd);

		int s = min(sOld, sNew);
		int e = max(eOld, eNew);

		int oldGoodNum = i08_sum_binvec( goodKnotBinVec+s-1, e-s+1L);				
		memset(goodKnotBinVec + sOld - 1, 1L, eOld - sOld + 1);
		memset(goodKnotBinVec + sNew - 1, 0L, eNew - sNew + 1);
		int newGoodNum = i08_sum_binvec(goodKnotBinVec+s - 1, e-s+1L);		
		goodKnotNumVec[vidx0] += newGoodNum - oldGoodNum;
	 
		// Replace the knot from knotList;	
		ReplaceElem(cptList, knotNum, oldIdx + 1L, knot, I32);

		bs->cptNumVec[vidx0]= bs->cptNumVec[vidx0];

		/* Update goodVarNum*/
		if (bs->goodVarBinVec[vidx0] == 0 && goodKnotNumVec[vidx0] > 0 && bs->cptNumVec[vidx0] < bs->maxCptNumVec[vidx0]) {
			bs->goodVarBinVec[vidx0] = 1L;
			bs->goodVarNum++;
		}
		else if (bs->goodVarBinVec[vidx0] == 1 && (goodKnotNumVec[vidx0] == 0 || bs->cptNumVec[vidx0] == bs->maxCptNumVec[vidx0])) {
			bs->goodVarBinVec[vidx0] = 0L;
			bs->goodVarNum--;
		}

		/* Update basisBirthable*/
		if (model->basisBirthable[bid] == 0 && bs->goodVarNum == 1) {  //goodVarNum may be always 1 if there is only 1 covariate
			model->basisBirthable[bid] = 1;
			model->goodBasisBirthableNum++;
		}
		else if (model->basisBirthable[bid] == 1L && bs->goodVarNum == 0) {
			model->basisBirthable[bid] = 0;
			model->goodBasisBirthableNum--;
		}
	}
	 
	/*
	if (bs->type == 2) {	
		I08 binvec[10000];
		__cvt_knot_binvec(knotList, bs->knotNumVec[vidx0], binvec, N, Npad, bs);
		if (memcmp(binvec, goodKnotBinVec, Npad)!=0) {
			char * errormsg="Error Occured!"
		}
	}
	*/
}
 
static void Update_Hinge2d_Hinge2dr_Hinge1dr(BVS_MODEL* model, int bid, int flag, NewPropTerm *  new, BVS_XINFO_PTR xinfo) {
	return;
}

typedef  void (*pfuncUpdate)(BVS_MODEL* model, int bid, int flag, NewPropTerm *  new, BVS_XINFO_PTR xinfo);

pfuncUpdate UpdateFuncs[] = {
	Update_LinearQudratic,
	Update_LinearQudratic,
	Update_StepStairPieceSeasonHingeChange, //step
	Update_StepStairPieceSeasonHingeChange, //stair
	Update_StepStairPieceSeasonHingeChange, //piece
	Update_StepStairPieceSeasonHingeChange, //season
	Update_StepStairPieceSeasonHingeChange, //hinge
	Update_HingePair,
	Update_StepStairPieceSeasonHingeChange,
	Update_Hinge2d_Hinge2dr_Hinge1dr,
	Update_Hinge2d_Hinge2dr_Hinge1dr,
	Update_Hinge2d_Hinge2dr_Hinge1dr,
};

 

 void Update_BasisState_All(NewPropTerm * new, BVS_MODEL *model, BVS_XINFO_PTR xinfo) {
	 
	int type = new->terms[0].type;
	int flag = new->jumpType;
	I32 bid  = model->type2id[type];
		
	((MemPointers*)GLOBAL_TMP)->verify_header(GLOBAL_TMP);

    /*********************************************/
	// Update basisState: use the old model.terms inside it
	/*********************************************/
	UpdateFuncs[type](model,bid,flag, new, xinfo);

	((MemPointers*)GLOBAL_TMP)->verify_header(GLOBAL_TMP);
	
	/*********************************************/
	// Update  model.terms
	/*********************************************/
	NEWCOLINFOv2* NewCol = new->xcols;
	swap_elem_bands(NewCol, model->terms, new->terms,sizeof(BVS_TERMS));
	model->Kterms = NewCol->Knew;

	((MemPointers*)GLOBAL_TMP)->verify_header(GLOBAL_TMP);

	// No need to update Kpositons at all
	if (NewCol->isEqualSwap) {
		return;
	}
	

	// NewCol->isEqualSwap == 0

	/*********************************************/
	// Update Kpositions: and basis->K
	/*********************************************/
	
	BSTATE_PIECE* b= &model->basisState[bid]; 

	// First the first swap band that is incosisent between X and Xnewterm
	int Kunequal0 =0;
	for (; Kunequal0 < NewCol->nbands && NewCol->kterms_x[Kunequal0] == NewCol->kterms_xnewterm[Kunequal0]; 
		  Kunequal0++	);	

	// if the unequal band ends at the last col of X, then no need to update other bases
	int UpdateOtherbasis   = NewCol->ks_x[Kunequal0] + NewCol->kterms_x[Kunequal0] - 1 == NewCol->K ? 0 : 1L;
	int UpdateCurrentBasis = 1;		

	((MemPointers*)GLOBAL_TMP)->verify_header(GLOBAL_TMP);
	// No need to update other basese and only update the current basis	
	if (UpdateOtherbasis == 0) {		
		if (NewCol->ks_x[Kunequal0] == NewCol->K + 1) {
			// If inserting after the X, no need to update other bases
			// Update the curent basis by simply addding the extra cols
			int    Kextra     = NewCol->Knew - NewCol->K;
			I16PTR Kposistion = b->Kposition;
			for (int i = 1; i <= Kextra; i++) {
				Kposistion[b->K++] = i+NewCol->K;
			}			
			UpdateCurrentBasis = 0;

			return;
		} else {
			//NewCol->ks_x[Kunequal0] <= K		
		}
	}

	((MemPointers*)GLOBAL_TMP)->verify_header(GLOBAL_TMP);

	// Must update the other bases
	if (UpdateOtherbasis == 1) {		
		
		// Use the current basis's Kposition as a temp buffer
		// This is buggy because some bases's Kpostion may be not long enough
		// I16PTR Ktmp = model->basisState[bid].Kposition;
		I16PTR Ktmp = new->KpositionBuf;

		// if a col index is >=Kthreahold, it need to be updated
		int Kthreshold = NewCol->ks_x[Kunequal0] + NewCol->kterms_x[Kunequal0]; 
		 

		for (int i = Kunequal0; i < NewCol->nbands; i++) {
			// i=0(x1, new1) i=1(x2, new2),...., i=bands-1(x, new), x_end
			// when i=0, take the x2 part
			int Kpart = 2 * i + 2;
			int Ksrc = NewCol->parts[Kpart].ks_src;
			int Kdst = NewCol->parts[Kpart].ks_dst;
			int Kterms = NewCol->parts[Kpart].kterms;
	 
			for (int j = Ksrc; j<= (Ksrc + Kterms - 1); j++) {
				Ktmp[j - 1] = Kdst++; 
			}
	 
		}

		((MemPointers*)GLOBAL_TMP)->verify_header(GLOBAL_TMP);

		// Update Kpositions for all the other bases
		for (int j = 0; j < model->NUMBASIS; j++) {
			if (model->basisType[j] == type) {
				continue; 
			}

			BSTATE_LINEAR* bstmp      = &model->basisState[j]; // only the first two feilds needed
			I16PTR         bKposition = model->basisState[j].Kposition;
			for (int i = 0; i < bstmp->K; i++) {
				int Kcur = bKposition[i];
				if (Kcur < Kthreshold) {
					continue;
				}
				bKposition[i] = Ktmp[Kcur - 1];
			}
			bstmp->K; // no not change it
		}
		((MemPointers*)GLOBAL_TMP)->verify_header(GLOBAL_TMP);

		UpdateCurrentBasis = 1;
	}

	((MemPointers*)GLOBAL_TMP)->verify_header(GLOBAL_TMP);

	if (UpdateCurrentBasis) {
		int         K = 0;
		for (int i = 0; i < NewCol->Knew; i++) {
			if (model->terms[i].type == type) {
				b->Kposition[K++] = i + 1;			
			}		
		}
		b->K = K;
	}
	((MemPointers*)GLOBAL_TMP)->verify_header(GLOBAL_TMP);
	
}
 
/*
* /////////////////
	for (int i = 0; i < b->P; i++) {
		I32PTR clist = b->cptMat + i * b->MAXCPTNUM;
		for (int j = 0; j < b->cptNumVec[i]; j++) {
			I32 knot = clist[j];
			int bingpo = 0;
			for (int k = 0; k < model->Kterms; k++) {
				if (model->terms[k].knotmid == knot && model->terms[k].vars[0] == i + 1) {
					bingpo = 1;
					break;
				}
			}
			assert(bingpo == 1);
		}

	}

	for (int k = 0; k < model->Kterms; k++) {
		if (model->terms[k].vars[0] != 3) {
			continue;
		}

		int bing = 0;
		int knot = model->terms[k].knotmid;
		if (knot == 1) continue;
		I32PTR clist = b->cptMat + 2 * b->MAXCPTNUM;
		for (int j = 0; j < b->cptNumVec[2]; j++) {
			if (knot == clist[j]) {
				bing = 1;
				break;
			}
		}
		assert(bing == 1);
	}

*/
 
#include "abc_000_warning.h"