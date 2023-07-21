#include "abc_000_warning.h"

#include "abc_001_config.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "abc_datatype.h"
#include "abc_blas_lapack_lib.h"
#include "abc_ide_util.h"  //printf
#include "abc_common.h"  //strcicmp
#include "abc_ts_func.h"
#include "abc_date.h"
#include "bvs_io.h"
#include "bvs_fun.h"

#define IsAlmostInteger(x)  ( fabs(x-round(x)) <1e-3 )

static int  GetArg_12_MetaData(VOIDPTR prhs[], int nrhs, BVS_METADATA_PTR meta) {

	// Check the number of inputs:  the first arg is the algorithm name.
	if (nrhs < 2L) {	
		r_error("ERROR: At least one input argument is needed!\n");
		return 0;
	}

	meta->nrhs = nrhs;
	meta->Yobj = prhs[1];
	meta->Xobj = prhs[2];

	if (!IsNumeric(meta->Xobj)  || !IsNumeric(meta->Yobj) ) {
		r_error("ERROR: metadata$time is not numeric. If times are strings, please use metadata$time$dateStr and metadata$time$strFmt to specify data observation times.\n");
		return 0;
	}
	/*
	 !(IsDouble(DATA) && numel > 2) && 	!(IsSingle(DATA) && numel > 2) &&
	   !(IsInt32(DATA)  && numel > 2) &&	!(IsInt16(DATA) && numel > 2) &&
	   !(IsStruct(DATA) && numel >= 1 ) &&  !(IsCell(DATA) )  
	*/
	meta->XobjDtype = GetDataType(meta->Xobj);
	meta->YobjDtype = GetDataType(meta->Yobj);

	/*************************************************/
	/* Get the dimension of the input              */
	/*************************************************/ 
	I32       ndims = GetNumOfDim(meta->Xobj);
	if (ndims == 0 || ndims == 1) {
		// the input is a vector: this branch is possible only for R vectors (not  Matlab)
		meta->ndimx = 1L;
		meta->Nobj = GetNumberOfElements(meta->Xobj);
		meta->Pobj = 1L; 
	} //ndims is impossible to be 1L
	else if (ndims == 2) {
		// Matlab: a vector or matrix;
		// R:      a matrix or a vector with a dim attribute.
		int N = GetDim1(meta->Xobj);
		int K = GetDim2(meta->Xobj);
		meta->ndimx = 2L;//  //PY is a matrix		

		if (min(N, K) == 1L) { //PY is a vector		
			N = max(N, K);	
			K = 1L;
			meta->ndimx = 1L;
		}

		meta->Nobj  = N;
		meta->Pobj  = K;
		//GetDimensions(Y, io->dims, 3L);
	}

	//BVS_fetch_next_timeSeries(&yInfo, i, MEMBUF, io);	
 
	//CopyStrideMEMToF32Arr(Y + i * N/*dst*/, io->pdata[i] /*src*/, N, stride, offset, io->dtype[i]);
	//f32_set_nan_by_value(Y + i * N, N, io->meta.missingValue);
 
	//	int  nMissing = f32_find_nans(Y, N, badRowIdx);
	//f32_mat_multirows_extract_set_by_scalar(Y, N, 1L, Ycopy, badRowIdx, nMissing, 0);

	//f32_mat_multirows_extract_set_by_scalar(Xtmp, N, K + 1, Xcopy, badRowIdx, nMissing, 0);
	//	linear_regression(Y, X, N, N, K, B, Yfit, NULL, XtX);
	//	f32_mat_multirows_set_by_submat(Xtmp, N, K + 1, Xcopy, badRowIdx, nMissing);
	//
	return 1;
}

static int i32_insert_noduplicate(I32PTR x, I32 N, I32PTR Xnew, I32 Nnew) {

 	for (int i = 0; i < Nnew; i++) {
		int newvalue = Xnew[i];
		int matched  = _False_;
		for (int j = 0; j < N; j++) {
			if (x[j] == newvalue) {
				matched = _True_;
				break;
			}
		}
		if (!matched) {
			x[N++] = newvalue;			
		}
	}

	return N;
}
 
static int i32_unique_inplace(I32PTR x, int N) {
	int Nuniq = 0;	
	int next  = 0;
	while (next < N) {
		int xcur = x[next];
		// first, move Inext foward to check if consectuve elements are the same
		while (next < (N-1) && x[next+1] == xcur) {
			next++;
		}

		int matched = 0;
		for (int i = 0; i < Nuniq; ++i) {
			if (x[i] == xcur) {
				matched = 1;
				break;
			}
		}
		if (!matched) {
			x[Nuniq++] = xcur;
		}
		next++;
	}

	return Nuniq;	 
}

static int i32_exclude_inplace(I32PTR x, int N, I32PTR excludeList, I32 Nexclude) {
	
	if (x == NULL || excludeList == NULL) return N;

	for (int i = 0; i < Nexclude && N>0; i++) {
		int xcur = excludeList[i];
		for (int j = 0; j < N; j++) {
			if (x[j] == xcur) {
			// find a match
				x[j] = x[N - 1];
				N--;
				break;
			}
		}
	}
	return N;
}

static void ParseArg_basisType(VOIDPTR S, BVS_PRIOR_PTR prior, BVS_METADATA_PTR meta) {

	VOIDPTR tmpField = GetField123Check(S, "basisType", 5);
	int     nBasis   = tmpField && IsChar(tmpField) ? GetNumberOfElements(tmpField) : 0;

	int basisIncluded[MAX_NUM_BASIS] = { 0 };
	int nGoodBasis                   = 0;
	prior->isBeast = 0;
	prior->isFFT   = 0;
	for (int i = 0; i < nBasis; i++) {
		char sbuf[100 + 1];
		int  len = GetCharVecElem(tmpField, i, sbuf, 100);
		if (len == 0) {
			continue;
		}

		if      (strcicmp_nfirst(sbuf, "beast", 3) == 0)
			prior->isBeast = 1, nGoodBasis++;
		else if (strcicmp_nfirst(sbuf, "fft", 3) == 0)
			prior->isFFT = 1,   nGoodBasis++;
		else if (strcicmp_nfirst(sbuf, "linear", 3) == 0)
			basisIncluded[LinearTerm] = 1, nGoodBasis++;
		else if (strcicmp_nfirst(sbuf, "quadratic", 3) == 0 ||
			     strcicmp_nfirst(sbuf, "square",    3) == 0 ) 
			basisIncluded[QuadraticTerm] = 1, nGoodBasis++;
		else if (strcicmp_nfirst(sbuf, "step", 3) == 0)
			basisIncluded[StepTerm] = 1, nGoodBasis++;
		else if (strcicmp_nfirst(sbuf, "stair", 3) == 0)
			basisIncluded[StairTerm] = 1, nGoodBasis++;
		else if (strcicmp_nfirst(sbuf, "piecewise", 3) == 0)
			basisIncluded[PieceTerm] = 1, nGoodBasis++;
		else if (strcicmp_nfirst(sbuf, "season", 3) == 0)
			basisIncluded[SeasonTerm] = 1, nGoodBasis++;
		else if (strcicmp(sbuf, "hinge") == 0)
			basisIncluded[HingeTerm] = 1, nGoodBasis++;
		else if (strcicmp(sbuf, "hingepair") == 0)
			basisIncluded[HingePairTerm] = 1, nGoodBasis++;
		else if (strcicmp_nfirst(sbuf, "change", 3) == 0)
			basisIncluded[ChangeTerm] = 1, nGoodBasis++;
		else if (strcicmp(sbuf, "hinge2d") == 0)
			basisIncluded[Hinge2DTerm] = 1, nGoodBasis++;
		else if (strcicmp(sbuf, "hinge2dr") == 0)
			basisIncluded[Hinge2DRTerm] = 1, nGoodBasis++;
		else if (strcicmp(sbuf, "hinge1dr") == 0)
			basisIncluded[Hinge1DRTerm] = 1, nGoodBasis++;
	}		 
	if (nGoodBasis == 0  ) {
		basisIncluded[LinearTerm] = 1;
	}

    #define pr (*prior)

	prior->beast_period         = (tmpField = GetField123Check(S, "period", 1))         ? GetScalar(tmpField) : getNaN();
	prior->beast_maxSeasonOrder = (tmpField = GetField123Check(S, "maxSeasonOrder", 0)) ? GetScalar(tmpField) : getNaN();
	prior->fft_maxOrder         = (tmpField = GetField123Check(S, "maxFFTOrder", 0))    ? GetScalar(tmpField) : getNaN();
	if (pr.beast_period <= 0)								             pr.beast_period = getNaN();
	if (pr.beast_maxSeasonOrder <= 0 || IsNaN(pr.beast_maxSeasonOrder))  pr.beast_maxSeasonOrder = 4;
	if (pr.fft_maxOrder <= 0         || IsNaN(pr.fft_maxOrder))          pr.fft_maxOrder = (int)meta->Nobj/2;
	
	int Pall = meta->Pobj;

    //////////////////////////////////////////////////////
	// Get the BEAST terms t- add
	//////////////////////////////////////////////////////
	pr.beast_Pbase   = Pall + 1;
	pr.beast_Pextra  = pr.beast_Ptime = -9999999999L;
	if (pr.isBeast) {
		VOID_PTR tmpBeast = GetField123Check(S, "beast", 0);
		pr.beast_Ptime = 1L;  // default value, updated if tmpBeast is not NULL
		if (tmpBeast) {
			pr.beast_Ptime = GetScalar(tmpBeast);
		}

		pr.beast_Pextra = 0;// is 0 for trend-only 
		if ( !IsNaN(pr.beast_period) ) {
			pr.beast_Pextra = pr.beast_maxSeasonOrder * 2L;
		}

		Pall = Pall + pr.beast_Pextra;	
	}

	//////////////////////////////////////////////////////
	// Get the  FFT terms t- add
	//////////////////////////////////////////////////////	
	pr.fft_Pbase = Pall + 1;
	pr.fft_Pnew  = pr.fft_Ptime = -1L;
	if (pr.isFFT) {
		VOID_PTR tmpFFT = GetField123Check(S, "fft", 0);
		pr.fft_Ptime = 1L;  // default value, updated if tmpBeast is not NULL
		if (tmpFFT) {
			pr.fft_Ptime = GetScalar(tmpFFT);
		}
		pr.fft_Pnew = pr.fft_maxOrder * 2;
		Pall = Pall + pr.fft_Pnew;
	}



	//////////////////////////////////////////////////////
	// Process individual basis Type
	//////////////////////////////////////////////////////
	memset(prior->Pindices0, 0, sizeof(F32PTR) * MAX_NUM_BASIS);
	memset(prior->Pbasis,    0, sizeof(I32)    * MAX_NUM_BASIS);


	int type;
	nGoodBasis = 0;

	type = LinearTerm;
	if (basisIncluded[type] || prior->isFFT ) {
	
		I32PTR varlist = pr.Pindices0[type] = malloc(Pall * sizeof(I32));
		I32    Pbasis  = 0;
		VOID_PTR tmpobj= GetField123Check(S, "linear", 0);
		if (tmpobj) {
			Pbasis = GetNumberOfElements(tmpobj);
			CopyNumericObjToI32Arr(varlist, tmpobj, Pbasis);
			Pbasis = i32_unique_inplace(varlist, Pbasis);
		}	else {
			Pbasis = meta->Pobj;
			i32_seq(varlist, 1, 1, Pbasis);     // if no linear field, insert all variables
		}
		

	    /**************************************************/
		// Process the fixedLiniarTerms
		/**************************************************/
		tmpobj = GetField123Check(S, "FixedLinearTerms", 0);
		if (tmpobj) {
			int Pfixed = GetNumberOfElements(tmpobj);
			prior->FixedLinearTerms0 = malloc(sizeof(I32) * Pfixed);
			CopyNumericObjToI32Arr(prior->FixedLinearTerms0, tmpobj, Pfixed);
			Pfixed                   = i32_unique_inplace(prior->FixedLinearTerms0, Pfixed);
			prior->nFixedLinearTerms = Pfixed;

		} else {
			prior->FixedLinearTerms0 = NULL;
			prior->nFixedLinearTerms = 0;
		}		
		Pbasis  = i32_insert_noduplicate(varlist, Pbasis, prior->FixedLinearTerms0, prior->nFixedLinearTerms);


		/**************************************************/
		// Add the fft terms, if any
		/**************************************************/
		if (pr.isFFT) {
			Pbasis=i32_insert_noduplicate(varlist, Pbasis, &pr.fft_Ptime, 1L);
			// pr.beast_Pextra>0 means taht there isa  seaosnal component
			for (int i = 0; i < pr.fft_Pnew; i++) {
				varlist[Pbasis++] = pr.fft_Pbase + i;
			}
		}

		pr.Pbasis[type]            = Pbasis;
		pr.basisType[nGoodBasis++] = type;
	}

	type = QuadraticTerm;
	if (basisIncluded[type] ) {
		VOID_PTR tmpobj  = GetField123Check(S, "square", 0);
		I32      Pbasis  = tmpobj ? GetNumberOfElements(tmpobj): meta->Pobj;
		I32PTR   varlist = pr.Pindices0[type] = malloc(Pbasis * sizeof(I32));
		if (tmpobj) {
			CopyNumericObjToI32Arr(varlist, tmpobj, Pbasis);
			Pbasis = i32_unique_inplace(varlist, Pbasis);
		}	else {		
			i32_seq(varlist, 1, 1, Pbasis);     // if no square field, insert all variables
		}		
		pr.Pbasis[type]           = Pbasis;
		pr.basisType[nGoodBasis++] = type;
	}


	type = StepTerm;
	if (basisIncluded[type] ) {
		VOID_PTR tmpobj = GetField123Check(S, "step", 0);
		I32      Pbasis = tmpobj ? GetNumberOfElements(tmpobj) : meta->Pobj;
		I32      Pmax   = Pbasis + pr.isBeast; // add one more term if isBeast=TRUE
		I32PTR   varlist = pr.Pindices0[type] = malloc(Pmax * sizeof(I32));
		if (tmpobj) {
			CopyNumericObjToI32Arr(varlist, tmpobj, Pbasis);
			Pbasis = i32_unique_inplace(varlist, Pbasis);
		}	else {
			i32_seq(varlist, 1, 1, Pbasis);     // if no square field, insert all variables
		}

		if (Pbasis > 0) {
			pr.Pbasis[type]           = Pbasis;
			pr.basisType[nGoodBasis++] = type;
		}	else {			
			free(pr.Pindices0[type]);
			pr.Pindices0[type] = NULL;
			pr.Pbasis[type]    = 0;
		}
	}

	type = StairTerm;
	if (basisIncluded[type] || pr.isBeast) {
		VOID_PTR tmpobj = GetField123Check(S,"stairterm", 3);
		I32      Pbasis = tmpobj ? GetNumberOfElements(tmpobj) : meta->Pobj;
		I32      Pmax   = Pbasis + pr.isBeast; // add one more term if isBeast=TRUE
		I32PTR   varlist = pr.Pindices0[type] = malloc(Pmax * sizeof(I32));
		if (tmpobj) {
			CopyNumericObjToI32Arr(varlist, tmpobj, Pbasis);
			Pbasis = i32_unique_inplace(varlist, Pbasis);
		} else {
			i32_seq(varlist, 1, 1, Pbasis);     // if no square field, insert all variables
		}

		if (pr.isBeast) {
			Pbasis = i32_insert_noduplicate(varlist, Pbasis, &pr.beast_Ptime, 1L);
		}
	
		if (Pbasis > 0) {
			pr.Pbasis[type]           = Pbasis;
			pr.basisType[nGoodBasis++] = type;
		}	else {			
			free(pr.Pindices0[type]);
			pr.Pindices0[type] = NULL;
			pr.Pbasis[type]    = 0;
		}

	}


	type = PieceTerm;
	if (basisIncluded[type] || prior->isBeast ) {
		VOID_PTR tmpobj = GetField123Check(S,"piecewise", 3);
		I32      Pbasis = tmpobj ? GetNumberOfElements(tmpobj) : meta->Pobj;
		I32      Pmax   = Pbasis + pr.isBeast; // add one more term if isBeast=TRUE
		I32PTR   varlist = pr.Pindices0[type] = malloc(Pmax * sizeof(I32));
		if (tmpobj) {
			CopyNumericObjToI32Arr(varlist, tmpobj, Pbasis);
			Pbasis = i32_unique_inplace(varlist, Pbasis);
		} else {
			i32_seq(varlist, 1, 1, Pbasis);     // if no square field, insert all variables
		}

		/**************************************************/
		// Add the beast terms, if any
		/**************************************************/
		if (pr.isBeast) {
			Pbasis = i32_insert_noduplicate(varlist, Pbasis, &pr.beast_Ptime, 1L);
			// pr.beast_Pextra>0 means taht there is a  seaosnal component
			//for (int i = 0; i < pr.beast_Pextra; i++) {
			//	varlist[Pbasis++] = pr.beast_Pbase + i;
			//}
		}

 		// Remove the lienar terms from the PieceWist terms
		Pbasis = i32_exclude_inplace(varlist, Pbasis, pr.Pindices0[LinearTerm], pr.Pbasis[LinearTerm]);
		if (Pbasis > 0) {
			pr.Pbasis[type]           = Pbasis;
			pr.basisType[nGoodBasis++] = type;
		}	else {			
			free(pr.Pindices0[type]);
			pr.Pindices0[type] = NULL;
			pr.Pbasis[type]    = 0;
		}

	}

	type = SeasonTerm;
	if (basisIncluded[type] || (pr.isBeast && pr.beast_Pextra > 0) ) {
		VOID_PTR tmpY   = GetField123Check(S, "seasonY", 0);
		VOID_PTR tmpX   = GetField123Check(S, "seasonX", 0);
		I32      Pbasis = tmpY && tmpX && GetNumberOfElements(tmpY) == GetNumberOfElements(tmpX) ?
			               GetNumberOfElements(tmpY) : 0;

		I32    Ptotal = Pbasis;
		if (pr.isBeast) {
			Ptotal = Ptotal + pr.beast_Pextra;// If there is a seaosnal component
		}	

		I32PTR   varlistX  = pr.Pindices0[type] = malloc(Ptotal * sizeof(I32));
		I32PTR   varlistY  = pr.PseasonY0       = malloc(Ptotal * sizeof(I32));
		if (Pbasis >0) {
			CopyNumericObjToI32Arr(varlistX,  tmpX, Pbasis);
			CopyNumericObjToI32Arr(varlistY, tmpY, Pbasis);
		}
        //TODO: need to check for any duppliation
		 
		if (pr.isBeast) {
			for (int i = 0; i < pr.beast_Pextra; i++) {
				varlistX[Pbasis]  = pr.beast_Ptime;  
				varlistY[Pbasis]  = pr.beast_Pbase + i;
				Pbasis++;
			}			 
		}
		if (Pbasis > 0) {
			pr.Pbasis[type]           = Pbasis;
			pr.basisType[nGoodBasis++] = type;
		}	else {			
			free(pr.Pindices0[type]);
			pr.Pindices0[type] = NULL;
			pr.Pbasis[type]    = 0;
		}
	}
 
	type = HingeTerm;
	if (basisIncluded[type] ) {
		VOID_PTR tmpobj = GetField123Check(S, "hinge", 0);
		I32      Pbasis = tmpobj ? GetNumberOfElements(tmpobj) : meta->Pobj;
		I32      Pmax   = Pbasis + pr.isBeast; // add one more term if isBeast=TRUE
		I32PTR   varlist = pr.Pindices0[type] = malloc(Pmax * sizeof(I32));
		if (tmpobj) {
			CopyNumericObjToI32Arr(varlist, tmpobj, Pbasis);
			Pbasis = i32_unique_inplace(varlist, Pbasis);
		} else {
			i32_seq(varlist, 1, 1, Pbasis);     // if no square field, insert all variables
		}

		if (Pbasis > 0) {
			pr.Pbasis[type]           = Pbasis;
			pr.basisType[nGoodBasis++] = type;
		}	else {			
			free(pr.Pindices0[type]);
			pr.Pindices0[type] = NULL;
			pr.Pbasis[type]    = 0;
		}
	}

	type = HingePairTerm;
	if (basisIncluded[type] ) {
		VOID_PTR tmpobj = GetField123Check(S, "hingepair", 0);
		I32      Pbasis = tmpobj ? GetNumberOfElements(tmpobj) : meta->Pobj;
		I32      Pmax    = Pbasis + pr.isBeast; // add one more term if isBeast=TRUE
		I32PTR   varlist = pr.Pindices0[type] = malloc(Pmax * sizeof(I32));
		if (tmpobj) {
			CopyNumericObjToI32Arr(varlist, tmpobj, Pbasis);
			Pbasis = i32_unique_inplace(varlist, Pbasis);
		}
		else {
			i32_seq(varlist, 1, 1, Pbasis);     // if no square field, insert all variables
		}

		if (Pbasis > 0) {
			pr.Pbasis[type]           = Pbasis;
			pr.basisType[nGoodBasis++] = type;
		}	else {			
			free(pr.Pindices0[type]);
			pr.Pindices0[type] = NULL;
			pr.Pbasis[type]    = 0;
		}
	}

	type = ChangeTerm;
	if (basisIncluded[type]  ) {
		VOID_PTR tmpY   = GetField123Check(S, "changeY", 0);
		VOID_PTR tmpX   = GetField123Check(S, "changeX", 0);
		I32      Pbasis = tmpY && tmpX && GetNumberOfElements(tmpY) == GetNumberOfElements(tmpX) ?
			               GetNumberOfElements(tmpY) : 0;
	 	 
		I32PTR   varlistX  = pr.Pindices0[type] = malloc(Pbasis * sizeof(I32));
		I32PTR   varlistY = pr.PchangeY0       = malloc(Pbasis * sizeof(I32));
		if (Pbasis >0) {
			CopyNumericObjToI32Arr(varlistX,  tmpX, Pbasis);
			CopyNumericObjToI32Arr(varlistY, tmpY, Pbasis);
		}
        //TODO: need to check for any duppliation
		 
		if (Pbasis > 0) {
			pr.Pbasis[type]           = Pbasis;
			pr.basisType[nGoodBasis++] = type;
		}	else {			
			free(pr.Pindices0[type]);
			pr.Pindices0[type] = NULL;
			pr.Pbasis[type]    = 0;
		}
	}
 
 
	type = Hinge2DTerm;
	if (basisIncluded[type]) {
		VOID_PTR tmp    = GetField123Check(S, "hinge2d", 0);
		int      Pbasis = 0;
		if (tmp==NULL || !IsStruct(tmp) ) {
			pr.nGrpHinge2d        = 1L;
			pr.sGrpHinge2d = malloc(1L * sizeof(I32));
			pr.eGrpHinge2d = malloc(1L * sizeof(I32));
			Pbasis = meta->Pobj;;
			pr.sGrpHinge2d[0] = 1;
			pr.eGrpHinge2d[0] = Pbasis;
			I32PTR varlist = pr.Pindices0[type] = malloc(Pbasis * sizeof(I32));
			i32_seq(varlist, 1, 1, Pbasis);     // if no square field, insert all variables
		} else {
			int nGrp = pr.nGrpHinge2d = GetNumberOfElements(tmp);
			pr.sGrpHinge2d = malloc(nGrp * sizeof(I32));
			pr.eGrpHinge2d = malloc(nGrp * sizeof(I32));
			Pbasis = 0;
			for (int i = 0; i < nGrp; i++) {
				VOID_PTR fld   = GetFieldByIdx(tmp, i);
				int      nelem = GetNumberOfElements(fld);
				pr.sGrpHinge2d[i] = Pbasis + 1;
				pr.eGrpHinge2d[i] = pr.sGrpHinge2d[i]+nelem-1;
				Pbasis = Pbasis + nelem;
			}
			I32PTR varlist = pr.Pindices0[type] = malloc(Pbasis * sizeof(I32));

			for (int i = 0; i < nGrp; i++) {
				VOID_PTR fld = GetFieldByIdx(tmp, i); 
				CopyNumericObjToI32Arr(varlist+ pr.sGrpHinge2d[i]-1, fld, pr.eGrpHinge2d[i]- pr.sGrpHinge2d[i]+1L);
			}
		}
		if (Pbasis > 0) {
			pr.Pbasis[type]           = Pbasis;
			pr.basisType[nGoodBasis++] = type;
		}	else {			
			free(pr.Pindices0[type]);
			pr.Pindices0[type] = NULL;
			pr.Pbasis[type]    = 0;
		}
	}


	type = Hinge2DRTerm;
	if (basisIncluded[type]) {
		VOID_PTR tmp    = GetField123Check(S, "hinge2dr", 0);
		int      Pbasis = 0;
		if (tmp==NULL || !IsStruct(tmp) ) {
			pr.nGrpHinge2dr        = 1L;
			pr.sGrpHinge2dr = malloc(1L * sizeof(I32));
			pr.eGrpHinge2dr = malloc(1L * sizeof(I32));
			Pbasis = meta->Pobj;;
			pr.sGrpHinge2dr[0] = 1;
			pr.eGrpHinge2dr[0] = Pbasis;
			I32PTR varlist = pr.Pindices0[type] = malloc(Pbasis * sizeof(I32));
			i32_seq(varlist, 1, 1, Pbasis);     // if no square field, insert all variables
		} else {
			int nGrp = pr.nGrpHinge2dr = GetNumberOfElements(tmp);
			pr.sGrpHinge2dr = malloc(nGrp * sizeof(I32));
			pr.eGrpHinge2dr = malloc(nGrp * sizeof(I32));
			Pbasis = 0;
			for (int i = 0; i < nGrp; i++) {
				VOID_PTR fld   = GetFieldByIdx(tmp, i);
				int      nelem = GetNumberOfElements(fld);
				pr.sGrpHinge2dr[i] = Pbasis + 1;
				pr.eGrpHinge2dr[i] = pr.sGrpHinge2dr[i]+nelem-1;
				Pbasis = Pbasis + nelem;
			}
			I32PTR varlist = pr.Pindices0[type] = malloc(Pbasis * sizeof(I32));

			for (int i = 0; i < nGrp; i++) {
				VOID_PTR fld = GetFieldByIdx(tmp, i); 
				CopyNumericObjToI32Arr(varlist+ pr.sGrpHinge2dr[i]-1, fld, pr.eGrpHinge2dr[i]- pr.sGrpHinge2dr[i]+1L);
			}
		}
		if (Pbasis > 0) {
			pr.Pbasis[type]           = Pbasis;
			pr.basisType[nGoodBasis++] = type;
		}	else {			
			free(pr.Pindices0[type]);
			pr.Pindices0[type] = NULL;
			pr.Pbasis[type]    = 0;
		}
	}


	type = Hinge1DRTerm;
	if (basisIncluded[type]) {
		VOID_PTR tmp    = GetField123Check(S, "hinge1dr", 0);
		int      Pbasis = 0;
		if (tmp==NULL || !IsStruct(tmp) ) {
			pr.nGrpHinge1dr        = 1L;
			pr.sGrpHinge1dr = malloc(1L * sizeof(I32));
			pr.eGrpHinge1dr = malloc(1L * sizeof(I32));
			Pbasis = meta->Pobj;;
			pr.sGrpHinge1dr[0] = 1;
			pr.eGrpHinge1dr[0] = Pbasis;
			I32PTR varlist = pr.Pindices0[type] = malloc(Pbasis * sizeof(I32));
			i32_seq(varlist, 1, 1, Pbasis);     // if no square field, insert all variables
		} else {
			int nGrp = pr.nGrpHinge1dr = GetNumberOfElements(tmp);
			pr.sGrpHinge1dr = malloc(nGrp * sizeof(I32));
			pr.eGrpHinge1dr = malloc(nGrp * sizeof(I32));
			Pbasis = 0;
			for (int i = 0; i < nGrp; i++) {
				VOID_PTR fld   = GetFieldByIdx(tmp, i);
				int      nelem = GetNumberOfElements(fld);
				pr.sGrpHinge1dr[i] = Pbasis + 1;
				pr.eGrpHinge1dr[i] = pr.sGrpHinge1dr[i]+nelem-1;
				Pbasis = Pbasis + nelem;
			}
			I32PTR varlist = pr.Pindices0[type] = malloc(Pbasis * sizeof(I32));

			for (int i = 0; i < nGrp; i++) {
				VOID_PTR fld = GetFieldByIdx(tmp, i); 
				CopyNumericObjToI32Arr(varlist+ pr.sGrpHinge1dr[i]-1, fld, pr.eGrpHinge1dr[i]- pr.sGrpHinge1dr[i]+1L);
			}
		}
		if (Pbasis > 0) {
			pr.Pbasis[type]           = Pbasis;
			pr.basisType[nGoodBasis++] = type;
		}	else {			
			free(pr.Pindices0[type]);
			pr.Pindices0[type] = NULL;
			pr.Pbasis[type]    = 0;
		}
	}


	// Pre-detemine the tyoe-to-id conversion array
	// re-used in the model struct
	prior->numBasis = min( nGoodBasis, MAX_NUM_BASIS);
	for (int i = 0; i < MAX_NUM_BASIS;   i++)     prior->type2id[i] = -1;
	for (int i = 0; i < prior->numBasis; i++)	  prior->type2id[prior->basisType[i]] = i;

	/****************************************************************/
	// The indices in the R input are one-based, and Pindices0 are zero-based
	// Here decrement the idnices by one
	/****************************************************************/
	for (int i = 0; i < prior->numBasis; i++) {
		int    type        = pr.basisType[i];
		I32PTR Pindices0   = pr.Pindices0[type];
		int    N           = pr.Pbasis[type];;
		i32_add_val_inplace(-1L, Pindices0, N);
	}
	// Here decrement the idnices by one
	// i32_add_val_inplace(value, Xarray, N)
	i32_add_val_inplace(-1L, prior->PseasonY0,          prior->Pbasis[SeasonTerm]);
	i32_add_val_inplace(-1L, prior->PchangeY0,			prior->Pbasis[ChangeTerm]);
	i32_add_val_inplace(-1L, prior->FixedLinearTerms0,  prior->nFixedLinearTerms);

	prior->Pall = Pall;
	prior->Pobj = meta->Pobj;
}

#define  hasTerm(x) (prior->type2id[x] >= 0)


void  ParseArg_minKnotNum_minSepDist(VOIDPTR S, BVS_PRIOR_PTR prior, BVS_IO_PTR io) {
	//maxSeasonOrder will be feteched if isBeast==1 or isFFt==1

	VOIDPTR tmp;
	memset(prior->maxKnotNum_Vec, 0, sizeof(I16PTR) * MAX_NUM_BASIS);
	memset(prior->minSepDist_Vec, 0, sizeof(I16PTR) * MAX_NUM_BASIS);
	memset(prior->maxKnotNum_Max, 0, sizeof(I32) * MAX_NUM_BASIS);

	// Get the default glboal maxKnotNum and minSepDist first
	int maxKnotNum_Global = 10L;
	int minSepDist_Global = 3L;
	if (S) {		
		maxKnotNum_Global = (tmp = GetField123Check(S, "maxKnotNum", 0)) ? GetScalar(tmp) : maxKnotNum_Global;
		minSepDist_Global = (tmp = GetField123Check(S, "minSepDist", 0)) ? GetScalar(tmp) : minSepDist_Global;
	}	

	// Get the basiss-specific and variable-specific maxKnotNum and minSpepDists
	char* fldnames[MAX_NUM_BASIS];
	char* fldnames1[MAX_NUM_BASIS];

	fldnames[StepTerm]       = "maxKnotNumStep";
	fldnames[StairTerm]      = "maxKnotNumStair";
	fldnames[PieceTerm]      = "maxKnotNumPiece";
	fldnames[SeasonTerm]     = "maxKnotNumSeason";
	fldnames[HingeTerm]      = "maxKnotNumHinge";
	fldnames[HingePairTerm]  = "maxKnotNumHingePair";
	fldnames[ChangeTerm]     = "maxKnotNumChange";

	fldnames1[StepTerm]      = "minSepDistStep";
	fldnames1[StairTerm]     = "minSepDistStair";
	fldnames1[PieceTerm]     = "minSepDistPiece";
	fldnames1[SeasonTerm]     = "minSepDistSeason";
	fldnames1[HingeTerm]     = "minSepDistHinge";
	fldnames1[HingePairTerm] = "minSepDistHingePair";
	fldnames1[ChangeTerm]    = "minSepDistChange";

		
	for (int type = StepTerm; type <= ChangeTerm; type++) {

		if (!hasTerm(type)) 	   
			continue;

		int P = prior->Pbasis[type];
		prior->maxKnotNum_Vec[type] = malloc(sizeof(I16) * P );
		prior->minSepDist_Vec[type] = malloc(sizeof(I16) * P );

		///////////////////////////////////////////////////////////////////////////////////
		tmp =  GetField123Check(S, fldnames[type], 0);
		if (tmp) {
			int Ptrue = GetNumberOfElements(tmp);
			for (int i = 0; i < P; i++)  
				prior->maxKnotNum_Vec[type][i] = GetNumericElement(tmp, i % Ptrue);
		}	else {
			for (int i = 0; i < P; i++)  
				prior->maxKnotNum_Vec[type][i] = maxKnotNum_Global;
		}

		///////////////////////////////////////////////////////////////////////////////////
		tmp =  GetField123Check(S, fldnames1[type], 0)  ;
		if (tmp) {
			int Ptrue = GetNumberOfElements(tmp);
			for (int i = 0; i < P; i++)
				prior->minSepDist_Vec[type][i] = GetNumericElement(tmp, i % Ptrue);
		} else {
			for (int i = 0; i < P; i++)  
				prior->minSepDist_Vec[type][i] = minSepDist_Global;
		}


		int maxKnotNum_Max = 0;
		for (int i = 0; i < P; ++i) {
			maxKnotNum_Max = max(maxKnotNum_Max, prior->maxKnotNum_Vec[type][i]);
		}
		prior->maxKnotNum_Max[type] = maxKnotNum_Max;

	}

}

void  BVS_DeallocatePriorPts(BVS_PRIOR_PTR prior) {
	// Deallocate the mem alloacted in ProcessArgBasisType 
	// Called in glue_code

	for (int i = 0; i < MAX_NUM_BASIS; i++) { 
			free0(prior->Pindices0[i]); 
			free0(prior->maxKnotNum_Vec[i]);
			free0(prior->minSepDist_Vec[i]);
		
	} 
	free0(prior->PchangeY0);
	free0(prior->PseasonY0);
	free0(prior->FixedLinearTerms0);

	free0(prior->sGrpHinge1dr);
	free0(prior->sGrpHinge2dr);
	free0(prior->sGrpHinge2d);
	free0(prior->eGrpHinge1dr);
	free0(prior->eGrpHinge2dr);
	free0(prior->eGrpHinge2d);
	free0(prior->con1d_Decrease);
	free0(prior->con1d_Increase);
	free0(prior->con1d_nShape);
	free0(prior->con1d_uShape);
	free0(prior->con1d_NoWiggle);
	memset(prior, 0, sizeof(BVS_PRIOR));

}
static int  GetArg_3_Prior__(VOIDPTR prhs[], int nrhs, BVS_PRIOR_PTR prior, BVS_IO_PTR io)
{	 
	// Before running this fuction,  period and N must have been determined

	#define o  (*prior)

	VOIDPTR S = NULL;
	if (nrhs >= 4) {
		S = prhs[3L];
		if (!IsStruct(S)) {
			r_warning("WARNING: The arg 'prior' is ignored because it is not a List/Struct variable.");
			S = NULL;
		}
	} 

	 ParseArg_basisType(S, prior, &io->meta);
	 ParseArg_minKnotNum_minSepDist(S, prior, io);

	 struct PRIOR_MISSING {
		 U08   modelPriorType;
 		 U08   Kmax;
		 U08   sigFactor;
		 U08   outlierSigFactor;
		 U08   sig2;
		 U08   precValue;
		 U08   alpha1, alpha2, delta1, delta2;
		 U08   con1d_Decrease, con1d_Increase, con1d_nShape, con1d_uShape, con1d_NoWiggle;
	 } m = { 0, };

	 if (nrhs < 4)
		 memset(&m, 1L, sizeof(struct PRIOR_MISSING));

	if (nrhs >= 4) {		 
	    S = prhs[3L];
		if (!IsStruct(S)) {
			r_warning("WARNING: The arg 'prior' is ignored because it is not a List/Struct variable.");
			memset(&m, 1L, sizeof(struct PRIOR_MISSING));
			S = NULL;
		}
		else {
			VOIDPTR tmp;	
			o.Kmax              = (tmp = GetField123Check(S, "Kmax",2)) ?			GetScalar(tmp) : (m.Kmax = 1);
			o.sigFactor         = (tmp = GetFieldCheck(S,    "sigFactor")) ?		GetScalar(tmp) : (m.sigFactor = 1);			
			o.sig2              = (tmp = GetField123Check(S, "sig2",2)) ?			GetScalar(tmp) : (m.sig2 = 1);
			o.precValue		    = (tmp = GetField123Check(S, "precValue",5)) ?		GetScalar(tmp) : (m.precValue = 1);
			o.alpha1			= (tmp = GetField123Check(S, "alpha1",0)) ?			GetScalar(tmp) : (m.alpha1 = 1);
			o.alpha2			= (tmp = GetField123Check(S, "alpha2",0)) ?			GetScalar(tmp) : (m.alpha2 = 1);
			o.delta1			= (tmp = GetField123Check(S, "delta1",0)) ?			GetScalar(tmp) : (m.delta1 = 1);
			o.delta2			= (tmp = GetField123Check(S, "delta2",0)) ?			GetScalar(tmp) : (m.delta2 = 1);						
			//o.precPriorType		    = (tmp = GetFieldCheck(S, "precPriorType")) ?	    GetScalar(tmp) : (m.precPriorType = 1); 			

			m.con1d_Decrease = 1;
			if (tmp = GetField123Check(S, "con1d_Decrease", 7)) {
				m.con1d_Decrease = 0;
				o.num_con1d_Decrease = GetNumberOfElements(tmp);
				o.con1d_Decrease = malloc(sizeof(I32) * o.num_con1d_Decrease);
				CopyNumericObjToI32Arr(o.con1d_Decrease, tmp,o.num_con1d_Decrease);
			}
			m.con1d_Increase = 1;
			if (tmp = GetField123Check(S, "con1d_Increase", 7)) {
				m.con1d_Increase = 0;
				o.num_con1d_Increase = GetNumberOfElements(tmp);
				o.con1d_Increase = malloc(sizeof(I32) * o.num_con1d_Increase);
				CopyNumericObjToI32Arr(o.con1d_Increase, tmp, o.num_con1d_Increase);
			}
			m.con1d_nShape = 1;
			if (tmp = GetField123Check(S, "con1d_nShape", 7)) {
				m.con1d_nShape = 0;
				o.num_con1d_nShape = GetNumberOfElements(tmp);
				o.con1d_nShape = malloc(sizeof(I32) * o.num_con1d_nShape);
				CopyNumericObjToI32Arr( o.con1d_nShape, tmp, o.num_con1d_nShape);
			}
			m.con1d_uShape = 1;
			if (tmp = GetField123Check(S, "con1d_uShape", 7)) {
				m.con1d_uShape = 0;
				o.num_con1d_uShape = GetNumberOfElements(tmp);
				o.con1d_uShape = malloc(sizeof(I32) * o.num_con1d_uShape);
				CopyNumericObjToI32Arr( o.con1d_uShape, tmp, o.num_con1d_uShape);
			}
			m.con1d_NoWiggle = 1;
			if (tmp = GetField123Check(S, "con1d_NoWiggle", 7)) {
				m.con1d_NoWiggle = 0;
				o.num_con1d_NoWiggle =  GetNumberOfElements(tmp);
				o.con1d_NoWiggle     = malloc(sizeof(I32) * o.num_con1d_NoWiggle);
				CopyNumericObjToI32Arr(o.con1d_NoWiggle, tmp, o.num_con1d_NoWiggle);
			}
		}

	} // if (nrhs >= 4)

	 
	if (m.Kmax)            o.Kmax       = 300;
	if (m.sigFactor)       o.sigFactor  = 1.8;      	
	
	if (m.sig2 )           o.sig2      = 0.2f;           
	if (m.precValue)       o.precValue = 1.5f;             
	if (m.alpha1)		   o.alpha1	 = 0.00000001f;
	if (m.alpha2)		   o.alpha2	 = 0.00000001f;
	if (m.delta1)		   o.delta1	 = 0.00000001f;
	if (m.delta2)		   o.delta2	 = 0.00000001f;
	if (m.con1d_Decrease) 	  o.num_con1d_Decrease = 0, o.con1d_Decrease = NULL;
	if (m.con1d_Increase) 	  o.num_con1d_Increase = 0, o.con1d_Increase = NULL;
	if (m.con1d_nShape) 	  o.num_con1d_nShape = 0,   o.con1d_nShape = NULL;
	if (m.con1d_uShape) 	  o.num_con1d_uShape = 0,   o.con1d_uShape = NULL;
	if (m.con1d_NoWiggle) o.num_con1d_NoWiggle = 0,  o.con1d_NoWiggle = NULL;

	o.sigFactor = min(o.sigFactor, 1.02);
	o.sig2		= max(o.sig2, 0.01);
	o.precValue = max(o.precValue, 0.01);

	return 1;

#undef o
}



static int  GetArg_4_MCMC___(VOIDPTR prhs[], int nrhs, BVS_MCMC_PTR mcmc,  BVS_OPTIONS_PTR opt)
{
	// Before running this function, trendMinSepDist must have been set.

	#define o (*mcmc)
	struct MCMC_MISSING {
		U08   seed;	                  // Unsigned long long seed;
		U08   credIntervalAlphaLevel;
		U08   trendResamplingOrderProb;
		U08   seasonResamplingOrderProb;
		U08   ridgeFactor;
		U08   burnin, samples, chainNumber;
		U08   maxMoveStepSize;
		U08   thinningFactor;
	} m = {0,};


	if (nrhs < 5) 
		memset(&m, 1L, sizeof(struct MCMC_MISSING));

	if (nrhs >= 5) {

		VOIDPTR S = prhs[4L];
		if (!IsStruct(S)) {
			r_warning("WARNING: The arg 'mcmc' is ignored because it is not a LIST variable.");
			memset(&m, 1L, sizeof(struct MCMC_MISSING));
		} else {
			VOIDPTR tmp;
			o.maxMoveStepSize = (tmp = GetField123Check(S, "maxMoveStepSize",2))? GetScalar(tmp) : (m.maxMoveStepSize = 1);
			o.samples         = (tmp = GetField123Check(S, "samples", 2)) ?        GetScalar(tmp) : (m.samples = 1);
			o.thinningFactor  = (tmp = GetField123Check(S, "thinningFactor", 2)) ? GetScalar(tmp) : (m.thinningFactor = 1);
			o.burnin          = (tmp = GetField123Check(S, "burnin", 2)) ?         GetScalar(tmp) : (m.burnin = 1);
			o.chainNumber     = (tmp = GetField123Check(S, "chainNumber", 2)) ?    GetScalar(tmp) : (m.chainNumber = 1);
			o.seed			  = (tmp = GetField123Check(S, "seed", 2)) ?			GetScalar(tmp) : (m.seed = 1);
			o.ridgeFactor	  = (tmp = GetField123Check(S, "ridgeFactor", 2)) ?	GetScalar(tmp) : (m.ridgeFactor = 1);

			o.trendResamplingOrderProb  = (tmp = GetField123Check(S, "trendResamplingOrderProb",  2)) ? GetScalar(tmp) : (m.trendResamplingOrderProb = 1);
			o.seasonResamplingOrderProb = (tmp = GetField123Check(S, "seasonResamplingOrderProb", 2)) ? GetScalar(tmp) : (m.seasonResamplingOrderProb = 1);
			o.credIntervalAlphaLevel    = (tmp = GetField123Check(S, "credIntervalAlphaLevel",    2)) ? GetScalar(tmp) : (m.credIntervalAlphaLevel = 1);
		}

	} // if (nrhs >= 5)

 
	//r_printf("move :%d  %d %d\n", o.maxMoveStepSize, opt->io.meta.hasSeasonCmpnt , opt->prior.trendMinSepDist  );	
	if (m.samples         || o.samples==0)		   o.samples		 = 3000;        o.samples		  = max(o.samples, 1);
	if (m.thinningFactor  || o.thinningFactor==0)  o.thinningFactor  = 1L;			o.thinningFactor = max(o.thinningFactor, 1L);
	if (m.burnin          || o.burnin==0)          o.burnin		     = 150L;		o.burnin		  = max(o.burnin, 150L);
	if (m.chainNumber || o.chainNumber==0)         o.chainNumber	 = 3;			o.chainNumber	  = max(o.chainNumber, 1L);
	if (m.seed)            o.seed			 = 0L;
	if (m.credIntervalAlphaLevel)      o.credIntervalAlphaLevel	 = .95;
	if (m.ridgeFactor)     o.ridgeFactor	 = 0.0001f;

	if (m.trendResamplingOrderProb)  o.trendResamplingOrderProb  = .1f;
	if (m.seasonResamplingOrderProb) o.seasonResamplingOrderProb = .17f;

	return 1;

#undef o
}
static int  GetArg_5_EXTRA__(VOIDPTR prhs[], int nrhs, BVS_EXTRA_PTR extra)
{
	// Before running this function, meta.whichDimIsTime must be first obtained

	#define o (*extra)
	struct OUTFLAGS_MISSING {
		I08   numThreadsPerCPU;
		I08   numParThreads;
		I08   numCPUCoresToUse;		

		I08   consoleWidth;
		I08   whichOutputDimIsTime;
		I08   removeSingletonDims;
		I08   dumpInputData;

		I08   ncpStatMethod;
		I08  smoothCpOccPrCurve;
		I08  useMeanOrRndBeta;
		I08  computeCredible;
		I08  fastCIComputation;

		I08  printOptions;
		I08  printProgressBar;
	} m = {0,};


	if (nrhs < 6) 
		memset(&m, 1L, sizeof(struct OUTFLAGS_MISSING));	

	if (nrhs >= 6) {
		VOIDPTR S = prhs[5L];
		if (!IsStruct(S)) {
			r_warning("WARNING: The arg 'extra' is ignored because it is not a LIST variable.");
			memset(&m, 1L, sizeof(struct OUTFLAGS_MISSING));
		}
		else {
			VOIDPTR tmp;
			o.whichOutputDimIsTime = (tmp = GetField123Check(S, "whichOutputDimIsTime",2)) ?	GetScalar(tmp) : (m.whichOutputDimIsTime = 1);
			o.removeSingletonDims  = (tmp = GetField123Check(S, "removeSingletonDims", 8)) ? GetScalar(tmp) : (m.removeSingletonDims = 1);			
			o.numThreadsPerCPU     = (tmp = GetField123Check(S, "numThreadsPerCPU", 4)) ? GetScalar(tmp) : (m.numThreadsPerCPU = 1);
			o.numParThreads        = (tmp = GetField123Check(S, "numParThreads", 4)) ?			GetScalar(tmp) : (m.numParThreads = 1);
			o.numCPUCoresToUse     = (tmp = GetField123Check(S, "numCPUCoresToUse", 4)) ?		GetScalar(tmp) : (m.numCPUCoresToUse = 1);
			o.consoleWidth         = (tmp = GetField123Check(S, "consoleWidth",2)) ?			GetScalar(tmp) : (m.consoleWidth = 1);
			o.dumpInputData        = (tmp = GetField123Check(S, "dumpInputData",2)) ? GetScalar(tmp) : (m.dumpInputData = 1);			
			#define _1(x)       o.x = (tmp=GetFieldCheck(S,#x))? GetScalar(tmp): (m.x=1)
			#define _2(x,y)     _1(x);_1(y)
			#define _3(x,y,z)   _1(x);_2(y,z)
			#define _4(x,y,z,w) _2(y,z);_2(y,z)
			#define _5(x,y,z,w) _2(y,z);_3(y,z)

			_2(printProgressBar, printOptions);
			_2(computeCredible,  fastCIComputation);
			_1(useMeanOrRndBeta);

		 
				 
		} // if (!IsStruct(S)) : S is a struct
	} // if (nrhs >= 5)

	if (m.removeSingletonDims)		o.removeSingletonDims  = 1;
 
	if (m.dumpInputData)		 o.dumpInputData       = 0;
	if (m.numThreadsPerCPU)      o.numThreadsPerCPU    = 2;
	if (m.numParThreads)         o.numParThreads		= 0;
	if (m.numCPUCoresToUse)      o.numCPUCoresToUse	= 0;	
	if (m.consoleWidth||o.consoleWidth<=0)  o.consoleWidth= GetConsoleWidth(); 	o.consoleWidth = max(o.consoleWidth, 40);
	if (m.printProgressBar)      o.printProgressBar	= 1;
	if (m.printOptions)          o.printOptions		= 1;

	if (m.computeCredible)       o.computeCredible		= 0L;
	if (m.fastCIComputation)     o.fastCIComputation	= 1L;

	if (m.useMeanOrRndBeta)      o.useMeanOrRndBeta = 0;
	
	return 1;
#undef o
}

void BVS_Print(BVS_OPTIONS_PTR opt){

	int           len   = 15;
	BVS_PRIOR_PTR prior = &opt->prior;
 
	r_printf("%-*s : Pobj=%d Pall=%d\n", len, "Dimension", prior->Pobj, prior->Pall);
	if (prior->isBeast) {
		r_printf("%-*s : Ptime=%d Pbase=%d Pnewterm=%d  maxSeasonOrder=%d\n", len, "BEAST", prior->beast_Ptime,  prior->beast_Pbase, prior->beast_Pextra, prior->beast_maxSeasonOrder); 
	}
	if (prior->isFFT) {
		r_printf("%-*s : Ptime=%d Pbase=%d Pnewterm=%d  maxSeasonOrder=%d\n", len, "FFT", prior->fft_Ptime,  prior->fft_Pbase, prior->fft_Pnew, prior->fft_maxOrder);
	}

	I32PTR Pbasis          = prior->Pbasis;
	I32PTR *Pindices       = prior->Pindices0;
	char*  names[]         = {"linear","FixedLinTerms", "quadratic","step", 
							   "stair", "piecewise", "season", "seasonY", "hinge", "hingepair","change(X)","change(Y)" };
	int    arrP[]          = { Pbasis[LinearTerm], prior->nFixedLinearTerms,  Pbasis[QuadraticTerm], Pbasis[StepTerm],
							   Pbasis[StairTerm], Pbasis[PieceTerm], Pbasis[SeasonTerm], Pbasis[SeasonTerm],  Pbasis[HingeTerm],  Pbasis[HingePairTerm], Pbasis[ChangeTerm],Pbasis[ChangeTerm] };
	I32PTR arrPindices0[] = { Pindices[LinearTerm], prior->FixedLinearTerms0, Pindices[QuadraticTerm], Pindices[StepTerm],
							Pindices[StairTerm],Pindices[PieceTerm],Pindices[SeasonTerm], prior->PseasonY0
		                    ,Pindices[HingeTerm], Pindices[HingePairTerm], Pindices[ChangeTerm],prior->PchangeY0,};
	
	/////////////////////////////////////////////////////////
	for (int j = 0; j < sizeof(arrP)/sizeof(int); j++) {  
		r_printf("%-*s :", len, names[j]);

		int   P  = arrP[j];		
		if (P == 0) {		
			r_printf(" none \n"); 
		} else {
			I32PTR Pindices0 = arrPindices0[j];
			r_printf(" [");
			for (int i = 0; i < P - 1; i++) {
				r_printf("%d,", Pindices0[i] + 1L);
			}		
		    r_printf("%d]\n", Pindices0[P - 1L] + 1L);
		}
	
	}
	/////////////////////////////////////////////////////////
 
	 
	/////////////////////////////////////////////////////////
	int    P		  = prior->Pbasis[Hinge2DTerm];
	I32PTR Pindices0  = prior->Pindices0[Hinge2DTerm];
	r_printf("%-*s :", len, "hinge2d");
	if (prior->nGrpHinge2d == 0) {
		r_printf(" none \n");
	} else {
		r_printf(" %d groups \n", prior->nGrpHinge2d);
		for (int j = 0; j < prior->nGrpHinge2d; j++) {
			int s = prior->sGrpHinge2d[j];
			int e = prior->eGrpHinge2d[j];
			r_printf("%-*s - grp#%-3d [", len, " ", j+1);
			for (int i = s; i < e; i++) {
				r_printf("%d,", Pindices0[i - 1] + 1L);
			}
			r_printf("%d]\n", Pindices0[e - 1L] + 1L); 
		}		
	}


	/////////////////////////////////////////////////////////
	P		  = prior->Pbasis[Hinge2DRTerm];
	Pindices0 = prior->Pindices0[Hinge2DRTerm];
	r_printf("%-*s :", len, "hinge2dr");
	if (prior->nGrpHinge2dr == 0) {
		r_printf(" none \n");
	} else {
		r_printf(" %d groups \n", prior->nGrpHinge2dr);
		for (int j = 0; j < prior->nGrpHinge2dr; j++) {
			int s = prior->sGrpHinge2dr[j];
			int e = prior->eGrpHinge2dr[j];
			r_printf("%-*s - grp#%-3d [", len, " ", j + 1);
			for (int i = s; i < e; i++) {
				r_printf("%d,", Pindices0[i - 1] + 1L);
			}
			r_printf("%d]\n", Pindices0[e - 1L] + 1L);
		}
	}
	;

	/////////////////////////////////////////////////////////
	P		  = prior->Pbasis[Hinge1DRTerm];
	Pindices0 = prior->Pindices0[Hinge1DRTerm];
	r_printf("%-*s :", len, "hinge1dr");
	if (prior->nGrpHinge1dr == 0) {
		r_printf(" none \n");
	}
	else {
		r_printf(" %d groups \n", prior->nGrpHinge1dr);
		for (int j = 0; j < prior->nGrpHinge1dr; j++) {
			int s = prior->sGrpHinge1dr[j];
			int e = prior->eGrpHinge1dr[j];
			r_printf("%-*s - grp#%-3d [", len, " ", j + 1);
			for (int i = s; i < e; i++) {
				r_printf("%d,", Pindices0[i - 1] + 1L);
			}
			r_printf("%d]\n", Pindices0[e - 1L] + 1L);
		}
	}

}


int BVS_GetArgs(VOIDPTR prhs[], int nrhs, A(OPTIONS_PTR) opt) {

	int  failed = !GetArg_12_MetaData(prhs, nrhs, &opt->io.meta)		 || 
			      !GetArg_3_Prior__(prhs,  nrhs, &opt->prior, &opt->io)   ||
			      !GetArg_4_MCMC___(prhs,  nrhs, &opt->mcmc,  opt)     ||
			      !GetArg_5_EXTRA__(prhs,  nrhs, &opt->extra ) ;
	int success = !failed;	
 
	BVS_Print(opt);
	return success;
}



#include "abc_000_warning.h"