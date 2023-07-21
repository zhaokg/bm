#include <math.h>
#include <string.h>
#include "assert.h"
#include "abc_000_warning.h"

#include "abc_ide_util.h"
#include "abc_common.h"
#include "abc_date.h"

 
// https:// stackoverflow.com/questions/26271154/how-can-i-make-a-mex-function-printf-while-its-running
extern void matlab_IOflush(void);


// C++ solution to define a keyword as a function name
// #pragma warning(suppress: 4483)
// extern void __identifier("?ioFlush@@YA_NXZ")(void);

I08 IDE_USER_INTERRUPT;

void printProgress(F32 pct, I32 width, char * buf, I32 firstTimeRun)
{//https:// stackoverflow.com/questions/2685435/cooler-ascii-spinners

	
 	static char spinnerChar[] =  "|/-\\";
	static I32  cnt = 1;
	cnt++;
	cnt = cnt == 4 ? 0 : cnt;
	width = max(width, 35); //cause errors when width is 0 or negative
	memset(buf, '*', width); // space = 20

	I32  len = 0;
	buf[len++] = spinnerChar[cnt];

	char prefix[] = "Progress:";
	I32 strLen    = sizeof(prefix)-1L; //the last byte is a Zero
	memcpy(buf+len, prefix, strLen);
	len += strLen;

	snprintf(buf + len, 15L, "%5.1f%% done", pct * 100);
	len += 5+1+5;
	buf[len++] = '[';

	I32 finishedLen = round((width - len - 1)*pct);
	memset(buf + len, '=', finishedLen);
	len += finishedLen;
	buf[len++] = '>';
	//memset(buf + len, ' ', width-len);
	buf[width - 1] = ']';
	buf[width] = 0;

#if R_INTERFACE==1
	Rprintf("\r%s", buf); 
#elif P_INTERFACE==1
	r_printf("\r%s", buf);
	//R doesnto allow io operationrs from external libraries
	//fflush(stdout);
#elif M_INTERFACE==1
	if (firstTimeRun == 1)
	{
		r_printf("\r\n");
		r_printf("%s", buf);		
		matlab_IOflush();
		//mexEvalString("drawnow");
		//mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
		//mexEvalString("pause(0.0000001);drawnow");
	}
	else
	{
		char * back = buf + width + 5;
		memset(back, '\b', width + 2);
		back[width + 2] = 0;

		r_printf(back);
		r_printf("%s\r\n", buf);
		matlab_IOflush();
		//mexEvalString("drawnow");
		//mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
		//mexEvalString("pause(0.0000001);drawnow");

	}


#endif	
}

void printProgress2(F32 pct, F64 time, I32 width, char * buf, I32 firstTimeRun)
{//https:// stackoverflow.com/questions/2685435/cooler-ascii-spinners
	
	static char spinnerChar[] = "|/-\\";
	static int  count         = 1;
	count = (++count) == 4 ? 0 : count;
	width = max(width, 40); //cause errors when width is 0 or negative
	memset(buf, '*', width); // space = 20

	I32  len = 0;
	buf[len++] = (pct<1.0) ? spinnerChar[count]: ' ';

	snprintf(buf + len, 10L, "%5.1f%%", pct * 100);
	len += 5 + 1;

	char prefix[] = "done";
	I32 strLen = sizeof(prefix)-1L; // the last byte is a Zero
	memcpy(buf + len, prefix, strLen);
	len += strLen;

	F64 SecsPerDay = 3600 * 24;
	I32 days = time / SecsPerDay;
	F64 tmp = time - days *SecsPerDay;
	I32 hrs = tmp / 3600;
	tmp = (tmp - hrs * 3600);
	I32 mins = tmp / 60;
	tmp = tmp - mins * 60;
	I32 secs = tmp;
	days = days >= 99 ? 99 : days;

	//https://stackoverflow.com/questions/71013186/warning-name0-directive-output-may-be-truncated-writing-10-bytes-into-a-regi
	if (time > SecsPerDay)
		snprintf(buf + len, 100L, "<Remaining%02dday%02dhrs%02dmin>", days, hrs, mins);
	else
		snprintf(buf + len, 100L,"<Remaining%02dhrs%02dmin%02dsec>", hrs, mins, secs);
	len += 26;

	buf[len++] = '[';

	I32 finishedLen = round((width - len - 1)*pct);
	memset(buf + len, '=', finishedLen);
	len += finishedLen;
	buf[len++] = '>';
	//memset(buf + len, ' ', width-len);
	buf[width - 1] = ']';
	buf[width] = 0;

#if R_INTERFACE==1
	r_printf("\r%s", buf);
	//R doesnto allow io operationrs from external libraries
	//fflush(stdout);
#elif P_INTERFACE==1
	r_printf("\r%s", buf);
#elif M_INTERFACE==1

	if (firstTimeRun == 1)
	{
		r_printf("\r\n");
		r_printf("%s", buf);
		//mexEvalString("drawnow");
		//mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
		matlab_IOflush();
	}
	else {
		char * back = buf + width + 5;
		memset(back, '\b', width + 2);
		back[width + 2] = 0;

		r_printf(back);
		r_printf("%s\r\n", buf);
		//mexEvalString("drawnow");
		//mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
		matlab_IOflush();
	}
#endif	
}

void RemoveField(FIELD_ITEM *fieldList, int nfields, char * fieldName)
{
	for (I64 i = 0; i < nfields && fieldList[i].name[0]!=0; i++) {
		if (strcmp(fieldList[i].name, fieldName) == 0)	{
			// The two strings are equal
			if (fieldList[i].ptr) {
				fieldList[i].ptr[0] = NULL; // the storage pointer is set to NULL
			}
			fieldList[i].ptr = NULL;        // the pointer to the storage pointer is set to NULL
			break;
		}
	}
}


void RemoveSingltonDims(FIELD_ITEM* flist, I32 nlist) {
	// 1 x1    -> 1
	// 2x1x1x2 -> 2x2
	// 1x10    -> 10

	// 10 : in R, this will be created as a vector, and in Matlab 
	// be created as a column vector.

	 #define  MAX_NUM_DIM   10 
	for (int i = 0; i < nlist && flist[i].name[0]!=0; i++) {

		if ( flist[i].ndim == 1) continue;		

		int goodN   = 0;
		int goodDims[MAX_NUM_DIM];
		// FInd all the non-one dims
		for (int j = 0; j < flist[i].ndim; j++) {
			if (flist[i].dims[j] != 1L) {
				goodDims[goodN++] = flist[i].dims[j];
			}
		}

		// Re-assign the non-one dims
		if (goodN == 0) {
			// it is 1x 1 x1 ..,
			flist[i].ndim    = 1;
			flist[i].dims[0] = 1;
		} else {
			flist[i].ndim = goodN;
			for (int j = 0; j < goodN; j++) {
				flist[i].dims[j] = goodDims[j];
			}
		}
	}

}

int CopyNumericObjToF32Arr(F32PTR outmem, VOID_PTR infield, int N) {

	VOID_PTR data = GetData(infield);

	if (IsSingle(infield))     		memcpy(outmem, data, sizeof(F32) * N);
	else if (IsDouble(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((double*)data + i);
	else if (IsInt32(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((int*)data + i);
	else if (IsInt64(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((I64*)data + i);
	else if (IsChar(infield))		return 0;
	else {	
		return 0;
	}
	return 1L;
}

int CopyNumericObjToI32Arr(I32PTR outmem, VOID_PTR infield, int N) {

	VOID_PTR data = GetData(infield);

	if      (IsInt32(infield))    	memcpy(outmem, data, sizeof(I32) * N);
	else if (IsDouble(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((double*)data + i);
	else if (IsSingle(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((int*)data + i);
	else if (IsInt64(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((I64*)data + i);
	else if (IsChar(infield))	  return 0;
	else {
		return 0;
	}
	return 1L;
}

int HaveEqualDimesions(const void* p1, const void* p2) {
	int dim1 = GetNumOfDim(p1);
	int dim2 = GetNumOfDim(p2);
	if (dim1 != dim2) return 0;

	I32 dims1[5], dims2[5];
	GetDimensions(p1, dims1, dim1);
	GetDimensions(p2, dims2, dim2);

	I32 equal = 1;
	for (int i = 0; i < dim1; ++i) {
		equal = equal & (dims1[i] == dims2[i]);
	}
	return equal;
}

int GetDataType(VOID_PTR Y) {

	if      (IsInt32(Y)) 						return DATA_INT32;
	else if (IsInt16(Y)) 						return DATA_INT16;
	else if (IsInt64(Y)) 						return DATA_INT64;
	else if (IsDouble(Y))  /* isReal(pY)*/		return DATA_DOUBLE;
	else if (IsSingle(Y))  /* isReal(pY)*/ 		return DATA_FLOAT;
	else                                        return DATA_UNKNOWN;
}

#if R_INTERFACE ==1 || M_INTERFACE ==1
F64 GetNumericElement(const void* Y, I32 idx) {

	if (!IsNumeric(Y)) {
		return getNaN();
	}

	I32 n = GetNumberOfElements(Y);
	if (n == 1) {
		// Y is a scalar
		if (idx==0)		return GetScalar(Y);
		else       		return getNaN();
	} else {
		// Y is a vector of more than one elements
		if (idx < n) {
			VOID_PTR y = GetData(Y);
			if (IsInt32(Y)) 						    return *((I32PTR)y + idx);
			else if (IsInt16(Y)) 						return *((I16PTR)y + idx);
			else if (IsDouble(Y))  /* isReal(pY)*/		return *((F64PTR)y + idx);
			else if (IsSingle(Y))  /* isReal(pY)*/ 		return *((F32PTR)y + idx);
			else                                        return getNaN();
		}	else {
			// idx is out of bounds
			return getNaN();
		}
	}

}

void* CvtToPyArray_NewRef(VOIDPTR Y) { return Y; }
#endif

void* GetField123Check(const void* structVar, char* fname, int nPartial) {
	// Check if the retured field is an empety or R_NilValue
	if (structVar == NULL) return NULL;

	VOID_PTR p = GetField123(structVar, fname, nPartial);
	if (p == NULL || IsEmpty(p))
		return NULL;
	else {
		return p;
	}

}

void* GetFieldCheck(const void* structVar, char* fname) {
	if (structVar == NULL) return NULL;
	// Check if the retured field is an empety or R_NilValue
	VOID_PTR p = GetField(structVar, fname);
	if (p == NULL || IsEmpty(p))
		return NULL;
	else {
		return p;
	}
}

int  GetNumElemTimeObject(  VOID_PTR timeObj ) {		
	// time should have been pre-allocated, with a length of Nraw
	if (timeObj == NULL) return -1L;
 
	/////////////////////////////////////////////////////
	// timeObj is not a struct variable
	/////////////////////////////////////////////////////	
	if (IsNumeric(timeObj) || IsChar(timeObj)) {
		return GetNumberOfElements(timeObj);
	} // 	if (IsNumeric(timeObj)==0) 

	// if it is not numeric or Struct.
	if ( IsStruct(timeObj) == 0 ) {
		return -1;
	} //if ( IsStruct(timeObj) == 0 )

	  /////////////////////////////////////////////////////
	// timeObj  is a struct variable
	/////////////////////////////////////////////////////	
	VOIDPTR yr  = GetField123Check(timeObj, "year",1);
	VOIDPTR mn  = GetField123Check(timeObj, "month",1);
	VOIDPTR day = GetField123Check(timeObj, "day",3);
	VOIDPTR doy = GetField123Check(timeObj, "doy",3);

	int isTimeProcessed = 0;
	if (!isTimeProcessed && yr && mn && IsNumeric(yr) && IsNumeric(mn) ) {
	
		int Nyr = GetNumberOfElements(yr);
		int Nmn = GetNumberOfElements(mn);
		int Ntime = Nyr;
		if (Nyr != Nmn) {		 
			return -1L;
		}	
   
		if (day && IsNumeric(day) && GetNumberOfElements(day) == Ntime) {	} 
		else {
			return -1L;
		}

		isTimeProcessed = 1; 
		return Ntime;
	}  

	if (!isTimeProcessed && yr && doy && IsNumeric(yr) && IsNumeric(doy) ) 	{
 
		int Nyr  = GetNumberOfElements(yr);
		int Ndoy = GetNumberOfElements(doy);
		int Ntime = Nyr;
		if (Nyr != Ndoy) {		 
			return -1;
		}	 
		isTimeProcessed = 1;	 
		return Ntime;
	}


	VOIDPTR datestr = GetField123Check(timeObj, "dateStr",3); 
	if (!isTimeProcessed && datestr && IsChar(datestr) ) 	{
		return GetNumberOfElements(datestr); 
	}
	return -1;
}
 
F32PTR  CvtTimeObjToF32Arr(VOID_PTR timeObj, int * Nactual ) {		

	//  Return an allocated mem buf

 	/////////////////////////////////////////////////////
	// timeObj is not a struct variable
	/////////////////////////////////////////////////////	
	if ( IsStruct(timeObj) == 0 ) {

		int Ntime = GetNumberOfElements(timeObj);
		
		if ( IsNumeric(timeObj)  ) {
			F32PTR ftime = malloc(sizeof(F32)*Ntime);
			if (IsClass(timeObj, "Date")) {
				// This branch is true only for R and Matlab's IsClass always return FALSE
				// in R, tmp is an integer array
				F64PTR days = GetData(timeObj);
				for (int i = 0; i < Ntime; ++i) {
					ftime[i] = fractional_civil_from_days((int)days[i]);
				}
			} else {
				int status = CopyNumericObjToF32Arr(ftime, timeObj, Ntime);
				if ( status == 0 ) {										
					free(ftime);
					r_error("ERROR: time has an unsupported data format or type.\n");
					return NULL;
				}
			}
			*Nactual = Ntime;
			return ftime;
			
		}
		else if (IsChar(timeObj)) {
			DynMemBuf        strarr  = { .max_len = 100 };
			DynAlignedBuf    start   = { .align = 4,.elem_size = 4,.p.raw = NULL, };
			DynAlignedBuf    nchars  = { .align = 4,.elem_size = 4,.p.raw = NULL, };

			CharObj2CharArr(timeObj, &strarr, &start, &nchars);
			float* ftime = strings_to_fyears(strarr.raw, start.p.i32, start.cur_len);
			// ftime is NULL upon return if failled

			Ntime=start.cur_len;
			dynbuf_kill(&strarr);
			adynbuf_kill(&start);
			adynbuf_kill(&nchars);
			*Nactual = Ntime;
			return ftime;
		}
		else {	
			r_error("ERROR: time is not numeric or strings. If times are strings, use time$dateStr and time$strFmt to specify data observation times.\n");
			return NULL;				
		} // 	if (IsNumeric(timeObj)==0) {
		
	} //if ( IsStruct(timeObj) == 0 )


	/////////////////////////////////////////////////////
	// timeObj  is a struct variable
	/////////////////////////////////////////////////////	

	VOIDPTR yr  = GetField123Check(timeObj, "year",1);
	VOIDPTR mn  = GetField123Check(timeObj, "month",1);
	VOIDPTR day = GetField123Check(timeObj, "day",3);
	VOIDPTR doy = GetField123Check(timeObj, "doy",3);

	int isTimeProcessed = 0;

	if (!isTimeProcessed && yr && mn && IsNumeric(yr) && IsNumeric(mn) ) {
	
		int Nyr = GetNumberOfElements(yr);
		int Nmn = GetNumberOfElements(mn);
		if (Nyr != Nmn) {
			q_warning("WARNING: time$year and time$month should have the same length.\n"); 
			goto __ENTRY_YEAR_DAY_LOC;
		}	
	
		int    Ntime = Nyr;
		F32PTR yr32  = malloc(sizeof(F32) * 3 * Ntime);
		F32PTR mn32  = yr32 + Ntime; // one for mn32 and another for day32
		F32PTR day32 = mn32 + Ntime;

		if (!CopyNumericObjToF32Arr(yr32, yr, Ntime)) {
			q_warning("WARNING: time$year has an unsupported data format or type.\n");
			free(yr32);
			goto __ENTRY_YEAR_DAY_LOC;
		}
		if (!CopyNumericObjToF32Arr(mn32, mn, Ntime)) {
			q_warning("WARNING: time$month has an unsupported data format or type.\n");
			free(yr32);
			goto __ENTRY_YEAR_DAY_LOC;
		}
 
		if (day && IsNumeric(day) && GetNumberOfElements(day) == Ntime) {
			/* time$day is present */
			if (!CopyNumericObjToF32Arr(day32,day, Ntime)) {
				q_warning("WARNING: time$day has an unsupported data format or type.\n");
				free(yr32);
				goto __ENTRY_YEAR_DAY_LOC;
			}
			
			for (int i = 0; i < Ntime; ++i) {
				yr32[i] = YMDtoF32time(yr32[i], mn32[i], day32[i]);
				if (yr32[i] < -1e9) {
					q_warning("WARNING: The (%d)-ith date (time$year=%d,time$month=%d, and time$day=%d) is not valid.\n",i+1, (int) yr32[i], (int)mn32[i], (int)day32[i] );
					free(yr32);
					goto __ENTRY_YEAR_DAY_LOC;
				}
			}
			
		} else {
			/*  time$day is not present */

			q_warning("WARNING: time$day is not specified, so only time$year and time$month are used!\n");
			for (int i = 0; i < Ntime; ++i) {
				yr32[i] =  yr32[i] + mn32[i]/12.0 -1./24.0;
				if (yr32[i] < -1e9) {
					q_warning("WARNING: The (%d)-ith date (metadata$time$year=%d,and metadata$time$month=%d) is not valid.\n", i+1, (int)yr32[i], (int)mn32[i]);
					free(yr32);
					goto __ENTRY_YEAR_DAY_LOC;
				}
			}			
		}

		isTimeProcessed = 1;
		*Nactual = Ntime;
		return yr32;
	}  


	__ENTRY_YEAR_DAY_LOC:
	if ( !isTimeProcessed && yr && doy && IsNumeric(yr) && IsNumeric(doy) ) {
 
		int Nyr  = GetNumberOfElements(yr);
		int Ndoy = GetNumberOfElements(doy);	
		if (Nyr != Ndoy) {
			q_warning("WARNING: time$year and time$doy should have the same length.\n");
			goto __ENTRY_DATESTR_LOC;
		}
		int    Ntime = Nyr;
		F32PTR yr32  = malloc(sizeof(F32) * 2 * Ntime);;
		F32PTR doy32 = yr32 + Ntime;

		if (!CopyNumericObjToF32Arr(yr32, yr, Ntime)) {
			q_warning("WARNING: time$year has an unsupported data format or type.\n");
			free(yr32);
			goto __ENTRY_DATESTR_LOC;
		}
		if (!CopyNumericObjToF32Arr(doy32, doy, Ntime)) {
			q_warning("WARNING: time$doy has an unsupported data format or type.\n");
			free(yr32);
			goto __ENTRY_DATESTR_LOC;
		}
 
		for (int i = 0; i < Ntime; ++i) {
			yr32[i] = YDOYtoF32time( yr32[i], doy32[i] );
			if (yr32[i] < -1e9) {
				q_warning("WARNING: The (%d)-ith date ( time$year=%d,and  time$doy=%d) is not valid.\n", i + 1, (int)yr32[i], (int)doy32[i] );				
				free(yr32);
				goto __ENTRY_DATESTR_LOC;
			}
		}

		isTimeProcessed = 1;
		*Nactual = Ntime;
		return yr32;
	}



	VOIDPTR datestr, strfmt;

__ENTRY_DATESTR_LOC:

	 datestr = GetField123Check(timeObj, "dateStr", 3);
	 strfmt = GetField123Check(timeObj, "strFmt", 3);
	if (!isTimeProcessed && datestr && strfmt && IsChar(datestr) && IsChar(strfmt)  ) 	{		

		int    Ntime = GetNumberOfElements(datestr);
		F32PTR ftime = malloc(sizeof(F32) * Ntime);

		char STRFmt[355 + 1];
		GetCharArray(strfmt, STRFmt, 355);
		
		DateFmtPattern1 fmt1;	
		if (GetStrPattern_fmt1(STRFmt, &fmt1)) {			
			for (int i = 0; i < Ntime; ++i) {
				char TMP[355 + 1];
				if (!GetCharVecElem(datestr, i, TMP, 355)) {			
					q_warning("WARNING: Unable to read the %d-ith date string from time$dateStr.\n", i + 1);
					free(ftime);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
				ftime[i] = Str2F32time_fmt1(TMP, &fmt1);
				if (ftime[i] < -1e9) {
					q_warning("WARNING:  The %d-th string (time$dateStr=\"%s\") is invalid, incompatiable with the specified "
							" time$strFmat=\"%s\".\n", i + 1, TMP, STRFmt);
					free(ftime);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
			}	
		}
	 
		DateFmtPattern2 fmt2;	
		if (GetStrPattern_fmt2(STRFmt, &fmt2)) {
			
			for (int i = 0; i < Ntime; ++i) {
				char TMP[255 + 1];
				if (!GetCharVecElem(datestr, i, TMP, 255)) {
					q_warning("WARNING: Unable to read the %d-ith date string from time$dateStr.\n", i + 1);
					free(ftime);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
				ftime[i] = Str2F32time_fmt2(TMP, &fmt2);
				if (ftime[i] < -1e9) {
					q_warning("WARNING: The %d-th string (time$dateStr=\"%s\") is invalid, incompatiable with the specified time$strFmat=\"%s\".\n", i + 1, TMP, STRFmt);
					free(ftime);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
			}	
		}


		DateFmtPattern3 fmt3;	
		if (GetStrPattern_fmt3(STRFmt, &fmt3)) {
			
			for (int i = 0; i < Ntime; ++i) {
				char TMP[255 + 1];
				if (!GetCharVecElem(datestr, i, TMP, 255)) {
					q_warning("WARNING: Unable to read the %d-ith date string from time$dateStr.\n", i + 1);
					free(ftime);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
				ftime[i] = Str2F32time_fmt3(TMP, &fmt3);
				if (ftime[i] < -1e9) {
					q_warning("WARNING: The %d-th string (time$dateStr=\"%s\") is invalid, incompatiable with the specified time$strFmat=\"%s\".\n", i + 1, TMP, STRFmt);
					free(ftime);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
			}	
		}

		isTimeProcessed = 1;
		*Nactual = Ntime;
		return ftime;

	}
	 
__ENTRY_DATESTR_SMART_LOC:
	if (!isTimeProcessed && datestr &&  IsChar(datestr)  ) {
		DynMemBuf        strarr = { .max_len = 100 };
		DynAlignedBuf    start  = { .align = 4,.elem_size = 4,.p.raw = NULL, };
		DynAlignedBuf    nchars = { .align = 4,.elem_size = 4,.p.raw = NULL, };

		CharObj2CharArr(timeObj, &strarr, &start, &nchars);
		float* ftime = strings_to_fyears(strarr.raw, start.p.i32, start.cur_len);
		// ftime is NULL upon return if failled

		int Ntime = start.cur_len;
		dynbuf_kill(&strarr);
		adynbuf_kill(&start);
		adynbuf_kill(&nchars);
		*Nactual = Ntime;
		return ftime;
	}

	return NULL;
}

void CharObj2CharArr(VOID_PTR o, DynMemBufPtr str , DynAlignedBufPtr charstart, DynAlignedBufPtr nchars){
	int n = GetNumberOfElements(o);
	dynbuf_init(str, n * 200);
	adynbuf_init(charstart, n);
	adynbuf_init(nchars, n);
	for (int i = 0; i < n; i++) {
		dynbuf_requestmore(str, 200);
		int  nlen = GetCharVecElem(o, i, str->raw+str->cur_len, 200);
		charstart->p.i32[i] = str->cur_len;
		nchars->p.i32[i]    = nlen;
		charstart->cur_len++;
		nchars->cur_len++;
		str->cur_len += (nlen+1); // 1L is the extra NULL.
	}	 
}

void obj_to_str(VOID_PTR o, DynMemBufPtr s, int leftMargin) {
	int nfld      = GetNumberOfFields(o);

	int maxKeyLen = 0;
	for (int i = 0; i < nfld; ++i) {
		char key[50];
		GetFieldNameByIdx(o, i, key, 50);
		int  keyLen = strlen(key);
		if (keyLen > maxKeyLen) {
			maxKeyLen = keyLen;
		}
	} 

	for (int i = 0; i < nfld; ++i) {
		char key[50];
		GetFieldNameByIdx(o, i, key, 50);
		int  keyLen = strlen(key);	    
		char tmp[200];
		snprintf(tmp, 199, "%*s%-*.*s : ", leftMargin, "", maxKeyLen, maxKeyLen, key);
		dynbuf_insert_bytes(s, tmp, strlen(tmp)+1);

		VOIDPTR e  = GetFieldByIdx(o, i);
		if (e == NULL) {
			snprintf(tmp, 199, "[]\n");
			dynbuf_insert_bytes(s, tmp, strlen(tmp) + 1);
		}
		else if (IsNumeric(e)) {
			int size  = GetNumberOfElements(e);
			int dtype = GetDataType(e);
			if (size == 1) {
				F64 value = GetScalar(e);				
				tmp[0] = 0;
				if (dtype == DATA_INT16 || dtype == DATA_INT32 || dtype == DATA_INT64) {
					snprintf(tmp,199, "%d\n",(int) value);
				} else if (dtype == DATA_FLOAT || dtype == DATA_DOUBLE) {
					snprintf(tmp,199, "%g\n", value);
				}	else {
					snprintf(tmp,199, "%g\n", value);
			
				}
				dynbuf_insert_bytes(s, tmp, strlen(tmp) + 1);
				continue;
			}

			char* etype;
			if      (dtype == DATA_INT16)  etype = "int16";
			else if (dtype == DATA_INT32)  etype = "int32";
			else if (dtype == DATA_INT64)  etype = "int64";
			else if (dtype == DATA_FLOAT)  etype = "float32";
			else if (dtype == DATA_DOUBLE)  etype = "float64";
			else etype = "others";

			int ndims = GetNumOfDim(e);
			int dims[10];			
			GetDimensions(e, dims, ndims);

			snprintf(tmp, 199, "[%d", dims[0]);
			dynbuf_insert_bytes(s, tmp, strlen(tmp) + 1);
			for (int j = 1; j < ndims; j++) {
				snprintf(tmp, 199, "x%d", dims[j]);
				dynbuf_insert_bytes(s, tmp, strlen(tmp) + 1);
			}
			snprintf(tmp, 199," %s] \n", etype);
			dynbuf_insert_bytes(s, tmp, strlen(tmp) + 1);
		}
		else if (IsChar(e)) {
			char tmpstr[30];
			GetCharArray(e, tmpstr, 30);
			snprintf(tmp, 199,"'%s'\n", tmpstr);
			dynbuf_insert_bytes(s, tmp, strlen(tmp) + 1);
		}
		else if (IsStruct(e)) {
			int len = GetNumberOfFields(e);
			snprintf(tmp, 199, "[ 1 object with %d fields] \n", len);
			dynbuf_insert_bytes(s, tmp, strlen(tmp) + 1);
			obj_to_str(e, s, leftMargin + maxKeyLen + 3);
		}
	}
}

int IDEPrintObject(VOID_PTR o) {
	if (!IsStruct(o)) {
		r_printf("Not an object, structure, or list.\n");
	}

	int nfld = GetNumberOfFields(o);
	r_printf("Object of %d field(s): \n\n", nfld);

	DynMemBuf s = { 0, };
	dynbuf_init(&s, 1000);
	int leftMargin = 1;
	obj_to_str(o, &s, leftMargin);
	s.raw[s.cur_len] = 0;
	r_printf("%s", s.raw);
	dynbuf_kill(&s);
	return 0;
}

#include "abc_000_warning.h"