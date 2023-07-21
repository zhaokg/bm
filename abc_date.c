#include<stdio.h>
#include<inttypes.h>
// include stdio.h library
#include <stdlib.h> // atof atoi
#include <string.h> //strchr strrchr strstr str

#include "abc_000_warning.h"

#include "abc_date.h"
#include "abc_common.h"   //ToUpper
#include "abc_ide_util.h" //ToUpper
#include "abc_vec.h" //ToUpper
#include "abc_sort.h" // insert_sort
 

//https://overiq.com/c-examples/c-program-to-calculate-the-day-of-year-from-the-date/
static int IsLeapYear(int year) { return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);}
static int GetNumDays(int year) { return IsLeapYear(year) ? 366 : 365; }
//https://stackoverflow.com/questions/19377396/c-get-day-of-year-from-date
static const int DAYS[2][13] = {
	{ 0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 },
	{ 0, 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 }
};

int Date2Doy(int year, int mon, int day) { 	return DAYS[IsLeapYear(year)][mon] + day; }
int Doy2Date(int doy,  int y, int* d, int* m){

		static int month[13] = { 0, 31, 28, 31, 30, 31, 30,	31, 31, 30, 31, 30, 31 };

		if (IsLeapYear(y))	month[2] = 29;

		int i;
		for (i = 1; i <= 12; i++) 		{
			if (doy <= month[i])
				break;
			doy = doy - month[i];
		}

		*d = doy;
		*m = i;

		return 0;
}


int IsValidDate(int year, int mon, int day )
{
	int is_valid = 1, is_leap = 0;

	if (year >= 1800 && year <= 9999)	{

		//  check whether year is a leap year
		is_leap = (year%4 == 0 && year%100 != 0) || (year % 400 == 0);
		
		// check whether mon is between 1 and 12
		if (mon >= 1 && mon <= 12) 		{
			// check for days in feb
			if (mon == 2) 	{
				if (is_leap && day == 29){
					is_valid = 1;
				} 
				else if (day > 28)
				{
					is_valid = 0;
				}
			}
			else if (mon == 4 || mon == 6 || mon == 9 || mon == 11)
			{	// check for days in April, June, September and November
				if (day > 30) {
					is_valid = 0;
				}
			}
			else if (day > 31)
			{  // check for days in rest of the months:Jan, Mar, May, July, Aug, Oct, Dec
				is_valid = 0;
			}
		}
		else
		{
			is_valid = 0;
		}

	}
	else
	{
		is_valid = 0;
	}

	return is_valid;

}

//codereview.stackexchange.com/questions/152017/simple-days-between-dates-calculator
int64_t CountLeapYears(int64_t year)
{
	// 438 - 17 + 4
	static int64_t fakeLeaps = (1753 / 4) - (1753 / 100) + (1753 / 400);
	// Leaps before 1753
	// We start at 0 in 1753
	int64_t leaps     = year / 4;
	int64_t badLeaps  = year / 100;
	int64_t extraLeaps = year / 400;

	return leaps - badLeaps + extraLeaps - fakeLeaps;
}

float YDOYtoF32time(int year, int doy) {
	return (float)year + ((float)doy - 0.5f) / (float) GetNumDays(year);
}


float YMDtoF32time(int year, int mon, int day) { 
	return YDOYtoF32time(year, Date2Doy(year, mon, day)); 
}

#include "math.h"

int F32time2YDOY(F32 fyear, int* doy) {
	int yr  = fyear;
	*doy    = (int) round( (fyear - yr) * GetNumDays(yr) + 0.5f);
	return yr;
}

int F32time2YMD(F32 fyear, int*mon, int *day) {
	int doy;
	int yr = F32time2YDOY(fyear, &doy);	 
	Doy2Date(doy, yr, day, mon);
	return yr;
}


int64_t datenum(int year, int mon, int day) {
	int64_t numYears   = year - 1753;
	int64_t leapYears  = CountLeapYears(year);
	int64_t daysInYear = DAYS[IsLeapYear(year)][mon];
	return  numYears * 365LL+ leapYears + daysInYear + (day - 1);
}
//https://overiq.com/c-examples/c-program-to-calculate-the-day-of-year-from-the-date/

//https://stackoverflow.com/questions/14218894/number-of-days-between-two-dates-c
//https://stackoverflow.com/questions/33712685/c-days-between-given-date-and-today%C2%B4s-date
// Civil calendar is the Gregorian calendar: the one we normally use
// Returns number of days since civil 1970-01-01.  Negative values indicate
//    days prior to 1970-01-01.
// Preconditions:  y-m-d represents a date in the civil (Gregorian) calendar
//                 m is in [1, 12]
//                 d is in [1, last_day_of_month(y, m)]
//                 y is "approximately" in
//                   [numeric_limits<Int>::min()/366, numeric_limits<Int>::max()/366]
//                 Exact range of validity is:
//                 [civil_from_days(numeric_limits<Int>::min()),
//                  civil_from_days(numeric_limits<Int>::max()-719468)]
int days_from_civil(int y, unsigned m, unsigned d) {	
	y -= m <= 2;
	const int      era = (y >= 0 ? y : y - 399) / 400;
	const unsigned yoe = (y - era * 400);      // [0, 399]
	const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;  // [0, 365]
	const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;         // [0, 146096]
	return era * 146097 +doe - 719468;
}
//http://howardhinnant.github.io/date_algorithms.html
int civil_from_days(int days, int * yr, int*mn, int* day) 
{ 
	days += 719468;
	const int era = (days >= 0 ? days : days - 146096) / 146097;
	const unsigned doe = (days - era * 146097);          // [0, 146096]
	const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;  // [0, 399]
	const int      y  = (yoe) + era * 400;
	const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);                // [0, 365]
	const unsigned mp  = (5 * doy + 2) / 153;                                   // [0, 11]
	const unsigned d  = doy - (153 * mp + 2) / 5 + 1;                             // [1, 31]
	const unsigned m  = mp < 10 ? mp + 3 : mp - 9;                            // [1, 12]
	*yr = y + (m <= 2);
	*mn = m;
	*day = d;
 

	return 0;
}


int F32time2DateNum(F32 fyear) {
	int mon, day;
	int yr = F32time2YMD(fyear, &mon, &day);
	return days_from_civil(yr, mon, day);
}


// //In R, dates are represented as the number of days since 1970-01-01
float fractional_civil_from_days(int days)
{
	days += 719468;
	const int era = (days >= 0 ? days : days - 146096) / 146097;
	const unsigned doe = (days - era * 146097);          // [0, 146096]
	const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;  // [0, 399]
	const int      y = (yoe)+era * 400;
	const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);                // [0, 365]
	const unsigned mp = (5 * doy + 2) / 153;                                   // [0, 11]
	const unsigned d = doy - (153 * mp + 2) / 5 + 1;                             // [1, 31]
	const unsigned m = mp < 10 ? mp + 3 : mp - 9;                            // [1, 12]
	int yr = y + (m <= 2);
	int mn = m;
	 int day = d;
	 
	 return YMDtoF32time(yr, mn, day);

}

void date_jump(int y, int m, int d, int jumpDays, int* y1, int* m1, int* d1) {
	int days=days_from_civil(y, m, d);
	civil_from_days(days + jumpDays, y1, m1, d1);
}


static int __FindPatternStart( char *str, char * token) {

	char * pchar = strstr(str, token);
	if (pchar)	{
		return (int)(pchar - str);
	}
	return -10000;
}

int    GetStrPattern_fmt1(char* fmtstr, DateFmtPattern1* pattern) {
	ToUpper(fmtstr);

	int yearIdx = __FindPatternStart(fmtstr, "YYYY");
	if (yearIdx <0) return 0;
	int monIdx  = __FindPatternStart(fmtstr, "MM");
	if (monIdx < 0) return 0;
	int dayIdx  = __FindPatternStart(fmtstr, "DD");
	if (dayIdx < 0) return 0;

	pattern->yearIdx = yearIdx;
	pattern->monIdx  = monIdx;
	pattern->dayIdx  = dayIdx;
	return 1;
}

float  Str2F32time_fmt1(char* datestr, DateFmtPattern1* pattern) {

	char s[5];
	memcpy(s, datestr + pattern->yearIdx, 4);	s[4] = 0;	int year = atoi(s);
	if (year == 0) { return -1e10; }

	memcpy(s, datestr + pattern->monIdx, 2); 	s[2] = 0;	int mon = atoi(s);
	if (mon < 1 || mon > 12) { return -1e10; }

	memcpy(s, datestr + pattern->dayIdx, 2);  	s[2] = 0;	int day = atoi(s);
	if (day < 1 || day >31) { return -1e10; }

	return YMDtoF32time(year, mon, day);
}

int    GetStrPattern_fmt2(char* fmtstr, DateFmtPattern2* pattern) {
	ToUpper(fmtstr);
	int yearIdx = __FindPatternStart(fmtstr, "YYYY");
	if (yearIdx <0) return 0;
	int doyIdx  = __FindPatternStart(fmtstr, "DOY");
	if (doyIdx < 0) return 0;
	
	pattern->yearIdx = yearIdx;
	pattern->doyIdx  = doyIdx;
	return 1;
}

float  Str2F32time_fmt2(char* datestr, DateFmtPattern2* pattern) {

	char s[5];
	memcpy(s, datestr + pattern->yearIdx, 4);	s[4] = 0;	int year = atoi(s);
	if (year == 0) { return -1e10; }

	memcpy(s, datestr + pattern->doyIdx, 3); 	s[3] = 0;	int doy = atoi(s);	
	if (doy < 0 || doy > 366) { return -1e10; }	 

	return YDOYtoF32time(year, doy);
}

 
static char* _FindCharOccurrence(char* s, char c, int* numTimes) {
	
	*numTimes = 0;	 	
	char *pLast =NULL;
	while (( s=strchr(s, c)) != NULL ) {
		pLast = s++;
		++*numTimes;				
	}	
	return pLast;
	
}

int    GetStrPattern_fmt3(char* fmtstr, DateFmtPattern3* pattern) {
    
	ToUpper(fmtstr);
	int nTimes;
	char* yearPt = _FindCharOccurrence(fmtstr, 'Y', &nTimes);
	if (nTimes == 0 || nTimes > 1) return 0;

	char* monPt = _FindCharOccurrence(fmtstr, 'M', &nTimes);
	if (nTimes == 0 || nTimes > 1) return 0;

	char* dayPt = _FindCharOccurrence(fmtstr, 'D', &nTimes);
	if (nTimes == 0 || nTimes > 1) return 0;

	pattern->order[0] = 'Y';
	pattern->order[1] = 'M';
	pattern->order[2] = 'D';

	char* pts[] = { yearPt, monPt, dayPt };
	InsertionSort_VOIDPTR(pts, pattern->order, 3);

	int64_t len;
	len =  (pts[1] - 1) - (pts[0] + 1) + 1;
	if (len <= 0) return 0;
	memcpy(pattern->sep1, pts[0] + 1, len); 	pattern->sep1[len] = 0;


	len =  (pts[2] - 1) - (pts[1] + 1) + 1 ;
	if (len <= 0) return 0;
	memcpy(pattern->sep2, pts[1] + 1, len);  	pattern->sep2[len] = 0;
	 
	return 1;
}


float  Str2F32time_fmt3(char* datestr, DateFmtPattern3* pattern) {

	int   N = (int) strlen(datestr);
	char  old;

	char* p0 = datestr;
	char *p1 = strstr(p0, pattern->sep1);
	if (p1 == NULL) return -1e10;
	old =p1[0]; p1[0] = 0;
	int n1=atoi(p0);
	p1[0] = old;

	p0 = p1 + strlen(pattern->sep1);
	p1 = strstr(p0, pattern->sep2);
	if (p1 == NULL) return -1e10;
	old = p1[0]; p1[0] = 0;
	int n2 = atoi(p0);
	p1[0] = old;

	p0 = p1 + strlen(pattern->sep2);
	if (p0 >= datestr+N) return -1e10;
	int n3 = atoi(p0);

	char* p = pattern->order;
	int year = p[0] == 'Y' ? n1 : (p[1] == 'Y' ? n2 : n3);
	int mon  = p[0] == 'M' ? n1 : (p[1] == 'M' ? n2 : n3);
	int day  = p[0] == 'D' ? n1 : (p[1] == 'D' ? n2 : n3);
	 
	return YMDtoF32time(year, mon, day);
}



INLINE static int is_dot(char c)          { return c == '.'; }
INLINE static int is_slash(char c)        { return c == '/'; }
INLINE static int is_digit(char c)        {  return c >= '0' && c <= '9';}
INLINE static int is_letter(char c)       {  return (c >= 'a' && c <= 'z' ) || (c >= 'A' && c <= 'Z');}
INLINE static int is_alphanumeric(char c) {  return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }


int get_word_size(char* s)         { int i = 0;	while (is_letter(s[i++])){};  	 return --i; }
int get_alphanumeric_size(char* s) { int i = 0;	while (is_alphanumeric(s[i++])) {}; return --i;}
int get_intger_size(char* s)       { int i = 0; while (is_digit(s[i++])) {};	 	 return --i; }
int get_slash_size(char* s)       { int i = 0; while (is_slash(s[i++])) {};	 	 return --i; }
int get_number_size(char* s, int * ndots) {
	int i = *ndots=0;
	while (is_digit(s[i]) || is_dot(s[i])) {
		*ndots += is_dot(s[i]);
		i++;
	};
	return i;
}
char* goto_validchar(char* s) { while (!is_alphanumeric(*s) && *s != 0) { s++; }	return s; }
char* goto_validchar_dot_slash(char* s) {while (!is_alphanumeric(*s) && !is_dot(*s) && !is_slash(*s)  && *s != 0) { s++; }	return s;}

static int split_numstr(char* s, int nPartMax, int* startidx, int* nchar, char* type) {

	char* s0 = s;
	int   nPart = 0;
	while (*s != 0 && nPart < nPartMax) {
		s = goto_validchar_dot_slash(s);
		if (is_digit(*s) || is_dot(*s) ) {
			int ndots;
			int nlen = get_number_size(s, &ndots);
			if (ndots >= 2) {
				return -1L;
			}
			nchar[nPart]    = nlen;
			startidx[nPart] = s - s0;
			type[nPart]     = 'N';
			nPart++;
			s = s + nlen;
		}
		else if (is_letter(*s)) {
			int nlen = get_word_size(s);
			nchar[nPart] = nlen;
			startidx[nPart] = s - s0;
			type[nPart] = 'W';
			nPart++;
			s = s + nlen;
		}
		else if (is_slash(*s)) {
			int nlen = get_slash_size(s);
			nchar[nPart] = nlen;
			startidx[nPart] = s - s0;
			type[nPart] = 'S';
			nPart++;
			s = s + nlen;
		}
	}
	return nPart;
}
INLINE static char char_toupper(char s) { return s >= 'a' && s <= 'z' ? s-32 : s; }

double extract_fyear(char* s) {

	int  nPartMax = 4;
	int  startIdx[4], nchar[4];
	char type[4] = { 0, };
	int  nPart = split_numstr(s, nPartMax, startIdx, nchar, type);

	double nan = (1.0 / 0.) * 0.;
	if (nPart <= 0 || type[0] != 'N' || type[nPart-1] !='W') {
		return nan;
	}

	double x[4];
	int   nNumber = 0;
	for (int i = 0; i < nPartMax; i++){
		if (type[i] != 'N') continue;

		char* ss = s + startIdx[i];
		int   len = nchar[i];
		char  old = ss[len]; ss[len] = 0; x[nNumber++] = atof(ss); ss[len] = old;
	}

	char* ss     = s + startIdx[nPart-1];
	char  letter = char_toupper(ss[0]);

	double y = nan;
	if (nPart==2 ) {       
		y = x[0];
	}
	else if (nPart == 4 && type[1] == 'S' && type[2] == 'N'){
		 y = x[0] / x[1];	
	}

#define _IsAlmostInteger(x)  ( fabs(x-round(x)) <1e-5 )

	double z = nan;
	if (letter == 'D') {
		double nyr1 = y / 366;
		double nyr2 = y / 365;
		if      (_IsAlmostInteger(nyr1) ) 	z = nyr1;
		else if (_IsAlmostInteger(nyr2)) 	z = nyr2;
		else  	z = y/365 ;
	}
	else if (letter == 'M') {
		z = y / 12;	  
	}
	else if (letter == 'Y') {
		z = y  ;
	}
	return z;
}
 

int split_datestr(char* s, int nPartMax, int* startidx, int* nchar, char* type) {

	char* s0    = s;
	int   nPart = 0;
	while (*s != 0 && nPart < nPartMax) {
		s = goto_validchar(s);
		if (is_digit(*s)) {
			int nlen         = get_intger_size(s);
			nchar[nPart]     = nlen;
			startidx[nPart]  = s - s0;
			char typePart    = 'N';
			if (startidx[nPart] > 0 && is_letter(s0[startidx[nPart]-1])) {
				typePart = 'A';
			}
			if (is_letter(s[nlen])) {
				typePart = 'A';
			}
			type[nPart] = typePart;
			nPart++;
			s = s + nlen;
		}
		else if (is_letter(*s)) {
			int nlen        = get_word_size(s);
			nchar[nPart]    = nlen;
			startidx[nPart] = s - s0;
			type[nPart]    = 'L';
			nPart++;
			s = s + nlen;
		}
	}
	return nPart;
}

static int cmp_months(char* s) {
	static char* months[] = { "Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct", "Nov","Dec"};
	for (int i = 0; i < 12; ++i) {
		int diff=strcicmp_nfirst(s, months[i], 3);
		if (diff == 0) {
			return i + 1;
		}
	}
	return -1;
}
float * strings_to_fyears(char *s, int * strstart, int n) {

    #define NPARTMAX 16
	 	
	int* startidx    = malloc(sizeof(int) * ( n * NPARTMAX * 2  ));
	int* partlength  = startidx + n * NPARTMAX;
		
	int  i     = 0; 
	char parttype[NPARTMAX] = { 0 };
	int  nPart = split_datestr(s + strstart[i], NPARTMAX, startidx + i * NPARTMAX, partlength + i * NPARTMAX,  parttype);
	for (i = 1; i < n; i++) {
		char newparttype[16] = { 0 };
		int  newnPart = split_datestr(s + strstart[i], NPARTMAX, startidx + i * NPARTMAX, partlength + i * NPARTMAX, newparttype);
		if (nPart != newnPart || memcmp(parttype, newparttype, NPARTMAX) != 0) {
			free(startidx);
			r_printf("ERROR: the input date strings have inconsisent formats and cann't be automatically parsed. Use time$datestr and time$strfmat to specify the format.\n");
			return NULL;
		}			
	}

	int nNumber=0, nWord = 0, nANumber = 0, nNumberANumber = 0, nNumber8=0, nNumber7=0;
	int idxNumber[NPARTMAX], idxWord[NPARTMAX], idxANumber[NPARTMAX], idxNumANumber[NPARTMAX], idxNumber8[NPARTMAX], idxNumber7[NPARTMAX];
	int partLenMin[NPARTMAX], partLenMax[NPARTMAX];
	for (int i = 0; i < nPart; i++) {
		char pattern = parttype[i];
		if (pattern == 'N') {
			idxNumber[nNumber++]            = i;
			idxNumANumber[nNumberANumber++] = i;
			
		}	else if (pattern == 'A') {
			idxANumber[nANumber++]           = i;
			idxNumANumber[nNumberANumber++ ] = i;
		
		}else if (pattern == 'L') {
			idxWord[nWord++] = i;
		}

		int* partLengthRow = partlength + i;
		partLenMin[i] = *partLengthRow;
		partLenMax[i] = *partLengthRow;
		for (int j = 0; j < n; j++) {
			partLenMin[i] = min(partLenMin[i], *partLengthRow);
			partLenMax[i] = max(partLenMax[i], *partLengthRow);	
			partLengthRow += NPARTMAX;
		}

		if (partLenMin[i] == 8 && partLenMax[i] == 8 && pattern !='L') {
			idxNumber8[nNumber8++] = i;	 
		}
		if (partLenMin[i] == 7 && partLenMax[i] == 7 && pattern != 'L') {
			idxNumber7[nNumber7++] = i;
		}
	}


	int* year  = malloc(sizeof(int) * n * 4);
	int* month = year + n;
	int* day   = month + n;
	int* tmp   = day + n;

	int DONE = 0;

	// 1991,12,1
	if (nNumberANumber ==3 || nNumber ==3) {
	
		int yearFound=0, monthFound = 0, dayFound = 0;			 
		for (int J = 0; J < 3; ++J) {
			int   minv  = 999999, maxv = -999999;
			F64   meanv = 0;
			int  idx    = nNumberANumber == 3 ? idxNumANumber[J] : idxNumber[J];
			for (int i = 0; i < n; ++i) {
				int   sidx = *(startidx   + i * NPARTMAX + idx);
				int   slen = *(partlength + i * NPARTMAX + idx);
				char* ss   = s + strstart[i] + sidx;
				char  old  = ss[slen]; ss[slen] = 0; int value = atoi(ss); ss[slen] = old;
				meanv += value;
				if (minv > value) minv = value;
				if (maxv < value) maxv = value;
				tmp[i] = value;
			}

			if (partLenMin[idx] == 4 && partLenMax[idx] == 4) {
				if (yearFound == 0) {
					yearFound = 1;
					memcpy(year, tmp, sizeof(int) * n);
					continue;
				}
			}
			if ( (maxv==12 || maxv == 11) && partLenMax[idx] < 4) {
				if (monthFound == 0) {
					monthFound = 1;
					memcpy(month, tmp, sizeof(int) * n);
					continue;
				}
			}
			if ( (maxv >=28 && maxv <=31 ) && partLenMax[idx] < 4) {
				if (dayFound == 0) {
					dayFound = 1;
					memcpy(day, tmp, sizeof(int) * n);
					continue;
				}
			}

			if (minv < 1 || maxv > 12) {
				// this must be Year or Month
				if (minv < 1 || maxv > 31) {
					// this must be year
					if (yearFound == 0) {
						yearFound = 1;
						memcpy(year, tmp, sizeof(int) * n);
					}
				}	else {
					//this should be days
					if (dayFound == 0) {
						dayFound = 1;
						memcpy(day, tmp, sizeof(int) * n);
					}
				}
			} else {
				// this shouldbe Month			 
				if (monthFound == 0) {
					monthFound = 1;
					memcpy(month, tmp, sizeof(int) * n);
				}
			}
		}

		if (yearFound == 1 && monthFound == 1 && dayFound == 1) {
			DONE = 1;
		}
	}
	if (nNumberANumber ==2 || nNumber ==2 ) {
	  // 1899,1
		int yearFound=0, monthFound = 0, dayFound = 0;			 
		for (int J = 0; J < 2; ++J) {
			int   minv  = 999999, maxv = -999999;
			F64   meanv = 0;
			int   idx    = nNumberANumber == 2 ? idxNumANumber[J] : idxNumber[J];
			for (int i = 0; i < n; ++i) {
				int   sidx = *(startidx   + i * NPARTMAX + idx);
				int   slen = *(partlength + i * NPARTMAX + idx);
				char* ss   = s + strstart[i] + sidx;
				char  old  = ss[slen]; ss[slen] = 0; int value = atoi(ss); ss[slen] = old;
				meanv += value;
				if (minv > value) minv = value;
				if (maxv < value) maxv = value;
				tmp[i] = value;
			}

			if (partLenMin[idx] == 4 && partLenMax[idx] == 4) {
				if (yearFound == 0) {
					yearFound = 1;
					memcpy(year, tmp, sizeof(int) * n);
					continue;
				}
			}
			if ((maxv == 12 || maxv == 11) && partLenMax[idx] <=2) {
				if (monthFound == 0) {
					monthFound = 1;
					memcpy(month, tmp, sizeof(int) * n);
					continue;
				}
			}

			if (maxv <=12 && minv >= 1 && partLenMax[idx] <= 2) {
				if (monthFound == 0) {
					monthFound = 1;
					memcpy(month, tmp, sizeof(int) * n);
				}
			} 
		}

		if (yearFound == 1 && monthFound == 1 ) {
			DONE = 3;
		}
	}


	
	if (!DONE && (nNumberANumber == 2 || nNumber == 2) && nWord == 1) {
		// 1991,Mar,1
		int yearFound=0, monthFound = 0, dayFound = 0;			 
		for (int J = 0; J < 2; ++J) {
			int   minv  = 999999, maxv = -999999;
			F64   meanv = 0;
			int   idx   = nNumberANumber == 2 ? idxNumANumber[J] : idxNumber[J];
			for (int i = 0; i < n; ++i) {
				int   sidx = *(startidx   + i * NPARTMAX + idx);
				int   slen = *(partlength + i * NPARTMAX + idx);
				char* ss   = s + strstart[i] + sidx;
				char  old  = ss[slen]; ss[slen] = 0; int value = atoi(ss); ss[slen] = old;
				meanv += value;
				if (minv > value) minv = value;
				if (maxv < value) maxv = value;
				tmp[i] = value;
			}

			if (minv >= 1 && maxv <= 31 && partLenMin[idx]<4) {
				//this should be days
				if (dayFound == 0) {
					dayFound = 1;
					memcpy(day, tmp, sizeof(int) * n);
				}

			} 
			else if (partLenMin[idx] >= 4) {
				// this should be year
				if (yearFound == 0) {
					yearFound = 1;
					memcpy(year, tmp, sizeof(int) * n);
				}
			}		 
		} // for (int J = 0; J < 2; ++J) 

		int allMatched = 1;
		int idx = idxWord[0];
		for (int i = 0; i < n; ++i) {
			int   sidx = *(startidx + i * NPARTMAX + idx);
			int   slen = *(partlength + i * NPARTMAX + idx);
			char* ss = s + strstart[i] + sidx;
			char  old = ss[slen]; ss[slen] = 0; int value = cmp_months(ss); ss[slen] = old;
			tmp[i] = value;
			if (value < 0) {
				allMatched = 0;
				break;
			}
			
		}
		if (allMatched) {
			monthFound = 1;
			memcpy(month, tmp, sizeof(int) * n);
		}

		if (yearFound == 1 && monthFound == 1 && dayFound == 1) {
			DONE = 1;
		}
 
	}
	if (!DONE && (nNumberANumber == 1) && nWord == 1) {
		// 1991, Marc
		int yearFound=0, monthFound = 0, dayFound = 0;			 
		for (int J = 0; J < 1; ++J) {
			int   minv  = 999999, maxv = -999999;
			F64   meanv = 0;
			int   idx   = idxNumANumber[J];
			for (int i = 0; i < n; ++i) {
				int   sidx = *(startidx   + i * NPARTMAX + idx);
				int   slen = *(partlength + i * NPARTMAX + idx);
				char* ss   = s + strstart[i] + sidx;
				char  old  = ss[slen]; ss[slen] = 0; int value = atoi(ss); ss[slen] = old;
				meanv += value;
				if (minv > value) minv = value;
				if (maxv < value) maxv = value;
				tmp[i] = value;
			}
			if (partLenMin[idx] == 4 && partLenMax[idx] == 4) {
				if (yearFound == 0) {
					yearFound = 1;
					memcpy(year, tmp, sizeof(int) * n);
					continue;
				}
			}
		} // for (int J = 0; J < 2; ++J) 

		int allMatched = 1;
		int idx = idxWord[0];
		for (int i = 0; i < n; ++i) {
			int   sidx = *(startidx + i * NPARTMAX + idx);
			int   slen = *(partlength + i * NPARTMAX + idx);
			char* ss = s + strstart[i] + sidx;
			char  old = ss[slen]; ss[slen] = 0; int value = cmp_months(ss); ss[slen] = old;
			tmp[i] = value;
			if (value < 0) {
				allMatched = 0;
				break;
			}
			
		}
		if (allMatched) {
			monthFound = 1;
			memcpy(month, tmp, sizeof(int) * n);
		}

		if (yearFound == 1 && monthFound == 1  ) {
			DONE = 3;
		}
 
	}
 
	if (!DONE && nNumber8 > 0) {
		// 19910301
		for (int J = 0; J < nNumber8; ++J) {		 
			int  idx    =   idxNumber8[J];
			int  zero   =  '0';
			for (int i = 0; i < n; ++i) {
				int   sidx = *(startidx   + i * NPARTMAX + idx);
				int   slen = *(partlength + i * NPARTMAX + idx);
				char* ss   = s + strstart[i] + sidx;
				year[i]  = (ss[0] - zero) * 1000 + (ss[1] - zero) * 100 + (ss[2] - zero) * 10 + (ss[3] - zero);
				month[i] = (ss[4] - zero) * 10 + (ss[5] - zero) * 1;
				day[i]   = (ss[6] - zero) * 10 + (ss[7] - zero) * 1;
				
			}
			int minv, maxv;
			i32_maxidx(year, n, &maxv);
			i32_minidx(year, n, &minv);
			if (maxv > 8000) continue;
			i32_maxidx(month, n, &maxv);
			i32_minidx(month, n, &minv);
			if (minv < 1 || maxv>12) {
				continue;
			}
			i32_maxidx(day, n, &maxv);
			i32_minidx(day, n, &minv);
			if (minv < 1 || maxv>31) {
				continue;
			}

			DONE = 1;
			break;		 
		}
	}
	if (!DONE && nNumber8 > 0) {
		// 03011991
		for (int J = 0; J < nNumber8; ++J) {
			int  idx = idxNumber8[J];
			int  zero = '0';
			for (int i = 0; i < n; ++i) {
				int   sidx = *(startidx + i * NPARTMAX + idx);
				int   slen = *(partlength + i * NPARTMAX + idx);
				char* ss = s + strstart[i] + sidx;
				year[i] = (ss[0] - zero) * 1000 + (ss[1] - zero) * 100 + (ss[2] - zero) * 10 + (ss[3] - zero);
				day[i] = (ss[4] - zero) * 10 + (ss[5] - zero) * 1;
				month[i] = (ss[6] - zero) * 10 + (ss[7] - zero) * 1;

			}
			int minv, maxv;
			i32_maxidx(year, n, &maxv);
			i32_minidx(year, n, &minv);
			if (maxv > 8000) continue;
			i32_maxidx(month, n, &maxv);
			i32_minidx(month, n, &minv);
			if (minv < 1 || maxv>12) {
				continue;
			}
			i32_maxidx(day, n, &maxv);
			i32_minidx(day, n, &minv);
			if (minv < 1 || maxv>31) {
				continue;
			}

			DONE = 1;
			break;
		}
	}
	if (!DONE && nNumber7 > 0) {
		// 1991231
		for (int J = 0; J < nNumber7; ++J) {
			int  idx    =   idxNumber7[J];
			int  zero   =  '0';
			for (int i = 0; i < n; ++i) {
				int   sidx = *(startidx   + i * NPARTMAX + idx);
				int   slen = *(partlength + i * NPARTMAX + idx);
				char* ss   = s + strstart[i] + sidx;
				year[i]  = (ss[0] - zero) * 1000 + (ss[1] - zero) * 100 + (ss[2] - zero) * 10 + (ss[3] - zero);
				day[i]   = (ss[4] - zero) * 100  + (ss[5] - zero) * 10  + (ss[6] - zero) * 1;
				
			}
			int minv, maxv;
			i32_maxidx(year, n, &maxv);
			i32_minidx(year, n, &minv);
			if (maxv > 4000) continue;
			
			i32_maxidx(day, n, &maxv);
			i32_minidx(day, n, &minv);
			if (minv < 1 || maxv>366) {
				continue;
			}

			DONE = 2;
			break;		 
		}
	}

	F32PTR out = NULL;

	if (DONE) {
		out = malloc(sizeof(F32) * n);
		if (DONE == 1) {
			r_printf("INFO: '%s' interpreted as %04d-%02d-%02d (Y-M-D)\n", s, year[0], month[0], day[0]);
			for (int i = 0; i < n; i++) {
				out[i] = YMDtoF32time(year[i], month[i], day[i]);
				//r_printf("%d %d %d  \n", year[i], month[i], day[i]);
			}
		} else if (DONE == 2)	 {
			r_printf("INFO: '%s' interpreted as %04d-%03d (Year-DOY)\n", s, year[0], day[0]);
			for (int i = 0; i < n; i++) {
				out[i] = YDOYtoF32time(year[i],  day[i]);
				//r_printf("%d %d  \n", year[i],  day[i]);
			}

		}
		else if (DONE == 3) {
			r_printf("INFO: '%s' interpreted as %04d-%02d (Year-Month)\n", s, year[0], month[0]);
			for (int i = 0; i < n; i++) {
				out[i] = year[i] + month[i]/12.0-1.0/24.0;
				//r_printf("%d %d  \n", year[i],  day[i]);
			}

		}
	}

	free(startidx);
	free(year);
	return out;
 
}
#include "abc_000_warning.h"