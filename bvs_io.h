#pragma once
#include "abc_datatype.h"
#include "abc_common.h"
#include "bvs_header.h"

/****************beastv2_in_args******************/
extern int  BVS_GetArgs(VOID_PTR prhs[], int nrhs, A(OPTIONS_PTR) opt);
extern void  BVS_DeallocatePriorPts(BVS_PRIOR_PTR prior);

 
/****************beastv2_out_allocmem******************/
extern void* Alloc_OutPut(BVS_OPTIONS_PTR  opt);
extern void  Init_Result(BVS_RESULT_PTR  result, BVS_OPTIONS_PTR  opt, const F32 nan);
extern void  Alloc_Result(BVS_RESULT_PTR  result, BVS_OPTIONS_PTR  opt, MemPointers* _restrict MEM);

/****************beastv2_out_print******************/
 

/****************beastv2_out_write******************/
void  BVS_WriteOutput(BVS_OPTIONS_PTR opt, BVS_RESULT_PTR result, BVS_XINFO* xinfo);

