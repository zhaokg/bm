#include "abc_000_warning.h"

#include "abc_001_config.h"


#include <stdio.h>

//#include "abc_common.h"
#include "abc_ide_util.h"
#include "abc_ts_func.h"
#include "abc_common.h"  //WriteF32ArrayToStrideMEM
#include "abc_vec.h"  //WriteF32ArrayToStrideMEM
#include "bvs_io.h"
 

void  BVS_WriteOutput(BVS_OPTIONS_PTR opt, BVS_RESULT_PTR result, BVS_XINFO *xinfo)
{

#define  _(x)          f32_to_strided_mem((F32PTR)result->x,  mat->x,  len,  stride, offset, datType)
#define  _2(x,y)        _(x), _(y)
#define  _3(x,y,z)      _2(x,y), _(z)
#define  _4(x,y,z,v)    _3(x,y,z), _(v)
#define  _5(x,y,z,v,v1) _4(x,y,z,v), _(v1)

	BVS_IO_PTR       io  =  &opt->io;
	BVS_RESULT_PTR   mat = io->out.result;

	DATA_TYPE datType = io->out.dtype;
 
	I64 len, offset, stride;

	len = 1, offset =0, stride = 1; //TODO: changed from 0 to 1
	f32_to_strided_mem(result->marg_lik, mat->marg_lik,	len, stride, offset, datType);
	f32_to_strided_mem(result->sig2,	 mat->sig2,		len, stride, offset, datType);
	f32_to_strided_mem(result->R2   ,    mat[0].R2, len, stride, offset, datType);
	f32_to_strided_mem(result->RMSE,     mat[0].RMSE, len, stride, offset, datType);
	f32_to_strided_mem(result->Kmean,    mat[0].Kmean, len, stride, offset, datType);

	len = xinfo->N;
	f32_to_strided_mem(result->Y, mat[0].Y, len, stride, offset, datType);
	f32_to_strided_mem(result->SD, mat[0].SD, len, stride, offset, datType);

	len = opt->prior.Kmax;
	f32_to_strided_mem(result->KProb, mat[0].KProb, len, stride, offset, datType);


	len = xinfo->Pall;
	f32_to_strided_mem(result->probAll, mat[0].probAll, len, stride, offset, datType);

	int type;

	type = LinearTerm;
	if (hasTerm(type)) {
		len = xinfo->xi_Pbasis[type] *1L ;
		f32_to_strided_mem(result->probLin, mat[0].probLin, len, stride, offset, datType);
		f32_to_strided_mem(result->Ylin, mat[0].Ylin, len, stride, offset, datType);
	}

	type = StepTerm;
	if ( hasTerm(type)  ) {
		len = xinfo->xi_Pbasis[type] * xinfo->N;
		f32_to_strided_mem(result->probStep, mat[0].probStep, len, stride, offset, datType);
		f32_to_strided_mem(result->Ystep,     mat[0].Ystep, len, stride, offset, datType);
	}
 
	type = StairTerm;
	if (hasTerm(type)) {
		len = xinfo->xi_Pbasis[type] * xinfo->N;
		f32_to_strided_mem(result->probStair, mat[0].probStair, len, stride, offset, datType);
		f32_to_strided_mem(result->Ystair,    mat[0].Ystair, len, stride, offset, datType);
	}

	type = PieceTerm;
	if (hasTerm(type)) {
		len = xinfo->xi_Pbasis[type] * xinfo->N;
		f32_to_strided_mem(result->probPiece, mat[0].probPiece, len, stride, offset, datType);
		f32_to_strided_mem(result->Ypiece,    mat[0].Ypiece, len, stride, offset, datType);
	}

	type = SeasonTerm;
	if (hasTerm(type)) {
		len = xinfo->xi_Pbasis[type] * xinfo->N;
		f32_to_strided_mem(result->probSeason, mat[0].probSeason, len, stride, offset, datType);
		f32_to_strided_mem(result->Yseason, mat[0].Yseason, len, stride, offset, datType);
	}

	type = HingeTerm;
	if (hasTerm(type)) {
		len = xinfo->xi_Pbasis[type] * xinfo->N;
		f32_to_strided_mem(result->probHinge, mat[0].probHinge, len, stride, offset, datType);
		f32_to_strided_mem(result->Yhinge, mat[0].Yhinge, len, stride, offset, datType);
	}

  
	type = HingePairTerm;
	if (hasTerm(type)) {
		len = xinfo->xi_Pbasis[type] * xinfo->N;
		f32_to_strided_mem(result->probHingePair, mat[0].probHingePair, len, stride, offset, datType);
		f32_to_strided_mem(result->YhingePair, mat[0].YhingePair, len, stride, offset, datType);
	}


	type = ChangeTerm;
	if (hasTerm(type)) {
		len = xinfo->xi_Pbasis[type] * xinfo->N;
		f32_to_strided_mem(result->probChange, mat[0].probChange, len, stride, offset, datType);
		f32_to_strided_mem(result->Ychange, mat[0].Ychange, len, stride, offset, datType);
	}

 
}
 
 
#include "abc_000_warning.h"