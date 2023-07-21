#pragma once
#include "abc_datatype.h"
#include "abc_common.h"
#include "bvs_header.h"

// bvs_basis_allocinit.c
void Alloc_BasisState_All(BVS_MODEL_PTR model, BVS_OPTIONS_PTR opt, MemPointers* MEM);
void Init_BasisState_All(BVS_MODEL_PTR model, BVS_OPTIONS_PTR opt, BVS_XINFO* xinfo);


// vcs_basis_gensegment.c
void GenTerms(BVS_TERMS_PTR terms, BVS_XINFO_PTR xinfo, F32PTR xcols, I32 K);
int GenOneTerm_hinge2d_GetGoodDataNpts(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol);
int GenOneTerm_hinge2dr_GetGoodDataNpts(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol);
int GenOneTerm_hinge1dr_GetGoodDataNpts(BVS_TERMS_PTR term, BVS_XINFO_PTR xinfo, F32PTR xcol);


//bvs_updategoodvec.c
void Update_BasisState_All(NewPropTerm * new, BVS_MODEL* model, BVS_XINFO_PTR xinfo);

//bvs_basis_proposeNew
void ProposeMove(BVS_MODEL_PTR model, NewPropTerm* NEW, NEWCOLINFOv2* newcol, PROP_DATA_PTR info);

//bvs_alloc_init
void Alloc_Xterms(F32PTR* Xmars, F32PTR* Xnewterm,BVS_OPTIONS_PTR opt, MemPointers* MEM);

//bvs-basis_io.c  
int ReadPreprocessInputData(BVS_OPTIONS_PTR opt, BVS_XINFO_PTR xinfo, BVS_YINFO_PTR yinfo, F32PTR MEMBUF);


//bvs_terms.c
void SaveOutput(BVS_XINFO* xinfo, BVS_TERMS_ENCODER* coder, BVS_OPTIONS_PTR opt);
void AllocInit_TermEncoder(BVS_TERMS_ENCODER* coder, MemPointers* MEM, BVS_OPTIONS_PTR opt); 
void InsertNewTermList(BVS_TERMS_PTR in, I32 K, BVS_TERMS_ENCODER* coder, F32PTR beta);
int ExtractNextTermFromList(BVS_TERMS_ENCODER* coder, BVS_TERMS_PTR out);
int  GetArg_Predict(VOIDPTR prhs[], int nrhs, BVS_XINFO_PTR xinfo, BVS_TERMS_ENCODER_PTR coder,BVS_PRED_XINFO* newxinfo, BVS_PRIOR_PTR* prior_ptr, BVS_PRED_OPT*newopt);

VOIDPTR Predict_Alloc_Output(BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, BVS_PRED_RESULT_PTR mat, BVS_PRIOR_PTR prior, int nModels, BVS_PRED_OPT* newopt);

//bvs_basis_genesegment_predict.c

void GenTermsPredict(BVS_TERMS_PTR terms, BVS_XINFO_PTR xinfo, BVS_PRED_XINFO_PTR newxinfo, F32PTR xcols, I32 K);


// predict_bvs.c
int bvs_predict_corev4(BVS_XINFO_PTR  xinfo, BVS_TERMS_ENCODER_PTR coder, BVS_PRED_XINFO* newxinfo, BVS_PRIOR_PTR* prior, BVS_PRED_OPT*);
