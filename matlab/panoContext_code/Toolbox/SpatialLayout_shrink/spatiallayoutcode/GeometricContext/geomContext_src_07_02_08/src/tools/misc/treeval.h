/*
 * MATLAB Compiler: 3.0
 * Date: Wed Aug 25 13:39:08 2004
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-W" "main" "-L" "C"
 * "-t" "-T" "link:lib" "-h" "test_track" "libmmfile.mlib" "-v" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __treeval_h
#define __treeval_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_treeval(void);
extern void TerminateModule_treeval(void);
extern _mexLocalFunctionTable _local_function_table_treeval;

extern mxArray * mlfNTreeval(int nargout,
                             mxArray * * nodes,
                             mxArray * * idname,
                             mxArray * Tree,
                             mxArray * X,
                             mxArray * subtrees);
extern mxArray * mlfTreeval(mxArray * * nodes,
                            mxArray * * idname,
                            mxArray * Tree,
                            mxArray * X,
                            mxArray * subtrees);
extern void mlfVTreeval(mxArray * Tree, mxArray * X, mxArray * subtrees);
extern void mlxTreeval(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
