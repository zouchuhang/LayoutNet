/*
 * MATLAB Compiler: 3.0
 * Date: Wed Nov 10 17:22:45 2004
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-t" "-W" "lib:libsd"
 * "-T" "link:lib" "block2fft" "compute_all_features_track"
 * "compute_data_representation_fft" "compute_integral_array"
 * "compute_windowed_fft" "makeblocks" "remove_close_matches" "test_track"
 * "treevalc.mexglx" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __treevalc_mex_interface_h
#define __treevalc_mex_interface_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_treevalc_mex_interface(void);
extern void TerminateModule_treevalc_mex_interface(void);
extern _mexLocalFunctionTable _local_function_table_treevalc;

extern mxArray * mlfNTreevalc(int nargout, mlfVarargoutList * varargout, ...);
extern mxArray * mlfTreevalc(mlfVarargoutList * varargout, ...);
extern void mlfVTreevalc(mxArray * synthetic_varargin_argument, ...);
extern void mlxTreevalc(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
