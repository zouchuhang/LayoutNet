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
#include "treevalc_mex_interface.h"

void InitializeModule_treevalc_mex_interface(void) {
}

void TerminateModule_treevalc_mex_interface(void) {
}

static mxArray * Mtreevalc(int nargout_, mxArray * varargin);

_mexLocalFunctionTable _local_function_table_treevalc
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfNTreevalc" contains the nargout interface for the
 * "treevalc" M-function from file
 * "/afs/cs.cmu.edu/user/dhoiem/solar_demo_matlab/treevalc.mexglx" (lines 0-0).
 * This interface is only produced if the M-function uses the special variable
 * "nargout". The nargout interface allows the number of requested outputs to
 * be specified via the nargout argument, as opposed to the normal interface
 * which dynamically calculates the number of outputs based on the number of
 * non-NULL inputs it receives. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
mxArray * mlfNTreevalc(int nargout, mlfVarargoutList * varargout, ...) {
    mxArray * varargin = NULL;
    mlfVarargin(&varargin, varargout, 0);
    mlfEnterNewContext(0, -1, varargin);
    nargout += mclNargout(varargout);
    *mlfGetVarargoutCellPtr(varargout) = Mtreevalc(nargout, varargin);
    mlfRestorePreviousContext(0, 0);
    mxDestroyArray(varargin);
    return mlfAssignOutputs(varargout);
}

/*
 * The function "mlfTreevalc" contains the normal interface for the "treevalc"
 * M-function from file
 * "/afs/cs.cmu.edu/user/dhoiem/solar_demo_matlab/treevalc.mexglx" (lines 0-0).
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
mxArray * mlfTreevalc(mlfVarargoutList * varargout, ...) {
    mxArray * varargin = NULL;
    int nargout = 0;
    mlfVarargin(&varargin, varargout, 0);
    mlfEnterNewContext(0, -1, varargin);
    nargout += mclNargout(varargout);
    *mlfGetVarargoutCellPtr(varargout) = Mtreevalc(nargout, varargin);
    mlfRestorePreviousContext(0, 0);
    mxDestroyArray(varargin);
    return mlfAssignOutputs(varargout);
}

/*
 * The function "mlfVTreevalc" contains the void interface for the "treevalc"
 * M-function from file
 * "/afs/cs.cmu.edu/user/dhoiem/solar_demo_matlab/treevalc.mexglx" (lines 0-0).
 * The void interface is only produced if the M-function uses the special
 * variable "nargout", and has at least one output. The void interface function
 * specifies zero output arguments to the implementation version of the
 * function, and in the event that the implementation version still returns an
 * output (which, in MATLAB, would be assigned to the "ans" variable), it
 * deallocates the output. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlfVTreevalc(mxArray * synthetic_varargin_argument, ...) {
    mxArray * varargin = NULL;
    mxArray * varargout = NULL;
    mlfVarargin(&varargin, synthetic_varargin_argument, 1);
    mlfEnterNewContext(0, -1, varargin);
    varargout = Mtreevalc(0, synthetic_varargin_argument);
    mlfRestorePreviousContext(0, 0);
    mxDestroyArray(varargin);
}

/*
 * The function "mlxTreevalc" contains the feval interface for the "treevalc"
 * M-function from file
 * "/afs/cs.cmu.edu/user/dhoiem/solar_demo_matlab/treevalc.mexglx" (lines 0-0).
 * The feval function calls the implementation version of treevalc through this
 * function. This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxTreevalc(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[1];
    int i;
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    mlfEnterNewContext(0, 0);
    mprhs[0] = NULL;
    mlfAssign(&mprhs[0], mclCreateVararginCell(nrhs, prhs));
    mplhs[0] = Mtreevalc(nlhs, mprhs[0]);
    mclAssignVarargoutCell(0, nlhs, plhs, mplhs[0]);
    mlfRestorePreviousContext(0, 0);
    mxDestroyArray(mprhs[0]);
}

/*
 * The function "Mtreevalc" is the implementation version of the "treevalc"
 * M-function from file
 * "/afs/cs.cmu.edu/user/dhoiem/solar_demo_matlab/treevalc.mexglx" (lines 1-1).
 * It contains the actual compiled code for that M-function. It is a static
 * function and must only be called from one of the interface functions,
 * appearing below.
 */
static mxArray * Mtreevalc(int nargout_, mxArray * varargin) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_treevalc);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return mclCExecMexFunction("treevalc", nargout_, varargin);
}
