/*
 * MATLAB Compiler: 3.0
 * Date: Wed Aug 25 13:39:08 2004
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-W" "main" "-L" "C"
 * "-t" "-T" "link:lib" "-h" "test_track" "libmmfile.mlib" "-v" 
 */
#include "treeval.h"
#include "libmatlbm.h"
#include "libmmfile.h"

static mxChar _array1_[6] = { 'm', 'e', 't', 'h', 'o', 'd' };
static mxArray * _mxarray0_;

static mxChar _array3_[43] = { 'T', 'h', 'e', ' ', 'f', 'i', 'r', 's', 't',
                               ' ', 'a', 'r', 'g', 'u', 'm', 'e', 'n', 't',
                               ' ', 'm', 'u', 's', 't', ' ', 'b', 'e', ' ',
                               'a', ' ', 'd', 'e', 'c', 'i', 's', 'i', 'o',
                               'n', ' ', 't', 'r', 'e', 'e', '.' };
static mxArray * _mxarray2_;

static mxChar _array5_[14] = { 'c', 'l', 'a', 's', 's', 'i', 'f',
                               'i', 'c', 'a', 't', 'i', 'o', 'n' };
static mxArray * _mxarray4_;

static mxChar _array7_[55] = { 'O', 'n', 'l', 'y', ' ', '2', ' ', 'o', 'u', 't',
                               'p', 'u', 't', ' ', 'a', 'r', 'g', 'u', 'm', 'e',
                               'n', 't', 's', ' ', 'a', 'v', 'a', 'i', 'l', 'a',
                               'b', 'l', 'e', ' ', 'f', 'o', 'r', ' ', 'r', 'e',
                               'g', 'r', 'e', 's', 's', 'i', 'o', 'n', ' ', 't',
                               'r', 'e', 'e', 's', '.' };
static mxArray * _mxarray6_;

static mxChar _array9_[34] = { 'T', 'h', 'e', ' ', 'X', ' ', 'm', 'a', 't',
                               'r', 'i', 'x', ' ', 'm', 'u', 's', 't', ' ',
                               'h', 'a', 'v', 'e', ' ', '%', 'd', ' ', 'c',
                               'o', 'l', 'u', 'm', 'n', 's', '.' };
static mxArray * _mxarray8_;
static mxArray * _mxarray10_;

static mxChar _array12_[26] = { 'S', 'U', 'B', 'T', 'R', 'E', 'E', 'S', ' ',
                                'm', 'u', 's', 't', ' ', 'b', 'e', ' ', 'a',
                                ' ', 'v', 'e', 'c', 't', 'o', 'r', '.' };
static mxArray * _mxarray11_;

static mxChar _array14_[24] = { 'S', 'U', 'B', 'T', 'R', 'E', 'E', 'S',
                                ' ', 'm', 'u', 's', 't', ' ', 'b', 'e',
                                ' ', 's', 'o', 'r', 't', 'e', 'd', '.' };
static mxArray * _mxarray13_;

static mxChar _array16_[9] = { 'p', 'r', 'u', 'n', 'e', 'l', 'i', 's', 't' };
static mxArray * _mxarray15_;

static mxChar _array18_[59] = { 'T', 'h', 'e', ' ', 'd', 'e', 'c', 'i', 's',
                                'i', 'o', 'n', ' ', 't', 'r', 'e', 'e', ' ',
                                'T', 'R', 'E', 'E', ' ', 'd', 'o', 'e', 's',
                                ' ', 'n', 'o', 't', ' ', 'i', 'n', 'c', 'l',
                                'u', 'd', 'e', ' ', 'a', ' ', 'p', 'r', 'u',
                                'n', 'i', 'n', 'g', ' ', 's', 'e', 'q', 'u',
                                'e', 'n', 'c', 'e', '.' };
static mxArray * _mxarray17_;
static double _ieee_plusinf_;
static mxArray * _mxarray19_;
static mxArray * _mxarray20_;
static mxArray * _mxarray21_;

void InitializeModule_treeval(void) {
    _mxarray0_ = mclInitializeString(6, _array1_);
    _mxarray2_ = mclInitializeString(43, _array3_);
    _mxarray4_ = mclInitializeString(14, _array5_);
    _mxarray6_ = mclInitializeString(55, _array7_);
    _mxarray8_ = mclInitializeString(34, _array9_);
    _mxarray10_ = mclInitializeDouble(0.0);
    _mxarray11_ = mclInitializeString(26, _array12_);
    _mxarray13_ = mclInitializeString(24, _array14_);
    _mxarray15_ = mclInitializeString(9, _array16_);
    _mxarray17_ = mclInitializeString(59, _array18_);
    _ieee_plusinf_ = mclGetInf();
    _mxarray19_ = mclInitializeDouble(_ieee_plusinf_);
    _mxarray20_ = mclInitializeDouble(1.0);
    _mxarray21_ = mclInitializeDouble(2.0);
}

void TerminateModule_treeval(void) {
    mxDestroyArray(_mxarray21_);
    mxDestroyArray(_mxarray20_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray17_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray10_);
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * mlfTreeval_doapply(mxArray * Tree,
                                    mxArray * X,
                                    mxArray * rows,
                                    mxArray * thisnode,
                                    mxArray * nodes_in,
                                    mxArray * subtrees,
                                    mxArray * prunelist,
                                    mxArray * endcol);
static void mlxTreeval_doapply(int nlhs,
                               mxArray * plhs[],
                               int nrhs,
                               mxArray * prhs[]);
static mxArray * Mtreeval(mxArray * * nodes,
                          mxArray * * idname,
                          int nargout_,
                          mxArray * Tree,
                          mxArray * X,
                          mxArray * subtrees);
static mxArray * Mtreeval_doapply(int nargout_,
                                  mxArray * Tree,
                                  mxArray * X,
                                  mxArray * rows,
                                  mxArray * thisnode,
                                  mxArray * nodes_in,
                                  mxArray * subtrees,
                                  mxArray * prunelist,
                                  mxArray * endcol);

static mexFunctionTableEntry local_function_table_[1]
  = { { "doapply", mlxTreeval_doapply, 8, 1, NULL } };

_mexLocalFunctionTable _local_function_table_treeval
  = { 1, local_function_table_ };

/*
 * The function "mlfNTreeval" contains the nargout interface for the "treeval"
 * M-function from file "/usr/local/lib/matlab6/toolbox/stats/treeval.m" (lines
 * 1-83). This interface is only produced if the M-function uses the special
 * variable "nargout". The nargout interface allows the number of requested
 * outputs to be specified via the nargout argument, as opposed to the normal
 * interface which dynamically calculates the number of outputs based on the
 * number of non-NULL inputs it receives. This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
mxArray * mlfNTreeval(int nargout,
                      mxArray * * nodes,
                      mxArray * * idname,
                      mxArray * Tree,
                      mxArray * X,
                      mxArray * subtrees) {
    mxArray * id = NULL;
    mxArray * nodes__ = NULL;
    mxArray * idname__ = NULL;
    mlfEnterNewContext(2, 3, nodes, idname, Tree, X, subtrees);
    id = Mtreeval(&nodes__, &idname__, nargout, Tree, X, subtrees);
    mlfRestorePreviousContext(2, 3, nodes, idname, Tree, X, subtrees);
    if (nodes != NULL) {
        mclCopyOutputArg(nodes, nodes__);
    } else {
        mxDestroyArray(nodes__);
    }
    if (idname != NULL) {
        mclCopyOutputArg(idname, idname__);
    } else {
        mxDestroyArray(idname__);
    }
    return mlfReturnValue(id);
}

/*
 * The function "mlfTreeval" contains the normal interface for the "treeval"
 * M-function from file "/usr/local/lib/matlab6/toolbox/stats/treeval.m" (lines
 * 1-83). This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
mxArray * mlfTreeval(mxArray * * nodes,
                     mxArray * * idname,
                     mxArray * Tree,
                     mxArray * X,
                     mxArray * subtrees) {
    int nargout = 1;
    mxArray * id = NULL;
    mxArray * nodes__ = NULL;
    mxArray * idname__ = NULL;
    mlfEnterNewContext(2, 3, nodes, idname, Tree, X, subtrees);
    if (nodes != NULL) {
        ++nargout;
    }
    if (idname != NULL) {
        ++nargout;
    }
    id = Mtreeval(&nodes__, &idname__, nargout, Tree, X, subtrees);
    mlfRestorePreviousContext(2, 3, nodes, idname, Tree, X, subtrees);
    if (nodes != NULL) {
        mclCopyOutputArg(nodes, nodes__);
    } else {
        mxDestroyArray(nodes__);
    }
    if (idname != NULL) {
        mclCopyOutputArg(idname, idname__);
    } else {
        mxDestroyArray(idname__);
    }
    return mlfReturnValue(id);
}

/*
 * The function "mlfVTreeval" contains the void interface for the "treeval"
 * M-function from file "/usr/local/lib/matlab6/toolbox/stats/treeval.m" (lines
 * 1-83). The void interface is only produced if the M-function uses the
 * special variable "nargout", and has at least one output. The void interface
 * function specifies zero output arguments to the implementation version of
 * the function, and in the event that the implementation version still returns
 * an output (which, in MATLAB, would be assigned to the "ans" variable), it
 * deallocates the output. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlfVTreeval(mxArray * Tree, mxArray * X, mxArray * subtrees) {
    mxArray * id = NULL;
    mxArray * nodes = NULL;
    mxArray * idname = NULL;
    mlfEnterNewContext(0, 3, Tree, X, subtrees);
    id = Mtreeval(&nodes, &idname, 0, Tree, X, subtrees);
    mlfRestorePreviousContext(0, 3, Tree, X, subtrees);
    mxDestroyArray(id);
    mxDestroyArray(nodes);
    mxDestroyArray(idname);
}

/*
 * The function "mlxTreeval" contains the feval interface for the "treeval"
 * M-function from file "/usr/local/lib/matlab6/toolbox/stats/treeval.m" (lines
 * 1-83). The feval function calls the implementation version of treeval
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxTreeval(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[3];
    int i;
    if (nlhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: treeval Line: 1 Column: "
            "1 The function \"treeval\" was called with mor"
            "e than the declared number of outputs (3)."),
          NULL);
    }
    if (nrhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: treeval Line: 1 Column:"
            " 1 The function \"treeval\" was called with m"
            "ore than the declared number of inputs (3)."),
          NULL);
    }
    for (i = 0; i < 3; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 3 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 3; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    mplhs[0]
      = Mtreeval(&mplhs[1], &mplhs[2], nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 3 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 3; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "mlfTreeval_doapply" contains the normal interface for the
 * "treeval/doapply" M-function from file
 * "/usr/local/lib/matlab6/toolbox/stats/treeval.m" (lines 83-145). This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
static mxArray * mlfTreeval_doapply(mxArray * Tree,
                                    mxArray * X,
                                    mxArray * rows,
                                    mxArray * thisnode,
                                    mxArray * nodes_in,
                                    mxArray * subtrees,
                                    mxArray * prunelist,
                                    mxArray * endcol) {
    int nargout = 1;
    mxArray * nodes = NULL;
    mlfEnterNewContext(
      0, 8, Tree, X, rows, thisnode, nodes_in, subtrees, prunelist, endcol);
    nodes
      = Mtreeval_doapply(
          nargout,
          Tree,
          X,
          rows,
          thisnode,
          nodes_in,
          subtrees,
          prunelist,
          endcol);
    mlfRestorePreviousContext(
      0, 8, Tree, X, rows, thisnode, nodes_in, subtrees, prunelist, endcol);
    return mlfReturnValue(nodes);
}

/*
 * The function "mlxTreeval_doapply" contains the feval interface for the
 * "treeval/doapply" M-function from file
 * "/usr/local/lib/matlab6/toolbox/stats/treeval.m" (lines 83-145). The feval
 * function calls the implementation version of treeval/doapply through this
 * function. This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
static void mlxTreeval_doapply(int nlhs,
                               mxArray * plhs[],
                               int nrhs,
                               mxArray * prhs[]) {
    mxArray * mprhs[8];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: treeval/doapply Line: 83 Colu"
            "mn: 1 The function \"treeval/doapply\" was called w"
            "ith more than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 8) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: treeval/doapply Line: 83 Colu"
            "mn: 1 The function \"treeval/doapply\" was called w"
            "ith more than the declared number of inputs (8)."),
          NULL);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 8 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 8; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(
      0,
      8,
      mprhs[0],
      mprhs[1],
      mprhs[2],
      mprhs[3],
      mprhs[4],
      mprhs[5],
      mprhs[6],
      mprhs[7]);
    mplhs[0]
      = Mtreeval_doapply(
          nlhs,
          mprhs[0],
          mprhs[1],
          mprhs[2],
          mprhs[3],
          mprhs[4],
          mprhs[5],
          mprhs[6],
          mprhs[7]);
    mlfRestorePreviousContext(
      0,
      8,
      mprhs[0],
      mprhs[1],
      mprhs[2],
      mprhs[3],
      mprhs[4],
      mprhs[5],
      mprhs[6],
      mprhs[7]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mtreeval" is the implementation version of the "treeval"
 * M-function from file "/usr/local/lib/matlab6/toolbox/stats/treeval.m" (lines
 * 1-83). It contains the actual compiled code for that M-function. It is a
 * static function and must only be called from one of the interface functions,
 * appearing below.
 */
/*
 * function [id,nodes,idname]=treeval(Tree,X,subtrees)
 */
static mxArray * Mtreeval(mxArray * * nodes,
                          mxArray * * idname,
                          int nargout_,
                          mxArray * Tree,
                          mxArray * X,
                          mxArray * subtrees) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_treeval);
    int nargin_ = mclNargin(3, Tree, X, subtrees, NULL);
    mxArray * id = NULL;
    mxArray * ntrees = NULL;
    mxArray * prunelist = NULL;
    mxArray * nc = NULL;
    mxArray * nr = NULL;
    mxArray * ans = NULL;
    mclCopyArray(&Tree);
    mclCopyArray(&X);
    mclCopyArray(&subtrees);
    /*
     * %TREEVAL Compute fitted value for decision tree applied to data.
     * %   YFIT = TREEVAL(TREE,X) takes a classification or regression tree TREE
     * %   as produced by the TREEFIT function, and a matrix X of predictor
     * %   values, and produces a vector YFIT of predicted response values.
     * %   For a regression tree, YFIT(j) is the fitted response value for a
     * %   point having the predictor values X(j,:).  For a classification tree,
     * %   YFIT(j) is the class number into which the tree would assign the point
     * %   with data X(j,:).  To convert the number into a class name, use the
     * %   third output argument (see below).
     * %
     * %   YFIT = TREEVAL(TREE,X,SUBTREES) takes an additional vector SUBTREES of
     * %   pruning levels, with 0 representing the full, unpruned tree.  TREE must
     * %   include a pruning sequence as created by the TREEFIT or TREEPRUNE function.
     * %   If SUBTREES has K elements and X has N rows, then the output YFIT is an
     * %   N-by-K matrix, with the Jth column containing the fitted values produced by
     * %   the SUBTREES(J) subtree.  SUBTREES must be sorted in ascending order.
     * %   (To compute fitted values for a tree that is not part of the optimal
     * %   pruning sequence, first use TREEPRUNE to prune the tree.)
     * %
     * %   [YFIT,NODE] = TREEVAL(...) also returns an array NODE of the same size
     * %   as YFIT containing the node number assigned to each row of X.  The
     * %   TREEDISP function can display the node numbers for any node you select.
     * %
     * %   [YFIT,NODE,CNAME] = TREEVAL(...) is valid only for classification trees.
     * %   It retuns a cell array CNAME containing the predicted class names.
     * %
     * %   NaN values in the X matrix are treated as missing.  If the TREEVAL
     * %   function encounters a missing value when it attempts to evaluate the
     * %   split rule at a branch node, it cannot determine whether to proceed to
     * %   the left or right child node.  Instead, it sets the corresponding fitted
     * %   value equal to the fitted value assigned to the branch node.
     * %
     * %   Example: Find predicted classifications for Fisher's iris data.
     * %      load fisheriris;
     * %      t = treefit(meas, species);  % create decision tree
     * %      sfit = treeval(t,meas);      % find assigned class numbers
     * %      sfit = t.classname(sfit);    % get class names
     * %      mean(strcmp(sfit,species))   % compute proportion correctly classified
     * %
     * %   See also TREEFIT, TREEPRUNE, TREEDISP, TREETEST.
     * 
     * %   Copyright 1993-2002 The MathWorks, Inc. 
     * %   $Revision: 1.2 $  $Date: 2002/03/21 20:36:31 $
     * 
     * if ~isstruct(Tree) | ~isfield(Tree,'method')
     */
    {
        mxArray * a_ = mclInitialize(mclNot(mlfIsstruct(mclVa(Tree, "Tree"))));
        if (mlfTobool(a_)
            || mlfTobool(
                 mclOr(
                   a_, mclNot(mlfIsfield(mclVa(Tree, "Tree"), _mxarray0_))))) {
            mxDestroyArray(a_);
            /*
             * error('The first argument must be a decision tree.');
             */
            mlfError(_mxarray2_, NULL);
        } else {
            mxDestroyArray(a_);
        }
    /*
     * end
     */
    }
    /*
     * if nargout>=3 & ~isequal(Tree.method,'classification')
     */
    {
        mxArray * a_ = mclInitialize(mclBoolToArray(nargout_ >= 3));
        if (mlfTobool(a_)
            && mlfTobool(
                 mclAnd(
                   a_,
                   mclNot(
                     mclFeval(
                       mclValueVarargout(),
                       mlxIsequal,
                       mlfIndexRef(mclVa(Tree, "Tree"), ".method"),
                       _mxarray4_,
                       NULL))))) {
            mxDestroyArray(a_);
            /*
             * error('Only 2 output arguments available for regression trees.');
             */
            mlfError(_mxarray6_, NULL);
        } else {
            mxDestroyArray(a_);
        }
    /*
     * end
     */
    }
    /*
     * [nr,nc] = size(X);
     */
    mlfSize(mlfVarargout(&nr, &nc, NULL), mclVa(X, "X"), NULL);
    /*
     * if nc~=Tree.npred
     */
    if (mlfTobool(
          mclFeval(
            mclValueVarargout(),
            mlxNe,
            mclVv(nc, "nc"),
            mlfIndexRef(mclVa(Tree, "Tree"), ".npred"),
            NULL))) {
        /*
         * error(sprintf('The X matrix must have %d columns.',Tree.npred));
         */
        mlfError(
          mlfSprintf(
            NULL, _mxarray8_, mlfIndexRef(mclVa(Tree, "Tree"), ".npred"), NULL),
          NULL);
    /*
     * end
     */
    }
    /*
     * 
     * if nargin<3
     */
    if (nargin_ < 3) {
        /*
         * subtrees = 0;
         */
        mlfAssign(&subtrees, _mxarray10_);
    /*
     * elseif prod(size(subtrees))>length(subtrees)
     */
    } else if (mclGtBool(
                 mlfProd(
                   mlfSize(
                     mclValueVarargout(), mclVa(subtrees, "subtrees"), NULL),
                   NULL),
                 mlfScalar(mclLengthInt(mclVa(subtrees, "subtrees"))))) {
        /*
         * error('SUBTREES must be a vector.');
         */
        mlfError(_mxarray11_, NULL);
    /*
     * elseif any(diff(subtrees)<0)
     */
    } else if (mlfTobool(
                 mlfAny(
                   mclLt(
                     mlfDiff(mclVa(subtrees, "subtrees"), NULL, NULL),
                     _mxarray10_),
                   NULL))) {
        /*
         * error('SUBTREES must be sorted.');
         */
        mlfError(_mxarray13_, NULL);
    /*
     * end
     */
    }
    /*
     * 
     * if isfield(Tree,'prunelist')
     */
    if (mlfTobool(mlfIsfield(mclVa(Tree, "Tree"), _mxarray15_))) {
        /*
         * prunelist = Tree.prunelist;
         */
        mlfAssign(&prunelist, mlfIndexRef(mclVa(Tree, "Tree"), ".prunelist"));
    /*
     * elseif ~isequal(subtrees,0)
     */
    } else if (mclNotBool(
                 mlfIsequal(mclVa(subtrees, "subtrees"), _mxarray10_, NULL))) {
        /*
         * error('The decision tree TREE does not include a pruning sequence.')
         */
        mlfError(_mxarray17_, NULL);
    /*
     * else
     */
    } else {
        /*
         * prunelist = repmat(Inf,size(Tree.node));
         */
        mlfAssign(
          &prunelist,
          mlfRepmat(
            _mxarray19_,
            mclFeval(
              mclValueVarargout(),
              mlxSize,
              mlfIndexRef(mclVa(Tree, "Tree"), ".node"),
              NULL),
            NULL));
    /*
     * end
     */
    }
    /*
     * 
     * ntrees = length(subtrees);
     */
    mlfAssign(&ntrees, mlfScalar(mclLengthInt(mclVa(subtrees, "subtrees"))));
    /*
     * nodes = doapply(Tree,X,1:nr,1,zeros(nr,ntrees),subtrees,prunelist,ntrees);
     */
    mlfAssign(
      nodes,
      mlfTreeval_doapply(
        mclVa(Tree, "Tree"),
        mclVa(X, "X"),
        mlfColon(_mxarray20_, mclVv(nr, "nr"), NULL),
        _mxarray20_,
        mlfZeros(mclVv(nr, "nr"), mclVv(ntrees, "ntrees"), NULL),
        mclVa(subtrees, "subtrees"),
        mclVv(prunelist, "prunelist"),
        mclVv(ntrees, "ntrees")));
    /*
     * id = Tree.class(nodes);
     */
    mlfAssign(
      &id,
      mlfIndexRef(mclVa(Tree, "Tree"), ".class(?)", mclVv(*nodes, "nodes")));
    /*
     * 
     * if nargout>=3
     */
    if (nargout_ >= 3) {
        /*
         * idname = Tree.classname(id);
         */
        mlfAssign(
          idname,
          mlfIndexRef(mclVa(Tree, "Tree"), ".classname(?)", mclVv(id, "id")));
    /*
     * end
     */
    }
    mclValidateOutput(id, 1, nargout_, "id", "treeval");
    mclValidateOutput(*nodes, 2, nargout_, "nodes", "treeval");
    mclValidateOutput(*idname, 3, nargout_, "idname", "treeval");
    mxDestroyArray(ans);
    mxDestroyArray(nr);
    mxDestroyArray(nc);
    mxDestroyArray(prunelist);
    mxDestroyArray(ntrees);
    mxDestroyArray(subtrees);
    mxDestroyArray(X);
    mxDestroyArray(Tree);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return id;
    /*
     * 
     * 
     * %------------------------------------------------
     */
}

/*
 * The function "Mtreeval_doapply" is the implementation version of the
 * "treeval/doapply" M-function from file
 * "/usr/local/lib/matlab6/toolbox/stats/treeval.m" (lines 83-145). It contains
 * the actual compiled code for that M-function. It is a static function and
 * must only be called from one of the interface functions, appearing below.
 */
/*
 * function nodes = doapply(Tree,X,rows,thisnode,nodes,subtrees,prunelist,endcol)
 */
static mxArray * Mtreeval_doapply(int nargout_,
                                  mxArray * Tree,
                                  mxArray * X,
                                  mxArray * rows,
                                  mxArray * thisnode,
                                  mxArray * nodes_in,
                                  mxArray * subtrees,
                                  mxArray * prunelist,
                                  mxArray * endcol) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_treeval);
    mxArray * nodes = NULL;
    mxArray * subrows = NULL;
    mxArray * ismissing = NULL;
    mxArray * isright = NULL;
    mxArray * isleft = NULL;
    mxArray * x = NULL;
    mxArray * ncols = NULL;
    mxArray * ntrees = NULL;
    mxArray * prunelevel = NULL;
    mxArray * catsplit = NULL;
    mxArray * kids = NULL;
    mxArray * assignedclass = NULL;
    mxArray * cutoff = NULL;
    mxArray * splitvar = NULL;
    mclCopyArray(&Tree);
    mclCopyArray(&X);
    mclCopyArray(&rows);
    mclCopyArray(&thisnode);
    mclCopyInputArg(&nodes, nodes_in);
    mclCopyArray(&subtrees);
    mclCopyArray(&prunelist);
    mclCopyArray(&endcol);
    /*
     * %DOAPPLY Apply classification rule to specified rows starting at a node.
     * %   This is a recursive function.  Starts at top node, then recurses over
     * %   child nodes.  THISNODE is the current node at each step.
     * %
     * %   NODES has one row per observation and one column per subtree.
     * %
     * %   X, NODES, PRUNELIST, and SUBTREES are the same in each recursive call
     * %   as they were in the top-level call.  ROWS describes the subset of X and
     * %   NODES to consider.  1:ENDCOL are colums of NODES and the elements of
     * %   SUBTREES to consider.
     * 
     * splitvar      = Tree.var(thisnode);
     */
    mlfAssign(
      &splitvar,
      mlfIndexRef(mclVa(Tree, "Tree"), ".var(?)", mclVa(thisnode, "thisnode")));
    /*
     * cutoff        = Tree.cut(thisnode);
     */
    mlfAssign(
      &cutoff,
      mlfIndexRef(mclVa(Tree, "Tree"), ".cut(?)", mclVa(thisnode, "thisnode")));
    /*
     * assignedclass = Tree.class(thisnode);
     */
    mlfAssign(
      &assignedclass,
      mlfIndexRef(
        mclVa(Tree, "Tree"), ".class(?)", mclVa(thisnode, "thisnode")));
    /*
     * kids          = Tree.children(thisnode,:);
     */
    mlfAssign(
      &kids,
      mlfIndexRef(
        mclVa(Tree, "Tree"),
        ".children(?,?)",
        mclVa(thisnode, "thisnode"),
        mlfCreateColonIndex()));
    /*
     * catsplit      = Tree.catsplit;
     */
    mlfAssign(&catsplit, mlfIndexRef(mclVa(Tree, "Tree"), ".catsplit"));
    /*
     * prunelevel    = prunelist(thisnode);
     */
    mlfAssign(
      &prunelevel,
      mclArrayRef1(mclVa(prunelist, "prunelist"), mclVa(thisnode, "thisnode")));
    /*
     * ntrees        = size(nodes,2);   % number of trees in sequence
     */
    mlfAssign(
      &ntrees,
      mlfSize(mclValueVarargout(), mclVa(nodes, "nodes"), _mxarray21_));
    /*
     * 
     * % For how many of the remaining trees is this a terminal node?
     * if splitvar==0      % all, if it's terminal on the unpruned tree
     */
    if (mclEqBool(mclVv(splitvar, "splitvar"), _mxarray10_)) {
        /*
         * ncols = endcol;
         */
        mlfAssign(&ncols, mclVa(endcol, "endcol"));
    /*
     * else                % some, if it's terminal only after pruning
     */
    } else {
        /*
         * ncols = sum(subtrees(1:endcol) >= prunelevel);
         */
        mlfAssign(
          &ncols,
          mlfSum(
            mclGe(
              mclArrayRef1(
                mclVa(subtrees, "subtrees"),
                mlfColon(_mxarray20_, mclVa(endcol, "endcol"), NULL)),
              mclVv(prunelevel, "prunelevel")),
            NULL));
    /*
     * end
     */
    }
    /*
     * if ncols>0          % for those trees, assign the node level now
     */
    if (mclGtBool(mclVv(ncols, "ncols"), _mxarray10_)) {
        /*
         * nodes(rows,(endcol-ncols+1:endcol)) = thisnode;
         */
        mclArrayAssign2(
          &nodes,
          mclVa(thisnode, "thisnode"),
          mclVa(rows, "rows"),
          mlfColon(
            mclPlus(
              mclMinus(mclVa(endcol, "endcol"), mclVv(ncols, "ncols")),
              _mxarray20_),
            mclVa(endcol, "endcol"),
            NULL));
        /*
         * endcol = endcol - ncols;
         */
        mlfAssign(
          &endcol, mclMinus(mclVa(endcol, "endcol"), mclVv(ncols, "ncols")));
    /*
     * end
     */
    }
    /*
     * 
     * % Now deal with non-terminal nodes
     * if endcol > 0
     */
    if (mclGtBool(mclVa(endcol, "endcol"), _mxarray10_)) {
        /*
         * % Determine if this point goes left, goes right, or stays here
         * x = X(rows,abs(splitvar));   
         */
        mlfAssign(
          &x,
          mclArrayRef2(
            mclVa(X, "X"),
            mclVa(rows, "rows"),
            mlfAbs(mclVv(splitvar, "splitvar"))));
        /*
         * if splitvar>0                % continuous variable
         */
        if (mclGtBool(mclVv(splitvar, "splitvar"), _mxarray10_)) {
            /*
             * isleft = (x < cutoff);
             */
            mlfAssign(&isleft, mclLt(mclVv(x, "x"), mclVv(cutoff, "cutoff")));
            /*
             * isright = ~isleft;
             */
            mlfAssign(&isright, mclNot(mclVv(isleft, "isleft")));
            /*
             * ismissing = isnan(x);
             */
            mlfAssign(&ismissing, mlfIsnan(mclVv(x, "x")));
        /*
         * else                         % categorical variable
         */
        } else {
            /*
             * isleft = ismember(x,catsplit{cutoff,1});
             */
            mlfAssign(
              &isleft,
              mclFeval(
                mclValueVarargout(),
                mlxIsmember,
                mclVv(x, "x"),
                mlfIndexRef(
                  mclVv(catsplit, "catsplit"),
                  "{?,?}",
                  mclVv(cutoff, "cutoff"),
                  _mxarray20_),
                NULL));
            /*
             * isright = ismember(x,catsplit{cutoff,2});
             */
            mlfAssign(
              &isright,
              mclFeval(
                mclValueVarargout(),
                mlxIsmember,
                mclVv(x, "x"),
                mlfIndexRef(
                  mclVv(catsplit, "catsplit"),
                  "{?,?}",
                  mclVv(cutoff, "cutoff"),
                  _mxarray21_),
                NULL));
            /*
             * ismissing = ~(isleft | isright);
             */
            mlfAssign(
              &ismissing,
              mclNot(
                mclOr(mclVv(isleft, "isleft"), mclVv(isright, "isright"))));
        /*
         * end
         */
        }
        /*
         * 
         * subrows = rows(isleft & ~ismissing);  % left child node
         */
        mlfAssign(
          &subrows,
          mclArrayRef1(
            mclVa(rows, "rows"),
            mclAnd(
              mclVv(isleft, "isleft"), mclNot(mclVv(ismissing, "ismissing")))));
        /*
         * if ~isempty(subrows)
         */
        if (mclNotBool(mlfIsempty(mclVv(subrows, "subrows")))) {
            /*
             * nodes = doapply(Tree,X,subrows,kids(1),nodes,subtrees,prunelist,endcol);
             */
            mlfAssign(
              &nodes,
              mlfTreeval_doapply(
                mclVa(Tree, "Tree"),
                mclVa(X, "X"),
                mclVv(subrows, "subrows"),
                mclIntArrayRef1(mclVv(kids, "kids"), 1),
                mclVa(nodes, "nodes"),
                mclVa(subtrees, "subtrees"),
                mclVa(prunelist, "prunelist"),
                mclVa(endcol, "endcol")));
        /*
         * end
         */
        }
        /*
         * 
         * subrows = rows(isright & ~ismissing); % right child node
         */
        mlfAssign(
          &subrows,
          mclArrayRef1(
            mclVa(rows, "rows"),
            mclAnd(
              mclVv(isright, "isright"),
              mclNot(mclVv(ismissing, "ismissing")))));
        /*
         * if ~isempty(subrows)
         */
        if (mclNotBool(mlfIsempty(mclVv(subrows, "subrows")))) {
            /*
             * nodes = doapply(Tree,X,subrows,kids(2),nodes,subtrees,prunelist,endcol);
             */
            mlfAssign(
              &nodes,
              mlfTreeval_doapply(
                mclVa(Tree, "Tree"),
                mclVa(X, "X"),
                mclVv(subrows, "subrows"),
                mclIntArrayRef1(mclVv(kids, "kids"), 2),
                mclVa(nodes, "nodes"),
                mclVa(subtrees, "subtrees"),
                mclVa(prunelist, "prunelist"),
                mclVa(endcol, "endcol")));
        /*
         * end
         */
        }
        /*
         * 
         * subrows = rows(ismissing);            % missing, treat as leaf
         */
        mlfAssign(
          &subrows,
          mclArrayRef1(mclVa(rows, "rows"), mclVv(ismissing, "ismissing")));
        /*
         * if ~isempty(subrows)
         */
        if (mclNotBool(mlfIsempty(mclVv(subrows, "subrows")))) {
            /*
             * nodes(subrows,1:endcol) = thisnode;
             */
            mclArrayAssign2(
              &nodes,
              mclVa(thisnode, "thisnode"),
              mclVv(subrows, "subrows"),
              mlfColon(_mxarray20_, mclVa(endcol, "endcol"), NULL));
        /*
         * end
         */
        }
    /*
     * end
     */
    }
    mclValidateOutput(nodes, 1, nargout_, "nodes", "treeval/doapply");
    mxDestroyArray(splitvar);
    mxDestroyArray(cutoff);
    mxDestroyArray(assignedclass);
    mxDestroyArray(kids);
    mxDestroyArray(catsplit);
    mxDestroyArray(prunelevel);
    mxDestroyArray(ntrees);
    mxDestroyArray(ncols);
    mxDestroyArray(x);
    mxDestroyArray(isleft);
    mxDestroyArray(isright);
    mxDestroyArray(ismissing);
    mxDestroyArray(subrows);
    mxDestroyArray(endcol);
    mxDestroyArray(prunelist);
    mxDestroyArray(subtrees);
    mxDestroyArray(thisnode);
    mxDestroyArray(rows);
    mxDestroyArray(X);
    mxDestroyArray(Tree);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return nodes;
    /*
     * 
     * 
     */
}
