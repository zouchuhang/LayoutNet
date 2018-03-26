//
// The source file for wrapping the ANN main functionality
// 
// It is based on Approximate Nearest Neighborhood ver 1.1
//
// Note: if you would like to mex this file, please incorporate the
// source files of ANN (v 1.1 or compatible versions)
//
// Author: Dahua Lin
// Date:   Apr 21, 2006
//


#include "mex.h"

#include "ANN/ANN.h"

#include <memory.h> // for memory copying

#pragma comment(lib, "ANN.lib")


// Auxiliary Functions

// allocate an ANN point array and copy points from mxArray to ANN point array
ANNpointArray CreateANNPointArray(int d, int n, const mxArray* arr)
{
    // allocate the ANN array
    ANNpointArray pts = annAllocPts(n, d);
    
    // copy the point coordinates
    ::memcpy(pts[0], mxGetPr(arr), sizeof(double) * d * n);    
    
    return pts;
}

// release the ANN point array allocated
void ReleaseANNPointArray(ANNpointArray pts)
{
    annDeallocPts(pts);
}

// create mxArray from double array
mxArray* GetMxArrayFromDoubleArray(int m, int n, const double* src)
{
    mxArray* M = mxCreateDoubleMatrix(m, n, mxREAL);
    int nelems = m * n;    
    ::memcpy(mxGetPr(M), src, sizeof(double) * nelems);  
    
    return M;
}


// create mxArray from integer array
mxArray* GetMxArrayFromIntArray(int m, int n, const int* src)
{
    mxArray* M = mxCreateDoubleMatrix(m, n, mxREAL);
    int nelems = m * n;
    
    double* pd = mxGetPr(M);
    for (int i = 0; i < nelems; ++i)
    {
        pd[i] = (double)src[i];
    }    
    
    return M;
}



// Search Functions

typedef void (*knnSearchFunc)(ANNkd_tree* tr, ANNpointArray pts, int n, int k, ANNidxArray nn_idx, ANNdistArray dists, double eb);

void KnnSearch_Normal(ANNkd_tree* tr, ANNpointArray pts, int n, int k, ANNidxArray nn_idx, ANNdistArray dists, double eb)
{
    ANNidx* cur_nn_idx = nn_idx;
    ANNdist* cur_dists = dists;
      
    for (int i = 0; i < n; ++i)
    {
        tr->annkSearch(pts[i], k, cur_nn_idx, cur_dists, eb);
        
        cur_nn_idx += k;
        cur_dists += k;
    }
}

void KnnSearch_PriorSearch(ANNkd_tree* tr, ANNpointArray pts, int n, int k, ANNidxArray nn_idx, ANNdistArray dists, double eb)
{
    ANNidx* cur_nn_idx = nn_idx;
    ANNdist* cur_dists = dists;
      
    for (int i = 0; i < n; ++i)
    {
        tr->annkPriSearch(pts[i], k, cur_nn_idx, cur_dists, eb);
        
        cur_nn_idx += k;
        cur_dists += k;
    }
}


knnSearchFunc SelectSearchFunc(int i)
{
    switch(i)
    {
        case 0:
            return KnnSearch_Normal;
        case 1:
            return KnnSearch_PriorSearch;
    }
    
    return 0;
}


///////////////////////////////////////////////////////////////////////////////
//
// the main entry for Matlab binary
// Input:
//  nrhs = 6
//      prhs[0]:  the points for building the KD Tree [d x n0 matrix]
//      prhs[1]:  the points whose neighbors are queried [d x n matrix]
//      prhs[2]:  the size of each neighborhood
//      prhs[3]:  the error bound
//      prhs[4]:  the index of the split rule to use (follow the ANN Library)
//          - 0: ANN_KD_STD (the optimized kd-splitting rule)
//          - 1: ANN_KD_MIDPT (midpoint split)
//          - 2: ANN_KD_FAIR (fair split)
//          - 3: ANN_KD_SL_MIDPT (sliding midpoint splitting)
//          - 4: ANN_KD_SL_FAIR (sliding fair splitting)
//          - 5: ANN_KD_SUGGEST (the authors' suggestion for best)
//      prhs[5]:  the index of search method
//          - 0: Normal approximate KNN search
//          - 1: Priority K near neighbor search
// Output:
//  nlhs <= 2
//      plhs[0]:  the array of indices for query points [k x n matrix]
//      plhs[1]:  the array of distance values [k x n matrix]
//
// Remarks:
//  No pre-condition checking is performed. If the condition is violated,
//  the behavior of the program is undefined. It is the responsible of the
//  invoker to guarantee that the pre-conditions are all satisfied.
//
///////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // name the input variables
    const mxArray *MP0         = prhs[0];     // the matrix of base points
    const mxArray *MP          = prhs[1];     // the matrix of query points
    const mxArray *MNNSize     = prhs[2];     // the scalar for K value
    const mxArray *MErrBound   = prhs[3];     // the scalar for error bound
    const mxArray *MISplitRule = prhs[4];     // the scalar for split rule
    const mxArray *MISearchWay = prhs[5];     // the scalar for search method
    
    // get basic information
    int d = mxGetM(MP0);        // point vector dimension
    int n0 = mxGetN(MP0);       // the number of base points
    int n = mxGetN(MP);         // the number of query points
    int k = (int)mxGetScalar(MNNSize);      // the K value (neighborhood size)
    double errbound = (double)mxGetScalar(MErrBound);   // the error bound
    int iSplitRule = (int)mxGetScalar(MISplitRule);     // the index of split rule
    int iSearchMethod = (int)mxGetScalar(MISearchWay);  // the index of search method
    
    // prepare the point arrays for ANN Library
    ANNpointArray pts0 = CreateANNPointArray(d, n0, MP0);
    ANNpointArray pts  = CreateANNPointArray(d, n, MP);
    
    // build the KD Tree (with bucket size using default value 1)
    ANNkd_tree *kdtr = new ANNkd_tree(pts0, n0, d, 1, (ANNsplitRule)iSplitRule);
    
    // prepare the ANN output
    ANNidxArray nn_idx = new ANNidx[k * n];
    ANNdistArray dists = new ANNdist[k * n];
    
    // conduct KNN search
    knnSearchFunc fs = SelectSearchFunc(iSearchMethod);
    fs(kdtr, pts, n, k, nn_idx, dists, errbound);
    
    // set Matlab outputs
    if (nlhs >= 1)  // output indices
    {
        plhs[0] = GetMxArrayFromIntArray(k, n, nn_idx);
    }
    if (nlhs >= 2) // output distances
    {
        plhs[1] = GetMxArrayFromDoubleArray(k, n, dists);
    }
        
    // release the KD Tree
    delete kdtr;
    
    // release ANN outputs
    delete[] nn_idx;
    delete[] dists;
    
    // release the point arrays
    ReleaseANNPointArray(pts0);
    ReleaseANNPointArray(pts);
    
    // finalize the ANN library
    annClose();
    
}


