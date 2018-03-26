/*******************************************************************************
* Sketch Token Toolbox     V0.95
* Copyright 2013 Joseph Lim [lim@csail.mit.edu]
* Please email me if you find bugs, or have suggestions or questions!
* Licensed under the Simplified BSD License [see bsd.txt]
*******************************************************************************/
#include "mex.h"
#include <math.h>
#include <omp.h>

typedef unsigned int uint32;

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  // get inputs
  float *chns = (float*) mxGetData(prhs[0]);
  float *chnsSs = (float*) mxGetData(prhs[1]);
  float *thrs = (float*) mxGetData(prhs[2]);
  uint32 *fids = (uint32*) mxGetData(prhs[3]);
  uint32 *child = (uint32*) mxGetData(prhs[4]);
  float *distr = (float*) mxGetData(prhs[5]);
  uint32 *cids1 =  (uint32*) mxGetData(prhs[6]);
  uint32 *cids2 =  (uint32*) mxGetData(prhs[7]);
  const int stride = (int) mxGetScalar(prhs[8]);
  const int rad = (int) mxGetScalar(prhs[9]);
  const int nChnFtrs = (int) mxGetScalar(prhs[10]);

  // get dimensions and constants
  const mwSize *chnsSize = mxGetDimensions(prhs[0]);
  const int height = (int) chnsSize[0];
  const int width = (int) chnsSize[1];
  const mwSize *distrSize = mxGetDimensions(prhs[5]);
  const int nTokens = (int) distrSize[0];
  const int nTreeNodes = (int) distrSize[1];
  const int nTrees = (int) distrSize[2];
  const int heightOut = (int) ceil((height-rad*2.0)/stride);
  const int widthOut = (int) ceil((width-rad*2.0)/stride);

  // create output
  const int outDims[3]={nTokens,heightOut,widthOut};
  float *S = (float*) mxCalloc(nTokens*heightOut*widthOut,sizeof(float));
  plhs[0] = mxCreateNumericMatrix(0,0,mxSINGLE_CLASS,mxREAL);
  mxSetData(plhs[0],S);
  mxSetDimensions(plhs[0],outDims,3);

  // apply forest to each patch
  #pragma omp parallel for
  for( int c=0; c<widthOut; c++ ) {
      for( int r=0; r<heightOut; r++ ) {
        // classify a single patch using all trees
        float *chns1 = chns + (r*stride) + (c*stride)*height;
        float *chnsSs1 = chnsSs + (r*stride) + (c*stride)*height;
        
        for( int t = 0; t < nTrees; t++ ) {
            uint32 k = t*nTreeNodes, res;
            
            while( child[k] ) {
                // compute feature (either lookup in channel or self-similarity feature)
                uint32 f = fids[k], cid1 = cids1[f], cid2 = cids2[f];
                float ftr;
                
                if( f<nChnFtrs )
                    ftr = chns1[cid1];
                else
                    ftr = chnsSs1[cid1] - chnsSs1[cid2];
                
                // compare ftr to threshold and move left or right accordingly
                if( ftr < thrs[k] )
                    k = child[k]-1;
                else
                    k = child[k];
                
                res = k;
                k += t*nTreeNodes;
            }
            
            // lookup probability and store results
            for( int i = 0; i < nTokens; i++ ) {
                S[r*nTokens + c*heightOut*nTokens + i] += distr[t*nTreeNodes*nTokens + res*nTokens + i];
            }
        }
      }
  }
}
