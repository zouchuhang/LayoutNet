#include <cstdio>
#include <cstdlib>
#include <image.h>
#include <misc.h>
#include <pnmfile.h>
#include "mex.h"
#include "segment-image.h"

void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[])
{
	// check arguments
	if (nrhs != 5) 
	{
		mexPrintf("Usage: [seg] = segmentGraphMex_edge(maxID, numEdge, edges, threshold, minSize);\n");
		return;
	}
	// convert edges memory from matlab to c++
    int maxID = (int)mxGetScalar(prhs[0]);
    int numEdge = (int)mxGetScalar(prhs[1]);
	double* edgeMat = (double*)mxGetData(prhs[2]);
	double c = mxGetScalar(prhs[3]);
	int min_size = (int)mxGetScalar(prhs[4]);
    printf("maxID: %d, numEdge: %d, c: %f, min_size: %d\n", maxID, numEdge, c, min_size);
    
	edge *edges = new edge[numEdge];
	for( int i = 0; i<numEdge; i++)
	{
		edges[i].a = edgeMat[i*3+0];
		edges[i].b = edgeMat[i*3+1];
		edges[i].w = edgeMat[i*3+2];
	}
    printf("a: %d, b: %d, w: %f\n", edges[0].a, edges[0].b, edges[0].w);
    printf("Loading finished!\n");
	universe *u = segment_graph( maxID, numEdge, edges, c);

    printf("get out of segment_graph\n");
	// post process
	for (int i = 0; i < numEdge; i++) 
	{
		int a = u->find(edges[i].a);
		int b = u->find(edges[i].b);
		if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
			u->join(a, b);
	}

    printf("finish post process\n");
	// pass result to matlab
	plhs[0] = mxCreateNumericMatrix((mwSize)maxID, 1, mxDOUBLE_CLASS, mxREAL);
	double* output = (double *)mxGetData(plhs[0]);
	for (int i = 0; i<maxID; i++)
    {
        output[i] = (double)(u->find(i+1));
    }
    
    printf("packed up output\n");
	delete[] edges;
    printf("delete edges\n");
	//delete u;
    printf("memory released\n");
}