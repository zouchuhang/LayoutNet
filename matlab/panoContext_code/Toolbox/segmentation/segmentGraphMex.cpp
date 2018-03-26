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
	if (nrhs != 6) 
	{
		mexPrintf("Usage: [seg] = segmentGraphMex(width, height, numEdge, edges, threshold, minSize);\n");
		return;
	}
	// convert edges memory from matlab to c++
	double* edgeMat = (double*)mxGetData(prhs[3]);
	int numEdge = (int)mxGetScalar(prhs[2]);
	int width = (int)mxGetScalar(prhs[0]);
	int height= (int)mxGetScalar(prhs[1]);
	double c = mxGetScalar(prhs[4]);
	int min_size = (int)mxGetScalar(prhs[5]);

	edge *edges = new edge[numEdge];
	for( int i = 0; i<numEdge; i++)
	{
		edges[i].a = edgeMat[i*3+0];
		edges[i].b = edgeMat[i*3+1];
		edges[i].w = edgeMat[i*3+2];
	}

	universe *u = segment_graph(width*height, numEdge, edges, c);

	// post process
	for (int i = 0; i < numEdge; i++) 
	{
		int a = u->find(edges[i].a);
		int b = u->find(edges[i].b);
		if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
			u->join(a, b);
	}


	// pass result to matlab
	plhs[0] = mxCreateNumericMatrix((mwSize)height, (mwSize)width, mxDOUBLE_CLASS, mxREAL);
	double* output = (double *)mxGetData(plhs[0]);
	
	for (int x=0; x<width; x++)
	{
		for (int y=0; y<height; y++)
		{
			int comp = u->find(y*width+x);
			int addr = x*height + y;
			output[addr] = (double)comp;
		}
	}
    


	delete[] edges;
	delete u;
}