#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
/** 
 * Return the decision tree node corresponding to the given value set
 *
 * var[n]: the attribute ids for node n
 * cut[n]: the threshold value for node n
 * left_child[n]: the node id of the left child of node n, 0 if node n is terminal
 * right_child[n]: the node id of the right child of node n, 0 if node n is terminal
 * ncatsplit[c]: the number of values resulting in a left branch
 * catsplit[c]: the values that would result in a left branch
 * attributes: the attribute (variable) values for each feature
 **/
void
treevalc(int* var, double* cut, int* left_child, int* right_child, 
	 int* ncatsplit, double** catsplit, 
	 double* attributes,
	 int* node_id) {
 
  int currnode = 0;

  int nextnode;
  int currvar;
  double currval;
  int cid, v;
  int numvals;
  double* vals;

  /*  printf("init nodes: %d  %d  \n", left_child[currnode], right_child[currnode]); */

  /* until reached terminal node */
  while ((left_child[currnode] != 0) && (right_child[currnode] != 0)) {
    
    /*printf("currnode: %d\n", currnode);*/

    nextnode = -1;
    
    currvar = abs(var[currnode])-1;
    currval = attributes[currvar];
    
    /* decision based on thresholded float value */
    if (var[currnode] > 0) {
  
      /*printf("currvar: %d\n", currvar);*/

      /* branch left */
      if (currval < cut[currnode]) {
	nextnode = left_child[currnode];
      }
      /* branch right */
      else {
	nextnode = right_child[currnode];
      }      
    }
    /* decision based on discrete value */
    else {
      numvals = ncatsplit[(int)cut[currnode]-1];
      vals = catsplit[(int)cut[currnode]-1];
      for (v = 0; v < numvals; v++) {
	if (currval == vals[v]) {
	  nextnode = left_child[currnode];
	  break;
	}
      }
      if (nextnode == -1) {
	nextnode = right_child[currnode];
      }
    }
    
    currnode = nextnode-1;
    /* printf("curr node: %d \n", currnode);*/
  }
	  	  
  *node_id = currnode+1;

}

/**
 * plhs = {var, cut, left_child, right_child, catsplit(cell array), attributes(numatt, numdata)}
 *
 */	 
void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])

{

	int *var, *left_child, *right_child, *left_chil, *ncatsplit, n, nsplits,numatt, numdata, tmp_id; 
	double *cut, **catsplit, *all_attributes,  *node_ids;
  if (nrhs != 6) {
    printf("Error: wrong number of input arguments: %d.\n", nlhs);
    printf("Syntax: node_ids = treevalc(var, cut, left_child, right_child, catsplit, attributes)\n");
  }
  
   var = (int*)mxGetPr(prhs[0]);
   cut = mxGetPr(prhs[1]);
   left_child = (int*)mxGetPr(prhs[2]);
   right_child = (int*)mxGetPr(prhs[3]);
  /* get catsplit variables */
  nsplits = mxGetNumberOfElements(prhs[4]);
  ncatsplit = malloc(sizeof(int) * nsplits);
  catsplit = malloc(sizeof(double*) * nsplits);

  n = 0;
  for (n = 0; n < nsplits; n++) {
    mxArray* catsplit_cell_mx = mxGetCell(prhs[4], n);
    if (catsplit_cell_mx == 0) {
      printf("null cell");
    }
    ncatsplit[n] = mxGetNumberOfElements(catsplit_cell_mx);
    catsplit[n] = (double*)mxGetPr(catsplit_cell_mx);
  }

   numatt = mxGetM(prhs[5]);
  numdata = mxGetN(prhs[5]);

  /*  printf("num data = %d   num att = %d\n", numdata, numatt);*/

 all_attributes = mxGetPr(prhs[5]);

  plhs[0] = mxCreateDoubleMatrix(numdata, 1, mxREAL);
  node_ids = mxGetPr(plhs[0]);
  

   tmp_id;
  for (n = 0; n < numdata; n++) {
    treevalc(var, cut, left_child, right_child, ncatsplit, catsplit,
	     &all_attributes[numatt*n], &tmp_id);
    node_ids[n] = (double)(tmp_id);
    /*    printf("final node id: %d\n", tmp_id); */
  }

  free(catsplit);
  free(ncatsplit);
  
}
