/*
 * pdistmex.c
 *
 * Calculates pairwise distances between observations.
 * Helper function to pdist.m
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1993-2004 The MathWorks, Inc.
 */

/* $Revision: 1.1.6.3 $  $Date: 2004/01/24 09:36:23 $ */

#include "mex.h"
#include <math.h>
#include <string.h>

/* Euclidean distance */
void eucdist(double *x, int m, int n, double *d)
{
    /*
     d = sqrt(sum((XI-XJ).^2,2));            % Euclidean
     */
    int i,j,k;
    double theSum,Y;
    double *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (*XI)-(*XJ);
                theSum += Y*Y;
            }
            *(d++) = sqrt(theSum);
        }
    }
}

/* Standardized Euclidean distance */
void seudist(double *x, int m, int n, double *arg, double *d)
{
    /*
     d = sqrt(((XI-XJ).^2) * arg);           % Standardized Euclidean
     */
    /* note that arg is a column vector of 1/var(X) */
    int i,j,k;
    double theSum,Y;
    double *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (*XI)-(*XJ);
                theSum += Y*Y*(*(arg+k));
            }
            *(d++) = sqrt(theSum);
        }
    }
}

/* City Block Distance */
void citdist(double *x, int m, int n, double *d)
{
    /*
     d = sum(abs((XI-XJ)),2);                % City Block
     */
    int i,j,k;
    double theSum,Y;
    double *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (*XI)-(*XJ);
                theSum += fabs(Y);
            }
            *(d++) = theSum;
        }
    }
}

/* Mahalanobis distance */
void mahdist(double *x, int m, int n, double *arg, double *d)
{
    /*
     Y = XI - XJ;
     d = sqrt(sum((Y*arg).*Y,2));            % Mahalanobis
     */
    /* arg is inv(cov(X) so we have to do the matrix multiplication */
    int i,j,k,l;
    double theSum, inner,Y,YY;
    double *XXI, *XXJ;
    double *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;

            for (l=0;l<n; l++,XI++,XJ++){
                XXJ = x+j*n;
                XXI = x+i*n;
                inner = 0;
                for (k=0;k<n;k++,XXI++,XXJ++){
                    YY = (*XXI)-(*XXJ);
                    inner += YY*arg[k+l*n];
                }
                Y = (*XI)-(*XJ);
                theSum += inner * Y;
            }
            *(d++) = sqrt(theSum);
        }
    }
}

/* Minkowski distance */
void mindist(double *x, int m, int n, double arg, double *d)
{
    /*
     d = sum(abs((XI-XJ)).^arg,2).^(1/arg);  % Minkowski
     */
    int i,j,k;
    double theSum,Y;
    double argRecip = 1/arg;
    double *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = fabs((*XI)-(*XJ));
                theSum += pow(Y,arg);
            }
            *(d++) = pow(theSum,argRecip);
        }
    }
}

/* Cosine and Correlation distances */
void coscordist(double *x, int m, int n, double *d)
{
    /*
     d = 1 - sum(XI.*XJ,2);   % Cosine  & Cor
     */
    int i,j,k;
    double theSum;
    double *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                theSum += (*XI)*(*XJ);
            }
            *(d++) = 1-theSum;
        }
    }
}

/* Hamming distance */
void hamdist(double *x, int m, int n, double *d)
{
    /*
     d = sum(XI ~= XJ,2) / size(XI,2);       % Hamming
     */
    int i,j,k;
    double theSum;
    double *XI, *XJ, *XI0;

    double *theNaN = mxCalloc(m,sizeof(double));

    /* Determine which rows of X have missing data.
     */
    XI = x;
    for (i=0; i<m; i++) {
        for (k=0;k<n;k++,XI++) if (mxIsNaN(*XI)) theNaN[i] = *XI;
    }

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            /* If XI or XJ have missing data, set their distance to NaN.
             */
            if (theNaN[i] || theNaN[j]) {
                XI += n;
                XJ += n;
                *(d++) = theNaN[i] + theNaN[j];
            } else {
                theSum = 0;
                for (k=0;k<n;k++,XI++,XJ++){
                    if ((*XI)!=(*XJ)) {
                        theSum++;
                    }
                }
                *(d++) = theSum/n;
            }
        }
    }
    mxFree(theNaN);
}

/* Jaccard distance */
void jacdist(double *x, int m, int n, double *d)
{
    /*
     nz = XI ~= 0 | XJ ~= 0;
     ne = XI ~= XJ;
     d = sum(ne&nz,2) ./ sum(nz,2);          % Jaccard
     */
    int i,j,k;
    double theSum,nz;
    double *XI, *XJ, *XI0;

    double *theNaN = mxCalloc(m,sizeof(double));

    /* Determine which rows of X have missing data.
     */
    XI = x;
    for (i=0; i<m; i++) {
        for (k=0;k<n;k++,XI++) if (mxIsNaN(*XI)) theNaN[i] = *XI;
    }

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            /* If XI or XJ have missing data, set their distance to NaN.
             */
            if (theNaN[i] || theNaN[j]) {
                XI += n;
                XJ += n;
                *(d++) = theNaN[i] + theNaN[j];
            } else {
                theSum = 0;
                nz = 0;
                for (k=0;k<n;k++,XI++,XJ++){
                    if ((*XI) || (*XJ)) {
                        nz++;
                        if ((*XI)!=(*XJ)) {
                            theSum++;
                        }
                    }
                }
                if (nz) {
                    *(d++) = theSum/nz;
                } else {
                    *(d++) = mxGetNaN();
                }
            }
        }
    }
    mxFree(theNaN);
}

/* Chebychev distance */
void chedist(double *x, int m, int n, double *d)
{
    /*
     d = max(abs(XI-XJ),[],2);
     */
    int i,j,k;
    double theMax,Y;
    double *XI, *XJ, *XI0;

    double *theNaN = mxCalloc(m,sizeof(double));

    /* Determine which rows of X have missing data.
     */
    XI = x;
    for (i=0; i<m; i++) {
        for (k=0;k<n;k++,XI++) if (mxIsNaN(*XI)) theNaN[i] = *XI;
    }

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            /* If XI or XJ have missing data, set their distance to NaN.
             */
            if (theNaN[i] || theNaN[j]) {
                XI += n;
                XJ += n;
                *(d++) = theNaN[i] + theNaN[j];
            } else {
                theMax = 0;
                for (k=0;k<n;k++,XI++,XJ++){
                    Y = fabs((*XI)-(*XJ));
                    if (Y>theMax) {
                        theMax = Y;
                    }
                }
                *(d++) = theMax;
            }
        }
    }
    mxFree(theNaN);
}

/************************************************************/
/* Euclidean distance */
void eucdistS(float *x, int m, int n, float *d)
{
    /*
     d = sqrt(sum((XI-XJ).^2,2));            % Euclidean
     */
    int i,j,k;
    float theSum,Y;
    float *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (*XI)-(*XJ);
                theSum += Y*Y;
            }
            *(d++) = (float)sqrt(theSum);
        }
    }
}

/* Standardized Euclidean distance */
void seudistS(float *x, int m, int n, float *arg, float *d)
{
    /*
     d = sqrt(((XI-XJ).^2) * arg);           % Standardized Euclidean
     */
    int i,j,k;
    float theSum,Y;
    float *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (*XI)-(*XJ);
                theSum += Y*Y*(*(arg+k));
            }
            *(d++) = (float)sqrt(theSum);
        }
    }
}

/* City Block Distance */
void citdistS(float *x, int m, int n, float *d)
{
    /*
     d = sum(abs((XI-XJ)),2);                % City Block
     */
    int i,j,k;
    float theSum,Y;
    float *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (*XI)-(*XJ);
                theSum += (float)fabs(Y);
            }
            *(d++) = theSum;
        }
    }
}

/* Mahalanobis distance */
void mahdistS(float *x, int m, int n, float *arg, float *d)
{
    /*
     Y = XI - XJ;
     d = sqrt(sum((Y*arg).*Y,2));            % Mahalanobis
     */

    int i,j,k,l;
    float theSum, inner,Y,YY;
    float *XXI, *XXJ;
    float *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;

            for (l=0;l<n; l++,XI++,XJ++){
                XXJ = x+j*n;
                XXI = x+i*n;
                inner = 0;
                for (k=0;k<n;k++,XXI++,XXJ++){
                    YY = (*XXI)-(*XXJ);
                    inner += YY*arg[k+l*n];
                }
                Y = (*XI)-(*XJ);
                theSum += inner * Y;
            }
            *(d++) = (float)sqrt(theSum);
        }
    }
}

/* Minkowski distance */
void mindistS(float *x, int m, int n, float arg, float *d)
{
    /*
     d = sum(abs((XI-XJ)).^arg,2).^(1/arg);  % Minkowski
     */
    int i,j,k;
    float theSum,Y;
    float argRecip = 1/arg;
    float *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (float)fabs((*XI)-(*XJ));
                theSum += (float)pow(Y,arg);
            }
            *(d++) = (float)pow(theSum,argRecip);
        }
    }
}

/* Cosine and Correlation distances */
void coscordistS(float *x, int m, int n, float *d)
{
    /*
     d = 1 - sum(XI.*XJ,2);   % Cosine  & Cor
     */
    int i,j,k;
    float theSum;
    float *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                theSum += (*XI)*(*XJ);
            }
            *(d++) = 1-theSum;
        }
    }
}

/* Hamming distance */
void hamdistS(float *x, int m, int n, float *d)
{
    /*
     d = sum(XI ~= XJ,2) / size(XI,2);       % Hamming
     */
    int i,j,k;
    float theSum;
    float *XI, *XJ, *XI0;

    float *theNaN = mxCalloc(m,sizeof(float));

    /* Determine which rows of X have missing data.
     */
    XI = x;
    for (i=0; i<m; i++) {
        for (k=0;k<n;k++,XI++) if (mxIsNaN(*XI)) theNaN[i] = *XI;
    }

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            /* If XI or XJ have missing data, set their distance to NaN.
             */
            if (theNaN[i] || theNaN[j]) {
                XI += n;
                XJ += n;
                *(d++) = theNaN[i] + theNaN[j];
            } else {
                theSum = 0;
                for (k=0;k<n;k++,XI++,XJ++){
                    if ((*XI)!=(*XJ)) {
                        theSum++;
                    }
                }
                *(d++) = theSum/n;
            }
        }
    }
    mxFree(theNaN);
}

/* Jaccard distance */
void jacdistS(float *x, int m, int n, float *d)
{
    /*
     nz = XI ~= 0 | XJ ~= 0;
     ne = XI ~= XJ;
     d = sum(ne&nz,2) ./ sum(nz,2);          % Jaccard
     */
    int i,j,k;
    float theSum,nz;
    float *XI, *XJ, *XI0;

    float *theNaN = mxCalloc(m,sizeof(float));

    /* Determine which rows of X have missing data.
     */
    XI = x;
    for (i=0; i<m; i++) {
        for (k=0;k<n;k++,XI++) if (mxIsNaN(*XI)) theNaN[i] = *XI;
    }

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            /* If XI or XJ have missing data, set their distance to NaN.
             */
            if (theNaN[i] || theNaN[j]) {
                XI += n;
                XJ += n;
                *(d++) = theNaN[i] + theNaN[j];
            } else {
                theSum = 0;
                nz = 0;
                for (k=0;k<n;k++,XI++,XJ++){
                    if ((*XI) || (*XJ)) {
                        nz++;
                        if ((*XI)!=(*XJ)) {
                            theSum++;
                        }
                    }
                }
                if (nz) {
                    *(d++) = theSum/nz;
                } else {
                    *(d++) = (float)mxGetNaN();
                }
            }
        }
    }
    mxFree(theNaN);
}

/* Chebychev distance */
void chedistS(float *x, int m, int n, float *d)
{
    /*
     d = max(abs(XI-XJ),[],2);
     */
    int i,j,k;
    float theMax,Y;
    float *XI, *XJ, *XI0;

    float *theNaN = mxCalloc(m,sizeof(float));

    /* Determine which rows of X have missing data.
     */
    XI = x;
    for (i=0; i<m; i++) {
        for (k=0;k<n;k++,XI++) if (mxIsNaN(*XI)) theNaN[i] = *XI;
    }

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            /* If XI or XJ have missing data, set their distance to NaN.
             */
            if (theNaN[i] || theNaN[j]) {
                XI += n;
                XJ += n;
                *(d++) = theNaN[i] + theNaN[j];
            } else {
                theMax = 0;
                for (k=0;k<n;k++,XI++,XJ++){
                    Y = (float)fabs((*XI)-(*XJ));
                    if (Y>theMax) {
                        theMax = Y;
                    }
                }
                *(d++) = theMax;
            }
        }
    }
    mxFree(theNaN);
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int     status,numCoords,numPoints;
    char	metric[4];
    double  *argD;

  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
    if (nrhs<2) {
        mexErrMsgIdAndTxt("stats:pdistmex:TooFewInputs",
                          "Two input arguments required.");
    } else if(nlhs>1) {
        mexErrMsgIdAndTxt("stats:pdistmex:TooManyOutputs",
                          "Too many output arguments.");
    }

    /*  get the metric */
    status = mxGetString(prhs[1],metric,4);

  /* Check the type of the input array */
  /* Currently only works with double or single(float) */
    if (mxIsDouble(prhs[0])) {
        double *x,*d,*arg,scalarArg;
  /*  create a pointer to the input matrix y */
        x = mxGetPr(prhs[0]);

  /*  get the dimensions of the matrix input y */
        numCoords = mxGetM(prhs[0]);
        numPoints = mxGetN(prhs[0]);
  /* get extra arg  */
        if (nrhs>2 && !mxIsEmpty(prhs[2])) {
            if (mxGetNumberOfElements(prhs[2]) == 1) {  /*scalar case */
                scalarArg = mxGetScalar(prhs[2]);
            } else if (mxIsDouble(prhs[2])) {
                arg = mxGetPr(prhs[2]);
            } else {
                mexErrMsgIdAndTxt("stats:pdistmex:MixedInputTypes",
                                  "Additional input arguments must be the same class as X.");
            }
        }
  /*  set the output pointer to the output matrix */
        plhs[0] = mxCreateDoubleMatrix(1,(numPoints * (numPoints-1))/2, mxREAL);

  /*  create a pointer to a copy of the output matrix */
        d = mxGetPr(plhs[0]);

  /*  call the appropriate distance subroutine */
        if (strcmp(metric,"euc") == 0)
            eucdist(x,numPoints,numCoords,d);
        else if(strcmp(metric,"seu") == 0)
            seudist(x,numPoints,numCoords,arg,d);
        else if(strcmp(metric,"cit") == 0)
            citdist(x,numPoints,numCoords,d);
        else if(strcmp(metric,"min") == 0)
            mindist(x,numPoints,numCoords,scalarArg,d);
        else if(strcmp(metric,"cos") == 0)
            coscordist(x,numPoints,numCoords,d);
        else if(strcmp(metric,"cor") == 0)
            coscordist(x,numPoints,numCoords,d);
        else if(strcmp(metric,"ham") == 0)
            hamdist(x,numPoints,numCoords,d);
        else if(strcmp(metric,"jac") == 0)
            jacdist(x,numPoints,numCoords,d);
        else if(strcmp(metric,"che") == 0)
            chedist(x,numPoints,numCoords,d);
        else if(strcmp(metric,"mah") == 0)
            mahdist(x,numPoints,numCoords,arg,d);

    /* deal with non double types */
    } else if (mxIsSingle(prhs[0])) {
        float *x,*d,*arg,scalarArg;

         /*  create a pointer to the input matrix y */
        x = (float*)mxGetData(prhs[0]);

  /*  get the dimensions of the matrix input y */
        numCoords = mxGetM(prhs[0]);
        numPoints = mxGetN(prhs[0]);

  /* get extra arg  */
        if (nrhs>2 && !mxIsEmpty(prhs[2])) {
            if (mxGetNumberOfElements(prhs[2]) == 1) {  /*scalar case */
                scalarArg = (float)mxGetScalar(prhs[2]);
            } else if (mxIsSingle(prhs[2])) {
                arg = mxGetData(prhs[2]);
            } else {
                mexErrMsgIdAndTxt("stats:pdistmex:MixedInputTypes",
                                  "Additional input arguments must be the same class as X.");
            }
        }

  /*  set the output pointer to the output matrix */
        plhs[0] = mxCreateNumericMatrix(1,(numPoints * (numPoints-1))/2, mxSINGLE_CLASS, mxREAL);

  /*  create a pointer to a copy of the output matrix */
        d = (float*)mxGetData(plhs[0]);

  /*  call the appropriate distance subroutine */
        if (strcmp(metric,"euc") == 0)
            eucdistS(x,numPoints,numCoords,d);
        else if(strcmp(metric,"seu") == 0)
            seudistS(x,numPoints,numCoords,arg,d);
        else if(strcmp(metric,"cit") == 0)
            citdistS(x,numPoints,numCoords,d);
        else if(strcmp(metric,"min") == 0)
            mindistS(x,numPoints,numCoords,scalarArg,d);
        else if(strcmp(metric,"cos") == 0)
            coscordistS(x,numPoints,numCoords,d);
        else if(strcmp(metric,"cor") == 0)
            coscordistS(x,numPoints,numCoords,d);
        else if(strcmp(metric,"ham") == 0)
            hamdistS(x,numPoints,numCoords,d);
        else if(strcmp(metric,"jac") == 0)
            jacdistS(x,numPoints,numCoords,d);
        else if(strcmp(metric,"che") == 0)
            chedistS(x,numPoints,numCoords,d);
        else if(strcmp(metric,"mah") == 0)
            mahdistS(x,numPoints,numCoords,arg,d);

    } else {
        mexErrMsgIdAndTxt("stats:pdistmex:BadInputType",
                          "PDISTMEX only supports real DOUBLE and SINGLE data.");
    }

}
/*
% ----------------------------------------------
function d = distcalc(XI,XJ,s,arg)
%DISTCALC Perform distance calculation for PDIST.
switch s
case 'euc'    d = sqrt(sum((XI-XJ).^2,2));            % Euclidean
case 'seu'    d = sqrt(((XI-XJ).^2) * arg);           % Standardized Euclidean
case 'cit'    d = sum(abs((XI-XJ)),2);                % City Block
case 'mah'    Y = XI - XJ;
              d = sqrt(sum((Y*arg).*Y,2));            % Mahalanobis
case 'min'    d = sum(abs((XI-XJ)).^arg,2).^(1/arg);  % Minkowski
case 'che'    d = max(abs(XI-XJ),[],2);               % Chebychev
case 'cos'    d = 1 - sum(XI.*XJ,2);                  % Cosine
case 'cor'    d = 1 - sum(XI.*XJ,2);                  % Correlation
case 'ham'    d = sum(XI ~= XJ,2) / size(XI,2);       % Hamming
case 'jac'    nz = XI ~= 0 | XJ ~= 0;
              ne = XI ~= XJ;
              d = sum(ne&nz,2) ./ sum(nz,2);          % Jaccard
end
*/
