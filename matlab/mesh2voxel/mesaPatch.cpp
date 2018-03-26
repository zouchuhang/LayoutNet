/*
Render patch, get depths. Can do min or max

Partly modified from Princeton vision toolbox
http://mit.edu/jxiao/
*/

#include "mex.h" 
#include <GL/osmesa.h>
#include <GL/glu.h>


// Input: 
//     arg0: 3x4 Projection matrix, 
//     arg1: image width, 
//     arg2: image height, 
//     arg3: 3xn double vertices matrix, 
//     arg4: 4xn uint32 face matrix, index from zero
// Output: you will need to transpose the result in Matlab manually
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // mexPrintf("RenderMex\n"); 

  double* proj_matrix = mxGetPr(prhs[0]); // 3x4 matrix
  double* mv_matrix   = mxGetPr(prhs[1]); // 3x4 matrix
  int m_width         = (int)mxGetScalar(prhs[2]);
  int m_height        = (int)mxGetScalar(prhs[3]);
  double*           vertex = mxGetPr(prhs[4]); // 3xn double vertices matrix
  unsigned int  num_vertex = mxGetN(prhs[4]);
  unsigned int*       face = (unsigned int*) mxGetData(prhs[5]); // 4xn uint32 face matrix
  unsigned int    num_face = mxGetN(prhs[5]);
  int flag_max             = (int)mxGetScalar(prhs[6]);

  // Step 1: setup off-screen mesa's binding 
  OSMesaContext ctx;
  ctx = OSMesaCreateContextExt(OSMESA_RGB, 32, 0, 0, NULL );
  unsigned char * pbuffer = new unsigned char [3 * m_width * m_height];
  // Bind the buffer to the context and make it current
  if (!OSMesaMakeCurrent(ctx, (void*)pbuffer, GL_UNSIGNED_BYTE, m_width, m_height)) {
    mexErrMsgTxt("OSMesaMakeCurrent failed!: ");
  }
  OSMesaPixelStore(OSMESA_Y_UP, 0);
  
  // Step 2: Setup basic OpenGL setting
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glDisable(GL_CULL_FACE);
  if (flag_max) {
    glDepthFunc(GL_GREATER);
    glClearDepth(0.0f);
  }
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, m_width, m_height);

  // Step 3: Set projection matrices
  // matrix is ready. use it
  glMatrixMode(GL_PROJECTION);
  glLoadMatrixd(proj_matrix);
  glMatrixMode(GL_MODELVIEW);
  glLoadMatrixd(mv_matrix);

  // Step 3: render the mesh with encoded color from their ID
  unsigned char colorBytes[3];
  unsigned int base_offset;

  // render face
  // printf("face begin %d\n", num_face);
  base_offset = 1;
  for (unsigned int i = 0; i < num_face; ++i) {
    glBegin(GL_TRIANGLES);
    glVertex3dv(vertex+3*(*face++));
    glVertex3dv(vertex+3*(*face++));
    glVertex3dv(vertex+3*(*face++));
    glEnd();
  }

  glFinish(); // done rendering

  unsigned int* pDepthBuffer;
  GLint outWidth, outHeight, bitPerDepth;
  OSMesaGetDepthBuffer(ctx, &outWidth, &outHeight, &bitPerDepth, (void**)&pDepthBuffer);
  // mexPrintf("w = %d, h = %d, bitPerDepth = %d\n", outWidth, outHeight, bitPerDepth);
  plhs[0] = mxCreateNumericMatrix((int)outWidth, (int)outHeight, mxUINT32_CLASS, mxREAL);
  unsigned int* result = (unsigned int*) mxGetData(plhs[0]);
  for(int i=0; i<outWidth*outHeight; i++){
    result[i] = pDepthBuffer[i];
  }

  OSMesaDestroyContext(ctx);
  delete [] pbuffer;
} 
