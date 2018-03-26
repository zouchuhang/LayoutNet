#include "mex.h"
#include "math.h"

#define min(X,Y) ((X) < (Y) ? (X) : (Y)) 
#define max(X,Y) ((X) > (Y) ? (X) : (Y))

int mindex3(int x, int y, int z, int sizx, int sizy, int sizz, int wrap) 
{ 
    int index;
 	if(wrap==1)
	{
        /* Positive modules */
		x=(x % sizx + sizx) % sizx;
        y=(y % sizy + sizy) % sizy;
        z=(z % sizz + sizz) % sizz;
    }
    else if(wrap>1)
    {
        /* Clamp */
        x=max(x,0);
        y=max(y,0); 
        z=max(z,0);
        x=min(x,sizx-1);
        y=min(y,sizy-1);
        z=min(z,sizz-1);
    }
    index=z*sizx*sizy+y*sizx+x;
    return index;
}

mxLogical *draw_or_split(mxLogical *Volume,double AX,double AY,double AZ,double BX,double BY,double BZ,double CX,double CY,double CZ, int *VolumeSize, int wrap)
{
    bool checkA, checkB, checkC;
    bool check1, check2, check3, check4, check5, check6;

    double dist1, dist2, dist3, maxdist;
    double DX,DY,DZ;
    /* Check if vertices outside */
    if(wrap==0)
	{
        checkA=(AX<0)||(AY<0)||(AZ<0)||(AX>(VolumeSize[0]-1))||(AY>(VolumeSize[1]-1))||(AZ>(VolumeSize[2]-1));
        checkB=(BX<0)||(BY<0)||(BZ<0)||(BX>(VolumeSize[0]-1))||(BY>(VolumeSize[1]-1))||(BZ>(VolumeSize[2]-1));
        checkC=(CX<0)||(CY<0)||(CZ<0)||(CX>(VolumeSize[0]-1))||(CY>(VolumeSize[1]-1))||(CZ>(VolumeSize[2]-1));

        check1=(AX<0)&&(BX<0)&&(CX<0);
		check2=(AY<0)&&(BY<0)&&(CY<0);
		check3=(AZ<0)&&(BZ<0)&&(CZ<0);
		check4=(AX>(VolumeSize[0]-1))&&(BX>(VolumeSize[0]-1))&&(CX>(VolumeSize[0]-1));
		check5=(AY>(VolumeSize[1]-1))&&(BY>(VolumeSize[1]-1))&&(CY>(VolumeSize[1]-1));
		check6=(AZ>(VolumeSize[2]-1))&&(BZ>(VolumeSize[2]-1))&&(CZ>(VolumeSize[2]-1));
		
		/* Return if all vertices outside, on the same side */
		if(check1||check2||check3||check4||check5||check6)
		{
			return Volume;
		}
	}
	
    dist1=(AX-BX)*(AX-BX)+(AY-BY)*(AY-BY)+(AZ-BZ)*(AZ-BZ);
    dist2=(CX-BX)*(CX-BX)+(CY-BY)*(CY-BY)+(CZ-BZ)*(CZ-BZ);
    dist3=(AX-CX)*(AX-CX)+(AY-CY)*(AY-CY)+(AZ-CZ)*(AZ-CZ);
    if(dist1>dist2)
    {
        if(dist1>dist3)
        {
            maxdist=dist1;
            if(maxdist>0.5)
            {
                DX=(AX+BX)/2; DY=(AY+BY)/2; DZ=(AZ+BZ)/2;
                Volume=draw_or_split(Volume,DX,DY,DZ,BX,BY,BZ,CX,CY,CZ,VolumeSize, wrap);
                Volume=draw_or_split(Volume,AX,AY,AZ,DX,DY,DZ,CX,CY,CZ,VolumeSize, wrap);
            }  
        }
        else
        {
            maxdist=dist3;
            if(maxdist>0.5)
            {
                DX=(AX+CX)/2; DY=(AY+CY)/2; DZ=(AZ+CZ)/2;
                Volume=draw_or_split(Volume,DX,DY,DZ,BX,BY,BZ,CX,CY,CZ,VolumeSize, wrap);
                Volume=draw_or_split(Volume,AX,AY,AZ,BX,BY,BZ,DX,DY,DZ,VolumeSize, wrap);
            }  

        }
    }
    else
    {
        if(dist2>dist3)
        {
            maxdist=dist2;
            DX=(CX+BX)/2; DY=(CY+BY)/2; DZ=(CZ+BZ)/2;
            if(maxdist>0.5)
            {
                Volume=draw_or_split(Volume,AX,AY,AZ,DX,DY,DZ,CX,CY,CZ,VolumeSize, wrap);
                Volume=draw_or_split(Volume,AX,AY,AZ,BX,BY,BZ,DX,DY,DZ,VolumeSize, wrap);
            }  
        }
        else
        {
            maxdist=dist3;
            if(maxdist>0.5)
            {
                DX=(AX+CX)/2; DY=(AY+CY)/2; DZ=(AZ+CZ)/2;
                Volume=draw_or_split(Volume,DX,DY,DZ,BX,BY,BZ,CX,CY,CZ,VolumeSize, wrap);
                Volume=draw_or_split(Volume,AX,AY,AZ,BX,BY,BZ,DX,DY,DZ,VolumeSize, wrap);
            }  

        }
    }
    if(wrap==0)
    {      
        if(checkA==false)
        {
            Volume[mindex3((int)(AX+0.5),(int)(AY+0.5), (int)(AZ+0.5), VolumeSize[0], VolumeSize[1], VolumeSize[2], wrap)]=1;
        }
        if(checkB==false)
        {
            Volume[mindex3((int)(BX+0.5),(int)(BY+0.5), (int)(BZ+0.5), VolumeSize[0], VolumeSize[1], VolumeSize[2], wrap)]=1;
        }
        if(checkC==false)
        {
            Volume[mindex3((int)(CX+0.5),(int)(CY+0.5), (int)(CZ+0.5), VolumeSize[0], VolumeSize[1], VolumeSize[2], wrap)]=1;
        }
    }
    else
    {
         Volume[mindex3((int)(AX+0.5),(int)(AY+0.5), (int)(AZ+0.5), VolumeSize[0], VolumeSize[1], VolumeSize[2], wrap)]=1;
         Volume[mindex3((int)(BX+0.5),(int)(BY+0.5), (int)(BZ+0.5), VolumeSize[0], VolumeSize[1], VolumeSize[2], wrap)]=1;
         Volume[mindex3((int)(CX+0.5),(int)(CY+0.5), (int)(CZ+0.5), VolumeSize[0], VolumeSize[1], VolumeSize[2], wrap)]=1;
    }
    return Volume;
}



/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
   double *FacesA, *FacesB, *FacesC, *VerticesX, *VerticesY, *VerticesZ, *VolumeSize;
   mxLogical *Volume;
   double *Wrapd;
   int wrap;
   double AX,AY,AZ;
   double BX,BY,BZ;
   double CX,CY,CZ;
   int i;
   int VolumeDims[3]={0,0,0};
   
   const mwSize *FacesDims;
   int FacesN=0;
    
   /* Check for proper number of arguments. */
   if(nrhs!=8) {
     mexErrMsgTxt("Eight inputs are required.");
   } else if(nlhs!=1) {
     mexErrMsgTxt("One output required");
   }
   
   /* Read all inputs */
   FacesA=mxGetPr(prhs[0]);
   FacesB=mxGetPr(prhs[1]);
   FacesC=mxGetPr(prhs[2]);
   VerticesX=mxGetPr(prhs[3]);
   VerticesY=mxGetPr(prhs[4]);
   VerticesZ=mxGetPr(prhs[5]);
   VolumeSize=mxGetPr(prhs[6]);
   Wrapd=mxGetPr(prhs[7]);
   wrap=(int)Wrapd[0];

   FacesDims = mxGetDimensions(prhs[0]);   
   FacesN=FacesDims[0]*FacesDims[1];
   
   /*  Create Output array */
   VolumeDims[0]=(int)VolumeSize[0]; 
   VolumeDims[1]=(int)VolumeSize[1]; 
   VolumeDims[2]=(int)VolumeSize[2];
   plhs[0] = mxCreateLogicalArray(3, VolumeDims);
   Volume = mxGetLogicals(plhs[0]);
  
     
   for (i=0; i<FacesN; i++)
   {
        AX=VerticesX[(int)FacesA[i]-1]-1.0;
        AY=VerticesY[(int)FacesA[i]-1]-1.0;
        AZ=VerticesZ[(int)FacesA[i]-1]-1.0;
        BX=VerticesX[(int)FacesB[i]-1]-1.0;
        BY=VerticesY[(int)FacesB[i]-1]-1.0;
        BZ=VerticesZ[(int)FacesB[i]-1]-1.0;
        CX=VerticesX[(int)FacesC[i]-1]-1.0;
        CY=VerticesY[(int)FacesC[i]-1]-1.0;
        CZ=VerticesZ[(int)FacesC[i]-1]-1.0;
        Volume=draw_or_split(Volume,AX,AY,AZ,BX,BY,BZ,CX,CY,CZ,VolumeDims,wrap);
   }
}
 