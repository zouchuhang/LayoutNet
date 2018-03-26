#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

//#include <future>
#include "mex.h"

#include <fstream>
#include <vector>
#include <algorithm> 
using namespace std;

#define INF 1000000
//#define MAXCOST 3.0f
#define NUMBTHREAD 4

float MAXCOST = 3.0f;
int verbose = 0;

struct OBJHYP
{
	float fObjSize[3];
	float fObjCenter[3];
	int iObjType;
	
	float fAspectRatioXY;
	float fAreaXY;
};

struct ROOMHYP
{
	int iRoomIdx;
	vector<OBJHYP> vObjects;
};

template <class T>
void ReadData(const mxArray *M, vector<ROOMHYP> &V)
{
	const mwSize *dims = mxGetDimensions(M);
	int hypRoomNum = (int) dims[0];

	V.assign( hypRoomNum, ROOMHYP());
	for (int i = 0; i < hypRoomNum; i++)
	{
		T *objcenter = (T *) mxGetData( mxGetField( M, i, "objcenter") );
		T *objsize = (T *) mxGetData( mxGetField( M, i, "objsize") );
		T *objtype = (T *) mxGetData( mxGetField( M, i, "objtype") );
		T *objaspectXY = (T *) mxGetData( mxGetField( M, i, "objaspectXY") );
		T *objregionXY = (T *) mxGetData( mxGetField( M, i, "objregionXY") );

		const mwSize *d = mxGetDimensions(mxGetField( M, i, "objtype"));
		int numObject = int(d[1]);

		V[i].vObjects.assign( numObject, OBJHYP());
		V[i].iRoomIdx = 0;

		for (int j = 0; j<numObject; j++)
		{
			OBJHYP &obj = V[i].vObjects[j];
			int iStartID = j*3;

			obj.iObjType = objtype[j];
			obj.fAspectRatioXY = (float) objaspectXY[j];
			obj.fAreaXY = (float) objregionXY[j];

			obj.fObjSize[0] = (float) objsize[iStartID];
			obj.fObjSize[1] = (float) objsize[iStartID+1];
			obj.fObjSize[2] = (float) objsize[iStartID+2];

			obj.fObjCenter[0] = (float) objcenter[iStartID];
			obj.fObjCenter[1] = (float) objcenter[iStartID+1];
			obj.fObjCenter[2] = (float) objcenter[iStartID+2];
		}
	}
}

template <typename T>
vector<int> sort_indexes(const vector<T> &v) {

	// initialize original index locations
	vector<int> idx(v.size());
	for (int i = 0; i != idx.size(); ++i) idx[i] = i;

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
		[&v](int i1, int i2) {return v[i1] < v[i2];});

	return idx;
}

void bipartiteMatching(vector<float> IN_costMat, int IN_hypnum, int IN_gndnum, 
					   float &OUT_cost, vector<int> &OUT_hypid, vector<int> &OUT_gndid)
{
	OUT_hypid.reserve(IN_hypnum);
	OUT_gndid.reserve(IN_gndnum);
	OUT_cost = 0;
	
	vector<int> idx = sort_indexes(IN_costMat);

	for ( int i=0; i<idx.size(); i++)
	{
		if (IN_costMat[idx[i]]<MAXCOST)
		{
			int hypid = idx[i]%IN_hypnum;
			int gndid = idx[i]/IN_hypnum;
			OUT_hypid.push_back(hypid);
			OUT_gndid.push_back(gndid);
			OUT_cost += IN_costMat[idx[i]]; //min( MAXCOST, IN_costMat[idx[i]]);

			for (int j=0; j<IN_hypnum; j++)
			{
				IN_costMat[gndid*IN_hypnum+j] = INF;
			}
			for (int j=0; j<IN_gndnum; j++)
			{
				IN_costMat[j*IN_hypnum+hypid] = INF;
			}
		}
	}
}

float objCtrDist( OBJHYP &ob1, OBJHYP &ob2, float *invnorm)
{
	float fDist[3];
	//float fSize[3];
	float buf = 0;
	for ( int i=0; i<3; i++)
	{
		buf = (ob1.fObjCenter[i] - ob2.fObjCenter[i])*invnorm[i];
		fDist[i] = buf*buf;
		//buf = (ob1.fObjSize[i]+ob2.fObjSize[i])/2;
		//fSize[i] = buf*buf;
	}
	return sqrt(fDist[0]+fDist[1]+fDist[2]);
	//return sqrt((fDist[0]+fDist[1]+fDist[2])/(fSize[0]+fSize[1]+fSize[2]));
}

float objVolInte( OBJHYP &ob1, OBJHYP &ob2)
{
	float fMinSize[3];
	for (int i=0; i<3; i++)
	{
		fMinSize[i] = min(ob1.fObjSize[i], ob2.fObjSize[i]);
	}
	return fMinSize[0]*fMinSize[1]*fMinSize[2];
}

float objVolume( OBJHYP &ob)
{
	return ob.fObjSize[0]*ob.fObjSize[1]*ob.fObjSize[2];
}

float objTypeCost( OBJHYP &ob1, OBJHYP &ob2)
{
	if ( ob1.iObjType == ob2.iObjType)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

void roomMatching( ROOMHYP &hyp, ROOMHYP &gnd, float &match_cost)
{
	float hypRoomAR = hyp.vObjects[hyp.iRoomIdx].fAspectRatioXY;
	float gndRoomAR = gnd.vObjects[gnd.iRoomIdx].fAspectRatioXY;
	float hypRoomRe = hyp.vObjects[hyp.iRoomIdx].fAreaXY;
	float gndRoomRe = gnd.vObjects[gnd.iRoomIdx].fAreaXY;

	float RegionProp = hypRoomRe/gndRoomRe;
	if ( RegionProp<0.5 || RegionProp>2)
	{
		match_cost = INF;
		return;
	}
	float AspectProp = hypRoomAR/gndRoomAR;
	if ( AspectProp<0.5 || AspectProp>2)
	{
		match_cost = INF;
		return;
	}

	int HypObjNum = hyp.vObjects.size();
	int GndObjNum = gnd.vObjects.size();
	float ctrDistNorm[3];
	ctrDistNorm[0] = 4.0/min(hyp.vObjects[hyp.iRoomIdx].fObjSize[0], gnd.vObjects[gnd.iRoomIdx].fObjSize[0]);
	ctrDistNorm[1] = 4.0/min(hyp.vObjects[hyp.iRoomIdx].fObjSize[1], gnd.vObjects[gnd.iRoomIdx].fObjSize[1]);
	ctrDistNorm[2] = 4.0/min(hyp.vObjects[hyp.iRoomIdx].fObjSize[2], gnd.vObjects[gnd.iRoomIdx].fObjSize[2]);
	float ctrDistThres = sqrt(2.0*(MAXCOST+0.125))-0.5;
	float volCostThres = 1.0/(1+MAXCOST);

	vector<float> costMat(HypObjNum*GndObjNum);
	float gndstart;
	float locCtrDist, locVolInte, locVolume1, locVolume2, locVolCost, locTypeCost;

	for ( int gid=0; gid<GndObjNum; gid++)
	{
		gndstart = gid*HypObjNum;
		for (int hid = 0; hid<HypObjNum; hid++)
		{
			locCtrDist = objCtrDist(hyp.vObjects[hid], gnd.vObjects[gid], ctrDistNorm);
			if ( locCtrDist>ctrDistThres)
			{
				costMat[gndstart+hid] = INF;
				continue;
			}
			locVolInte = objVolInte(hyp.vObjects[hid], gnd.vObjects[gid]);
			locVolume1 = objVolume(hyp.vObjects[hid]);
			locVolume2 = objVolume(gnd.vObjects[gid]);		
			locVolCost = locVolInte/(locVolume1+locVolume2-locVolInte);
			if ( locVolCost<volCostThres)
			{
				costMat[gndstart+hid] = INF;
				continue;
			}
			locTypeCost = objTypeCost(hyp.vObjects[hid], gnd.vObjects[gid]);
			costMat[gndstart+hid] = 0.5*(locCtrDist+0.5)*(locCtrDist+0.5)-0.125 + 1.0/(locVolCost+0.000001)-1;

			/*if (verbose>2)
			{
				printf("hypid: %d, gndid: %d, CtrDist: %f, VolInte: %f, Vol1: %f, Vol2: %f, Volc: %f, Type: %f\n", hid+1, gid+1, locCtrDist, locVolInte, locVolume1, locVolume2, locVolInte/(locVolume1+locVolume2-locVolInte), locTypeCost);
			}*/

			if (locTypeCost>0.5)
			{
				costMat[gndstart+hid] += 1;
			}
			if (hyp.vObjects[hid].iObjType==29 || gnd.vObjects[gid].iObjType==29)
			{
				costMat[gndstart+hid] *= 2;
			}
			
		}
	}

	vector<int> match_hypid;
	vector<int> match_gndid;
	bipartiteMatching(costMat, HypObjNum, GndObjNum, match_cost, match_hypid, match_gndid);

	if (verbose>1)
	{
		printf("hypID: ");
		for (int i=0; i<match_hypid.size(); i++)
		{
			printf("%d ", match_hypid[i]+1);
		}
		printf("\n");
		printf("gndID: ");
		for (int i=0; i<match_gndid.size(); i++)
		{
			printf("%d ", match_gndid[i]+1);
		}
		printf("\n");
		printf("Cost: ");
		for (int i=0; i<match_gndid.size(); i++)
		{
			printf("%f ", costMat[match_gndid[i]*HypObjNum+match_hypid[i]]);
		}
		printf("\n");
	}
	
	vector<bool> unmatchedGnd;
	unmatchedGnd.assign(GndObjNum, true);
	vector<bool> unmatchedHyp;
	unmatchedHyp.assign(HypObjNum, true);
	for ( int i=0; i<match_gndid.size(); i++)
	{
		unmatchedGnd[match_gndid[i]] = false;
	}
	for ( int i=0; i<match_hypid.size(); i++)
	{
		unmatchedHyp[match_hypid[i]] = false;
	}

	float normVolume = 1.0/ctrDistNorm[0]/ctrDistNorm[1]/ctrDistNorm[2];
	float missVolume, missCost;
	float hypUnmatchCost = MAXCOST/2/2;
	float gndUnmatchCost = MAXCOST/2;
	for (int i=0; i<GndObjNum; i++)
	{
		if ( unmatchedGnd[i] )
		{
			missVolume = gnd.vObjects[i].fObjSize[0]*gnd.vObjects[i].fObjSize[1]*gnd.vObjects[i].fObjSize[2];
			missCost = gndUnmatchCost * missVolume/normVolume * 2;
			match_cost += min(max(missCost, gndUnmatchCost), gndUnmatchCost*4);

			if (verbose>1)
			{
				printf("%d unmatched gnd adds %f cost\n", i+1, min(max(missCost, gndUnmatchCost), gndUnmatchCost*4));
			}
		}
	}
	for (int i=0; i<HypObjNum; i++)
	{
		if ( unmatchedHyp[i] )
		{
			missVolume = hyp.vObjects[i].fObjSize[0]*hyp.vObjects[i].fObjSize[1]*hyp.vObjects[i].fObjSize[2];
			missCost = hypUnmatchCost * missVolume/normVolume * 2;
			match_cost += min(max(missCost, hypUnmatchCost), hypUnmatchCost*4);

			if (verbose>1)
			{
				printf("%d unmatched hyp adds %f cost\n", i+1, min(max(missCost, hypUnmatchCost), hypUnmatchCost*4));
			}
		}
	}

	//match_cost += MAXCOST*(GndObjNum-match_gndid.size()) + 0.5*MAXCOST*(HypObjNum-match_hypid.size());
	match_cost /= GndObjNum;

	if (verbose>0)
	{
		printf("Total: %f\n", match_cost);
	}
	
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	vector<ROOMHYP> vHypRooms;
	ReadData<float>( prhs[0], vHypRooms);

	vector<ROOMHYP> vGndRooms;
	ReadData<float>( prhs[1], vGndRooms);

	double *c = (double *) mxGetData(prhs[2]);
	MAXCOST = *c;

	if ( nrhs==3)
	{
		verbose = 0;
	}
	else if (nrhs==4)
	{
		double *c = (double *) mxGetData(prhs[3]);
		verbose = (int)(*c);
	}
	else
	{
		printf("Wrong number of input.\n");
		return;
	}
	

	mwSize dims[2];
	dims[0] = vHypRooms.size();
	dims[1] = vGndRooms.size();
	mxArray *out1 = mxCreateNumericArray( 2, dims, mxSINGLE_CLASS, mxREAL);
	float* p = (float *)mxGetData(out1);

	//omp_set_dynamic(0);
	//omp_set_num_threads(NUMBTHREAD);
	//#pragma omp parallel for

	//vector<future<void>> f(16);
	//int count = 0;


	//omp_set_dynamic(0);
	//omp_set_num_threads(NUMBTHREAD);
	#pragma omp parallel for
	for (int i=0; i<vGndRooms.size(); i++)
	{
		for (int j=0; j<vHypRooms.size(); j++)
		{
			roomMatching( vHypRooms[j], vGndRooms[i], p[i*vHypRooms.size()+j]);

			//f[count++] = async([&](){
			//	roomMatching( vHypRooms[j], vGndRooms[i], p[i*vHypRooms.size()+j]);
			//}, std::launch::async);
			//
			//if (count == 16) {
			//	for(int k = 0; k != 16; ++k) {
			//		f[k].get();
			//	}
			//	count = 0;
			//}
		}		
	}
	
	nlhs = 1;
	plhs[0] = out1;
	
}