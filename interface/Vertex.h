#ifndef VERTEX_H 
#define VERTEX_H 

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/format.h"

class Vertex{
 	public:
		Vertex(){}
		Vertex( VertexInfoBranches& VtxInfo, int idx ){
			Fill( VtxInfo, idx );
		}
		void Fill( VertexInfoBranches& VtxInfo, int idx ){
			isValid = VtxInfo.isValid[idx];
			isFake = VtxInfo.isFake[idx]; 
			Type = VtxInfo.Type[idx];    //0 - Offline Primary Vertices, 1 - Offline Primary Vertices with beam spot constraint, 2 - Pixel Vertices
			Ndof = VtxInfo.Ndof[idx];
			NormalizedChi2 = VtxInfo.NormalizedChi2[idx];
			Pt_Sum = VtxInfo.Pt_Sum[idx];
			Pt_Sum2 = VtxInfo.Pt_Sum2[idx];
			x = VtxInfo.x[idx];
			y = VtxInfo.y[idx];
			z = VtxInfo.z[idx];
			Rho = VtxInfo.Rho[idx];
		}

		int     isValid;
		bool    isFake; //Uly 2011-04-04
		int     Type;   //0 - Offline Primary Vertices, 1 - Offline Primary Vertices with beam spot constraint, 2 - Pixel Vertices
		float   Ndof;
		float   NormalizedChi2;
		float   Pt_Sum;
		float   Pt_Sum2;
		float   x;
		float   y;
		float   z;
		float   Rho;
};

#endif
