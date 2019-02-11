/*CWBoon 2015 */
#ifdef YADE_POTENTIAL_BLOCKS
#pragma once

#include<vector>
#include<core/Shape.hpp>
#include<Eigen/Core>
#include <Eigen/LU> 
#include <Eigen/QR>
#include<lib/base/openmp-accu.hpp>
namespace yade{
class PotentialBlock : public Shape
{

	public:
		
		
		struct Planes{
			vector<int> vertexID;
		};
		struct Vertices{
			vector<int> edgeID;
			vector<int> planeID;
		};
		struct Edges{
			vector<int> vertexID;
		};
		
		void addPlaneStruct();
		void addVertexStruct();
		void addEdgeStruct();
		
		vector<Planes> planeStruct;	
		vector<Vertices> vertexStruct;	
		vector<Edges> edgeStruct;
	
		Eigen::MatrixXd Amatrix;
		Eigen::MatrixXd Dmatrix;
		virtual ~PotentialBlock ();
		void postLoad(PotentialBlock&);
	
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(PotentialBlock,Shape,"Geometry of PotentialBlock.",
		((bool, isLining, false,, "false"))
		((double, liningStiffness, pow(10.0,8),, "lining stiffness"))
		((double, liningFriction, 20.0,, "lining friction"))
		((double, liningLength, 0.0,, "lining spacing between nodes"))
		((double, liningTensionGap, 0.0,, "numerical gap between lining and block to allowing tension to be calculated"))
		((Vector3r, liningNormalPressure, Vector3r(0,0,0),, "normal pressure acting on lining"))
		((Vector3r, liningTotalPressure, Vector3r(0,0,0),, "total pressure acting on lining"))
		((bool, isBoundary, false,, "boolean whether a block is part of a boundary block"))
		((bool, isEastBoundary, false,, "boolean used for a case history"))
		((bool, isBolt, false,, "boolean whether a block is part of a bolt"))
		((bool, fixedNormal, false,, "boolean whether to fix the normal contact at a boundary"))
		((Vector3r, boundaryNormal, Vector3r::Zero(),,"normal direction of boundary"))
		((bool, AabbMinMax, false,, "aabb"))
		((Vector3r, minAabb, Vector3r::Zero(),,"min from box centre"))
		((Vector3r, maxAabb, Vector3r::Zero(),,"max from box centre"))
		((Vector3r, minAabbRotated, Vector3r::Zero(),,"updated aabb from box centre"))
		((Vector3r, maxAabbRotated, Vector3r::Zero(),,"updated aabb from box centre"))
		((Vector3r, halfSize, Vector3r::Zero(),,"max from box centre"))
		((vector<Vector3r> , node, ,, "not used "))
		((Real, gridVol, ,, "not used "))
		((Quaternionr , oriAabb, Quaternionr::Identity(),, "original aabb "))
		((Real , r, 0.1,, "r in Potential Particles "))
		((Real , R, 1.0,, "R in Potential Particles"))
		((Real , k, 0.01,, "k in Potential Particles"))
		((Real , volume, 0.1,, "volume"))
		((int, id, -1,, " for graphics"))
		((bool, erase, false,, " for excavation"))
		((vector<bool>, intactRock, false,, " property for plane"))
		((vector<bool>, isBoundaryPlane, ,, "property for plane"))
		((vector<double>, hwater, ,, "property for plane (not used), height of pore water"))
		((vector<double>, JRC, ,, "property for plane (not used),rock joint"))
		((vector<double>, JCS, ,, "property for plane (not used),rock joint"))
		((vector<double>, asperity, ,, "property for plane (not used),rock joint"))
		((vector<double>, sigmaC, ,,"property for plane (not used),rock joint"))
		((vector<double>, phi_b, ,, "property for plane ,rock joint"))
		((vector<double>, phi_r, ,, "property for plane ,rock joint"))
		((vector<double>, cohesion, ,, "property for plane ,rock joint"))
		((vector<double>, tension, ,, "property for plane ,rock joint"))
		((vector<double>, lambda0, ,, "property for plane (not used),heat"))
		((vector<double>, kn, ,, "property for plane ,rock joint"))
		((vector<double>, ks, ,, "property for plane ,rock joint"))
		((vector<double>, heatCapacity, ,, "property for plane ,rock joint"))
		((vector<double>, rFactor, ,, "not used, individual factor for r"))
		((vector<Vector3r>, vertices,,,"vertices"))
		//((Eigen::MatrixXd , Amatrix, ,, "a "))
		//((Eigen::MatrixXd , Dmatrix, ,, "b "))
		((double, waterVolume, ,, "volume of body submerged in water"))
		((vector<Vector3r> ,verticesCD, ,, "vertices "))
		((vector<double> , a, ,, "list of a of plane normal"))
		((vector<double> , b, ,, "list of b of plane normal"))
		((vector<double> , c, ,, "list of c of plane normal"))
		((vector<double> , d, ,, "list of d of plane normal"))
		((vector<int> , jointType, ,, "jointType"))
		,
		createIndex(); /*ctor*/
		#if 0
		for (int i=0; i<a.size(); i++){
			Amatrix(i,0) = a[i]; Amatrix(i,1)=b[i]; Amatrix(i,2)=c[i];
			Dmatrix(i,0) = d[i] + r;
	 	}
		#endif
	
	);
	//#endif

	REGISTER_CLASS_INDEX(PotentialBlock,Shape);
	
};
}
using namespace yade; 

REGISTER_SERIALIZABLE(PotentialBlock);

#ifdef __cplusplus
extern "C" {
#endif
	 void dgesv_(const int *N, const int *nrhs, double *Hessian, const int *lda, int *ipiv, double *gradient, const int *ldb, int *info);
	void dsyev_(const char *jobz, const char *uplo, const int *N, double *A, const int *lda, double *W, double *work, int *lwork, int *info);
#ifdef __cplusplus
};
#endif
#endif // YADE_POTENTIAL_BLOCKS
