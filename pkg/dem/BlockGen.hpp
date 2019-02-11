/*CWBoon 2016 */
/* Please cite: */
/* CW Boon, GT Houlsby, S Utili (2015).  A new rock slicing method based on linear programming.  Computers and Geotechnics 65, 12-29. */
/* The numerical library is changed from CPLEX to CLP because subscription to the academic initiative is required to use CPLEX for free */
#ifdef YADE_POTENTIAL_BLOCKS
#pragma once
#include<core/FileGenerator.hpp>

#include<Python.h>


#include<lib/base/Math.hpp>
#include<lib/serialization/Serializable.hpp>
#include<lib/multimethods/Indexable.hpp>

#include <ClpSimplexDual.hpp>
#include <ClpSimplex.hpp>
#include <CoinHelperFunctions.hpp>
#include <CoinTime.hpp>
#include <CoinBuild.hpp>
#include <CoinModel.hpp>
#include <iomanip>
#include <cassert>

class GlobalStiffnessTimeStepper;

class BlockGen : public FileGenerator
{
	private	:
				
		
		
		void createActors(shared_ptr<Scene>& scene);
		void positionRootBody(shared_ptr<Scene>& scene);
		
		shared_ptr<GlobalStiffnessTimeStepper> globalStiffnessTimeStepper;
		

	protected:
//		std::ofstream output2;        // it was always creating files "BlkGen" "BlockGenFindExtreme.txt", but they are not used in the code, so I commented this out, Janek
//		std::string myfile;
		std::string Key;
//		static std::ofstream output;
		
	public : 
		~BlockGen ();
		bool generate(string&);
		template<class T>
		double gen_normal_3(T &generator){ return generator();}
		

		struct Discontinuity{
			Vector3r centre;
			Discontinuity(Vector3r pos){centre = pos; persistence = false; phi_b = 30.0; phi_r= 30.0; JRC=15; JCS = pow(10,6); asperity = 5; sigmaC = JCS;isBoundary= false;sliceBoundaries=false; lambda0=0.0; heatCapacity=0.0;hwater=-1.0;intactRock=false; throughGoing = false; constructionJoints=false;jointType=0;}
			Real a;
			Real b;
			Real c;
			Real d;
			Real a_p;
			Real b_p;
			Real c_p;
			Real d_p;
			bool persistence;
			bool isBoundary;
			bool sliceBoundaries;
			bool constructionJoints;
			vector<double> persistence_a;
			vector<double> persistence_b;
			vector<double> persistence_c;
			vector<double> persistence_d;
			/* Joint properties */
			double phi_b;
			double phi_r;
			double JRC;
			double JCS;
			double asperity;
			double sigmaC;
			double cohesion;
			double tension;
			double lambda0;
			double heatCapacity;
			double hwater;
			bool intactRock;
			bool throughGoing;
			int jointType;
		};
		struct Planes{
			vector<int> vertexID;
		};
		struct Block{
			Vector3r tempCentre;
			Vector3r centre;
			Block(Vector3r pos, Real kPP, Real rPP, Real RPP){centre = pos; k=kPP; r=rPP; R=RPP; tooSmall=false;isBoundary=false;tempCentre=pos;}
			vector<double> a;
			vector<double> b;
			vector<double> c;
			vector<double> d;
			vector<bool> redundant;	
			vector<bool> isBoundaryPlane;
			bool isBoundary;	
			vector<struct Block> subMembers;
			vector<Vector3r> falseVertex;
			vector<Vector3r> node; Real gridVol;
			Real r;
			Real R;
			Real k;
			bool tooSmall;
			/* Joint properties */
			vector<double> phi_b;
			vector<double> phi_r;
			vector<double> JRC;
			vector<double> JCS;
			vector<double> asperity;
			vector<double> sigmaC;
			vector<double> cohesion;
			vector<double> tension;
			vector<double> lambda0;
			vector<double> heatCapacity;
			vector<double> hwater;
			vector<bool> intactRock;
			vector<int> jointType;
			vector<struct Planes> planes;
		};
		

	

		

	double getSignedArea(const Vector3r pt1,const Vector3r pt2, const Vector3r pt3);
	double getDet(const Eigen::MatrixXd A);
	double getCentroidTetrahedron(const Eigen::MatrixXd A);


		bool createBlock(shared_ptr<Body>& body, struct BlockGen::Block block, int no);
		bool contactDetection(struct BlockGen::Discontinuity joint, struct BlockGen::Block block,  Vector3r& touchingPt);
		
		bool contactDetectionLPCLP(struct BlockGen::Discontinuity joint, struct BlockGen::Block block,  Vector3r& touchingPt);
		bool contactDetectionLPCLPglobal(struct BlockGen::Discontinuity joint, struct BlockGen::Block block,  Vector3r& touchingPt);
	
		bool checkRedundancyLPCLP(struct BlockGen::Discontinuity joint, struct BlockGen::Block block,  Vector3r& touchingPt);
		bool startingPointFeasibility(struct BlockGen::Block block,  Vector3r& initialPoint);
	
		double inscribedSphereCLP(struct BlockGen::Block block,  Vector3r& initialPoint, bool twoDimension);
		
		
		bool contactBoundaryLPCLP(struct BlockGen::Discontinuity joint, struct BlockGen::Block block,  Vector3r& touchingPt);
		bool contactBoundaryLPCLPslack(struct BlockGen::Discontinuity joint, struct BlockGen::Block block,  Vector3r& touchingPt);
		
		

	
		
		
		
		void calculateInertia(struct Block block, Real& Ixx, Real& Iyy, Real& Izz,Real& Ixy, Real& Ixz, Real& Iyz);
		
		Vector3r calCentroid(struct Block block, double & blockVol);
		
		bool checkCentroid(struct Block block, Vector3r presentTrial);
		Real evaluateFNoSphere(struct Block block, Vector3r presentTrial);
		
	

		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(
		BlockGen,FileGenerator,"Prepare a scene for triaxial tests. See full documentation at http://yade-dem.org/wiki/TriaxialTest."
		,
   		/* public */
		((Real,dampingMomentum,0.2,,"Coefficient of global damping"))
		((Real, maxClosure, 0.0002, ,"not used"))
		((Real, peakDisplacement, 0.02, ,"not used"))
		((double, brittleLength, 2.0, ,"not used"))
		((Real,damp3DEC,0.8,,"not used"))
		((Real,unitWidth2D,1.0,,"out of plane distance"))
		((Real,density,2600,,"density of blocks"))
		((Real,Kn,pow(10,-5),,"contact normal stiffness"))
		((Real,Ks,pow(10,-5),,"contact shear stiffness"))
		((Real,frictionDeg,18.0,,"friction angle [°]"))
		((Vector3r,globalOrigin,Vector3r(0.0,0.0,0.0),,"friction angle [°]"))
		((Real,inertiaFactor,1.0,,"scaling of inertia"))
		((Real,rForPP,0.1,,"r in Potential Particles"))
		((Real,kForPP,0.01,,"k in Potential Particles"))
		((Real,RForPP,3.5,,"R in Potential Particles"))
		((int,numberOfGrids,1,,"not used"))
		((bool,probabilisticOrientation,false,,"boolean to generate rock joints randomly"))
		((bool,Talesnick,false,,"boolean for contact law, used for validating code previously against model test"))
		((bool,neverErase,false,,"boolean whether to erase non interacting contacts"))
		((bool,calJointLength,false,,"boolean whether to calculate jointLength"))
		((bool,twoDimension,false,,"boolean whether it is 2D"))
                ((Real,shrinkFactor,1.0,,"ratio to shrink r"))
		((Real,viscousDamping,0.8,,"viscous damping"))
		((bool, intactRockDegradation, false, ,"boolean to active degradation of parameters for contact"))
		((Real,initialOverlap,0.0,,"initial overlap between blocks"))
		((Vector3r,gravity,Vector3r(0.0,-9.81,0.0),,"gravity"))
		((Real,defaultDt,-1,,"Max time-step. Used as initial value if defined. Latter adjusted by the time stepper."))
		((int,timeStepUpdateInterval,50,,"interval for :yref:`GlobalStiffnessTimeStepper`"))
		((bool,traceEnergy,true,,"boolean to traceEnergy - not used"))
		((bool,exactRotation,true,,"boolean to use more accurate quarternions"))
		((Real,minSize,50.0,,"minimum size for all blocks"))
		((Real,minSize2,50.0,,"minimum size for blocks with joint Type=2, minSize2 is smaller than minSize"))
		((Real,maxRatio,3.0,,"minimum ratio for all blocks"))
		((Real,maxRatio2,1000.0,,"not used"))
		((Real,boundarySizeXmax,1.0,,"max X of domain"))
		((Real,boundarySizeYmax,1.0,,"max Y of domain"))
		((Real,boundarySizeZmax,1.0,,"max Z of domain"))
		((Real,boundarySizeXmin,1.0,,"min X of domain"))
		((Real,boundarySizeYmin,1.0,,"min Y of domain"))
		((Real,boundarySizeZmin,1.0,,"min Z of domain"))
		((Vector3r,directionA,Vector3r(1,0,0),,"local x-direction to check minSize"))
		((Vector3r,directionB,Vector3r(0,1,0),,"local y-direction to check minSize"))
		((Vector3r,directionC,Vector3r(0,0,1),,"local z-direction to check minSize"))
		((Real,calAreaStep,10.0,,"length Z of domain"))
		((Real,extremeDist,0.5,,"boundary to base calculation of octree algorithms"))  //FIXME remove when not needed anymore
		((Real,subdivisionRatio,0.1,,"smallest size/boundary of octree algorithms"))
		((vector<double>,joint_a,,,"list of a of plane normal"))
		((vector<double>,joint_b,,,"list of b of plane normal"))
		((vector<double>,joint_c,,,"list of c of plane normal"))
		((vector<double>,joint_d,,,"list of d in plane equation"))
		/* check for different joint types */
		((bool, persistentPlanes,false,,"boolean to check persistence"))
		((bool, jointProbabilistic,false,,"boolean to check for filename jointProbabilistic"))	
		((bool, jointProbabilisticRockBridge,false,,"boolean to check for filename jointProbabilisticRockBridge"))	
		((bool, opening,false,,"boolean to check for filename opening"))	
		((bool, boundaries,false,,"boolean to check for filename boundaries"))
		((bool, slopeFace,false,,"boolean to check for filename "))
		((bool, sliceBoundaries,false,,"boolean to check for filename "))
		((bool, useGlobalStiffnessTimeStepper,false,,"boolean "))
		/* which contact law to use */
		((bool, useBartonBandis,false,,"not used"))
		((bool, useFaceProperties, false, , "boolean to use face properties"))
		((bool, useOverlapVol,false,,"not used"))
		((std::string,filenamePersistentPlanes,"/home/booncw/yade/BranchA/scripts/Tunnel/jointPersistent.csv",,"filename to look for joint properties"))
		((std::string,filenameProbabilisticRockBridge,"/home/booncw/yade/BranchA/scripts/boon/Tunnel/jointProbabilistic.csv",,"filename to look for joint with probabilistic models"))
		((std::string,filenameProbabilistic,"/home/booncw/yade/BranchA/scripts/boon/Tunnel/jointProbabilistic.csv",,"filename to look for joint with probabilistic models"))
		((std::string,filenameBoundaries,"/home/booncw/yade/BranchA/scripts/boon/Tunnel/boundaries.csv",,"filename to look for joint with probabilistic models"))
		((std::string,filenameOpening,"/home/booncw/yade/BranchA/scripts/boon/Tunnel/opening.csv",,"filename to look for joint outline of joints"))
		((std::string,filenameSlopeFace,"/home/booncw/yade/BranchA/scripts/boon/Tunnel/opening.csv",,"filename to look for joint outline of joints"))
		((std::string,filenameSliceBoundaries,"/home/booncw/yade/BranchA/scripts/boon/Tunnel/sliceBoundaries.csv",,"filename to look for joint outline of joints"))
		//((double*, array_a,,,"a"))	
		,
		/* init */
		,
		/* constructor for private */
		Key			="";
//		myfile = "./BlkGen"+Key;
//		output2.open(myfile.c_str(), fstream::app);

		,
		//.def("setContactProperties",&TriaxialCompressionEngine::setContactProperties,"Assign a new friction angle (degrees) to dynamic bodies and relative interactions")
		 );
	DECLARE_LOGGER;
};

#ifdef __cplusplus
extern "C" {
#endif
	 void dgesv_(const int *N, const int *nrhs, double *Hessian, const int *lda, int *ipiv, double *gradient, const int *ldb, int *info);
	void dsyev_(const char *jobz, const char *uplo, const int *N, double *A, const int *lda, double *W, double *work, int *lwork, int *info);
#ifdef __cplusplus
};
#endif

REGISTER_SERIALIZABLE(BlockGen);
#endif // YADE_POTENTIAL_BLOCKS
