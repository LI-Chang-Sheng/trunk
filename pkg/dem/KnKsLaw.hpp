#pragma once
#ifdef YADE_POTENTIAL_PARTICLES
#include <pkg/common/ElastMat.hpp>
#include <pkg/common/Dispatching.hpp>
#include <pkg/common/NormShearPhys.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <set>
#include <boost/tuple/tuple.hpp>



class KnKsPhys: public FrictPhys {
	public:
		virtual ~KnKsPhys();
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(KnKsPhys,FrictPhys,"EXPERIMENTAL. IPhys originally for potential particles",
			((vector<double>,lambdaIPOPT,0.0,,"not used, Lagrange multiplier for equality constraints"))
			((vector<int>,cstatCPLEX,,,"not used"))
			((vector<int>,rstatCPLEX,,,"not used"))
			((Real,frictionAngle,0.0,,"fric angle"))
			((Real,tanFrictionAngle,0.0,,"tangent of fric angle"))
			((Vector3r,contactDetectionPt,Vector3r(0,0,0),,"contact detection result"))
			((Real, viscousDamping, 0.8, ,"viscousDamping"))
			((Real, unitWidth2D, 1.0, ,"viscousDamping"))
			((Real, maxClosure, 0.0002, ,"not used, vmi"))
			((Real, u_peak, 0.05, ,"peak shear displacement, not fully in use"))
			((Real, u_elastic, 0.0, ,"elastic shear displacement, not fully in use"))
			((double, brittleLength, 5.0, ,"shear length where strength degrades, not fully in use"))
			((double, kn_i, 5.0, ,"currently assumes kn_i and Knormal are adopting the same value in Ip2 initialisation"))
			((double, ks_i, 5.0, ,"currently assumes kn_i and Knormal are adopting the same value in Ip2 initialisation"))
			((Vector3r, normalViscous, Vector3r(0,0,0), ,"viscousDamping"))
			((Vector3r, shearViscous, Vector3r(0,0,0), ,"viscousDamping"))
			((double, hwater, 0.0, ,"not used, height of pore water"))
			((bool, rockJointContact, false, ," rock joint"))
			((bool, intactRock, false, ,"brittle rock"))
			((int, jointType, 0, ,"jointType"))
			((Real,Kshear_area,0.0,,"kshear area"))
			((Real,Knormal_area, 0.0,,"knormal area"))
			((Vector3r, shearDir, Vector3r(0,0,0), ,"shear direction"))
			((vector<Vector3r>, shearForces, , ,"shear force"))
			((Vector3r, shear, Vector3r(0,0,0),, "shear displacement"))
			((Vector3r, prevNormal, Vector3r(0.0,0.0,0.0),, "previous Normal"))
			((Vector3r, normal, Vector3r(0.0,0.0,0.0),, " normalVector"))
			((vector<Vector3r>, pointsArea, ,, "not used, intermediate contact points"))
			((vector<Vector3r>, pointsShear, ,, "not used, points to calculate shear"))
			((vector<double>, areaShear, ,, "area to attribute shear"))
			((vector<double>, overlapDistances, ,, "not used, overlap distance"))
			((Real, finalSize, 0.0,, "not used, finalgridsize"))
			((int, finalGridNo, 0,, "not used, final number of grids"))
			((vector<double>, dualityGap, , ,"not used, duality gap for SOCP"))
			((bool, warmstart, false, ,"not used, warmstart for SOCP"))
			((int, generation, 0,, "not used, number of subdivisions"))
			((int, triNoMain, 24,, "not used, number of subdivisions"))
			((int, triNoSub, 6,, "not used, number of subdivisions"))
			((Vector3r, initial1, Vector3r(0.0,0.0,0.0),, "midpoint"))
			((Vector3r, ptOnP1, Vector3r(0.0,0.0,0.0),, "pt on particle"))
			((Vector3r, ptOnP2, Vector3r(0.0,0.0,0.0),, " pt on particle 2"))
			((vector<bool>, redundantA, , ,"not used, activePlanes for interaction.id1"))
			((vector<bool>, redundantB, , ,"not used, activePlanes for interaction.id1"))
			((vector<bool>, activePlanes1, , ,"not used, activePlanes for interaction.id1"))
			((vector<bool>, activePlanes2, , ,"not used, activePlanes for interaction.id2"))
			((vector<Real>, activeA1, , ,"not used, activePlanes for interaction.id2"))
			((vector<Real>, activeB1, , ,"not used, activePlanes for interaction.id2"))
			((vector<Real>, activeC1, , ,"not used, activePlanes for interaction.id2"))
			((vector<Real>, activeD1, , ,"not used, activePlanes for interaction.id2"))
			((vector<Real>, activeA2, , ,"not used, activePlanes for interaction.id2"))
			((vector<Real>, activeB2, , ,"not used, activePlanes for interaction.id2"))
			((vector<Real>, activeC2, , ,"not used, activePlanes for interaction.id2"))
			((vector<Real>, activeD2, , ,"not used, activePlanes for interaction.id2"))
			((int, noActive1, 0, ,"not used, activePlanes for interaction.id1"))
			((int, noActive2,0 , ,"not used, activePlanes for interaction.id2"))
			((int, smallerID,1 , ,"id of particle with smaller plane"))
			((Real, cumulative_us, 0.0,, "cumulative translation"))
			((Real, cumulativeRotation, 0.0,, "cumulative rotation"))
			((Real, mobilizedShear, , ,"mobilizedShear"))
			((Real, contactArea, 0.0, ,"contactArea"))
			((Real, radCurvFace, , ,"not used, face"))
			((double, prevJointLength, 0.0, ,"previous joint length"))
			((Real, radCurvCorner, , ,"not used, corners"))
			((Real, prevSigma,0.0 , ,"previous normal stress"))
			((vector<Real>, prevSigmaList,0.0 , ,"previous normal stress"))
			((bool, calJointLength, false,, "calculate joint length"))
			((bool, useOverlapVol, false,, "calculate overlap volume"))
			((bool, calContactArea, false,, "calculate contact area"))
			((double, jointLength, 1.0,, "approximatedJointLength"))
			((double, shearIncrementForCD, 0.0,, "toSeeWhether it is necessary to update contactArea"))
			((Real, overlappingVol,0.0 , ,"not used, overlapping vol"))
			((Real, overlappingVolMulti,0.0 , ,"not used, overlapping vol"))
			((double, gap_normalized, 0.0,, "not used, distance between particles normalized by particle size. Estimated using Taubin Distance "))
			((double, gap, 0.0,, "not used, distance between particles normalized by particle size. Estimated using Taubin Distance "))
			((bool, findCurv, false,, "not used, to get radius of curvature"))
			((bool, useFaceProperties, false,,"boolean to get face properites"))
			((Real, cohesion, 0.0, ,"cohesion"))
			((Real, tension, 0.0, ,"tension"))
			((bool, cohesionBroken, true, ,"cohesion already broken"))
			((bool, tensionBroken, true, ,"tension already broken"))
			((bool, twoDimension, false, ,"tension already broken"))
			((Real, phi_b, 0.0, ,"basic friction angle (degrees)"))
			((Real, phi_r, 0.0, ,"residual friction angle (degrees)"))
			((Real, asperity, 3.0, ,"not used, asperity height"))
			((Real, JRC, 0.0, ,"not used, Joint Roughness Coefficient"))
			((Real, JRCmobilized, 0.0, ,"not used, Joint Roughness Coefficient"))
			((Real, JCS, 0.0, ,"not used, Joint Roughness Coefficient"))
			((Real, sigmaC,0.0 , ,"not used, Joint Roughness Coefficient"))
			((Real, u_dilate,0.0 , ,"not used, dilation distance"))
			((Real, dilation_angle,0.0 , ,"not used, dilation distance"))
			/* pore water pressure */
			((Real, lambda0,0.0 , ,"not used, initial pore water pressure to stress ratio"))
			((Real, lambda_present,0.0 , ,"not used, Voight&Faust (1992) Sitar et al. (2005)"))
			((Real, u_cumulative, 0.0,, "not used, cumulative translation"))
			((Vector3r, prevShearDir, Vector3r(0.0,0.0,0.0),, "previous shear direction"))
			((Vector3r, initialShearDir, Vector3r(0.0,0.0,0.0),, "initial shear direction"))
			((double, delta_porePressure, 0.0,, "not used, change in pore water pressure"))
			((double, porePressure, 0.0,, "not used, pore water pressure"))
			((double, bandThickness, 0.1,, "not used, clay layer thickness"))
			((double, heatCapacities, 0.0,, "not used, clay layer thickness"))
			((double, effective_phi, 0.0,, "friction angle in clay after displacement"))
			((double, prevOverlap, 0.0,, "previous overlap"))
			((Real, h, 0.0,,"not used, cd")),

			//((Real, cumulativeRotation, 0.0,, "cumulative rotation"))
			//((Quaternionr, initialOrientation1, Quaternionr(1.0,0.0,0.0,0.0),, "orientation1"))
			//((Quaternionr, initialOrientation2, Quaternionr(1.0,0.0,0.0,0.0),, "orientation2")),
			createIndex();

		);

		REGISTER_CLASS_INDEX(KnKsPhys,FrictPhys);
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(KnKsPhys);




class Ip2_FrictMat_FrictMat_KnKsPhys: public IPhysFunctor {
	public:
		virtual void go(const shared_ptr<Material>& pp1, const shared_ptr<Material>& pp2, const shared_ptr<Interaction>& interaction);
		YADE_CLASS_BASE_DOC_ATTRS(Ip2_FrictMat_FrictMat_KnKsPhys,IPhysFunctor,"EXPERIMENTAL. Ip2 functor for :yref:`KnKsPhys`",
			((Real,Knormal,0.0,,"allows user to input values directly from python scripts"))
			((Real,Kshear,0.0,,"allows user to input values directly from python scripts"))
			((Real, unitWidth2D, 1.0, ,"width in 2D"))
			((double, brittleLength, 3.0, ,"shear length for degradation"))
			((double, kn_i, 0.002, ,"assumes kn_i and Knormal are adopting the same value in Ip2 initialisation"))
			((double, ks_i, 0.002, ,"assumes ks_i and Knormal are adopting the same value in Ig2 initialisation"))
			((double, u_peak, -1.0, ,"not used"))
			((double, maxClosure, 0.002, ,"not used"))
			((Real, viscousDamping, 0.8, ,"viscousDamping"))
			((Real, cohesion, 0.0, ,"cohesion"))
			((Real, tension, 0.0, ,"tension"))
			((bool, cohesionBroken, true, ,"cohesion"))
			((bool, tensionBroken, true, ,"tension"))
			((Real, phi_b, 0.0, ,"friction angle"))
			((bool, useOverlapVol, false,, "not used, calculate overlap volume"))
			((bool, useFaceProperties, false,,"boolean to get face properites"))
			((bool, calJointLength, false,, "calculate joint length"))
			((bool, twoDimension, false, ,"boolean for 2D"))
		);
		FUNCTOR2D(FrictMat,FrictMat);
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ip2_FrictMat_FrictMat_KnKsPhys);



class Law2_SCG_KnKsPhys_KnKsLaw: public LawFunctor {
	public:

		static Real Real0;
		//OpenMPAccumulator<Real,&Law2_SCG_KnKsPhys_KnKsLaw::Real0> plasticDissipation;
		virtual bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I);
		FUNCTOR2D(ScGeom,KnKsPhys);
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_SCG_KnKsPhys_KnKsLaw,LawFunctor,"Law for linear compression, without cohesion and Mohr-Coulomb plasticity surface.\n\n.. note::\n This law uses :yref:`ScGeom`; there is also functionally equivalent :yref:`Law2_Dem3DofGeom_FrictPhys_Basic`, which uses :yref:`Dem3DofGeom` (sphere-box interactions are not implemented for the latest).",
			((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))
			((bool,preventGranularRatcheting,false,,"bool to avoid granular ratcheting"))
			((bool,traceEnergy,false,,"Define the total energy dissipated in plastic slips at all contacts."))
			((bool,Talesnick,false,,"use contact law developed for validation against model test"))
			((double,waterLevel,0.0,,"not used"))
			((bool, allowBreakage, false, ,"allow cohesion to break"))
			((double,initialOverlapDistance,0.0,,"initial over lap distance"))
			,,

		);

		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_SCG_KnKsPhys_KnKsLaw);

#endif // YADE_POTENTIAL_PARTICLES
