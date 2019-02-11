// 2016 © William Chèvremont <william.chevremont@univ-grenoble-alpes.fr>

#pragma once

#include<pkg/dem/CohesiveFrictionalContactLaw.hpp>
#include<pkg/common/Dispatching.hpp>
#include<pkg/common/Sphere.hpp>
#include<pkg/common/PeriodicEngines.hpp>
#include<pkg/common/NormShearPhys.hpp>
#include<pkg/dem/DemXDofGeom.hpp>
#include<pkg/dem/ScGeom.hpp>
#include<pkg/dem/FrictPhys.hpp>
#include<pkg/dem/ElasticContactLaw.hpp>
#include<pkg/dem/Lubrication.hpp>

namespace py=boost::python;

class ElectrostaticMat: public FrictMat {
	public:
                YADE_CLASS_BASE_DOC_ATTRS_CTOR(ElectrostaticMat,FrictMat,"Electrostatic material, used in :yref:`Ip2_ElectrostaticMat_ElectrostaticMat_ElectrostaticPhys` and :yref:`Law2_ScGeom_ElectrostaticPhys`.",
					((Real, A, 1.e-19,,"Hamaker constant for this material[J]"))
			,
                        createIndex();
		);
                REGISTER_CLASS_INDEX(ElectrostaticMat,FrictMat);
};
REGISTER_SERIALIZABLE(ElectrostaticMat);

 
class ElectrostaticPhys: public LubricationPhys {
        public:
				explicit ElectrostaticPhys(LubricationPhys const&); // Inheritance constructor
                YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(ElectrostaticPhys,LubricationPhys,"IPhys class containing parameters of DLVO interaction Inherits from LubricationPhys. Used by :yref:`Law2_ScGeom_ElectrostaticPhys`.",
                        ((Real,DebyeLength,1e-6,Attr::readonly,"Debye Length $\\kappa^{-1}$[m]"))
                        ((Real,Z,1e-12,Attr::readonly,"Double layer interaction constant $Z$ [N]"))
                        ((Real,A,1e-19,Attr::readonly,"Hamaker constant $A = \\sqrt{A_1A_2}$ [J]"))
						((Real,vdw_cutoff, 1.e-4,,"Van der Waals cutoff ratio []"))
						((Vector3r,normalDLForce,Vector3r::Zero(),Attr::readonly,"Normal force due to Double Layer electrostatic repulsion"))
						((Vector3r,normalVdWForce,Vector3r::Zero(),Attr::readonly,"Normal force due to Van-der-Waals attraction"))
			, // ctors
                        createIndex();,
		);
		DECLARE_LOGGER;
                REGISTER_CLASS_INDEX(ElectrostaticPhys,LubricationPhys);
};
REGISTER_SERIALIZABLE(ElectrostaticPhys);

class Ip2_ElectrostaticMat_ElectrostaticMat_ElectrostaticPhys: public Ip2_FrictMat_FrictMat_LubricationPhys{
	public:
		virtual void go(const shared_ptr<Material>& material1, const shared_ptr<Material>& material2, const shared_ptr<Interaction>& interaction);
		static Real getInteractionConstant(Real const& epsr=78, Real const& T=293, Real const& z=1, Real const& phi0=0.050);
		//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getInteractionConstant_overloads, getInteractionConstant, 0, 4)
                FUNCTOR2D(ElectrostaticMat,ElectrostaticMat);
		DECLARE_LOGGER;
                YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Ip2_ElectrostaticMat_ElectrostaticMat_ElectrostaticPhys, Ip2_FrictMat_FrictMat_LubricationPhys,"Ip2 creating Electrostatic_Phys from two ElectrostaticMat instances.",
                        ((Real,DebyeLength,1.e-6,,"Debye length [m]."))
                        ((Real,Z,1.e-12,,"Interaction constant [N]."))
						((Real,vdw_cutoff, 1.e-4,,"Van der Waals cutoff ratio []"))
                    ,,
					.def("getInteractionConstant", &getInteractionConstant, (py::args("epsr")=78,py::args("T")=293,py::args("z")=1,py::args("phi0")=0.050),"Get the interaction constant from thermal properties").staticmethod("getInteractionConstant")
		);
};
REGISTER_SERIALIZABLE(Ip2_ElectrostaticMat_ElectrostaticMat_ElectrostaticPhys);



class Law2_ScGeom_ElectrostaticPhys: public Law2_ScGeom_ImplicitLubricationPhys{
	public:
		Real normalForce_DLVO_Adim(ElectrostaticPhys *phys, ScGeom* geom, Real const& undot, bool isNew, bool dichotomie);
		Real DLVO_NRAdimExp_integrate_u(Real const& un, Real const& eps, Real const& alpha, Real const& A, Real const& vdwcut, Real const& Z, Real const& K, Real & prevDotU, Real const& dt, Real const& prev_d, Real const& undot, int depth=1);
		
		Real DLVO_DichoAdimExp_integrate_u(Real const& un, Real const& eps, Real const& alpha, Real const& A, Real const& vdwcut, Real const& Z, Real const& K, Real & prevDotU, Real const& dt, Real const& prev_d, Real const& undot);
		Real ObjF(Real const& un, Real const& eps, Real const& alpha, Real const& prevDotU, Real const& dt, Real const& prev_d, Real const& undot, Real const& A, Real const& vdwcut, Real const& Z, Real const& K, Real const& d);
		
		static void getStressForEachBody(vector<Matrix3r>& DLStresses, vector<Matrix3r>& VdWStresses);
		static py::tuple PyGetStressForEachBody();
		static void getTotalStresses(Matrix3r& DLStresses, Matrix3r& VdWStresses);
		static py::tuple PyGetTotalStresses();
		
		bool go(shared_ptr<IGeom>& iGeom, shared_ptr<IPhys>& iPhys, Interaction* interaction);
		FUNCTOR2D(GenericSpheresContact,ElectrostaticPhys);
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_ScGeom_ElectrostaticPhys,Law2_ScGeom_ImplicitLubricationPhys,"Material law for lubricated spheres with DLVO interaction between 2 spheres. Solved either by Newton-Rafson method (only Double-Layer repulsion) or dichotomy method (whole interaction). No exact resolution method available",,,
			.def("getStressForEachBody",&Law2_ScGeom_ElectrostaticPhys::PyGetStressForEachBody,"Get stresses tensors for each bodies: normal contact stress, shear contact stress, normal lubrication stress, shear lubrication stress, DL stress, VdWStress.")
			.staticmethod("getStressForEachBody")
			.def("getTotalStresses",&Law2_ScGeom_ElectrostaticPhys::PyGetTotalStresses,"Get total stresses tensors: normal contact stress, shear contact stress, normal lubrication stress, shear lubrication stress, DL stress, VdW Stress")
			.staticmethod("getTotalStresses")
		);
		DECLARE_LOGGER;

};
REGISTER_SERIALIZABLE(Law2_ScGeom_ElectrostaticPhys);
