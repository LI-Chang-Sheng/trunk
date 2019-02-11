/*************************************************************************
*  Copyright (C) 2018 by Robert Caulk <rob.caulk@gmail.com>  		 *
*  Copyright (C) 2018 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*									 *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
/* This engine is under active development. Experimental only */

//#define THERMAL
#ifdef THERMAL
#pragma once 

#include<core/PartialEngine.hpp>
#include<core/State.hpp>
#include<pkg/dem/JointedCohesiveFrictionalPM.hpp>
#include<pkg/dem/ScGeom.hpp>
#include<pkg/common/Dispatching.hpp>

#ifdef FLOW_ENGINE
//#include<pkg/pfv/FlowEngine.hpp>
#include<lib/triangulation/Tesselation.h>
#include<lib/triangulation/FlowBoundingSphere.hpp>
#include "FlowEngine_FlowEngineT.hpp"
#include<pkg/dem/TesselationWrapper.hpp>
#include<lib/triangulation/Network.hpp>
#endif

class ThermalEngine : public PartialEngine
{
	public:
		typedef TemplateFlowEngine_FlowEngineT<FlowCellInfo_FlowEngineT,FlowVertexInfo_FlowEngineT> FlowEngineT;
		typedef FlowEngineT::Tesselation					Tesselation;
		typedef FlowEngineT::RTriangulation					RTriangulation;
		typedef FlowEngineT::FiniteCellsIterator				FiniteCellsIterator;
		typedef FlowEngineT::CellHandle						CellHandle;
		typedef FlowEngineT::VertexHandle	VertexHandle;
		typedef std::vector<CellHandle>		VectorCell;
		typedef typename VectorCell::iterator		VCellIterator;

	public:
		Scene* scene;
		bool energySet; //initialize the internal energy of the particles
		FlowEngineT* flow;
        	bool timeStepEstimated;
       		Real thermalDT;
        	bool conductionIter;
        	int elapsedIters;
        	int conductionIterPeriod;
        	Real elapsedTime;
        	bool first;
        	bool runConduction;
        	Real maxTimeStep;
		
		virtual ~ThermalEngine();
		virtual void action();
		void setInitialValues();
		void resetBoundaryFluxSums();
		void setConductionBoundary();
		void thermalExpansion();
		void initializeInternalEnergy();
        	void computeNewPoreTemperatures();
        	void computeNewParticleTemperatures();
		void computeSolidFluidFluxes();
		void computeVertexSphericalArea();
		void computeFlux(CellHandle& cell, const shared_ptr<Body>& b, const double surfaceArea);
		void computeSolidSolidFluxes();
        	void timeStepEstimate();
        	const Real computeCellPressureChangeFromDeltaTemp(CellHandle& cell);
        	const Real computeCellPressureChangeFromDeltaVolume(CellHandle& cell);
        	Real getThermalDT() {return thermalDT;}
        	int getConductionIterPeriod() {return conductionIterPeriod;}
        	Real getMaxTimeStep() {return maxTimeStep;}
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ThermalEngine,PartialEngine,"preliminary",
		/*attributes*/
		((bool,advection,true,,"Activates advection"))
		((bool,debug,false,,"debugging flags"))
		((bool,conduction,true,,"Activates conduction"))
		((bool,thermoMech,true,,"Activates thermoMech"))
        	((bool,fluidThermoMech,true,,"Activates thermoMech"))
        	((bool,solidThermoMech,true,,"Activates thermoMech"))
		((vector<bool>, bndCondIsTemperature, vector<bool>(6,false),,"defines the type of boundary condition for each side. True if temperature is imposed, False for no heat-flux. Indices can be retrieved with :yref:`FlowEngine::xmin` and friends."))
		((vector<double>, thermalBndCondValue, vector<double>(6,0),,"Imposed value of a boundary condition."))
		((vector<double>, thermalBndFlux, vector<double>(6,0),,"Flux through thermal boundary."))
		((bool,boundarySet,false,,"set false after changing boundary conditions"))
		((bool,useKernMethod,false,,"flag to use Kern method for thermal conductivity area calc"))
        	((Real,fluidBeta,0.0002,,"volumetric temperature coefficient m^3/m^3C, default water, <= 0 deactivates"))
        	((Real,particleT0,0,,"Initial temperature of particles"))
		((Real,particleK,3.0,,"Particle thermal conductivity (W/(mK)"))
		((Real,particleCp,750.,,"Particle thermal heat capacity (J/(kgK)"))
		((Real,particleAlpha,11.6e-6,,"Particle volumetric thermal expansion coeffcient"))
        	((double, fluidK, 0.650,,"Thermal conductivity of the fluid."))
        	((double, tsSafetyFactor, 0.8,,"Allow user to control the timstep estimate with a safety factor. Default 0.8"))
		,
		/* extra initializers */
		,
		/* ctor */
		energySet=false;timeStepEstimated=false;thermalDT=0;elapsedTime=0;first=true;runConduction=false;maxTimeStep=10000;
		,
		/* py */
        	.def("getThermalDT",&ThermalEngine::getThermalDT,"let user check estimated thermalDT .")
        	.def("getConductionIterPeriod",&ThermalEngine::getConductionIterPeriod,"let user check estimated conductionIterPeriod .")
        	.def("getMaxTimeStep",&ThermalEngine::getMaxTimeStep,"let user check estimated maxTimeStep.")
	)
	DECLARE_LOGGER;

	
};
REGISTER_SERIALIZABLE(ThermalEngine);


#endif//THERMAL
	
