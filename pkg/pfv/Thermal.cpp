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
#include<pkg/pfv/Thermal.hpp>
#include<core/Scene.hpp>
#include<core/Omega.hpp>
#include<pkg/common/Sphere.hpp>
#include<pkg/dem/Shop.hpp>
 
//YADE_PLUGIN((ThermalState));
//ThermalState::~ThermalState(){};

CREATE_LOGGER(ThermalEngine);
YADE_PLUGIN((ThermalEngine));

ThermalEngine::~ThermalEngine(){}

void ThermalEngine::action()
{
	scene = Omega::instance().getScene().get();
	if (first) setInitialValues();	
	elapsedTime += scene->dt;
	if (elapsedIters > conductionIterPeriod || first){
	        runConduction = true;
	        timeStepEstimated = false;
	        if (first) thermalDT = scene->dt;
	        else thermalDT = elapsedTime;
	        first = false;
	}
	
	FOREACH(const shared_ptr<Engine> e, Omega::instance().getScene()->engines) {
		if (e->getClassName() == "FlowEngine") {
			flow = dynamic_cast<FlowEngineT*>(e.get());
		}
	}
	resetBoundaryFluxSums();
	//resetStepFluxSums();
	if (!boundarySet) setConductionBoundary();
	if (!conduction) thermoMech = false; //don't allow thermoMech if conduction is not activated
	if (advection) {
		flow->solver->initializeInternalEnergy(); // internal energy of cells
		flow->solver->augmentConductivityMatrix(scene->dt);
	}
	if (conduction && runConduction) {
		// if (!energySet) initializeInternalEnergy();  // internal energy of particles // dont need if we stick to tracking temps as we are now
		computeSolidSolidFluxes(); 
		if (advection) computeSolidFluidFluxes();
	}
	
	if (conduction && runConduction) computeNewParticleTemperatures();
	if (advection) computeNewPoreTemperatures();
	if (thermoMech && runConduction) thermalExpansion();
	if (!timeStepEstimated) timeStepEstimate();
	elapsedIters += 1;
}

void ThermalEngine::setInitialValues() {
	const shared_ptr<BodyContainer>& bodies=scene->bodies;
	const long size=bodies->size();
	for (long i=0; i<size; i++){
		const shared_ptr<Body>& b=(*bodies)[i];
		if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
		auto* state =b->state.get();
		state->temp = particleT0;
		state->k = particleK;
		state->Cp = particleCp;
		state->alpha = particleAlpha;
	}
}

void ThermalEngine::timeStepEstimate() {
	const shared_ptr<BodyContainer>& bodies=scene->bodies;
	const long size=bodies->size();
//	#pragma omp parallel for
	for (long i=0; i<size; i++){
		const shared_ptr<Body>& b=(*bodies)[i];
		if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
		auto* thState = b->state.get();
		const Real bodyTimeStep = thState->mass * thState->Cp / thState->stabilityCoefficient;
        thState->stabilityCoefficient = 0; // reset the stability coefficient
        if (!maxTimeStep) maxTimeStep = bodyTimeStep;
        if (bodyTimeStep < maxTimeStep) maxTimeStep = bodyTimeStep;
	}
	timeStepEstimated = true;
    	// estimate the conduction iterperiod based on current mechanical/fluid timestep
   	conductionIterPeriod = int(tsSafetyFactor * maxTimeStep/scene->dt);
	elapsedIters = 0;
	elapsedTime = 0;
	timeStepEstimated = true;
	runConduction = false;
}

void ThermalEngine::resetBoundaryFluxSums() {
	for (int i=0; i<6; i++) thermalBndFlux[i]=0;
}

void ThermalEngine::setConductionBoundary() {

	for (int k=0;k<6;k++)	{
		flow->solver->conductionBoundary(flow->wallIds[k]).fluxCondition=!bndCondIsTemperature[k];
        flow->solver->conductionBoundary(flow->wallIds[k]).value=thermalBndCondValue[k];
	}

        RTriangulation& Tri = flow->solver->T[flow->solver->currentTes].Triangulation();
	const shared_ptr<BodyContainer>& bodies=scene->bodies;
        for (int bound=0; bound<6;bound++) {
                int& id = *flow->solver->boundsIds[bound];
		flow->solver->conductionBoundingCells[bound].clear();
		if (id<0) continue;
                CGT::ThermalBoundary& bi = flow->solver->conductionBoundary(id);
	
                if (!bi.fluxCondition) {
                        VectorCell tmpCells;
                        tmpCells.resize(10000);
                        VCellIterator cells_it = tmpCells.begin();
                        VCellIterator cells_end = Tri.incident_cells(flow->solver->T[flow->solver->currentTes].vertexHandles[id],cells_it);
			
                        for (VCellIterator it = tmpCells.begin(); it != cells_end; it++){
				CellHandle& cell = *it;
				for (int v=0;v<4;v++) {
					if (!cell->vertex(v)->info().isFictious){
						const long int id = cell->vertex(v)->info().id();
						const shared_ptr<Body>& b =(*bodies)[id];
						if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
						auto* thState = b->state.get();
						thState->Tcondition=true;
						thState->temp=bi.value;
						thState->boundaryId=bound;
					}
				}
				flow->solver->conductionBoundingCells[bound].push_back(cell);
			}
                }
        }
	boundarySet = true;
}

void ThermalEngine::initializeInternalEnergy() {
	const shared_ptr<BodyContainer>& bodies=scene->bodies;
	const long size=bodies->size();
//	#pragma omp parallel for
	for (long i=0; i<size; i++){
		const shared_ptr<Body>& b=(*bodies)[i];
		if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
		auto* thState = b->state.get();
		thState->U = thState->Cp*thState->mass*thState->temp;
	}
	energySet = true;
}

void ThermalEngine::computeSolidFluidFluxes() {		
	if (!flow->solver->sphericalVertexAreaCalculated) computeVertexSphericalArea();
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	Tesselation& Tes = flow->solver->T[flow->solver->currentTes];
//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
//	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
    	for (long i=0; i<size; i++){
		CellHandle& cell = Tes.cellHandles[i];
//	#else	
		if (cell->info().isGhost) continue; // Do we need special cases for fictious cells?
		for (int v=0;v<4;v++){
			//if (cell->vertex(v)->info().isFictious) continue;
			const long int id = cell->vertex(v)->info().id();
			const shared_ptr<Body>& b =(*bodies)[id];
			if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
			const double surfaceArea = cell->info().sphericalVertexSurface[v];
			computeFlux(cell,b,surfaceArea);
		}
	}
}


void ThermalEngine::computeVertexSphericalArea(){
	Tesselation& Tes = flow->solver->T[flow->solver->currentTes];
//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
//	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
    	for (long i=0; i<size; i++){
		CellHandle& cell = Tes.cellHandles[i];
//	#else		
		if (cell->info().isFictious) continue;
		VertexHandle W[4];
		for (int k=0;k<4;k++) W[k] = cell->vertex(k);
        	cell->info().sphericalVertexSurface[0]=flow->solver->fastSphericalTriangleArea(W[0]->point(),W[1]->point().point(),W[2]->point().point(),W[3]->point().point());
       		cell->info().sphericalVertexSurface[1]=flow->solver->fastSphericalTriangleArea(W[1]->point(),W[0]->point().point(),W[2]->point().point(),W[3]->point().point());
        	cell->info().sphericalVertexSurface[2]=flow->solver->fastSphericalTriangleArea(W[2]->point(),W[1]->point().point(),W[0]->point().point(),W[3]->point().point());
        	cell->info().sphericalVertexSurface[3]=flow->solver->fastSphericalTriangleArea(W[3]->point(),W[1]->point().point(),W[2]->point().point(),W[0]->point().point());
   	 }
	flow->solver->sphericalVertexAreaCalculated=true;
}


void ThermalEngine::computeFlux(CellHandle& cell,const shared_ptr<Body>& b, const double surfaceArea) {
	const Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
	auto* thState = b->state.get(); 
	const double h = 2. * fluidK / (2.*sphere->radius); // heat transfer coeff assuming Re<<1 (stokes flow)
	const double flux = h*surfaceArea*(cell->info().temp() - thState->temp);
   	if (runConduction) thState->stabilityCoefficient += h*surfaceArea; // for auto time step estimation
	if (!cell->info().Tcondition) cell->info().internalEnergy -= flux*thermalDT; 
	if (!thState->Tcondition) thState->stepFlux += flux;  
}


void ThermalEngine::computeSolidSolidFluxes() {
//	#ifdef YADE_OPENMP
	const shared_ptr<InteractionContainer>& interactions=scene->interactions;
	const long size=interactions->size();
//	#pragma omp parallel for
	for (long i=0; i<size; i++){
		const shared_ptr<Interaction>& I=(*interactions)[i];
//	#else
//	for (const auto & I : *scene->interactions){
//	#endif
		const ScGeom* geom;
		if (!I || !I->geom.get() || !I->phys.get() || !I->isReal()) continue;
		if (I->geom.get())
		{
		geom = YADE_CAST<ScGeom*> (I->geom.get());
		if (!geom) continue;
		const double pd = geom->penetrationDepth;
		const shared_ptr<Body>& b1_=Body::byId(I->getId1(),scene);
		const shared_ptr<Body>& b2_=Body::byId(I->getId2(),scene);
		if (b1_->shape->getClassIndex()!=Sphere::getClassIndexStatic() || b2_->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b1_ || !b2_) continue;
		auto* thState1 = b1_->state.get(); 
		auto* thState2 = b2_->state.get(); 
		Sphere* sphere1 = dynamic_cast<Sphere*>(b1_->shape.get());
		Sphere* sphere2 = dynamic_cast<Sphere*>(b2_->shape.get());
		
		//double& U1 = thState1->U;
		//double& U2 = thState2->U;	
		const double k1 = thState1->k;
		const double k2 = thState2->k;
		const double r1 = sphere1->radius;
		const double r2 = sphere2->radius;
		const double T1 = thState1->temp;
		const double T2 = thState2->temp;
		const double d = r1 + r2 - pd;
		if (d==0) continue;
		double R = 0;
		double r = 0;
		// for equation: 
		if (r1 >= r2) {R = r1;r = r2;}
		else if (r1 < r2) {R = r2;r = r1;} 
		// The radius of the intersection found by: Kern, W. F. and Bland, J. R. Solid Mensuration with Proofs, 2nd ed. New York: Wiley, p. 97, 1948.	http://mathworld.wolfram.com/Sphere-SphereIntersection.html	
        //const double area = M_PI*pow(rc,2);
		
		//const double dt = scene->dt;		
		//const double fluxij = 4.*rc*(T1-T2) / (1./k1 + 1./k2);
		
		// compute the overlapping volume for thermodynamic considerations
//		const double capHeight1 = (r1-r2+d)*(r1+r2-d)/2*d;
//		const double capHeight2 = (r2-r1+d)*(r2+r1-d)/2*d;
//		thState1->capVol += (1./3.)*M_PI*pow(capHeight1,2)*(3.*r1-capHeight1);
//		thState2->capVol += (1./3.)*M_PI*pow(capHeight2,2)*(3.*r2-capHeight2);

		// compute the thermal resistance of the pair and the associated flux
		double thermalResist;
		if (useKernMethod) {
            		const double numerator = pow((-d+r-R)*(-d-r+R)*(-d+r+R)*(d+r+R),0.5);	
            		const double rc = numerator / (2.*d);
            		thermalResist = 4.*rc / (1./k1 + 1./k2);
		}//thermalResist = ((k1+k2)/2.)*area/(r1+r2-pd);}
		else {thermalResist = 2.*(k1+k2)*r1*r2 / (r1+r2-pd);}
		const double fluxij = thermalResist * (T1-T2);

		//cout << "Flux b/w "<< b1_->id << " & "<< b2_->id << " fluxij " << fluxij << endl;
        	if (runConduction){
           		thState1->stabilityCoefficient += thermalResist;
            		thState2->stabilityCoefficient += thermalResist;
        	}
		if (!thState1->Tcondition) thState1->stepFlux -= fluxij; //U1 -= fluxij*dt;
		else thermalBndFlux[thState1->boundaryId] -= fluxij;
		if (!thState2->Tcondition) thState2->stepFlux += fluxij; // U2 += fluxij*dt;
		else thermalBndFlux[thState2->boundaryId] += fluxij;
		}
	}
}

void ThermalEngine::computeNewPoreTemperatures() {
    flow->solver->setNewCellTemps();
}

void ThermalEngine::computeNewParticleTemperatures() {
    const shared_ptr<BodyContainer>& bodies=scene->bodies;
    const long size=bodies->size();
//	#pragma omp parallel for
    for (long i=0; i<size; i++){
        const shared_ptr<Body>& b=(*bodies)[i];
        if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
        auto* thState = b->state.get();
        Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
        const double density = b->material->density;
        const double volume = 4./3. * M_PI *pow(sphere->radius,3); // - thState->capVol;
        if (thState->Tcondition) continue;
//		if (!thState->oldTempSet){
        thState->oldTemp = thState->temp; 
        thState->temp = thState->stepFlux*thermalDT/(thState->Cp*density*volume) + thState->oldTemp; // first order forward difference
        thState->stepFlux=0;
    }
}


void ThermalEngine::thermalExpansion() {
	const shared_ptr<BodyContainer>& bodies=scene->bodies;
	const long size=bodies->size();
	// adjust particle size
	for (long i=0; i<size; i++){
		const shared_ptr<Body> b=(*bodies)[i];
		if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
		Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());		
		auto* thState = b->state.get();
       		if (!thState->Tcondition) {
			thState->delRadius = thState->alpha * sphere->radius * (thState->temp - thState->oldTemp);
			sphere->radius += thState->delRadius;
            //sphere->radius += thState->alpha * sphere->radius * (thState->temp - thState->oldTemp);
        	}
	}
	
	// adjust cell pressure
	if (fluidBeta <= 0) return;
   	Tesselation& Tes = flow->solver->T[flow->solver->currentTes];
//	#ifdef YADE_OPENMP
	const long sizeCells = Tes.cellHandles.size();
//	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
    	for (long i=0; i<sizeCells; i++){
		CellHandle& cell = Tes.cellHandles[i];
//	#else	
		if (cell->info().isGhost || cell->info().Pcondition || cell->info().Tcondition) continue; // Do we need special cases for fictious cells?
		Real term1 = 0;
        	Real term2 = 0;
		if (solidThermoMech)  term1 = computeCellPressureChangeFromDeltaVolume(cell);
        	if (fluidThermoMech)  term2 = computeCellPressureChangeFromDeltaTemp(cell);
		cell->info().p() += term1 + term2;
//        cell->info().p() += fluidK * ( (1./(1+fluidBeta*(cell->info().dtemp()))) - 1.);
	}
	
	
}


const Real ThermalEngine::computeCellPressureChangeFromDeltaVolume(CellHandle& cell){
    const shared_ptr<BodyContainer>& bodies=scene->bodies;
    const Real K = flow->fluidBulkModulus;
    const Real Vo = 1./cell->info().invVoidVolume(); // FIXME: depends on volumeSolidPore(), which only uses CGAL triangulations so this value is only updated with remeshes. We need our own function for volumeSolidPore...
    Real sum = 0;
    for (int v=0;v<4;v++){
        const long int id = cell->vertex(v)->info().id();
        const shared_ptr<Body>& b =(*bodies)[id];
        if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
        Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());	
        auto* thState = b->state.get();
        const double surfaceArea = cell->info().sphericalVertexSurface[v]; // from last fluid engine remesh. We could store this in cell->info() if it becomes a burden
        const Real rCubedDiff = pow(sphere->radius,3) - pow(sphere->radius - thState->delRadius,3);
        const double areaTot = 4. * M_PI * pow(sphere->radius,2); // We could store this in cell->info() if it becomes a burden
        sum += surfaceArea/areaTot * rCubedDiff * M_PI * 4./3. ;
    }
    const Real delP = K*Vo*sum;
    return delP;
}

const Real ThermalEngine::computeCellPressureChangeFromDeltaTemp(CellHandle& cell){
    const Real K = flow->fluidBulkModulus;
    const Real delP = K * ( (1./(1+fluidBeta*(cell->info().dtemp()))) - 1.); //dtemp() here is only the most recent step, which includes the change of temp due to our large conductionStep by design. Specifically does NOT include temp changes due to advection for conductionIterPeriod.
    return delP;
}

// void ThermalEngine::computeVolumeSolidPore(CellHandle& Cell){
//     const shared_ptr<BodyContainer>& bodies=scene->bodies;
// 	Tesselation& Tes = flow->solver->T[flow->solver->currentTes];
// //	#ifdef YADE_OPENMP
// 	const long size = Tes.cellHandles.size();
// //	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
//     for (long i=0; i<size; i++){
// 		CellHandle& cell = Tes.cellHandles[i];
//         double vSolid; 
//         VertexHandle W[4];
//         Body
// 		for (int k=0;k<4;k++){
//             
//             const long int id = cell->vertex(k)->info().id();
// 			const shared_ptr<Body>& b =(*bodies)[id];
//             for (int l=0;l<4;k++) W[l] = cell->vertex(l);
// 			if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
//             Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());	
//             W[k] = cell->vertex(k);
//             vSolid += sphericalTriangleVolume(*sphere,W[1]->point().point,W[2]->point().point(),W[3]->point().point())
//             
// }
// 
// double ThermalEngine::sphericalTriangleVolume(const Sphere& ST1, const Point& PT1, const Point& PT2, const Point& PT3)
// {
//         double rayon = sqrt(ST1.weight());
//         if (rayon == 0.0) return 0.0;
//         return ((ONE_THIRD * rayon) * (fastSphericalTriangleArea(ST1, PT1, PT2, PT3))) ;
// }
// 
// Real ThermalEngine::fastSolidAngle(const Point& STA1, const Point& PTA1, const Point& PTA2, const Point& PTA3)
// {
//         //! This function needs to be fast because it is used heavily. Avoid using vector operations which require constructing vectors (~50% of cpu time in the non-fast version), and compute angle using the 3x faster formula of Oosterom and StrackeeVan Oosterom, A; Strackee, J (1983). "The Solid Angle of a Plane Triangle". IEEE Trans. Biom. Eng. BME-30 (2): 125-126. (or check http://en.wikipedia.org/wiki/Solid_angle)
//         using namespace CGAL;
//         double M[3][3];
//         M[0][0] = PTA1.x() - STA1.x();
//         M[0][1] = PTA2.x() - STA1.x();
//         M[0][2] = PTA3.x() - STA1.x();
//         M[1][0] = PTA1.y() - STA1.y();
//         M[1][1] = PTA2.y() - STA1.y();
//         M[1][2] = PTA3.y() - STA1.y();
//         M[2][0] = PTA1.z() - STA1.z();
//         M[2][1] = PTA2.z() - STA1.z();
//         M[2][2] = PTA3.z() - STA1.z();
// 
//         double detM = M[0][0]* (M[1][1]*M[2][2]-M[2][1]*M[1][2]) +
//                       M[1][0]* (M[2][1]*M[0][2]-M[0][1]*M[2][2]) +
//                       M[2][0]* (M[0][1]*M[1][2]-M[1][1]*M[0][2]);
// 
//         double pv12N2 = pow(M[0][0],2) +pow(M[1][0],2) +pow(M[2][0],2);
//         double pv13N2 = pow(M[0][1],2) +pow(M[1][1],2) +pow(M[2][1],2);
//         double pv14N2 = pow(M[0][2],2) +pow(M[1][2],2) +pow(M[2][2],2);
// 
//         double pv12N = sqrt(pv12N2);
//         double pv13N = sqrt(pv13N2);
//         double pv14N = sqrt(pv14N2);
// 
//         double cp12 = (M[0][0]*M[0][1]+M[1][0]*M[1][1]+M[2][0]*M[2][1]);
//         double cp13 = (M[0][0]*M[0][2]+M[1][0]*M[1][2]+M[2][0]*M[2][2]);
//         double cp23 = (M[0][1]*M[0][2]+M[1][1]*M[1][2]+M[2][1]*M[2][2]);
// 
//         double ratio = detM/ (pv12N*pv13N*pv14N+cp12*pv14N+cp13*pv13N+cp23*pv12N);
//         return abs(2*atan(ratio));
// }


#endif//THERMAL
