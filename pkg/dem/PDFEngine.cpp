#include "PDFEngine.hpp"

void PDFEngine::getSpectrums(vector<PDFEngine::PDF> & pdfs)
{
	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	
	vector<int> nTheta, nPhi;
	vector<Real> dTheta, dPhi;
	Real V;
	int N = scene->bodies->size();
	
	if(!scene->isPeriodic) {
		LOG_WARN("Volume is based on periodic cell volume. Set V = 1");
		V = 1.;
	}
	else
		V = scene->cell->getVolume();
	
	for(uint i(0);i<pdfs.size();i++)
	{
		nTheta.push_back(pdfs[i].shape()[0]);
		nPhi.push_back(pdfs[i].shape()[1]);
		dTheta.push_back(Mathr::PI / nTheta[i]);
		dPhi.push_back(2.*Mathr::PI / nPhi[i]);
	}
	
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions) {
		if(!I->isReal()) continue;
		GenericSpheresContact* geom=dynamic_cast<GenericSpheresContact*>(I->geom.get());
		
		if(geom)
		{
			Real theta 	= acos(geom->normal.y()); //[0;pi]
 			Real phi 	= atan2(geom->normal.z(), geom->normal.x()) + Mathr::PI; //[-pi;pi] => [0; 2pi] for index calculation
 			
 			for(uint i(0);i<pdfs.size();i++)
			{
				Real dS = sin(theta) * dTheta[i] * dPhi[i];
				
				// Calculate indexes
				int idT1 = ((int)round((theta) / dTheta[i])) % nTheta[i];
				int idP1 = ((int)round((phi) / dPhi[i])) % nPhi[i];
				int idT2 = ((int)round((Mathr::PI - theta) / dTheta[i])) % nTheta[i];
				int idP2 = ((int)round((Mathr::PI + phi) / dPhi[i])) % nPhi[i];
				
				if(pdfs[i][idT1][idP1]) pdfs[i][idT1][idP1]->addData(I, dS, V, N);
				if(pdfs[i][idT2][idP2]) pdfs[i][idT2][idP2]->addData(I, dS, V, N);
			}
		}
	}
}

void PDFEngine::writeToFile(vector<PDFEngine::PDF> const& pdfs)
{
	FILE* fid = fopen(filename.c_str(), (firstRun) ? "w" : "a");
	
	if(fid) {
		if(firstRun) {
			
			fprintf(fid, "# ");
			for(uint i(0);i<pdfs.size();i++) {
				uint nTheta = pdfs[i].shape()[0];
				uint nPhi = pdfs[i].shape()[1];
				Real dTheta = (Mathr::PI / nTheta);
				Real dPhi = (2.*Mathr::PI / nPhi);
				
				for(uint t(0);t<nTheta;t++) for(uint p(0);p<nPhi;p++) if(pdfs[i][t][p]) {
					vector<string> ss = pdfs[i][t][p]->getSuffixes();
					
					if(ss.size() > 1)
						for(uint j(0);j<ss.size();j++)
							fprintf(fid, "%s_%s(%f,%f)\t",pdfs[i][t][p]->name.c_str(),  ss[j].c_str(), t*dTheta, p*dPhi - Mathr::PI);
					else
						fprintf(fid, "%s(%f,%f)\t", pdfs[i][t][p]->name.c_str(), t*dTheta, p*dPhi - Mathr::PI);
				}
			}
			firstRun = false;
			fprintf(fid, "\n");
		}
		
		for(uint i(0);i<pdfs.size();i++) for(uint t(0);t<pdfs[i].shape()[0];t++) for(uint p(0);p<pdfs[i].shape()[1];p++)
			if(pdfs[i][t][p]) {
				vector<string> dat = pdfs[i][t][p]->getDatas();

				for(uint j(0);j<dat.size();j++)
					fprintf(fid, "%s\t",dat[j].c_str());
			}
		fprintf(fid, "\n");
		fclose(fid);
	}
	else {
		if(!warnedOnce) LOG_ERROR("Unable to open " << filename << " for PDF writing");
		warnedOnce = true;
	}
}

void PDFEngine::action()
{
	vector<PDFEngine::PDF> pdfs;
	pdfs.resize(5);
	
	for(uint i(0);i<pdfs.size();i++) {
		pdfs[i].resize(boost::extents[numDiscretizeAngleTheta][numDiscretizeAnglePhi]);
	}
	
	// Hint: If you want data on particular points, allocate only those pointers.
	for(uint t(0);t<numDiscretizeAngleTheta;t++) for(uint p(0);p<numDiscretizeAnglePhi;p++) {
		pdfs[0][t][p] = shared_ptr<PDFCalculator>(new PDFSpheresStressCalculator<NormPhys>(&NormPhys::normalForce, "normalStress"));
		pdfs[1][t][p] = shared_ptr<PDFCalculator>(new PDFSpheresStressCalculator<NormShearPhys>(&NormShearPhys::shearForce, "shearStress"));
		pdfs[2][t][p] = shared_ptr<PDFCalculator>(new PDFSpheresDistanceCalculator("h"));
		pdfs[3][t][p] = shared_ptr<PDFCalculator>(new PDFSpheresVelocityCalculator("v"));
		pdfs[4][t][p] = shared_ptr<PDFCalculator>(new PDFSpheresIntrsCalculator("P"));
	}
	
	getSpectrums(pdfs); // Where the magic happen :)
	writeToFile(pdfs);
}


CREATE_LOGGER(PDFEngine);

PDFSpheresDistanceCalculator::PDFSpheresDistanceCalculator(string name) :
	PDFEngine::PDFCalculator(name),
	m_h(0.),
	m_N(0)
{
	
}

vector<string> PDFSpheresDistanceCalculator::getDatas() const
{
	return vector<string>({std::to_string(m_h/m_N)});
}

void PDFSpheresDistanceCalculator::cleanData()
{
	m_h = 0.;
	m_N = 0;
}

bool PDFSpheresDistanceCalculator::addData(const shared_ptr<Interaction>& I, Real const& dS ,Real const& V, int const& N)
{
	if(!I->isReal()) return false;
	ScGeom* geom=dynamic_cast<ScGeom*>(I->geom.get());
    Real a((geom->radius1+geom->radius2)/2.);
	
	if(!geom)
		return false;
	
	m_N++;
	m_h -= geom->penetrationDepth/a;
	
	return true;
}

PDFSpheresVelocityCalculator::PDFSpheresVelocityCalculator(string name) :
	PDFEngine::PDFCalculator(name),
	m_vel(Vector3r::Zero()),
	m_N(0)
{
	
}

vector<string> PDFSpheresVelocityCalculator::getSuffixes() const
{
	return vector<string>({"x","y","z"});
}

vector<string> PDFSpheresVelocityCalculator::getDatas() const
{
	vector<string> ret;
	for(int i(0);i<3;i++) ret.push_back(std::to_string(m_vel(i)/m_N));
	return ret;
}

void PDFSpheresVelocityCalculator::cleanData()
{
	m_vel = Vector3r::Zero();
	m_N = 0;
}

bool PDFSpheresVelocityCalculator::addData(const shared_ptr<Interaction>& I, Real const& dS ,Real const& V, int const& N)
{
	if(!I->isReal()) return false;
	
	// Geometry
    ScGeom* geom=static_cast<ScGeom*>(I->geom.get());
	if(!geom) return false;

    // geometric parameters
   // Real a((geom->radius1+geom->radius2)/2.);
    Vector3r relV = geom->getIncidentVel_py(I, false);
	
	m_N++;
	m_vel += relV;
	
	return true;
}

PDFSpheresIntrsCalculator::PDFSpheresIntrsCalculator(string name) :
	PDFEngine::PDFCalculator(name),
	m_P(0.)
{
	
}

vector<string> PDFSpheresIntrsCalculator::getDatas() const
{
	return vector<string>({std::to_string(m_P)});
}

void PDFSpheresIntrsCalculator::cleanData()
{
	m_P = 0.;
}

bool PDFSpheresIntrsCalculator::addData(const shared_ptr<Interaction>& I, Real const& dS ,Real const& V, int const& N)
{
	if(!I->isReal()) return false;
	
	m_P += 1./(dS*N);
	
	return true;
}
