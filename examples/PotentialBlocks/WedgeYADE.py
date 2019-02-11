#CW BOON 2018
# Use the following algorithms:
# CW Boon, GT Houlsby, S Utili (2012).  A new algorithm for contact detection between convex polygonal and polyhedral particles in the discrete element method.  Computers and Geotechnics 44, 73-82. 
# CW Boon, GT Houlsby, S Utili (2015).  A new rock slicing method based on linear programming.  Computers and Geotechnics 65, 12-29. 

#Users need to update westBodyId, midBodyId, eastBodyID by displaying the IDs of the body and re-run. Otherwise will receive error. 
#The code runs some steps to stabilise.  After that a vertical cut is carried out.  And the simulations shows the failure mechanism of the slope.
#Display is saved to a vtk file in the "vtk folder" and the user is required to load it using ParaView.  Control the frequency of printing a vtk file using vtkRecorder.iterPeriod in python

#Disclaimer: This script is just for illustration to demonstrate the function of Block Gen and PotentialBlock Code, and not for other purpose outside for this intended use

#To use this script:
#Compile with
#ENABLE_POTENTIAL_BLOCKS=ON, and add  sudo apt-get install
#coinor-clp, 
#coinor-libclp-dev, 
#coinor-libclp1, 
#coinor-libosi1v5

import os
import errno
try:
   os.mkdir('./vtk/')
except OSError as exc:
   if exc.errno != errno.EEXIST:
      raise
   pass

name ='BlockGeneration' # PLEASE CHECK THIS  

p=BlockGen()
p.maxRatio = 10000.0;
p.minSize =7.0;
p.density = 2700 # 23.6/9.81*1000.0 	#kg/m3
p.dampingMomentum = 0.8
p.damp3DEC = 0.0
p.viscousDamping = 0.0
p.Kn = 2.0e8  	#1GPa
p.Ks= 0.1e8 	#1GPa
p.frictionDeg = 50.0 #degrees
p.traceEnergy = False
p.defaultDt = 1e-4
p.rForPP = 0.4 #0.04
p.kForPP = 0.55
p.RForPP = 1800.0
p.gravity = [0,0,9.81]
p.inertiaFactor = 1.0
p.initialOverlap = 1e-6 
p.numberOfGrids = 100
p.exactRotation = False
p.shrinkFactor = 1.0
p.boundarySizeXmin = 20.0; #South
p.boundarySizeXmax = 20.0; #North
p.boundarySizeYmin = 20.0; #West
p.boundarySizeYmax = 20; #East 4100
p.boundarySizeZmin = 0.0; #Up
p.boundarySizeZmax = 40.0; #Down
p.extremeDist = 50.5;
p.subdivisionRatio = 0.05;
p.persistentPlanes = False
p.jointProbabilistic = True
p.opening = False
p.boundaries = True
p.slopeFace = False
p.calJointLength = False
p.twoDimension = False
p.unitWidth2D = 9.0
p.intactRockDegradation = True
p.calAreaStep = 20.0
p.useFaceProperties = False
p.neverErase = False # Must be used when tension is on
p.peakDisplacement = 0.1
p.brittleLength = 1.0
p.maxClosure = 0.0003 
p.sliceBoundaries = True
p.filenamePersistentPlanes = ''
p.filenameProbabilistic = './joints/jointC.csv'  #  PLEASE CHEK THIS  
p.filenameOpening = '' #'./joints/opening.csv'
p.filenameBoundaries = ''
p.filenameSlopeFace =''
p.filenameSliceBoundaries = ''
p.directionA=Vector3(1,0,0)
p.directionB=Vector3(0,1,0)
p.directionC=Vector3(0,0,1)
p.load()
O.engines[2].lawDispatcher.functors[0].initialOverlapDistance = p.initialOverlap - 1e-6
O.engines[2].lawDispatcher.functors[0].allowBreakage = False
#Gl1_PotentialBlock.sizeX = 50
#Gl1_PotentialBlock.sizeY = 50
#Gl1_PotentialBlock.sizeZ = 50

O.engines[2].physDispatcher.functors[0].ks_i = 0.1e9
O.engines[2].physDispatcher.functors[0].kn_i = 2e9

from yade import plot

rockFriction = 40.0
boundaryFriction = 40.0
targetFriction = 40.0
waterHeight = 460.0  # PLEASE CHEK THIS 460.0 # 
startCountingBrokenBonds = False
minTimeStep = 1000000.0
westBodyId=10 #PLEASE CHEK THIS  
midBodyId =10 #PLEASE CHEK THIS  
eastBodyId =10 #PLEASE CHEK THIS  
originalPositionW = O.bodies[westBodyId].state.pos
originalPositionE = O.bodies[eastBodyId].state.pos
originalPositionM = O.bodies[midBodyId].state.pos
velocityDependency = False

O.engines[3].label='integration'
O.dt = 10.0e-5 #10e-4
O.engines = O.engines + [PotentialBlockVTKRecorderTunnel(fileName='./vtk/Wedge',iterPeriod=5000,twoDimension=False,sampleX=70,sampleY=70,sampleZ=70,maxDimension=100.0,REC_INTERACTION=True,REC_VELOCITY=True,label='vtkRecorder')]
O.engines = O.engines+ [PyRunner(iterPeriod=5000,command='calTimeStep()')] 
O.engines=O.engines+[PyRunner(iterPeriod=200,command='myAddPlotData()')]

O.engines = O.engines+ [PyRunner(iterPeriod=200,command='goToNextStage2()',label='dispChecker')] 


O.step()


idCountBegin = len(O.bodies)
count = idCountBegin

# material
young = 600.0e6 # [N/m^2]
poisson = 0.6 
frictionAngle=radians(60.0) # [rad]
# append geometry and material
O.materials.append(FrictMat(young=young,poisson=poisson,density=p.density,frictionAngle=frictionAngle, label='frictionless'))

#VERTICAL BOUNDARIES
bBottom = Body()
thickness = 0.1*p.boundarySizeZmax

wire=False
color=[0,0,1]
highlight=False
kPP = p.kForPP
rPP = p.rForPP
length = 2.0*thickness+2*rPP
RPP = 0.5*p.boundarySizeXmax
aPP = [1,-1,0,0,0,0]
bPP = [0,0,1,-1,0,0]
cPP = [0,0,0,0,1,-1]
dPP = [p.boundarySizeXmax , p.boundarySizeXmax ,p.boundarySizeYmax,p.boundarySizeYmax,thickness,thickness]
minmaxAabb = Vector3(1.05*(dPP[0]+rPP),1.05*(dPP[2]-rPP),1.05*(dPP[4]+rPP))
bBottom.shape=PotentialBlock(k=kPP, r=rPP , R=RPP,a=aPP, b=bPP, c=cPP, d=dPP,id=count,isBoundary=True,color=color ,wire=wire,highlight=highlight,AabbMinMax=True,minAabb=minmaxAabb ,maxAabb=minmaxAabb,minAabbRotated=minmaxAabb ,maxAabbRotated=minmaxAabb,fixedNormal=True,boundaryNormal=Vector3(0,0,-1))	
V = (2*(thickness+rPP))*(p.boundarySizeXmax+rPP)*(p.boundarySizeYmax)
geomInert = 1/12*(p.boundarySizeYmax)*(length**3)
geomInertX = 1/12*V*p.density*(p.boundarySizeYmax**2 + (2.0*thickness+2.0*rPP)**2 ) 
geomInertY = 1/12*V*p.density*((2.0*p.boundarySizeXmax+2.0*rPP)**2 + (2.0*thickness+2.0*rPP)**2 ) 
geomInertZ = 1/12*V*p.density*((2.0*p.boundarySizeXmax+2.0*rPP)**2 + p.boundarySizeYmax**2) 
utils._commonBodySetup(bBottom,V,Vector3(geomInertX,geomInertY,geomInertZ), material='frictionless',pos= [0.0,0.0,p.boundarySizeZmax+rPP+thickness+rPP], dynamic=True, fixed=True) 
bBottom.state.pos = [0.0,0.0,p.boundarySizeZmax+rPP+thickness+rPP]
bBottom.shape.isBoundary=True
bBottom.state.mass = p.density*V
bBottom.dynamic=False
O.bodies.append(bBottom)
count=count+1

bSideA = Body()
thickness = 0.1*p.boundarySizeZmax

wire=False
color=[0,0,1]
highlight=False
kPP = p.kForPP
rPP = p.rForPP
length = 2.0*thickness+2*rPP
RPP = 0.5*p.boundarySizeXmax
aPP = [1,-1,0,0,0,0]
bPP = [0,0,1,-1,0,0]
cPP = [0,0,0,0,1,-1]
dPP = [p.boundarySizeXmax , p.boundarySizeXmax ,thickness,thickness,0.5*p.boundarySizeZmax,0.5*p.boundarySizeZmax]
minmaxAabb = Vector3(1.05*(dPP[0]+rPP),1.05*(dPP[2]-rPP),1.05*(dPP[4]+rPP))
bSideA.shape=PotentialBlock(k=kPP, r=rPP , R=RPP,a=aPP, b=bPP, c=cPP, d=dPP,id=count,isBoundary=True,color=color ,wire=wire,highlight=highlight,AabbMinMax=True,minAabb=minmaxAabb ,maxAabb=minmaxAabb,minAabbRotated=minmaxAabb ,maxAabbRotated=minmaxAabb,fixedNormal=True,boundaryNormal=Vector3(0,-1.0,0))	
V = (2*(thickness+rPP))*(p.boundarySizeXmax+rPP)*(p.boundarySizeYmax)
geomInert = 1/12*(p.boundarySizeYmax)*(length**3)
geomInertX = 1/12*V*p.density*(p.boundarySizeYmax**2 + (2.0*thickness+2.0*rPP)**2 ) 
geomInertY = 1/12*V*p.density*((2.0*p.boundarySizeXmax+2.0*rPP)**2 + (2.0*thickness+2.0*rPP)**2 ) 
geomInertZ = 1/12*V*p.density*((2.0*p.boundarySizeXmax+2.0*rPP)**2 + p.boundarySizeYmax**2) 
utils._commonBodySetup(bSideA,V,Vector3(geomInertX,geomInertY,geomInertZ), material='frictionless',pos= [0.0,0.0,p.boundarySizeZmax+rPP+thickness+rPP], dynamic=True, fixed=True) 
bSideA.state.pos = [0.0,p.boundarySizeYmax+rPP+thickness+rPP,0.5*p.boundarySizeZmax]
bSideA.shape.isBoundary=True
bSideA.state.mass = p.density*V
bSideA.dynamic=False
O.bodies.append(bSideA)
count=count+1

bSideB = Body()


wire=False
color=[0,0,1]
highlight=False
kPP = p.kForPP
rPP = p.rForPP
length = 2.0*thickness+2*rPP
RPP = 0.5*p.boundarySizeXmax
aPP = [1,-1,0,0,0,0]
bPP = [0,0,1,-1,0,0]
cPP = [0,0,0,0,1,-1]
dPP = [p.boundarySizeXmax , p.boundarySizeXmax ,thickness,thickness,0.5*p.boundarySizeZmax,0.5*p.boundarySizeZmax]
minmaxAabb = Vector3(1.05*(dPP[0]+rPP),1.05*(dPP[2]-rPP),1.05*(dPP[4]+rPP))
bSideB.shape=PotentialBlock(k=kPP, r=rPP , R=RPP,a=aPP, b=bPP, c=cPP, d=dPP,id=count,isBoundary=True,color=color ,wire=wire,highlight=highlight,AabbMinMax=True,minAabb=minmaxAabb ,maxAabb=minmaxAabb,minAabbRotated=minmaxAabb ,maxAabbRotated=minmaxAabb,fixedNormal=True,boundaryNormal=Vector3(0,1.0,0.0))	
V = (2*(thickness+rPP))*(p.boundarySizeXmax+rPP)*(p.boundarySizeYmax)
geomInert = 1/12*(p.boundarySizeYmax)*(length**3)
geomInertX = 1/12*V*p.density*(p.boundarySizeYmax**2 + (2.0*thickness+2.0*rPP)**2 ) 
geomInertY = 1/12*V*p.density*((2.0*p.boundarySizeXmax+2.0*rPP)**2 + (2.0*thickness+2.0*rPP)**2 ) 
geomInertZ = 1/12*V*p.density*((2.0*p.boundarySizeXmax+2.0*rPP)**2 + p.boundarySizeYmax**2) 
utils._commonBodySetup(bSideB,V,Vector3(geomInertX,geomInertY,geomInertZ), material='frictionless',pos= [0.0,0.0,p.boundarySizeZmax+rPP+thickness+rPP], dynamic=True, fixed=True) 
bSideB.state.pos = [0.0,-p.boundarySizeYmax-rPP-thickness-rPP,0.5*p.boundarySizeZmax]
bSideB.shape.isBoundary=True
bSideB.state.mass = p.density*V
bSideB.dynamic=False
O.bodies.append(bSideB)
count=count+1

bSideC = Body()


wire=False
color=[0,0,1]
highlight=False
kPP = p.kForPP
rPP = p.rForPP
length = 2.0*thickness+2*rPP
RPP = 0.5*p.boundarySizeXmax
aPP = [1,-1,0,0,0,0]
bPP = [0,0,1,-1,0,0]
cPP = [0,0,0,0,1,-1]
dPP = [thickness , thickness ,p.boundarySizeYmax,p.boundarySizeYmax,0.5*p.boundarySizeZmax,0.5*p.boundarySizeZmax]
minmaxAabb = Vector3(1.05*(dPP[0]+rPP),1.05*(dPP[2]-rPP),1.05*(dPP[4]+rPP))
bSideC.shape=PotentialBlock(k=kPP, r=rPP , R=RPP,a=aPP, b=bPP, c=cPP, d=dPP,id=count,isBoundary=True,color=color ,wire=wire,highlight=highlight,AabbMinMax=True,minAabb=minmaxAabb ,maxAabb=minmaxAabb,minAabbRotated=minmaxAabb ,maxAabbRotated=minmaxAabb,fixedNormal=True,boundaryNormal=Vector3(1.0,0,0.0))	
V = (2*(thickness+rPP))*(p.boundarySizeXmax+rPP)*(p.boundarySizeYmax)
geomInert = 1/12*(p.boundarySizeYmax)*(length**3)
geomInertX = 1/12*V*p.density*(p.boundarySizeYmax**2 + (2.0*thickness+2.0*rPP)**2 ) 
geomInertY = 1/12*V*p.density*((2.0*p.boundarySizeXmax+2.0*rPP)**2 + (2.0*thickness+2.0*rPP)**2 ) 
geomInertZ = 1/12*V*p.density*((2.0*p.boundarySizeXmax+2.0*rPP)**2 + p.boundarySizeYmax**2) 
utils._commonBodySetup(bSideC,V,Vector3(geomInertX,geomInertY,geomInertZ), material='frictionless',pos= [0.0,0.0,p.boundarySizeZmax+rPP+thickness+rPP], dynamic=True, fixed=True) 
bSideC.state.pos = [-p.boundarySizeXmax-rPP-thickness-rPP,0.0,0.5*p.boundarySizeZmax]
bSideC.shape.isBoundary=True
bSideC.state.mass = p.density*V
bSideC.dynamic=False
O.bodies.append(bSideC)
count=count+1

bSideD = Body()


wire=False
color=[0,0,1]
highlight=False
kPP = p.kForPP
rPP = p.rForPP
length = 2.0*thickness+2*rPP
RPP = 0.5*p.boundarySizeXmax
aPP = [1,-1,0,0,0,0]
bPP = [0,0,1,-1,0,0]
cPP = [0,0,0,0,1,-1]
dPP = [thickness , thickness ,p.boundarySizeYmax,p.boundarySizeYmax,0.5*p.boundarySizeZmax,0.5*p.boundarySizeZmax]
minmaxAabb = Vector3(1.05*(dPP[0]+rPP),1.05*(dPP[2]-rPP),1.05*(dPP[4]+rPP))
bSideD.shape=PotentialBlock(k=kPP, r=rPP , R=RPP,a=aPP, b=bPP, c=cPP, d=dPP,id=count,isBoundary=True,color=color ,wire=wire,highlight=highlight,AabbMinMax=True,minAabb=minmaxAabb ,maxAabb=minmaxAabb,minAabbRotated=minmaxAabb ,maxAabbRotated=minmaxAabb,fixedNormal=True,boundaryNormal=Vector3(-1.0,0,0.0))	
V = (2*(thickness+rPP))*(p.boundarySizeXmax+rPP)*(p.boundarySizeYmax)
geomInert = 1/12*(p.boundarySizeYmax)*(length**3)
geomInertX = 1/12*V*p.density*(p.boundarySizeYmax**2 + (2.0*thickness+2.0*rPP)**2 ) 
geomInertY = 1/12*V*p.density*((2.0*p.boundarySizeXmax+2.0*rPP)**2 + (2.0*thickness+2.0*rPP)**2 ) 
geomInertZ = 1/12*V*p.density*((2.0*p.boundarySizeXmax+2.0*rPP)**2 + p.boundarySizeYmax**2) 
utils._commonBodySetup(bSideD,V,Vector3(geomInertX,geomInertY,geomInertZ), material='frictionless',pos= [0.0,0.0,p.boundarySizeZmax+rPP+thickness+rPP], dynamic=True, fixed=True) 
bSideD.state.pos = [p.boundarySizeXmax+rPP+thickness+rPP,0.0,0.5*p.boundarySizeZmax]
bSideD.shape.isBoundary=True
bSideD.state.mass = p.density*V
bSideD.dynamic=False
O.bodies.append(bSideD)





def calTimeStep():
	global minTimeStep
	mkratio = 99999999.9
	maxK = 0.0
	minMass = 1.0e15
	for i in O.interactions:
		if i.isReal==True:
			dt1 = O.bodies[i.id1].state.mass/i.phys.Knormal_area
			dt2 = O.bodies[i.id2].state.mass/i.phys.Knormal_area
			if dt1 < dt2:
				presentDt = 0.15*sqrt(dt1)
				if minTimeStep > presentDt:
					minTimeStep = presentDt
					O.dt = minTimeStep
			else:
				presentDt = 0.15*sqrt(dt2)
				if minTimeStep > presentDt:
					minTimeStep = presentDt
					O.dt = minTimeStep
			
def excavate():
	for b in O.bodies:
		if NorthWall(b.state.pos[0],b.state.pos[1],b.state.pos[2]) <0.001:
			if b.isClumpMember == True and b.shape.isBoundary==False:
				O.bodies.deleteClumpMember(O.bodies[b.clumpId], b)
			elif b.isClump == False and b.shape.isBoundary==False:
				O.bodies.erase(b.id)
				continue
	O.bodies.erase(bSideA.id)
	O.bodies.erase(bSideC.id)



prevDistance = O.bodies[westBodyId].state.pos[0]
tolDistance = 0.003 #0.1
tolDistance2 = 0.05
SRinProgress = False
SRcounter = 0
checkIter = 0
prevKE = 0.0
currentKE = 0.0
tolKE =10e10
initBondedContacts = 0
initDispRate = -1.0
prevDispRate = 0

def changeKE(newKE):
	global tolKE
	tolKE = newKE

def changeTolDist(newTol):
	global tolDistance
	tolDistance = newTol


removeDampingBool = False

def goToNextStage2():
	global startCountingBrokenBonds
	global velocityDependency
	global waterHeight
	global boundaryFriction
	global rockFriction
	global targetFriction
	global prevDistance
	global originalPositionW
	global tolDistance
	global tolDistance2
	global checkIter
	global SRinProgress
	global SRcounter
	global prevKE
	global currentKE
	global tolKE
	global initDispRate
	global prevDispRate
	global removeDampingBool
	prevKE = currentKE
	KE = utils.kineticEnergy()
	currentKE = KE
	uf = utils.unbalancedForce()
	if O.iter>2000 and removeDampingBool == False:
		removeDamping()
		removeDampingBool=True	
	if O.iter>5000 and SRcounter == 0: # and uf<0.005:
		print O.iter
		O.pause()
		vtkRecorder.iterPeriod=1
		for b in O.bodies:
			b.state.vel = Vector3(0.0,0.0,0.0)
			b.state.angVel = Vector3(0.0,0.0,0.0)
		calTimeStep()
		excavate()
		dispChecker.iterPeriod=1
		SRcounter = 1		
		O.step()	
		vtkRecorder.iterPeriod=2000		
		O.run(20000)		
		return

def SouthWall(x,y,z):
	Xcentre1 = 0.0 	
	Ycentre1 = 0.0
	Zcentre1 = -0.0
	dip = 90.0
	dipdir = 315.0
	dipdirN = 0.0
	dipN = 90.0-dip
	if dipdir > 180.0:
		dipdirN = dipdir - 180.0
	else:
		dipdirN = dipdir + 180.0
	dipRad = dipN/180.0*pi
	dipdirRad = dipdirN/180.0*pi
	a = cos(dipdirRad)*cos(dipRad)
	b = sin(dipdirRad)*cos(dipRad)
	c = sin(dipRad)
	l = sqrt(a*a + b*b +c*c)
	a = a/l
	b = b/l
	c = c/l
	d = a*Xcentre1 + b*Ycentre1 + c*Zcentre1
	plane = a*x+ b*y +c*z - d
	if plane < 0.0:
		plane = 0.0
	return plane

def NorthWall(x,y,z):
	Xcentre1 = 0.0 	
	Ycentre1 = 0.0
	Zcentre1 = -0.0
	dip = 90.0
	dipdir = 315.0
	dipdirN = 0.0
	dipN = 90.0-dip
	if dipdir > 180.0:
		dipdirN = dipdir - 180.0
	else:
		dipdirN = dipdir + 180.0
	dipRad = dipN/180.0*pi
	dipdirRad = dipdirN/180.0*pi
	a = cos(dipdirRad)*cos(dipRad)
	b = sin(dipdirRad)*cos(dipRad)
	c = sin(dipRad)
	l = sqrt(a*a + b*b +c*c)
	a = -a/l
	b = -b/l
	c = -c/l
	d = a*Xcentre1 + b*Ycentre1 + c*Zcentre1
	plane = a*x+ b*y +c*z - d
	if plane < 0.0:
		plane = 0.0
	return plane



def removeDamping():
	for i in O.interactions:
		i.phys.viscousDamping = 0.5
		i.phys.cumulative_us = 0.0
	O.engines[2].physDispatcher.functors[0].viscousDamping = 0.5
	integration.damping= 0.0


def myAddPlotData():
	global westBodyId
	global midBodyId
	global eastBodyId
	global originalPositionW
	global originalPositionM
	global originalPositionE
	global boundaryFriction
	KE = utils.kineticEnergy()
	uf = utils.unbalancedForce()
	displacementWx = O.bodies[westBodyId].state.pos[0] - originalPositionW[0]
	displacementW = (O.bodies[westBodyId].state.pos - originalPositionW).norm()
	displacementMx = O.bodies[midBodyId].state.pos[0] - originalPositionM[0]
	displacementM = (O.bodies[midBodyId].state.pos - originalPositionM).norm()
	displacementEx = O.bodies[eastBodyId].state.pos[0] - originalPositionE[0]
	displacementE = (O.bodies[eastBodyId].state.pos - originalPositionE).norm()
	plot.addData(timeStep=O.iter,timeStep1=O.iter,timeStep2=O.iter,timeStep3=O.iter,timeStep4=O.iter,timeStep5=O.iter,kineticEn=KE,unbalancedForce=uf,waterLevel=waterHeight,boundary_phi=boundaryFriction,displacement=displacementW,displacementWest=displacementW,dispWx=displacementWx,displacementMid=displacementM,dispMx=displacementMx,displacementEast=displacementE,dispEx=displacementEx)
	

plot.plots={'timeStep2':('kineticEn'),'timeStep3':(('displacementWest','ro-'),('dispWx','go-')),'timeStep1':(('displacementMid','ro-'),('dispMx','go-')),'timeStep5':(('displacementEast','ro-'),('dispEx','go-')),'timeStep4':('unbalancedForce')} #PLEASE CHECK





O.step()
O.step()
#excavate()
O.step()
calTimeStep()
O.run(20000)

