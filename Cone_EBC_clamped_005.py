# --------------------------------------------------------------
# Python script for the stability analysis of cylindrical shells with 
# geometric imperfections
#
#
# Ronald Wagner
# DLR Braunschweig
# 05.04.2018
# --------------------------------------------------------------

#cmd_line ='abaqus cae script=SPDA_clamped_001.py'
#os.system(cmd_line)
#print(cmd_line)

#---------------------------------------------------------------

# Define the navigation path

#---------------------------------------------------------------

#import os
#myScriptPath = 'C:/Users/Harald/Desktop/Paper_SPA/Abaqus Temp'
#myErgPath = '/res/'
#myTempPath = '/tmp/'

#os.chdir(myScriptPath+myTempPath)

#---------------------------------------------------------------

# Import the required libraries

#---------------------------------------------------------------


from abaqus import *
from abaqusConstants import *
import regionToolset
import __main__
import section
import regionToolset
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import connectorBehavior
import odbAccess
from operator import add

import numpy as np


def createConeGeometry(radiusR,height,angle):
    s = myModel.ConstrainedSketch(name='__profile__', sheetSize=500.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.ConstructionLine(point1=(0.0, -250.0), point2=(0.0, 250.0))
    s.FixedConstraint(entity=g[2])
    
    s.Line(point1=(-radiusR, 0.0), point2=(-radiusR+height*tan(angle*np.pi/180), height))
    p = myModel.Part(name='Part-1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = myModel.parts['Part-1']
    p.BaseShellRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)
    s.unsetPrimaryObject()
    a = myModel.rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = myModel.parts['Part-1']
    a.Instance(name='Part-1-1', part=p, dependent=ON)
    a = myModel.rootAssembly
    a.rotate(instanceList=('Part-1-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0), angle=90.0)
    p = myModel.parts['Part-1']
    a.Instance(name='Part-1-2', part=p, dependent=ON)
    a.rotate(instanceList=('Part-1-2', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0), angle=90.0)
    a.InstanceFromBooleanMerge(name='Part-2', instances=(a.instances['Part-1-1'], 
    a.instances['Part-1-2'], ), originalInstances=SUPPRESS, domain=GEOMETRY)
    del myModel.parts['Part-1']
    myModel.parts.changeKey(fromName='Part-2', toName='Part-1')
    p = myModel.parts['Part-1']
    myDatum = p.DatumCsysByThreePoints(name='Zylinder_KOS', coordSysType=CYLINDRICAL, origin=(0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))
        
    myID = myDatum.id
    return myID
    

def createCylinderGeometry(radius,length):
    """
    
    Die Funktion erstellt die Zylinder Geometrie Set 
    Input:
    Radius and Länge des Zylinders (radius,length - Float
    Ouput:
    Nummer der Zylinderkoordinatensystems
    -    

    """
    s1 = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(radius, 0.0))
    p = myModel.Part(name='Part-1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = myModel.parts['Part-1']
    p.BaseShellExtrude(sketch=s1, depth=length)
    s1.unsetPrimaryObject()
    
    
    # Create the reference coordinate system for the cylinder
    
    p = myModel.parts['Part-1']
    myDatum = p.DatumCsysByThreePoints(name='Zylinder_KOS', coordSysType=CYLINDRICAL, origin=(0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))
    
    myID = myDatum.id
    return myID


      
#------------------------------------------------------------------------------    
def createGeometrySet_Edge(x,y,z,setname):
    """
    
    Die Funktion erstellt ein Geometry Set (Typ - Edge)
    
    Input:
    Koordinaten der "Edge" (x,y,z - Float)
    Name "Edge" (setname - String)
    Ouput:
    -    

    """
    edge = ()    
    p = myModel.parts['Part-1']
    e = p.edges
    myedge = e.findAt((x,y,z),)
    edge = edge + (e[myedge.index:myedge.index+1], )
    p.Set(edges=edge, name=setname)   

#------------------------------------------------------------------------------
    
def createMaterial(e11,e22,nu12,g12,g13,g23,materialname):
    """
    
    Die Funktion erstellt eine Materialmodell
    
    Input:
    Materialparameter (e11,e22,nu12,g12,g13,g23 - Float)
    Name "Material" (materialname - String)
    Ouput:
    -    

    """
    myMaterial = myModel.Material(name=materialname)
    myMaterial.Density(table=((myLaminateDensity, ), ))
    myMaterial.Elastic(type=LAMINA, table=((e11,e22,nu12,g12,g13,g23), ))
    return

#------------------------------------------------------------------------------
def createReferencePoint(x,y,z):
    """
    
    Die Funktion erstellt eine Referenzpunkt (Typ - reference point)
    
    Input:
    Koordinaten des Referenzpunkes (x,y,z - Float)
    Ouput:
    Name und Position des Referenzpunktes    

    """
    refpointname = a.ReferencePoint(point=(x, y, z))
    refpointposition = r.findAt((x, y, z), )
    return refpointname,refpointposition

#------------------------------------------------------------------------------

def createMesh(ele_size):
    """
    
    Die Funktion vernetzt das FE Model mit S4R Schalenelementen
    
    Input:
    Elementgröße (ele_size - Float)
    Ouput:
    -   

    """
    p = myModel.parts['Part-1']
    p.seedPart(size=ele_size, deviationFactor=0.1, minSizeFactor=0.1)
    p = myModel.parts['Part-1']
    elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=S3R, elemLibrary=STANDARD)
    p = myModel.parts['Part-1']
    f = p.faces[:]
    faces = f
    pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    p.generateMesh()
    a.regenerate()
    p = myModel.parts['Part-1']
    regions = p.sets['Set-Cone']
    p.flipNormal(regions=regions)
    
#------------------------------------------------------------------------------

def createRigidBodyInteraction(RefPoint,RefEdge,RefName):
    """
    
    Die Funktion erstelle eine "Rigid Body Interaction" zwischen 
    einen Referenzpunkt und einer Edge
    
    Input:
    Referenzpunt und Edge für die Interaction (RefPoint,RefEdge - Float)
    Name der Interaction (RefName - String)
    Ouput:
    -   

    """
    myRefPoint = (RefPoint, )
    myRegion = regionToolset.Region(referencePoints=myRefPoint)
    myIndex = RefEdge.index
    myRefRegion = regionToolset.Region(edges=e[myIndex:myIndex+1])
    myRigidbody = myModel.RigidBody(name=RefName, refPointRegion=myRegion, tieRegion=myRefRegion)
    
#------------------------------------------------------------------------------

def createStaticStep(Stepname,Pre_Stepname,numInc,Initial_inc,Max_inc,Min_inc):
    """
    
    Die Funktion erstellt einen "Static Simulation" Step
    
    Input:
    Anzahl der Inkremente, Start-Inkrement, max. Inkrement und min. Inkrement (numInc,Initial_inc,Max_inc,Min_inc - Float)
    Name des "aktuellen" und des "vorherigen" Steps (Stepname,Pre_Stepname - String)
    Ouput:
    -   

    """
    myModel.StaticStep(name=Stepname, previous=Pre_Stepname, maxNumInc=numInc, initialInc=Initial_inc, maxInc=Max_inc,minInc=Min_inc, extrapolation=NONE, nlgeom=ON)
    #myModel.steps[Stepname].setValues(maxNumInc=numInc, stabilizationMagnitude=1e-09, stabilizationMethod=DAMPING_FACTOR, continueDampingFactors=False, adaptiveDampingRatio=None, minInc=1e-15)
    


#------------------------------------------------------------------------------
def createDispBC_Edge_RP(RP,RPname,BCname,Step_name,trans1,trans2,trans3,rot1,rot2,rot3):
    """
    
    Die Funktion erstellt eine Verschiebungsrandbedingung an einem Referenzpunkt
    
    Input:
    Name des Punktes (RP - interne Variable)
    Name des Referenz-Punktes (RPname - string)
    Name der Randbedingung (BCname - string)
    Name des Steps (Step_name - string)
    Verschiebungs und Rotationsfreiheitsgrade (trans1,trans2,trans3,rot1,rot2,rot3 - Float - oder UNSET wenn Freiheitsgrad uneingeschränkt ist)
    
    Ouput:
    -   

    """
    a = myModel.rootAssembly
    p = a.instances['Part-1-1']
    r = a.referencePoints
    refPoints1=(RP, )
    a.Set(referencePoints=refPoints1, name=RPname)
    myRefPoint = (RP, )
    myRegion = regionToolset.Region(referencePoints=myRefPoint)
    myBC = myModel.DisplacementBC(name=BCname, createStepName=Step_name, region=myRegion, u1=trans1, u2=trans2, u3=trans3, ur1=rot1, ur2=rot2, ur3=rot3)

#------------------------------------------------------------------------------
def createPlane_by_Prinicpal(Type,Offset):
    """
    
    Die Funktion erstellt eine Principal-Plane
    
    Input:
    Typ der Prinicpal-Plane (XY,XZ,YZ - interne Varaible)
    Offset der Plane in Bezug auf das globale Koordinatensystem (Offset - Float)
    
    Ouput:
    Nummer der Principal Plane

    """
    p = myModel.parts['Part-1']
    myPlane = p.DatumPlaneByPrincipalPlane(principalPlane=Type, offset=Offset)
    myID = myPlane.id
    return myID

#------------------------------------------------------------------------------
def createVertice_Set(x,y,z,VerticeSet_name):
    """
    
    Die Funktion erstellt ein Vertice Set (Typ - Vertice)
    
    Input:
    Koordinaten des "Vertice" (x,y,z - Float)
    Name "Vertice" (setname - String)
    Ouput:
    -    

    """
    a = myModel.rootAssembly
    v = a.instances['Part-1-1'].vertices
    myVertice = v.findAt((x,y,z),)
    myVertice = (v[myVertice.index:myVertice.index+1],)
    a.Set(vertices=myVertice, name=VerticeSet_name)

#------------------------------------------------------------------------------

def createDispBC_Edge_Vertice(VerticeSet,BCname,Step_name,trans1,trans2,trans3,rot1,rot2,rot3):
    """
    
    Die Funktion erstellt eine Verschiebungsrandbedingung an einem Vertice
    
    Input:
    Name des Vertice (VerticeSet - interne Variable)
    Name der Randbedingung (BCname - string)
    Name des Steps (Step_name - string)
    Verschiebungs und Rotationsfreiheitsgrade (trans1,trans2,trans3,rot1,rot2,rot3 - Float - oder UNSET wenn Freiheitsgrad uneingeschränkt ist)
    
    Ouput:
    -   

    """
    a = myModel.rootAssembly
    myModel.DisplacementBC(name=BCname, createStepName=Step_name, region=region, u1=UNSET, u2=-myPerturbation, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
#------------------------------------------------------------------------------

        
            

#------------------------------------------------------------------------------
def createEdge_Set(x,y,z,EdgeSet_name):
    """
    
    Die Funktion erstellt ein Edge Set (Typ - Edge)
    
    Input:
    Koordinaten der "Edge" (x,y,z - Float)
    Name "Edge" (EdgeSet_name - String)
    Ouput:
    -    

    """
    e = p.edges
    myEdge = ()
    Edge = e.findAt((x,y,z),)
    myEdge = myEdge + (e[Edge.index:Edge.index+1], )
    p.Set(edges=myEdge, name=EdgeSet_name)    


def createFace_Partition_byPlane(ID_Plane):
    """
    
    Die Funktion erstellt eine Geometrie Partition mit einer Ebene(Plane)
    
    Input:
    Nummer der Ebene (ID_Plane - integer)
    
    Ouput:
    -   

    """
    p = myModel.parts['Part-1']
    f = p.faces[:]
    pickedFaces = f
    d = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d[ID_Plane], faces=pickedFaces)        

#------------------------------------------------------------------------------      

# Layup and number of plies


myLaminate6 = [45,-45,0,90,90,0,-45,45]
myLayerNumber_v = [len(myLaminate6)]
myShell = [myLaminate6]

# name of the numerical model

myName = ['S1']

# Material parameter

# Material parameter

myE1_v =  [210000]
myE2_v =  [210000]
myG12_v = [210000/(2*(1+0.3))]
myNu12_v = [0.3]
myG23_v = [210000/(2*(1+0.3))]

# Cylinder length and radius

myHeight_v = [600.0]                   # Height
myRadius_v = [400.0]  

# ply thickness

myLaminateThickness_v = [1.2/8]
myCore_v = [12.7]





LRSM_Factor = 1000


# element length
# can be estimated with = 0.5*np.sqrt(myRadius*myLayerNumber*myLaminateThickness)

myMesh_Size_v = [0.5*np.sqrt(myRadius_v[0]*myLayerNumber_v[0]*myLaminateThickness_v[0])]
#myMesh_Size_v = [0.91]

# axial shortening of the cylinder (simulation is displacement controlled)

my_disp_v = [0.25]


# number of cores for the simulation

myCpu = 8

# start and end of iterations for SPDA (20 increments)
    
my_START = [1]
my_END = [11]

myLoadKDF_v = [0.01]

# Limit for the outer loop (number of shells investigated)

Limit = 1
 
# outer loop for different shells with vectors for result extraction
myAngle = -2
mySemi_Vertex_Angle = 0

for jc in range(0,Limit,1):
  somefloats_f = []
  somefloats_f2 = []
  somefloats_f3 = []#
  somefloats_u = []
  pert_L_v = []
  
  # inner loop for SPDA iteration
  
  for ic in range(my_START[0],my_END[0],1):
   # Define the material variables   
    myHeight = myHeight_v[0]    # Cylinder length [mm]                                                                
    myRadius = myRadius_v[0]    # Cylinder radius [mm]                                                              
    myLaminateThickness = myLaminateThickness_v[0]   # Laminate Thickness [mm]
    myCore = myCore_v[0]
    myE1 = myE1_v[0]   # Laminate E1 [MPa]                                                                     
    myE2 = myE2_v[0] # Laminate E2 [MPa]                                                                           
    myG12 = myG12_v[0]    # Laminate G12 [MPa]                                                                    
    myNu12 = myNu12_v[0]    # Laminate nu12 [-]                                                                         
    myG13 = myG12   # Outplane shear modulus [MPa]                                                          
    myG23 = myG12    # Outplane shear modulus[MPa]                                                           
    myNu13 = myNu12   # Poissons coefficient [-]                                                            
    myNu23 = myNu13   # Poissons coefficient [-] 
    myLaminateDensity = 1.6e-9    # Laminate density [t/mm³]
    myLayerNumber = myLayerNumber_v[0]   # Number of lamina layers                                                          
    myLayerOrientation = myShell[jc]    # Orientation lamina layers 
    disp = my_disp_v[0]  # axial shortening
    Y = 280
    Mesh_Size = myMesh_Size_v[0] # mesh
    SlantLength = myHeight/(np.cos(mySemi_Vertex_Angle*np.pi/180))
    my_Small_Radius = myRadius-np.sin(mySemi_Vertex_Angle*np.pi/180)*SlantLength
    Average_Radius = (my_Small_Radius+myRadius)/(2*np.cos(mySemi_Vertex_Angle*np.pi/180))
    L_R = int(SlantLength/Average_Radius)
    R_t = int(Average_Radius/(myLayerNumber*myLaminateThickness))
    Zb = int((SlantLength**2*(np.sqrt(1-0.3**2)))/(Average_Radius*myLayerNumber*myLaminateThickness))
    An = int(mySemi_Vertex_Angle)
    myString = str(myName[jc])+'_Z_'+str(Zb)+'_angle_'+str(An)+'_EBC_Loop-'+str(ic)
    myString_result = str(myName[jc])+'_Z_'+str(Zb)+'_angle_'+str(An)
    #myString_result = str(myName[jc])+'_Z_'+str(Zb)+'_ß_'+str(An)
    myModel = mdb.Model(name=myString)
    Ncl = (2*np.pi*myE1*(myLayerNumber*myLaminateThickness)**2)/(np.sqrt(3*(1-myNu12**2)))*(np.cos(mySemi_Vertex_Angle*np.pi/180.0))**2
    Ny = 2*3.141*(myLayerNumber*myLaminateThickness)*myRadius*Y
    #-----------------------------------------------------------
    
    myPerturbation = Ny/my_END[0]*ic
    #myPerturbation = 460528+(ic*7675.4)
    # myPerturbation = 2400*ic
    
    myID = createConeGeometry(myRadius,myHeight,mySemi_Vertex_Angle)
    
    # Create Cylinder Material
         
    createMaterial(myE1, myE2, myNu12, myG12, myG13, myG23,'CFRP_Cone')
    myModel.materials['CFRP_Cone'].Plastic(table=((Y, 0.0), ))
    
    p = myModel.parts['Part-1']
    faces = p.faces
    p.Set(faces=faces, name='Set-Cone')
    
    
    
    # Create Cylinder Bottom and Top Edge for Boundary Condition and Load Application
    
    createGeometrySet_Edge(myRadius,0.0,0.0,'Set-Top_EDGE')
    createGeometrySet_Edge(my_Small_Radius,0.0,myHeight,'Set-Bottom_EDGE')
    
    #-------------------------------------------------------------------------- 
    
    # Define material orientation
    
    #-------------------------------------------------------------------------- 
    
#    def createMaterial_CompositeLayup():
    
    region = regionToolset.Region(faces=faces)
    orientation = myModel.parts['Part-1'].datums[myID]
    myModel.parts['Part-1'].MaterialOrientation(region=region, orientationType=SYSTEM, axis=AXIS_2, localCsys=orientation, fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, additionalRotationField='')
    #: Specified material orientation has been assigned to the selected regions.
    
    
    #-------------------------------------------------------------------------- 
    
    # Lagenaufbau definieren
    
    #-------------------------------------------------------------------------- 
    
    layupOrientation = myModel.parts['Part-1'].datums[myID]
    
    p = myModel.parts['Part-1']
    
    myRegion = p.sets['Set-Cone']
        
    compositeLayup = myModel.parts['Part-1'].CompositeLayup(name='CompositeLayup-1', description='', elementType=SHELL, offsetType=MIDDLE_SURFACE, symmetric=False, thicknessAssignment=FROM_SECTION)
    
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, useDensity=OFF)
    
    compositeLayup.ReferenceOrientation(orientationType=SYSTEM, localCsys=layupOrientation, fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, additionalRotationField='', axis=AXIS_2)    
    
    for j in range(0,myLayerNumber):
      compositeLayup.CompositePly(suppressed=False, plyName='myPly-'+str(j), region=myRegion, material='CFRP_Cone', thicknessType=SPECIFY_THICKNESS, thickness=myLaminateThickness, orientationType=SPECIFY_ORIENT, orientationValue=myLayerOrientation[j])
    
      
    #-------------------------------------------------------------------------- 
    
    # Assembly 
    
    #-------------------------------------------------------------------------- 
    
    a = myModel.rootAssembly
    p = myModel.parts['Part-1']
    a.Instance(name='Part-1-1', part=p, dependent=ON)
        
    #--------------------------------------------------------------------------    
    
    # Reference Punkte erstellen
    
    #-------------------------------------------------------------------------- 
    
    a = myModel.rootAssembly
    r = a.referencePoints
    e = a.instances['Part-1-1'].edges
                   
    myRefpoint1,myVertice1 = createReferencePoint(0.0,0.0,0.0)
    myRefpoint2,myVertice2 = createReferencePoint(0.0,0.0,myHeight)
    #-------------------------------------------------------------------------- 
    
    # Create the rigid body for the top of the cylinder
    
    #-------------------------------------------------------------------------- 
    myEdgeTop = e.findAt((myRadius,0,0),)
    createRigidBodyInteraction(myVertice1,myEdgeTop,'Constraint-Top')
    #-------------------------------------------------------------------------- 
    
    # Create the rigid body for the bottom of the cylinder
    
    #-------------------------------------------------------------------------- 
    myEdgeBottom = e.findAt((my_Small_Radius,0,myHeight),)
    createRigidBodyInteraction(myVertice2,myEdgeBottom,'Constraint-Bottom')
    #-------------------------------------------------------------------------- 
    
    # Create the simulation step for the SPDA 
    
    #-------------------------------------------------------------------------- 
    createStaticStep('Step-1','Initial',1000,1,1,1E-05)
    # Step 2 for compression of the cylinder
    createStaticStep('Step-2','Step-1',75,0.01,0.01,1E-05)
    #-------------------------------------------------------------------------- 
    
    # Create the bondary conditions at the top of the cylinder and at the bottom of the cylinder (clamped)
    
    #-------------------------------------------------------------------------- 
    
    createDispBC_Edge_RP(myVertice1,'RP-1','BC Top','Initial',0,0,UNSET,0,0,0)
    createDispBC_Edge_RP(myVertice2,'RP-2','BC Bottom','Initial',0,0,0,0,0,0)
    
    #-------------------------------------------------------------------------- 
    
    # Apply the loading (displacement controlled)
    
    #-------------------------------------------------------------------------- 
    
    #createDispBC_Edge_RP(myVertice1,'RP-1','DISP','Step-2',UNSET,UNSET,disp,UNSET,UNSET,UNSET)
    
    #-------------------------------------------------------------------------- 
    
    # Define Datum plane and Edge for Perturbation application
    myID_YZ = createPlane_by_Prinicpal(YZPLANE,0.0)
    myID_XY = createPlane_by_Prinicpal(XYPLANE,myHeight/2.0)
    myID_XZ = createPlane_by_Prinicpal(XZPLANE,0)
    createFace_Partition_byPlane(myID_YZ)
    createEdge_Set(0,(myRadius+my_Small_Radius)/2.0,myHeight/2.0,'TOP_EDGE')
    createFace_Partition_byPlane(myID_XY)
    createVertice_Set(0,-(myRadius+my_Small_Radius)/2.0,myHeight/2.0,'SPDA_Point')
    a = myModel.rootAssembly
    region = a.sets['SPDA_Point']
    myModel.DisplacementBC(name='BC-4', createStepName='Step-2', region=region, u1=UNSET, u2=myRadius/R_t*15.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    #createDispBC_Edge_Vertice(region,'SPDA_Point','Step-2',UNSET,myRadius/R_t*15.0,UNSET,UNSET,UNSET,UNSET)
    #myModel.boundaryConditions['Pert_BC'].deactivate('Step-2')
    createFace_Partition_byPlane(myID_XZ)
    a = myModel.rootAssembly
    region = a.sets['RP-1']
    myModel.ConcentratedForce(name='Load-3', createStepName='Step-1', region=region, cf3=myPerturbation, distributionType=UNIFORM, field='', localCsys=None)
    a = myModel.rootAssembly
    a.features['Part-2-1'].suppress()
    createMesh(Mesh_Size)    
    
    #-------------------------------------------------------------------------- 
            
    # Generate Job name and submit Job
    
    #--------------------------------------------------------------------------   
    myJobName = myString
    myJob = mdb.Job(name=myJobName, model=myJobName, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', multiprocessingMode=DEFAULT, numCpus=myCpu, numDomains=myCpu, numGPUs=0)
    myJob.submit()
    myJob.waitForCompletion()
    
    myPerturbation = [myPerturbation]
    pert_L_v = pert_L_v + myPerturbation
    
    # open ODB file and extract values for reaction forces
    odb = session.openOdb(str(myString)+'.odb')
    time_v = []
    Force_v = []
    AxialDisp_v = []
    lastStep=odb.steps['Step-2']
    
    for x in range(len(lastStep.frames)):
        lastFrame = lastStep.frames[x]
        Time1=lastFrame.frameValue
        rForce = lastFrame.fieldOutputs['RF']
        center = odb.rootAssembly.nodeSets['SPDA_POINT']
        centerRForce = rForce.getSubset(region=center)
        time_vr = [Time1]
        time_v = time_v + time_vr
        rDisp = lastFrame.fieldOutputs['U']
        center = odb.rootAssembly.nodeSets['SPDA_POINT']
        centerRDisp = rDisp.getSubset(region=center)
    
        for i in centerRDisp.values:
            AxialDisp_vr = [i.data[1]]
            AxialDisp_v = AxialDisp_v + AxialDisp_vr
    
        for i in centerRForce.values:
            Force_vr = [i.data[1]]
            Force_v = Force_v + Force_vr
    
   
  
    Force = np.max(Force_v) 
    somefloats_fr = [Force]
    somefloats_f = somefloats_f + somefloats_fr
    Index = np.argmax(Force_v)  
    Disp = AxialDisp_v[Index]
    somefloats_fr3 = [Disp]
    
    # write input to text files
    
    somefloats_f3 = somefloats_f3 + somefloats_fr3
    np.savetxt('N_'+str(myString_result)+'.txt',somefloats_f)
    np.savetxt('u_'+str(myString_result)+'.txt',somefloats_f3)
    np.savetxt('v_'+str(myString_result)+'.txt',pert_L_v)
    np.savetxt('Rf_'+str(myString)+'.txt',Force_v)
    np.savetxt('u_'+str(myString)+'.txt',AxialDisp_v)
  

    
    