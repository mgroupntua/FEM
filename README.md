![alt text](http://mgroup.ntua.gr/wp-content/uploads/2018/05/MGroup52.png "MGroup")

# FEM
MSolve library that support the spatial discretization method of Finite Elements.

## Features

- **Elements:** Both continuum and structural elements are supported. The continuum elementa can be utilized for both structural and thermal analysis. Specifically:
  * Continuum Elements 2D
    1. Tri3
    2. Tri6
    3. Quad4
    4. Quad8
    5. Quad9
   
  * Continuum Elements 3D
    1. Tet4
    2. Tet10
    3. Hexa8
    4. Hexa20
    5. Hexa27
    6. Wedge6
    7. Wedge15
    8. Wedge18
    9. Pyra5
    10. Pyra13
    11. Pyra14 
    
  * Structural Elements
    1. Beam2DCorotational
    2. Beam3DCorotationalQuaternion
    3. EulerBeam2D
    4. EulerBeam3D
    5. Rod2D
    6. Shell8 (Geometrically non-linear formulation)
    
  * Special elements
    1. CohesiveShellToHexa8
    2. ConcentratedMass3D
    3. SpringDamper3D
    
- **Loading conditions:**    
  * Nodal loads
  * MassAccelerationLoads
  

## Installation instructions
You can choose either to clone the solution or downloads it as a zip file.

### Clone solution
1. Under the repository name, click **Clone or Download** option.

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/CloneOrDownload.png "1")

2. In the popup appearing choose the **Use HTTPS** option.

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/2.png "2")

3. Use the ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/3.png "3") to copy the link provided.

4. Open Visual Studio. In Team Explorer window appearing in your screen under Local Git Repositories click the **Clone** option. If Team Explorer window is not visible you can enable in View -> Team Explorer

  ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/4.png "4")
  
5. In the text box appearing paste the link.

 ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/5.png "5")

6. Click clone and Visual Studio will automatically download and import **MGroup.FEM**


### Download as ZIP
1. Under the repository name, click **Clone or Download** option

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/CloneOrDownload.png "1")

2. Click **Download ZIP** option. **MGroup.FEM** will be downloaded as a ZIP file.

3. Extract the ZIP file to the folder of choice.

4. Double click on **MSolve.IGA.sln** file to open the code with Visual Studio
