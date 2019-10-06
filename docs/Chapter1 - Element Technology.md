# Element technology

A summary of the basic elements og **MSolve.FEM** is provided.

## Continuum elements
This category contains structural elements with transitional olny degrees of freedom. These elements are split into two categories, two-dimensional and three-dimensioanl respectively depending on dimensionality of the model.



### Two-dimensional elements
The first  element category contains the two-dimensional elements. Below a full list of the elements can be found including their respective local node numbering, gauss integrations for stiffness and mass. The code snippet shows a sample 2D dimensional element creation. Initially, a Continuum element factory is generated which is responsible for the element definition. Its parameters are the thickness, the material and the dynamic properties which are considered constant throughout the model.
Then, the specific type of the two-dimensional element is defined by picking an element cell type and providing a list with element nodes.

```csharp
var factory = new ContinuumElement2DFactory(thickness, material, dynamicMaterial);
var quad4 = factory.CreateElement(CellType.Quad4, nodeSet);
```

#### Tri3
<p align="center">
  <img src="../docs/Images/Tri3.png" width="400"/>
</p>

For the stiffness integration of the Tri3 element an one-point integration is sufficient.

**GaussPoints for stiffness:**
| X   | Y   | Z | W   |
|-----|-----|---|-----|
| 1/3 | 1/3 | 0 | 0.5 | 
For the mass integration of the Tri3 element an three-point integration is needed.

**GaussPoints for mass:**
| X                | Y                | Z | W                      |
|------------------|------------------|---|------------------------|
| 0.66666666666667 | 0.16666666666667 | 0 | 0.5 * 0.33333333333333 |   
| 0.16666666666667 | 0.66666666666667 | 0 | 0.5 * 0.33333333333333 |   
| 0.16666666666667 | 0.16666666666667 | 0 | 0.5 * 0.33333333333333 |   

#### Tri6
<p align="center">
  <img src="../docs/Images/Tri6.png" width="400"/>
</p>

**GaussPoints for stiffness:**

For the stiffness integration of the Tri3 element an three-point integration is needed.
| X                | Y                | Z | W                      |
|------------------|------------------|---|------------------------|
| 0.66666666666667 | 0.16666666666667 | 0 | 0.5 * 0.33333333333333 |   
| 0.16666666666667 | 0.66666666666667 | 0 | 0.5 * 0.33333333333333 |   
| 0.16666666666667 | 0.16666666666667 | 0 | 0.5 * 0.33333333333333 |  

For the mass integration of the Tri3 element an three-point integration is needed.

**GaussPoints for mass:**

| X                | Y                | Z | W                      |
|------------------|------------------|---|------------------------|
| 0.44594849091597 | 0.44594849091597 | 0 | 0.5 * 0.22338158967801 |
| 0.44594849091597 | 0.10810301816807 | 0 | 0.5 * 0.22338158967801 |
| 0.10810301816807 | 0.44594849091597 | 0 | 0.5 * 0.22338158967801 |
| 0.09157621350977 | 0.09157621350977 | 0 | 0.5 * 0.10995174365532 |
| 0.09157621350977 | 0.81684757298046 | 0 | 0.5 * 0.10995174365532 |
| 0.81684757298046 | 0.09157621350977 | 0 | 0.5 * 0.10995174365532 | 

#### Quad4
<p align="center">
  <img src="../docs/Images/Quad4.png" width="400"/>
</p>
In case of quadrilateral elements the same Gauss points are used for both stiffness and mass.

**GaussPoints for stiffness and mass:**

| X                            | Y                            | Z | W |
|------------------------------|------------------------------|---|---|
| -0.5773502691896257645091488 | -0.5773502691896257645091488 | 0 | 1 |
| 0.5773502691896257645091488  | -0.5773502691896257645091488 | 0 | 1 |
| -0.5773502691896257645091488 | 0.5773502691896257645091488  | 0 | 1 |
| 0.5773502691896257645091488  | 0.5773502691896257645091488  | 0 | 1 |

#### Quad8
<p align="center">
  <img src="../docs/Images/Quad8.png" width="400"/>
</p>

**GaussPoints for stiffness and mass:**

| X                           | Y                            | Z | W             |
|-----------------------------|------------------------------|---|---------------|
| -0.774596669241483377035853 | -0.7745966692414833770358531 | 0 | 0.3086419753  |
| -0.774596669241483377035853 | 0.0                          | 0 | 0.49382716049 |
| -0.774596669241483377035853 | 0.7745966692414833770358531  | 0 | 0.3086419753  |
| 0.0                         | -0.7745966692414833770358531 | 0 | 0.3086419753  |
| 0.0                         | 0.0                          | 0 | 0.79012345679 |
| 0.0                         | 0.7745966692414833770358531  | 0 | 0.3086419753  |
| 0.7745966692414833770358531 | -0.7745966692414833770358531 | 0 | 0.3086419753  |
| 0.7745966692414833770358531 | 0.0                          | 0 | 0.49382716049 |
| 0.7745966692414833770358531 | 0.7745966692414833770358531  | 0 | 0.3086419753  |


#### Quad9
<p align="center">
  <img src="../docs/Images/Quad9.png" width="400"/>
</p>

**GaussPoints for stiffness and mass:**

| X                           | Y                            | Z | W             |
|-----------------------------|------------------------------|---|---------------|
| -0.774596669241483377035853 | -0.7745966692414833770358531 | 0 | 0.3086419753  |
| -0.774596669241483377035853 | 0.0                          | 0 | 0.49382716049 |
| -0.774596669241483377035853 | 0.7745966692414833770358531  | 0 | 0.3086419753  |
| 0.0                         | -0.7745966692414833770358531 | 0 | 0.3086419753  |
| 0.0                         | 0.0                          | 0 | 0.79012345679 |
| 0.0                         | 0.7745966692414833770358531  | 0 | 0.3086419753  |
| 0.7745966692414833770358531 | -0.7745966692414833770358531 | 0 | 0.3086419753  |
| 0.7745966692414833770358531 | 0.0                          | 0 | 0.49382716049 |
| 0.7745966692414833770358531 | 0.7745966692414833770358531  | 0 | 0.3086419753  |

### Three-dimensional elements
The first  element category contains the three-dimensional elements. Below a full list of the elements can be found including their respective local node numbering, gauss integrations for stiffness and mass. The code snippet shows a sample 3D dimensional element creation. 

```csharp
var factory = new ContinuumElement3DFactory(Material0, DynamicMaterial0);
var tet4 = factory.CreateElement(CellType.Tet4, NodeSet0);
```

#### Tet4
<p align="center">
  <img src="../docs/Images/Tet4.png" height="400"/>
</p>

#### Tet10
<p align="center">
  <img src="../docs/Images/Tet10.png" height="400"/>
</p>


#### Hexa8
<p align="center">
  <img src="../docs/Images/Hexa8.png" width="400"/>
</p>

#### Hexa20
<p align="center">
  <img src="../docs/Images/Hexa20.png" height="400"/>
</p>

#### Hexa27
<p align="center">
  <img src="../docs/Images/Hexa27.png" height="400"/>
</p>

#### Wedge6
<p align="center">
  <img src="../docs/Images/Wedge6.png" height="400"/>
</p>

#### Wedge15
<p align="center">
  <img src="../docs/Images/Wedge15.png" height="400"/>
</p>


#### Wedge18
<p align="center">
  <img src="../docs/Images/Wedge18.png" height="400"/>
</p>


#### Pyra5
<p align="center">
  <img src="../docs/Images/Pyra5.png" width="400"/>
</p>


#### Pyra13
<p align="center">
  <img src="../docs/Images/Pyra13.png" width="400"/>
</p>


#### Pyra14
<p align="center">
  <img src="../docs/Images/Pyra13.png" width="400"/>
</p>





## Structural elements




## Special elements