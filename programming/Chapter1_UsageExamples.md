## Hexa8NonLinear
```c#
public class Hexa8NonLinear : IStructuralFiniteElement, IEmbeddedHostElement
```

A Hexa8NonLinear element can be defined by use of a constructor.
A material and an appropriate quadrature scheme should be chosen. 
```c#
ElasticMaterial3D material1 = new ElasticMaterial3D()
{
    YoungModulus = E_disp/totalDuplicates,
    PoissonRatio = ni_disp,
};

Element e1 = new Element()
{
    ID = hexaID,
    ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))
};
``` 
Elemeant nodes should be given in the right order in view of the interpolation functions given in theory manual
```c#
e1.NodesDictionary.Add(node.ID,node);
``` 


## Hexa8NonLinearDefGrad
```c#
public class Hexa8NonLinearDefGrad : IStructuralFiniteElement, IEmbeddedHostElement
```
A Hexa8NonLinearDefGrad element can be defined by use of a constructor.
A material and an appropriate quadrature scheme should be chosen. 
```c#
e1 = new Element()
{
    ID = nElement + 1,
    ElementType = new Hexa8NonLinearDefGrad(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))
};
``` 
Elemeant nodes shouuld be given in the right order in view of the interpolation functions given in theory manual
```c#
e1.NodesDictionary.Add(node.ID,node);
``` 

## Shell8NonLinear

```c#
public class Shell8NonLinear : IStructuralFiniteElement
```
A Shell8NonLinear element can be defined by use of a constructor.
A material, an appropriate quadrature scheme, the direction vectors and the shell thickness should be chosen. 
```c#
e1 = new Element()
{
    ID = nElement + 1,
    ElementType = new Shell8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 2))
    {
        oVn_i = new double[][] { new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },},
        tk = new double[] { tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate },
    }
};
                
``` 
Elemeant nodes shouuld be given in the right order in view of the interpolation functions given in theory manual
```c#
for (int j = 0; j < 8; j++)
{
    e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
}
``` 

## CohesiveShell8ToHexa20

```c#
public class CohesiveShell8ToHexa20 : IStructuralFiniteElement, IEmbeddedElement
```
A CohesiveShell8ToHexa20 element can be defined by use of a constructor.
A material, an appropriate quadrature scheme, the direction vectors, the shell thickness and the side of the element should be chosen. 
```c#
e2 = new Element()
{
    ID = ElementID,
    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
    {
        oVn_i = new double[][] { new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },
                                 new double[] { 0,0,1 },},
        tk = Tk_vec,
        ShellElementSide = 1,
    }
};
``` 
Element nodes shouuld be given in the right order in view of the interpolation functions given in theory manual
```c#
for (int j = 0; j < 8; j++)
{
    e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
}
for (int j = 0; j < 8; j++)
{
    e1.NodesDictionary.Add(elementData[nElement, j + 1+8], model.NodesDictionary[elementData[nElement, j + 1+8]]);
}
``` 

Shell element side denotes which nodes are given first.

```c#
(ShellElementSide == 0) 
``` 

means that the sheel nodes are given first and then the hexa 20 nodes.
