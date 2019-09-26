using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Jacobians;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Geometry.Coordinates;


//TODO: In XFEM I used dictionaries with nodes as keys, but it is less efficient and offers little extra safety, since the
//      shape functions and derivatives will be used by classes that have direct access to the nodes.
namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Stores the shape functions, 1st order derivatives with respect to the global cartesian coordinates and the Jacobian
	/// of an interpolation, evaluated at a certain natural point of a finite element. These quantities are needed in many 
	/// places, thus passing an instance of this class is less verbose and error prone.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class EvalInterpolation2D
	{
		private readonly IReadOnlyList<Node> elementNodes;

		public EvalInterpolation2D(IReadOnlyList<Node> elementNodes, double[] shapeFunctions, Matrix shapeGradientsNatural,
			IsoparametricJacobian2D jacobian)
		{
			int numNodes = elementNodes.Count;
#if DEBUG
			if ((shapeFunctions.Length != numNodes) || (shapeGradientsNatural.NumRows != numNodes))
			{
				throw new ArgumentException($"There are {numNodes} nodes, but {ShapeFunctions.Length} shape functions"
					+ $" and {shapeGradientsNatural.NumRows} natural shape derivatives.");
			}
#endif
			this.elementNodes = elementNodes;
			this.ShapeFunctions = shapeFunctions;
			this.ShapeGradientsNatural = shapeGradientsNatural;
			this.Jacobian = jacobian;
			this.ShapeGradientsCartesian = jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural);
		}

		/// <summary>
		/// The inverse Jacobian matrix of the interpolation and its determinant.
		/// </summary>
		public IsoparametricJacobian2D Jacobian { get; }

		/// <summary>
		/// A vector that contains the shape functions in the same order as the nodes of the interpolation.
		/// </summary>
		public double[] ShapeFunctions { get; }

		/// <summary>
		/// A matrix that contains the 1st order shape function derivatives with respect to the global cartesian coordinate 
		/// system at the integration points defined by a given quadrature. Each row corresponds to the gradient of a single 
		/// shape function. Each column corresponds to the derivatives of all shape functions with respect to a single 
		/// coordinate.
		/// </summary>
		public Matrix ShapeGradientsCartesian { get; }

		/// <summary>
		/// A matrix that contains the 1st order shape function derivatives with respect to the natural coordinate 
		/// system at the integration points defined by a given quadrature. Each row corresponds to the gradient of a single 
		/// shape function. Each column corresponds to the derivatives of all shape functions with respect to a single 
		/// coordinate.
		/// </summary>
		public Matrix ShapeGradientsNatural { get; }

		public CartesianPoint TransformPointNaturalToGlobalCartesian()
		{

			double x = 0, y = 0;
			for (int i = 0; i < ShapeFunctions.Length; ++i)
			{
				Node node = elementNodes[i];
				double val = ShapeFunctions[i];
				x += val * node.X;
				y += val * node.Y;
			}
			return new CartesianPoint(x, y);
		}
	}
}
