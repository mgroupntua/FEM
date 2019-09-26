using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Mesh;
using MGroup.MSolve.Geometry.Coordinates;

//TODO: perhaps I should return lists instead of Dictionaries with Gauss points as keys. It would be faster. The order of GPs is
//      defined by the interpolation anyway.
namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Interpolation (or shape or basis) functions used by isoparametric finite elements. Instances are able to compute the 
	/// values of the shape functions, their derivatives, the Jacobian of the matrix, etc.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public interface IIsoparametricInterpolation2D
	{
		/// <summary>
		/// The shape of a cell. Useful for interacting with other modules and software.
		/// </summary>
		CellType CellType { get; }

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
		/// nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// The number of shape functions that define this interpolation.
		/// </summary>
		int NumFunctions { get; }

		/// <summary>
		/// Checks the number and possibly the order of an element's nodes.
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		void CheckElementNodes(IReadOnlyList<Node> nodes);

		/// <summary>
		/// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		IInverseInterpolation2D CreateInverseMappingFor(IReadOnlyList<Node> nodes);

		/// <summary>
		/// Evaluate the shape functions, the 1st order derivatives with respect to the global cartesian coordinate system and 
		/// the jacobian at a given natural point.
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
		///     element.</param>
		EvalInterpolation2D EvaluateAllAt(IReadOnlyList<Node> nodes, NaturalPoint naturalPoint);

		/// <summary>
		/// Evaluate the shape functions, the 1st order derivatives with respect to the global cartesian coordinate system and
		/// the jacobian at the integration points defined by a given quadrature. 
		/// This method caches all possible quantities from previous calls. Use it instead of 
		/// <see cref="EvaluateAllAt(IReadOnlyList{Node}, NaturalPoint)"/> when the integration points of an element are
		/// the same across multiple elements or multiple iterations of a non linear or dynamic analysis (in the latter cases
		/// we could also cache the <see cref="EvalInterpolation2D"/> at the integration points of each element). 
		/// The <see cref="EvalInterpolation2D"/>s per integration point are returned in the same order as the integration 
		/// points in <paramref name="quadrature"/>.<see cref="IQuadrature2D.IntegrationPoints"/>.
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <param name="quadrature">The integration rule that defines integration points where shape functions and derivatives  
		///     are calculated. The integration points of this instance of <see cref="IQuadrature2D"/> are always the 
		///     same.</param>
		IReadOnlyList<EvalInterpolation2D> EvaluateAllAtGaussPoints(IReadOnlyList<Node> nodes, IQuadrature2D quadrature);

		/// <summary>
		/// Evaluate the shape functions at a given natural point and returns them in a vector in the same order as the nodes
		/// of the interpolation.
		/// </summary>
		/// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
		///     element.</param>
		double[] EvaluateFunctionsAt(NaturalPoint naturalPoint);

		/// <summary>
		/// Evaluate the shape functions at the integration points defined by a given quadrature. 
		/// This method caches all possible quantities from previous calls. Use it instead of 
		/// <see cref="EvaluateFunctionsAt(NaturalPoint)"/> when the integration points of an element are the same across 
		/// multiple elements or multiple iterations of a non linear or dynamic analysis.
		/// The shape functions vectors per integration point are returned in the same order as the integration 
		/// points in <paramref name="quadrature"/>.<see cref="IQuadrature2D.IntegrationPoints"/>. Each vector contains the
		/// shape functions in the same order as the nodes of the interpolation.
		/// </summary>
		/// <param name="quadrature">The integration rule that defines integration points where shape functions are calculated. 
		///     The integration points of this instance of <see cref="IQuadrature2D"/> are always the same.</param>
		IReadOnlyList<double[]> EvaluateFunctionsAtGaussPoints(IQuadrature2D quadrature);

		/// <summary>
		/// Evaluate the 1st order shape function derivatives with respect to the natural coordinate system 
		/// at a given natural point. 
		/// Each row of the returned matrix corresponds to the gradient of a single shape function. Each column corresponds to 
		/// the derivatives of all shape functions with respect to a single coordinate.
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
		///     element.</param>
		Matrix EvaluateNaturalGradientsAt(NaturalPoint naturalPoint);

		/// <summary>
		/// Evaluate the 1st order shape function derivatives with respect to the natural coordinate system 
		/// at the integration points defined by a given quadrature. 
		/// This method caches all possible quantities from previous calls. Use it instead of 
		/// <see cref="EvaluateNaturalGradientsAt(NaturalPoint)"/>, when the integration points of an element 
		/// are the same across multiple elements or multiple iterations of a non linear or dynamic analysis (in the latter 
		/// cases we could also cache the cartesian gradients at the integration points of each element). 
		/// The shape gradients matrices per integration point are returned in the same order as the integration 
		/// points in <paramref name="quadrature"/>.<see cref="IQuadrature2D.IntegrationPoints"/>. Each row of a matrix  
		/// corresponds to the gradient of a single shape function. Each column corresponds to the derivatives of all shape
		/// functions with respect to a single coordinate.
		/// </summary>
		/// <param name="quadrature">The integration rule that defines integration points where shape function derivatives are 
		///     calculated. The integration points of this instance of <see cref="IQuadrature2D"/> are always the same.</param>
		IReadOnlyList<Matrix> EvaluateNaturalGradientsAtGaussPoints(IQuadrature2D quadrature);

		/// <summary>
		/// Transforms the coordinates from the natural (element local) coordinate system to the the global
		/// coordinate system of a point that is internal to the finite element.
		/// </summary>
		/// <param name="nodes">The coordinates of the finite element's nodes in the global cartesian system.</param>
		/// <param name="naturalPoint">The coordinates in the natural system of a point that is internal to the finite 
		///     element.</param>
		CartesianPoint TransformNaturalToCartesian(IReadOnlyList<Node> nodes, NaturalPoint naturalPoint);
	}
}