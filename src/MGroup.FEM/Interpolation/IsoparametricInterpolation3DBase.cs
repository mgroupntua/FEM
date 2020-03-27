using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.FEM.Interpolation.Jacobians;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Integration;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Mesh;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Basic implementation of <see cref="IIsoparametricInterpolation3D"/>
	/// Authors: Dimitris Tsapetis, Serafeim Bakalakos
	/// </summary>
	public abstract class IsoparametricInterpolation3DBase : IIsoparametricInterpolation3D
	{
		private readonly Dictionary<IQuadrature3D, IReadOnlyList<double[]>> cachedFunctionsAtGPs;
		private readonly Dictionary<IQuadrature3D, IReadOnlyList<Matrix>> cachedNaturalGradientsAtGPs;

		public IsoparametricInterpolation3DBase(int numFunctions)
		{
			this.NumFunctions = numFunctions;
			this.cachedFunctionsAtGPs = new Dictionary<IQuadrature3D, IReadOnlyList<double[]>>();
			this.cachedNaturalGradientsAtGPs = new Dictionary<IQuadrature3D, IReadOnlyList<Matrix>>();
		}

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation3D.CellType"/>.
		/// </summary>
		public virtual CellType CellType { get; }

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation3D.NumFunctions"/>.
		/// </summary>
		public int NumFunctions { get; }

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation3D.NodalNaturalCoordinates"/>.
		/// </summary>
		public abstract IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation3D.CreateInverseMappingFor(IReadOnlyList{Node})"/>.
		/// </summary>
		public abstract IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> nodes);

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation3D.EvaluateAllAt(IReadOnlyList{Node}, NaturalPoint)"/>.
		/// </summary>
		public EvalInterpolation3D EvaluateAllAt(IReadOnlyList<Node> nodes, NaturalPoint naturalPoint)
		{
			double xi = naturalPoint.Xi;
			double eta = naturalPoint.Eta;
			double zeta = naturalPoint.Zeta;
			var shapeFunctions = EvaluateAt(xi, eta, zeta);
			var naturalShapeDerivatives = EvaluateGradientsAt(xi, eta, zeta);
			return new EvalInterpolation3D(nodes, shapeFunctions, naturalShapeDerivatives,
				new IsoparametricJacobian3D(nodes, naturalShapeDerivatives));
		}

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation3D.EvaluateAllAtGaussPoints(IReadOnlyList{Node}, IQuadrature3D)"/>
		/// </summary>
		public IReadOnlyList<EvalInterpolation3D> EvaluateAllAtGaussPoints(IReadOnlyList<Node> nodes, IQuadrature3D quadrature)
		{
			// The shape functions and natural derivatives at each Gauss point are probably cached from previous calls
			IReadOnlyList<double[]> shapeFunctionsAtGPs = EvaluateFunctionsAtGaussPoints(quadrature);
			IReadOnlyList<Matrix> naturalShapeDerivativesAtGPs = EvaluateNaturalGradientsAtGaussPoints(quadrature);

			// Calculate the Jacobians and shape derivatives w.r.t. global cartesian coordinates at each Gauss point
			int numGPs = quadrature.IntegrationPoints.Count;
			var interpolationsAtGPs = new EvalInterpolation3D[numGPs];
			for (int gp = 0; gp < numGPs; ++gp)
			{
				interpolationsAtGPs[gp] = new EvalInterpolation3D(nodes, shapeFunctionsAtGPs[gp],
					naturalShapeDerivativesAtGPs[gp], new IsoparametricJacobian3D(nodes, naturalShapeDerivativesAtGPs[gp]));
			}
			return interpolationsAtGPs;
		}

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation3D.EvaluateFunctionsAt(NaturalPoint)"/>.
		/// </summary>
		public double[] EvaluateFunctionsAt(NaturalPoint naturalPoint)
			=> EvaluateAt(naturalPoint.Xi, naturalPoint.Eta, naturalPoint.Zeta);

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation3D.EvaluateFunctionsAtGaussPoints(IQuadrature3D)"/>.
		/// </summary>
		public IReadOnlyList<double[]> EvaluateFunctionsAtGaussPoints(IQuadrature3D quadrature)
		{
			bool isCached = cachedFunctionsAtGPs.TryGetValue(quadrature,
				out IReadOnlyList<double[]> shapeFunctionsAtGPs);
			if (isCached) return shapeFunctionsAtGPs;
			else
			{
				int numGPs = quadrature.IntegrationPoints.Count;
				var shapeFunctionsAtGPsArray = new double[numGPs][];
				for (int gp = 0; gp < numGPs; ++gp)
				{
					GaussPoint gaussPoint = quadrature.IntegrationPoints[gp];
					shapeFunctionsAtGPsArray[gp] = EvaluateAt(gaussPoint.Xi, gaussPoint.Eta, gaussPoint.Zeta);
				}
				cachedFunctionsAtGPs.Add(quadrature, shapeFunctionsAtGPsArray);
				return shapeFunctionsAtGPsArray;
			}
		}

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation3D.EvaluateNaturalGradientsAt(NaturalPoint)".
		/// </summary>
		public Matrix EvaluateNaturalGradientsAt(NaturalPoint naturalPoint)
			=> EvaluateGradientsAt(naturalPoint.Xi, naturalPoint.Eta, naturalPoint.Zeta);

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation3D.EvaluateNaturalGradientsAtGaussPoints(IQuadrature3D)"/>.
		/// </summary>
		/// <param name="quadrature"></param>
		public IReadOnlyList<Matrix> EvaluateNaturalGradientsAtGaussPoints(IQuadrature3D quadrature)
		{
			bool isCached = cachedNaturalGradientsAtGPs.TryGetValue(quadrature,
				out IReadOnlyList<Matrix> naturalGradientsAtGPs);
			if (isCached) return naturalGradientsAtGPs;
			else
			{
				int numGPs = quadrature.IntegrationPoints.Count;
				var naturalGradientsAtGPsArray = new Matrix[numGPs];
				for (int gp = 0; gp < numGPs; ++gp)
				{
					GaussPoint gaussPoint = quadrature.IntegrationPoints[gp];
					naturalGradientsAtGPsArray[gp] = EvaluateGradientsAt(gaussPoint.Xi, gaussPoint.Eta, gaussPoint.Zeta);
				}
				cachedNaturalGradientsAtGPs.Add(quadrature, naturalGradientsAtGPsArray);
				return naturalGradientsAtGPsArray;
			}
		}

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation3D.TransformNaturalToCartesian(IReadOnlyList{Node}, NaturalPoint)"/>.
		/// </summary>
		public CartesianPoint TransformNaturalToCartesian(IReadOnlyList<Node> nodes, NaturalPoint naturalPoint)
		{
			double[] shapeFunctionValues = EvaluateAt(naturalPoint.Xi, naturalPoint.Eta, naturalPoint.Zeta);
			double x = 0, y = 0, z = 0;
			for (int i = 0; i < nodes.Count; i++)
			{
				x += shapeFunctionValues[i] * nodes[i].X;
				y += shapeFunctionValues[i] * nodes[i].Y;
				z += shapeFunctionValues[i] * nodes[i].Z;
			}
			return new CartesianPoint(x, y, z);
		}

		public abstract void CheckElementNodes(IReadOnlyList<Node> nodes);

		/// <summary>
		/// Evaluate shape function at a given point expressed in the natural coordinate system. Each entry corresponds to a
		/// different shape function.
		/// </summary>
		/// <param name="xi">The coordinate of the point along local axis Xi.</param>
		/// <param name="eta">The coordinate of the point along local axis Eta.</param>
		/// <param name="zeta">The coordinate of the point along local axis Zeta.</param>
		/// <returns></returns>
		protected abstract double[] EvaluateAt(double xi, double eta, double zeta);

		/// <summary>
		/// Evaluate derivatives of shape functions with respect to natural coordinates at a given point expressed in the 
		/// natural coordinate system. Each row corresponds to a different shape function, column 0 corresponds to derivatives
		/// with respect to Xi coordinate, column 1 corresponds to derivatives with respect to Eta coordinate, etc.
		/// </summary>
		/// <param name="xi">The coordinate of the point along local axis Xi.</param>
		/// <param name="eta">The coordinate of the point along local axis Eta.</param>
		/// <param name="zeta">The coordinate of the point along local axis Zeta.</param>
		/// <returns></returns>
		protected abstract Matrix EvaluateGradientsAt(double xi, double eta, double zeta);
	}
}
