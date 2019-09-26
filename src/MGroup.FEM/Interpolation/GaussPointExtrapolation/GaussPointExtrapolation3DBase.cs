using System.Collections.Generic;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation.GaussPointExtrapolation
{
	/// <summary>
	/// Basic implementation of <see cref="IGaussPointExtrapolation2D"/>.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public abstract class GaussPointExtrapolation3DBase : IGaussPointExtrapolation3D
	{
		/// <summary>
		/// Each <see cref="IIsoparametricInterpolation3D"/> is mapped to a 2D array that contains the values of the 
		/// extrapolation functions calculated at the nodes of the finite element using that interpolation. Each row corresponds
		/// to a different node. Each columns corresponds to a different extrapolation function.
		/// </summary>
		private Dictionary<IIsoparametricInterpolation3D, double[][]> cachedExtrapolationFunctionsAtNodes;

		protected GaussPointExtrapolation3DBase(IQuadrature3D quadrature)
		{
			this.cachedExtrapolationFunctionsAtNodes = new Dictionary<IIsoparametricInterpolation3D, double[][]>();
			this.Quadrature = quadrature;
		}

		/// <summary>
		/// The integration rule which defines the integration points used for extrapolating values and defining an auxiliary 
		/// coordinate system.
		/// </summary>
		public IQuadrature3D Quadrature { get; }

		/// <summary>
		/// Calculates a scalar quantity at a given point by extrapolating (or interpolating) its known values at 
		/// the integration points.
		/// </summary>
		/// <param name="scalarAtGaussPoints">Their order must be the same as the order of integration points defined by 
		///     <see cref="Quadrature"/>.</param>
		/// <param name="naturalPoint">The point where the scalar will be computed. Its coordinates are expressed in the natural
		///     (element local) system, instead of the coordinate system defined by the integration points.</param>
		/// <returns></returns>
		public double ExtrapolateScalarFromGaussPoints(IReadOnlyList<double> scalarAtGaussPoints,
			NaturalPoint naturalPoint)
		{
			double[] extrapolationFunctions = EvaluateExtrapolationFunctionsAt(naturalPoint);
			double scalar = 0;
			for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; gp++)
			{
				scalar += extrapolationFunctions[gp] * scalarAtGaussPoints[gp];
			}

			return scalar;
		}

		/// <summary>
		/// Calculates a scalar quantity at the nodes of a finite element by extrapolating its known values at the integration 
		/// points.
		/// </summary>
		/// <param name="scalarsAtGaussPoints">Their order must be the same as the order of integration points defined by 
		///     <see cref="Quadrature"/>.</param>
		/// <param name="interpolation">Defines the natural coordinates of the finite element's nodes.</param>
		/// <returns></returns>
		public IReadOnlyList<double> ExtrapolateScalarFromGaussPointsToNodes(IReadOnlyList<double> scalarsAtGaussPoints,
			IIsoparametricInterpolation3D interpolation)
		{
			double[][] nodalExtrapolationFunctions = EvaluateExtrapolationFunctionsAtNodes(interpolation);
			IReadOnlyList<NaturalPoint> nodes = interpolation.NodalNaturalCoordinates;
			var nodalScalars = new double[nodes.Count];
			for (int i = 0; i < nodes.Count; i++)
			{
				double scalar = 0;
				for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; gp++)
				{
					scalar += nodalExtrapolationFunctions[i][gp] * scalarsAtGaussPoints[gp];
				}

				nodalScalars[i] = scalar;
			}

			return nodalScalars;
		}

		/// <summary>
		/// Calculates a tensor quantity at a given point by extrapolating (or interpolating) its known values at 
		/// the integration points.
		/// </summary>
		/// <param name="tensorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
		///     <see cref="Quadrature"/>.</param>
		/// <param name="naturalPoint">The point where the tensor will be computed. Its coordinates are expressed in the natural
		///     (element local) system, instead of the coordinate system defined by the integration points.</param>
		/// <returns></returns>
		public double[] ExtrapolateTensorFromGaussPoints(IReadOnlyList<double[]> tensorsAtGaussPoints,
			NaturalPoint naturalPoint)
		{
			double[] extrapolationFunctions = EvaluateExtrapolationFunctionsAt(naturalPoint);
			var tensor = new double[6];
			for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; gp++)
			{
				tensor[0] += extrapolationFunctions[gp] * tensorsAtGaussPoints[gp][0];
				tensor[1] += extrapolationFunctions[gp] * tensorsAtGaussPoints[gp][1];
				tensor[2] += extrapolationFunctions[gp] * tensorsAtGaussPoints[gp][2];
				tensor[3] += extrapolationFunctions[gp] * tensorsAtGaussPoints[gp][3];
				tensor[4] += extrapolationFunctions[gp] * tensorsAtGaussPoints[gp][4];
				tensor[5] += extrapolationFunctions[gp] * tensorsAtGaussPoints[gp][5];
			}

			return tensor;
		}

		/// <summary>
		/// Calculates a tensor quantity at the nodes of a finite element by extrapolating its known values at the integration 
		/// points.
		/// </summary>
		/// <param name="tensorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
		///     <see cref="Quadrature"/>.</param>
		/// <param name="interpolation">Defines the natural coordinates of the finite element's nodes.</param>
		/// <returns></returns>
		public IReadOnlyList<double[]> ExtrapolateTensorFromGaussPointsToNodes(
			IReadOnlyList<double[]> tensorsAtGaussPoints, IIsoparametricInterpolation3D interpolation)
		{
			double[][] nodalExtrapolationFunctions = EvaluateExtrapolationFunctionsAtNodes(interpolation);
			IReadOnlyList<NaturalPoint> nodes = interpolation.NodalNaturalCoordinates;
			var nodalTensors = new double[nodes.Count][];
			for (int i = 0; i < nodes.Count; i++)
			{
				var tensor = new double[6];
				for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; gp++)
				{
					tensor[0] += nodalExtrapolationFunctions[i][gp] * tensorsAtGaussPoints[gp][0];
					tensor[1] += nodalExtrapolationFunctions[i][gp] * tensorsAtGaussPoints[gp][1];
					tensor[2] += nodalExtrapolationFunctions[i][gp] * tensorsAtGaussPoints[gp][2];
					tensor[3] += nodalExtrapolationFunctions[i][gp] * tensorsAtGaussPoints[gp][3];
					tensor[4] += nodalExtrapolationFunctions[i][gp] * tensorsAtGaussPoints[gp][4];
					tensor[5] += nodalExtrapolationFunctions[i][gp] * tensorsAtGaussPoints[gp][5];
				}

				nodalTensors[i] = tensor;
			}

			return nodalTensors;
		}

		/// <summary>
		/// Calculates a vector quantity at a given point by extrapolating (or interpolating) its known values at 
		/// the integration points.
		/// </summary>
		/// <param name="vectorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
		///     <see cref="Quadrature"/>.</param>
		/// <param name="naturalPoint">The point where the tensor will be computed. Its coordinates are expressed in the natural
		///     (element local) system, instead of the coordinate system defined by the integration points.</param>
		/// <returns></returns>
		public double[] ExtrapolateVectorFromGaussPoints(IReadOnlyList<double[]> vectorsAtGaussPoints,
			NaturalPoint naturalPoint)
		{
			double[] extrapolationFunctions = EvaluateExtrapolationFunctionsAt(naturalPoint);
			var vector = new double[3];
			for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; gp++)
			{
				vector[0] += extrapolationFunctions[gp] * vectorsAtGaussPoints[gp][0];
				vector[1] += extrapolationFunctions[gp] * vectorsAtGaussPoints[gp][1];
				vector[2] += extrapolationFunctions[gp] * vectorsAtGaussPoints[gp][1];
			}

			return vector;
		}

		/// <summary>
		/// Calculates a vector quantity at the nodes of a finite element by extrapolating its known values at the integration 
		/// points.
		/// </summary>
		/// <param name="vectorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
		///     <see cref="Quadrature"/>.</param>
		/// <param name="interpolation">Defines the natural coordinates of the finite element's nodes.</param>
		/// <returns></returns>
		public IReadOnlyList<double[]> ExtrapolateVectorFromGaussPointsToNodes(
			IReadOnlyList<double[]> vectorsAtGaussPoints, IIsoparametricInterpolation3D interpolation)
		{
			double[][] nodalExtrapolationFunctions = EvaluateExtrapolationFunctionsAtNodes(interpolation);
			IReadOnlyList<NaturalPoint> nodes = interpolation.NodalNaturalCoordinates;
			var nodalVectors = new double[nodes.Count][];
			for (int i = 0; i < nodes.Count; i++)
			{
				var vector = new double[2];
				for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; gp++)
				{
					vector[0] += nodalExtrapolationFunctions[i][gp] * vectorsAtGaussPoints[gp][0];
					vector[1] += nodalExtrapolationFunctions[i][gp] * vectorsAtGaussPoints[gp][1];
					vector[2] += nodalExtrapolationFunctions[i][gp] * vectorsAtGaussPoints[gp][2];
				}

				nodalVectors[i] = vector;
			}

			return nodalVectors;
		}

		/// <summary>
		/// Calculates the functions used for extrapolating quantities from the integration points to a given point, at the 
		/// given point.
		/// </summary>
		/// <param name="naturalPoint">The coordinates of the point where the extrapolation functions will be calculated, in the 
		///     natural (element local) system.</param>
		/// <returns></returns>
		protected abstract double[] EvaluateExtrapolationFunctionsAt(NaturalPoint naturalPoint);

		private double[][] EvaluateExtrapolationFunctionsAtNodes(IIsoparametricInterpolation3D interpolation)
		{
			bool isCached =
				cachedExtrapolationFunctionsAtNodes.TryGetValue(interpolation, out double[][] nodalExtrapolationFunctions);
			if (!isCached)
			{
				IReadOnlyList<NaturalPoint> nodes = interpolation.NodalNaturalCoordinates;
				nodalExtrapolationFunctions = new double[nodes.Count][];
				for (int i = 0; i < nodes.Count; i++)
				{
					nodalExtrapolationFunctions[i] = EvaluateExtrapolationFunctionsAt(nodes[i]);
				}
				cachedExtrapolationFunctionsAtNodes.Add(interpolation, nodalExtrapolationFunctions);
			}

			return nodalExtrapolationFunctions;
		}
	}
}
