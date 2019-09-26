using System.Collections.Generic;
using MGroup.MSolve.Discretization.Commons;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation.GaussPointExtrapolation
{
	/// Calculates extrapolations of scalar, vector and tensor fields from the integration points of symmetric Gauss quadrature
	/// for triangles with 1 Gauss point. This can be done at any point, but utility methods for directly outputting the 
	/// extrapolated fields at the nodes of finite elements are also provided. Note that since there is only 1 Gauss point,
	/// the scalar, vector and tensor fields are constant at all points and equal to their values at the Gauss point.
	/// Implements Singleton pattern.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class ExtrapolationGaussTriangular1Point : IGaussPointExtrapolation2D
	{
		private static readonly ExtrapolationGaussTriangular1Point uniqueInstance = new ExtrapolationGaussTriangular1Point();

		/// <summary>
		/// See <see cref="IGaussPointExtrapolation2D.Quadrature"/>
		/// </summary>
		public IQuadrature2D Quadrature { get { return TriangleQuadratureSymmetricGaussian.Order1Point1; } }

		/// <summary>
		/// Get the unique <see cref="ExtrapolationGaussTriangular1Point"/> object for the whole program. Thread safe.
		/// </summary>
		public static ExtrapolationGaussTriangular1Point UniqueInstance { get { return uniqueInstance; } }

		/// <summary>
		/// See <see cref="IGaussPointExtrapolation2D.ExtrapolateScalarFromGaussPoints(IReadOnlyList{double}, NaturalPoint)"/>.
		/// </summary>
		public double ExtrapolateScalarFromGaussPoints(IReadOnlyList<double> scalarsAtGaussPoints, NaturalPoint point)
			=> scalarsAtGaussPoints[0];

		/// <summary>
		/// See <see cref="IGaussPointExtrapolation2D.ExtrapolateScalarFromGaussPointsToNodes(
		/// IReadOnlyList{double}, IIsoparametricInterpolation2D)"/>
		/// </summary>
		public IReadOnlyList<double> ExtrapolateScalarFromGaussPointsToNodes(IReadOnlyList<double> scalarsAtGaussPoints,
			IIsoparametricInterpolation2D interpolation)
		{
			var nodalScalars = new double[interpolation.NumFunctions];
			for (int i = 0; i < nodalScalars.Length; ++i) nodalScalars[i] = scalarsAtGaussPoints[0];
			return nodalScalars;
		}

		/// <summary>
		/// See <see cref="IGaussPointExtrapolation2D.ExtrapolateTensorFromGaussPoints(IReadOnlyList{double[]}, NaturalPoint)"/>.
		/// </summary>
		public double[] ExtrapolateTensorFromGaussPoints(IReadOnlyList<double[]> tensorsAtGaussPoints, NaturalPoint point)
			=> tensorsAtGaussPoints[0];

		/// <summary>
		/// See <see cref="IGaussPointExtrapolation2D.ExtrapolateTensorFromGaussPoints(IReadOnlyList{Tensor2D}, NaturalPoint)"/>.
		/// </summary>
		public Tensor2D ExtrapolateTensorFromGaussPoints(IReadOnlyList<Tensor2D> tensorsAtGaussPoints, NaturalPoint point)
			=> tensorsAtGaussPoints[0];

		/// <summary>
		/// See <see cref="IGaussPointExtrapolation2D.ExtrapolateTensorFromGaussPointsToNodes(
		/// IReadOnlyList{double[]}, IIsoparametricInterpolation2D)"/>
		/// </summary>
		public IReadOnlyList<double[]> ExtrapolateTensorFromGaussPointsToNodes(IReadOnlyList<double[]> tensorsAtGaussPoints,
			IIsoparametricInterpolation2D interpolation)
		{
			var nodalTensors = new double[interpolation.NumFunctions][];
			for (int i = 0; i < nodalTensors.Length; ++i) nodalTensors[i] = tensorsAtGaussPoints[0];
			return nodalTensors;
		}

		/// <summary>
		/// See <see cref="IGaussPointExtrapolation2D.ExtrapolateTensorFromGaussPointsToNodes(
		/// IReadOnlyList{Tensor2D}, IIsoparametricInterpolation2D)"/>
		/// </summary>
		public IReadOnlyList<Tensor2D> ExtrapolateTensorFromGaussPointsToNodes(IReadOnlyList<Tensor2D> tensorsAtGaussPoints,
			IIsoparametricInterpolation2D interpolation)
		{
			var nodalTensors = new Tensor2D[interpolation.NumFunctions];
			for (int i = 0; i < nodalTensors.Length; ++i) nodalTensors[i] = tensorsAtGaussPoints[0];
			return nodalTensors;
		}

		/// <summary>
		/// See <see cref="IGaussPointExtrapolation2D.ExtrapolateVectorFromGaussPoints(IReadOnlyList{double[]}, NaturalPoint)"/>.
		/// </summary>
		public double[] ExtrapolateVectorFromGaussPoints(IReadOnlyList<double[]> vectorsAtGaussPoints, NaturalPoint point)
			=> vectorsAtGaussPoints[0];

		/// <summary>
		/// See <see cref="IGaussPointExtrapolation2D.ExtrapolateVectorFromGaussPointsToNodes(
		/// IReadOnlyList{double[]}, IIsoparametricInterpolation2D)"/>
		/// </summary>
		public IReadOnlyList<double[]> ExtrapolateVectorFromGaussPointsToNodes(IReadOnlyList<double[]> vectorsAtGaussPoints,
			IIsoparametricInterpolation2D interpolation)
		{
			var nodalVectors = new double[interpolation.NumFunctions][];
			for (int i = 0; i < nodalVectors.Length; ++i) nodalVectors[i] = vectorsAtGaussPoints[0];
			return nodalVectors;
		}
	}
}
