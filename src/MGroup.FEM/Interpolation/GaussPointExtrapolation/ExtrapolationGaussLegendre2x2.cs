using System;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Geometry.Coordinates;

//TODO: the order of the shape functions must be the same as the order of Gauss points in GaussLegendre2D.Order2x2 integration.  
//      This is extremely error prone. Find a way to avoid this necessity (e.g. GaussLegendre2D provides the shape functions) or
//      order the shape functions automatically (e.g. static fields that are ordered in a static constructor, which derives
//      the order by looking at GaussLegendre2D.Order2x2. When this is done, remove the warning from GaussLegendre2D
namespace MGroup.FEM.Interpolation.GaussPointExtrapolation
{
	/// <summary>
	/// Calculates extrapolations of scalar, vector and tensor fields from the integration points of 2-by-2 Gauss-Legendre 
	/// quadrature. This can be done at any point, but utility methods for directly outputting the extrapolated fields at the
	/// nodes of finite elements are also provided.
	/// Implements Singleton pattern.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class ExtrapolationGaussLegendre2x2 : GaussPointExtrapolation2DBase
	{
		private static readonly double sqrt3 = Math.Sqrt(3.0);
		private static readonly ExtrapolationGaussLegendre2x2 uniqueInstance = new ExtrapolationGaussLegendre2x2();

		private ExtrapolationGaussLegendre2x2() : base(GaussLegendre2D.GetQuadratureWithOrder(2, 2))
		{ }

		/// <summary>
		/// Get the unique <see cref="ExtrapolationGaussLegendre2x2"/> object for the whole program. Thread safe.
		/// </summary>
		public static ExtrapolationGaussLegendre2x2 UniqueInstance => uniqueInstance;

		protected override double[] EvaluateExtrapolationFunctionsAt(NaturalPoint point)
		{
			// Coordinates of the point in the auxiliary coordinate system of an imaginary "Gauss element" that has the Gauss 
			// points as its nodes.
			double r = sqrt3 * point.Xi;
			double s = sqrt3 * point.Eta;

			// Shape functions of the imaginary "Gauss element". 
			// Each shape function corresponds to an integration point of Gauss-Legendre 2x2. Therefore their order is the same
			// as the one defined by GaussLegendre2D.Order2x2, namely Xi minor/Eta major, instead of the usual Quad4 order.
			var shapeFunctions = new double[4];             // The usual Quad4 shape function would be:
			shapeFunctions[0] = 0.25 * (1 - r) * (1 - s);   // N0
			shapeFunctions[1] = 0.25 * (1 + r) * (1 - s);   // N1

			shapeFunctions[2] = 0.25 * (1 - r) * (1 + s);   // N3
			shapeFunctions[3] = 0.25 * (1 + r) * (1 + s);   // N2
			return shapeFunctions;
		}
	}
}
