using System;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Geometry.Coordinates;

// Quad9 nodes:
// 3 -- 6 -- 2
// |    |    |
// 7 -- 8 -- 5
// |    |    |
// 0 -- 4 -- 1

//TODO: the order of the shape functions must be the same as the order of Gauss points in GaussLegendre2D.Order3x3 integration.  
//      This is extremely error prone. Find a way to avoid this necessity (e.g. GaussLegendre2D provides the shape functions) or
//      order the shape functions automatically (e.g. static fields that are ordered in a static constructor, which derives
//      the order by looking at GaussLegendre2D.Order3x3. When this is done, remove the warning from GaussLegendre2D
namespace MGroup.FEM.Interpolation.GaussPointExtrapolation
{
	/// <summary>
	/// Calculates extrapolations of scalar, vector and tensor fields from the integration points of 3-by-3 Gauss-Legendre 
	/// quadrature. This can be done at any point, but utility methods for directly outputting the extrapolated fields at the
	/// nodes of finite elements are also provided.
	/// Implements Singleton pattern.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class ExtrapolationGaussLegendre3x3 : GaussPointExtrapolation2DBase
	{
		private static readonly double sqrt5over3 = Math.Sqrt(5.0 / 3.0);
		private static readonly ExtrapolationGaussLegendre3x3 uniqueInstance = new ExtrapolationGaussLegendre3x3();

		private ExtrapolationGaussLegendre3x3() : base(GaussLegendre2D.GetQuadratureWithOrder(3, 3))
		{ }

		/// <summary>
		/// Get the unique <see cref="ExtrapolationGaussLegendre3x3"/> object for the whole program. Thread safe.
		/// </summary>
		public static ExtrapolationGaussLegendre3x3 UniqueInstance => uniqueInstance;

		protected override double[] EvaluateExtrapolationFunctionsAt(NaturalPoint point)
		{
			// Coordinates of the point in the auxiliary coordinate system of an imaginary "Gauss element" that has the Gauss 
			// points as its nodes.
			double r = sqrt5over3 * point.Xi;
			double s = sqrt5over3 * point.Eta;

			// Cache some values. TODO: also cache the 1-r, 1+s, etc.
			double rsOver4 = 0.25 * r * s;
			double r_2 = r * r;
			double s_2 = s * s;

			// Shape functions of the imaginary "Gauss element". 
			// Each shape function corresponds to an integration point of Gauss-Legendre 3x3. Therefore their order is the same
			// as the one defined by GaussLegendre2D.Order3x3, namely Xi minor/Eta major, instead of the usual Quad9 order.
			var shapeFunctions = new double[9];                 // The usual Quad9 shape function would be:
			shapeFunctions[0] = rsOver4 * (1 - r) * (1 - s);    // N0
			shapeFunctions[1] = -0.5 * s * (1 - r_2) * (1 - s); // N4
			shapeFunctions[2] = -rsOver4 * (1 + r) * (1 - s);   // N1

			shapeFunctions[3] = -0.5 * r * (1 - r) * (1 - s_2); // N7
			shapeFunctions[4] = (1 - r_2) * (1 - s_2);          // N8
			shapeFunctions[5] = 0.5 * r * (1 + r) * (1 - s_2);  // N5

			shapeFunctions[6] = -rsOver4 * (1 - r) * (1 + s);   // N3
			shapeFunctions[7] = 0.5 * s * (1 - r_2) * (1 + s);  // N6
			shapeFunctions[8] = rsOver4 * (1 + r) * (1 + s);    // N2
			return shapeFunctions;
		}
	}
}