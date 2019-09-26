using System;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation.GaussPointExtrapolation
{
	/// <summary>
	/// Calculates extrapolations of scalar , vector and tensor fields from the integration points of 2-by-2-by-2 Gauss-Legendre
	/// quadrature. this can be done at any point , but utility methods for directly outputting the extrapolated fields at the
	/// nodes of finite elements are also provided.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class ExtrapolationGaussLegendre2x2x2 : GaussPointExtrapolation3DBase
	{
		private static readonly double sqrt3 = Math.Sqrt(3.0);
		private static readonly ExtrapolationGaussLegendre2x2x2 uniqueInstance = new ExtrapolationGaussLegendre2x2x2();

		private ExtrapolationGaussLegendre2x2x2() : base(GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2)) { }

		/// <summary>
		/// Get the unique <see cref="ExtrapolationGaussLegendre2x2x2"/> object for the whole program. Thread safe.
		/// </summary>
		public static ExtrapolationGaussLegendre2x2x2 UniqueInstance => uniqueInstance;

		protected override double[] EvaluateExtrapolationFunctionsAt(NaturalPoint naturalPoint)
		{
			// Coordinates of the point in the auxiliary coordinate system of an imaginary "Gauss element" that has the Gauss 
			// points as its nodes.
			double ksi = sqrt3 * naturalPoint.Xi;
			double heta = sqrt3 * naturalPoint.Eta;
			double zeta = sqrt3 * naturalPoint.Zeta;

			// Shape functions of the imaginary "Gauss element". 
			// Each shape function corresponds to an integration point of Gauss-Legendre 2x2x2. Therefore their order is the same
			// as the one defined by GaussLegendre2D.Order2x2, namely Xi major/Eta minor, instead of the usual Quad4 order.
			var values = new double[8];                                 // The usual Hexa8 shape function would be:
			values[0] = 0.125 * (1 - ksi) * (1 - heta) * (1 - zeta);    // N0
			values[1] = 0.125 * (1 + ksi) * (1 - heta) * (1 - zeta);    // N1
			values[2] = 0.125 * (1 - ksi) * (1 + heta) * (1 - zeta);    // N3
			values[3] = 0.125 * (1 + ksi) * (1 + heta) * (1 - zeta);    // N2

			values[4] = 0.125 * (1 - ksi) * (1 - heta) * (1 + zeta);    // N4
			values[5] = 0.125 * (1 + ksi) * (1 - heta) * (1 + zeta);    // N5
			values[6] = 0.125 * (1 - ksi) * (1 + heta) * (1 + zeta);    // N7
			values[7] = 0.125 * (1 + ksi) * (1 + heta) * (1 + zeta);    // N6
			return values;
		}
	}
}
