


// This extrapolation is used for Tri6 interpolation with 3 Gauss points, where the node and GP numbering is:
//
// eta
// ^
// | s
// | ^
//   |
// 1
// | \
// | 1 \
// |     \
// 4       3
// |         \
// | 2      0  \    --> r
// |             \
// 2 ---- 5 ----- 0   --> xi

using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation.GaussPointExtrapolation
{
	/// <summary>
	/// Calculates extrapolations of scalar, vector and tensor fields from the integration points of 
	/// <see cref="TriangleSymmetricGaussianQuadrature.Order2Points3"/>. This can be done at any point, but utility methods for  
	/// directly outputting the extrapolated fields at the nodes of finite elements are also provided.
	/// Implements Singleton pattern.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class ExtrapolationGaussTriangular3Points : GaussPointExtrapolation2DBase
	{
		private const double oneOverThree = 1.0 / 3.0;
		private static readonly ExtrapolationGaussTriangular3Points uniqueInstance = new ExtrapolationGaussTriangular3Points();

		private ExtrapolationGaussTriangular3Points() : base(TriangleQuadratureSymmetricGaussian.Order2Points3)
		{ }

		/// <summary>
		/// Get the unique <see cref="ExtrapolationGaussTriangular3Points"/> object for the whole program. Thread safe.
		/// </summary>
		public static ExtrapolationGaussTriangular3Points UniqueInstance => uniqueInstance;

		protected override double[] EvaluateExtrapolationFunctionsAt(NaturalPoint point)
		{
			// Coordinates of the point in the auxiliary coordinate system of an imaginary "Gauss element" that has the Gauss 
			// points as its nodes.
			double r = 2.0 * point.Xi - oneOverThree;
			double s = 2.0 * point.Eta - oneOverThree;

			// Shape functions of the imaginary "Gauss element".
			// Each shape function corresponds to an integration point of GaussQuadratureForTrianglesSymmetric.Order2Points3.  
			// Therefore their order must be the same: point on Xi, point on Eta, point at right angle. This might differ from 
			// the node order in InterpolationTri3.
			var shapeFunctions = new double[3];
			shapeFunctions[0] = r;
			shapeFunctions[1] = s;
			shapeFunctions[2] = 1 - r - s;
			return shapeFunctions;
		}
	}
}
