using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization.Integration;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Geometry.Coordinates;

// Truss nodes:
// 0 -- 1

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a quadrilateral finite element with 4 nodes. Linear shape functions.
	/// Implements Singleton pattern.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class InterpolationTruss1D
	{
		private static readonly InterpolationTruss1D uniqueInstance = new InterpolationTruss1D();
		private readonly Dictionary<IQuadrature1D, IReadOnlyList<Matrix>> cachedNaturalGradientsAtGPs;
		private readonly Dictionary<IQuadrature1D, IReadOnlyList<Vector>> cachedFunctionsAtGPs;

		private InterpolationTruss1D()
		{
			NodalNaturalCoordinates = new NaturalPoint[]
			{
				new NaturalPoint(-1.0),
				new NaturalPoint(+1.0)
			};
		}

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
		/// nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique <see cref="InterpolationTruss1D"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationTruss1D UniqueInstance => uniqueInstance;

		/// <summary>
		/// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		//public override IInverseInterpolation1D CreateInverseMappingFor(IReadOnlyList<Node> nodes)
		//    => new InverseInterpolationTruss1D(nodes);

		public double[] EvaluateAt(double xi)
		{
			var values = new double[2];
			values[0] = 0.50 * (1 - xi);
			values[1] = 0.50 * (1 + xi);
			return values;
		}

		public Matrix EvaluateGradientsAt()
		{
			var derivatives = Matrix.CreateZero(1, 2);
			derivatives[0, 0] = -0.50; // N1,ksi
			derivatives[0, 1] = +0.50; // N2,ksi
			return derivatives;
		}

		public IReadOnlyList<EvalInterpolation1D> EvaluateAllAtGaussPoints(IReadOnlyList<Node> nodes, IQuadrature1D quadrature)
		{
			// The shape functions and natural derivatives at each Gauss point are probably cached from previous calls
			IReadOnlyList<Vector> shapeFunctionsAtGPs = EvaluateFunctionsAtGaussPoints(quadrature);
			IReadOnlyList<Matrix> naturalShapeDerivativesAtGPs = EvaluateNaturalGradientsAtGaussPoints(quadrature);

			// Calculate the Jacobians and shape derivatives w.r.t. global cartesian coordinates at each Gauss point
			int numGPs = quadrature.IntegrationPoints.Count;
			var interpolationsAtGPs = new EvalInterpolation1D[numGPs];
			//for (int gp = 0; gp < numGPs; ++gp)
			//{
			//    interpolationsAtGPs[gp] = new EvalInterpolation2D(shapeFunctionsAtGPs[gp],
			//        naturalShapeDerivativesAtGPs[gp], new IsoparametricJacobian2D(nodes, naturalShapeDerivativesAtGPs[gp]));
			//}
			return interpolationsAtGPs;
		}

		private IReadOnlyList<Matrix> EvaluateNaturalGradientsAtGaussPoints(IQuadrature1D quadrature)
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
					naturalGradientsAtGPsArray[gp] = EvaluateGradientsAt();
				}
				cachedNaturalGradientsAtGPs.Add(quadrature, naturalGradientsAtGPsArray);
				return naturalGradientsAtGPsArray;
			}
		}

		public IReadOnlyList<Vector> EvaluateFunctionsAtGaussPoints(IQuadrature1D quadrature)
		{
			bool isCached = cachedFunctionsAtGPs.TryGetValue(quadrature,
				out IReadOnlyList<Vector> shapeFunctionsAtGPs);
			if (isCached) return shapeFunctionsAtGPs;
			else
			{
				int numGPs = quadrature.IntegrationPoints.Count;
				var shapeFunctionsAtGPsArray = new Vector[numGPs];
				for (int gp = 0; gp < numGPs; ++gp)
				{
					GaussPoint gaussPoint = quadrature.IntegrationPoints[gp];
					shapeFunctionsAtGPsArray[gp] = Vector.CreateFromArray(EvaluateAt(gaussPoint.Xi));
				}
				cachedFunctionsAtGPs.Add(quadrature, shapeFunctionsAtGPsArray);
				return shapeFunctionsAtGPsArray;
			}
		}

		public class EvalInterpolation1D
		{
		}
	}


}
