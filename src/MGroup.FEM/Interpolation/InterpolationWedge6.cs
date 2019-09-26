using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a wedge finite element with 6 nodes. Linear shape functions.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InterpolationWedge6 : IsoparametricInterpolation3DBase
	{
		private static readonly InterpolationWedge6 uniqueInstance = new InterpolationWedge6();

		private InterpolationWedge6() : base(6)
		{
			NodalNaturalCoordinates = new NaturalPoint[]
			{
				new NaturalPoint(-1, 1, 0),
				new NaturalPoint(-1, 0, 1),
				new NaturalPoint(-1, 0, 0),
				new NaturalPoint(1, 1, 0),
				new NaturalPoint(1, 0, 1),
				new NaturalPoint(1, 0, 0),
			};
		}

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order
		/// of these nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique instance <see cref="InterpolationWedge6"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationWedge6 UniqueInstance => uniqueInstance;

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation2D.CheckElementNodes(IReadOnlyList{Node})"/>
		/// </summary>
		public override void CheckElementNodes(IReadOnlyList<Node> nodes)
		{
			if (nodes.Count != 6) throw new ArgumentException(
				$"A Wedge6 finite element has 6 nodes, but {nodes.Count} nodes were provided.");
			// TODO: Also check the order of the nodes too and perhaps even the shape
		}

		/// <summary>
		/// The reverse mapping for this interpolation, namely from global cartesian coordinates to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> node)
			=> throw new NotImplementedException("Implementation pending");

		protected sealed override double[] EvaluateAt(double xi, double eta, double zeta)
		{
			var values = new double[6];

			values[0] = 0.5 * eta * (1 - xi);
			values[1] = 0.5 * zeta * (1 - xi);
			values[2] = 0.5 * (1 - eta - zeta) * (1 - xi);
			values[3] = 0.5 * eta * (xi + 1);
			values[4] = 0.5 * zeta * (xi + 1);
			values[5] = 0.5 * (1 - eta - zeta) * (xi + 1);

			return values;
		}

		// Evaluation based on: https://www.code-aster.org/V2/doc/v11/en/man_r/r3/r3.01.01.pdf
		protected sealed override Matrix EvaluateGradientsAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;

			var derivatives = Matrix.CreateZero(6, 3);

			derivatives[0, 0] = -y / 2;
			derivatives[1, 0] = -z / 2;
			derivatives[2, 0] = (y + z - 1) / 2;
			derivatives[3, 0] = y / 2;
			derivatives[4, 0] = z / 2;
			derivatives[5, 0] = (1 - z - y) / 2;

			derivatives[0, 1] = (1 - x) / 2;
			derivatives[1, 1] = 0.0;
			derivatives[2, 1] = (x - 1) / 2;
			derivatives[3, 1] = (x + 1) / 2;
			derivatives[4, 1] = 0.0;
			derivatives[5, 1] = (-x - 1) / 2;

			derivatives[0, 2] = 0.0;
			derivatives[1, 2] = (1 - x) / 2;
			derivatives[2, 2] = (x - 1) / 2;
			derivatives[3, 2] = 0.0;
			derivatives[4, 2] = (x + 1) / 2;
			derivatives[5, 2] = (-x - 1) / 2;

			return derivatives;
		}
	}
}
