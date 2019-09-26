using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Isoparamteric interpolation of a wedge finite element with 15 nodes. Quadratic shape functions.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InterpolationWedge15 : IsoparametricInterpolation3DBase
	{
		private static readonly InterpolationWedge15 uniqueInstance = new InterpolationWedge15();

		private InterpolationWedge15() : base(15)
		{
			NodalNaturalCoordinates = new NaturalPoint[]
			{
				new NaturalPoint(-1, 1, 0),
				new NaturalPoint(-1, 0, 1),
				new NaturalPoint(-1, 0, 0),
				new NaturalPoint(1, 1, 0),
				new NaturalPoint(1, 0, 1),
				new NaturalPoint(1, 0, 0),

				new NaturalPoint(-1, 0.5, 0.5),
				new NaturalPoint(-1, 0.5, 0),
				new NaturalPoint(0, 1, 0),
				new NaturalPoint(-1, 0, 0.5),

				new NaturalPoint(0, 0, 1),
				new NaturalPoint(0, 0, 0),
				new NaturalPoint(1, 0.5, 0.5),
				new NaturalPoint(1, 0.5, 0),
				new NaturalPoint(1, 0, 0.5)
			};
		}

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order
		/// of these nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique instance <see cref="InterpolationWedge15"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationWedge15 UniqueInstance => uniqueInstance;

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation2D.CheckElementNodes(IReadOnlyList{Node})"/>
		/// </summary>
		public override void CheckElementNodes(IReadOnlyList<Node> nodes)
		{
			if (nodes.Count != 15) throw new ArgumentException(
				$"A Wedge15 finite element has 15 nodes, but {nodes.Count} nodes were provided.");
			// TODO: Also check the order of the nodes too and perhaps even the shape
		}

		/// <summary>
		/// The reverse mapping for this interpolation, namely from global cartesian coordinates to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> node) =>
			throw new NotImplementedException("Iterative procedure needed");

		// Evaluated according to https://www.code-aster.org/V2/doc/v11/en/man_r/r3/r3.01.01.pdf
		protected sealed override double[] EvaluateAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;

			var values = new double[15];
			values[0] = y * (1 - x) * (2 * y - 2 - x) / 2;
			values[1] = z * (1 - x) * (2 * z - 2 - x) / 2;
			values[2] = (x - 1) * (1 - y - z) * (x + 2 * y + 2 * z) / 2;
			values[3] = y * (1 + x) * (2 * y - 2 + x) / 2;
			values[4] = z * (1 + x) * (2 * z - 2 + x) / 2;
			values[5] = (-x - 1) * (1 - y - z) * (-x + 2 * y + 2 * z) / 2;
			values[6] = 2 * y * z * (1 - x);
			values[7] = 2 * y * (1 - y - z) * (1 - x);
			values[8] = y * (1 - x * x);
			values[9] = 2 * z * (1 - y - z) * (1 - x);
			values[10] = z * (1 - x * x);
			values[11] = (1 - y - z) * (1 - x * x);
			values[12] = 2 * y * z * (1 + x);
			values[13] = 2 * y * (1 - y - z) * (1 + x);
			values[14] = 2 * z * (1 - y - z) * (1 + x);

			return values;
		}

		protected sealed override Matrix EvaluateGradientsAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;

			var derivatives = Matrix.CreateZero(15, 3);

			derivatives[0, 0] = (y * (x - 2 * y + 2)) / 2 + (y * (x - 1)) / 2;
			derivatives[1, 0] = (z * (x - 2 * z + 2)) / 2 + (z * (x - 1)) / 2;
			derivatives[2, 0] = -((y + z - 1) * (x + 2 * y + 2 * z)) / 2 - ((x - 1) * (y + z - 1)) / 2;
			derivatives[3, 0] = (y * (x + 2 * y - 2)) / 2 + (y * (x + 1)) / 2;
			derivatives[4, 0] = (z * (x + 2 * z - 2)) / 2 + (z * (x + 1)) / 2;
			derivatives[5, 0] = ((y + z - 1) * (2 * y - x + 2 * z)) / 2 - ((x + 1) * (y + z - 1)) / 2;
			derivatives[6, 0] = -2 * y * z;
			derivatives[7, 0] = 2 * y * (y + z - 1);
			derivatives[8, 0] = -2 * x * y;
			derivatives[9, 0] = 2 * z * (y + z - 1);
			derivatives[10, 0] = -2 * x * z;
			derivatives[11, 0] = 2 * x * (y + z - 1);
			derivatives[12, 0] = 2 * y * z;
			derivatives[13, 0] = -2 * y * (y + z - 1);
			derivatives[14, 0] = -2 * z * (y + z - 1);

			derivatives[0, 1] = ((x - 1) * (x - 2 * y + 2)) / 2 - y * (x - 1);
			derivatives[1, 1] = 0.0;
			derivatives[2, 1] = -(x - 1) * (y + z - 1) - ((x - 1) * (x + 2 * y + 2 * z)) / 2;
			derivatives[3, 1] = y * (x + 1) + ((x + 1) * (x + 2 * y - 2)) / 2;
			derivatives[4, 1] = 0.0;
			derivatives[5, 1] = ((x + 1) * (2 * y - x + 2 * z)) / 2 + (x + 1) * (y + z - 1);
			derivatives[6, 1] = -2 * z * (x - 1);
			derivatives[7, 1] = 2 * (x - 1) * (y + z - 1) + 2 * y * (x - 1);
			derivatives[8, 1] = 1 - x * x;
			derivatives[9, 1] = 2 * z * (x - 1);
			derivatives[10, 1] = 0.0;
			derivatives[11, 1] = x * x - 1;
			derivatives[12, 1] = 2 * z * (x + 1);
			derivatives[13, 1] = -2 * (x + 1) * (y + z - 1) - 2 * y * (x + 1);
			derivatives[14, 1] = -2 * z * (x + 1);

			derivatives[0, 2] = 0.0;
			derivatives[1, 2] = ((x - 1) * (x - 2 * z + 2)) / 2 - z * (x - 1);
			derivatives[2, 2] = -(x - 1) * (y + z - 1) - ((x - 1) * (x + 2 * y + 2 * z)) / 2;
			derivatives[3, 2] = 0.0;
			derivatives[4, 2] = z * (x + 1) + ((x + 1) * (x + 2 * z - 2)) / 2;
			derivatives[5, 2] = ((x + 1) * (2 * y - x + 2 * z)) / 2 + (x + 1) * (y + z - 1);
			derivatives[6, 2] = -2 * y * (x - 1);
			derivatives[7, 2] = 2 * y * (x - 1);
			derivatives[8, 2] = 0.0;
			derivatives[9, 2] = 2 * (x - 1) * (y + z - 1) + 2 * z * (x - 1);
			derivatives[10, 2] = 1 - x * x;
			derivatives[11, 2] = x * x - 1;
			derivatives[12, 2] = 2 * y * (x + 1);
			derivatives[13, 2] = -2 * y * (x + 1);
			derivatives[14, 2] = -2 * (x + 1) * (y + z - 1) - 2 * z * (x + 1);

			return derivatives;
		}
	}
}