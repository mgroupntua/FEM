using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a wedge with 18 nodes. Quadratic shape functions.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InterpolationWedge18 : IsoparametricInterpolation3DBase
	{
		private static readonly InterpolationWedge18 uniqueInstance = new InterpolationWedge18();

		private InterpolationWedge18() : base(18)
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
				new NaturalPoint(1, 0, 0.5),
				new NaturalPoint(0,0.5,0.5),
				new NaturalPoint(0,0.5,0),
				new NaturalPoint(0,0,0.5),
			};
		}

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order
		/// of these nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique instance <see cref="InterpolationWedge18"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationWedge18 UniqueInstance => uniqueInstance;

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation2D.CheckElementNodes(IReadOnlyList{Node})"/>
		/// </summary>
		public override void CheckElementNodes(IReadOnlyList<Node> nodes)
		{
			if (nodes.Count != 18) throw new ArgumentException(
				$"A Wedge18 finite element has 18 nodes, but {nodes.Count} nodes were provided.");
			// TODO: Also check the order of the nodes too and perhaps even the shape
		}

		/// <summary>
		/// The reverse mapping for this interpolation, namely from global cartesian coordinates to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> node)
			=> throw new NotImplementedException("Iterative procedure needed");

		protected sealed override double[] EvaluateAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;

			var values = new double[18];

			values[0] = x * y * (x - 1) * (2 * y - 1) / 2;
			values[1] = x * z * (x - 1) * (2 * z - 1) / 2;
			values[2] = x * (x - 1) * (z + y - 1) * (2 * z + 2 * y - 1) / 2;
			values[3] = x * y * (x + 1) * (2 * y - 1) / 2;
			values[4] = x * z * (x + 1) * (2 * z - 1) / 2;
			values[5] = x * (x + 1) * (z + y - 1) * (2 * z + 2 * y - 1) / 2;
			values[6] = 2 * x * y * z * (x - 1);
			values[7] = -2 * x * y * (x - 1) * (z + y - 1);
			values[8] = y * (1 - x * x) * (2 * y - 1);
			values[9] = -2 * x * z * (x - 1) * (z + y - 1);
			values[10] = z * (1 - x * x) * (2 * z - 1);
			values[11] = (1 - x * x) * (z + y - 1) * (2 * z + 2 * y - 1);
			values[12] = 2 * x * y * z * (x + 1);
			values[13] = -2 * x * y * (x + 1) * (z + y - 1);
			values[14] = -2 * x * z * (x + 1) * (z + y - 1);
			values[15] = 4 * y * z * (1 - x * x);
			values[16] = 4 * y * (x * x - 1) * (z + y - 1);
			values[17] = 4 * z * (x * x - 1) * (z + y - 1);

			return values;
		}

		protected sealed override Matrix EvaluateGradientsAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;

			var derivatives = Matrix.CreateZero(18, 3);

			derivatives[0, 0] = (x * y * (2 * y - 1)) / 2 + (y * (2 * y - 1) * (x - 1)) / 2;
			derivatives[1, 0] = (x * z * (2 * z - 1)) / 2 + (z * (2 * z - 1) * (x - 1)) / 2;
			derivatives[2, 0] = (x * (2 * y + 2 * z - 1) * (y + z - 1)) / 2 + ((x - 1) * (2 * y + 2 * z - 1) * (y + z - 1)) / 2;
			derivatives[3, 0] = (x * y * (2 * y - 1)) / 2 + (y * (2 * y - 1) * (x + 1)) / 2;
			derivatives[4, 0] = (x * z * (2 * z - 1)) / 2 + (z * (2 * z - 1) * (x + 1)) / 2;
			derivatives[5, 0] = (x * (2 * y + 2 * z - 1) * (y + z - 1)) / 2 + ((x + 1) * (2 * y + 2 * z - 1) * (y + z - 1)) / 2;
			derivatives[6, 0] = 2 * y * z * (x - 1) + 2 * x * y * z;
			derivatives[7, 0] = -2 * y * (x - 1) * (y + z - 1) - 2 * x * y * (y + z - 1);
			derivatives[8, 0] = -2 * x * y * (2 * y - 1);
			derivatives[9, 0] = -2 * z * (x - 1) * (y + z - 1) - 2 * x * z * (y + z - 1);
			derivatives[10, 0] = -2 * x * z * (2 * z - 1);
			derivatives[11, 0] = -2 * x * (2 * y + 2 * z - 1) * (y + z - 1);
			derivatives[12, 0] = 2 * y * z * (x + 1) + 2 * x * y * z;
			derivatives[13, 0] = -2 * y * (x + 1) * (y + z - 1) - 2 * x * y * (y + z - 1);
			derivatives[14, 0] = -2 * z * (x + 1) * (y + z - 1) - 2 * x * z * (y + z - 1);
			derivatives[15, 0] = -8 * x * y * z;
			derivatives[16, 0] = 8 * x * y * (y + z - 1);
			derivatives[17, 0] = 8 * x * z * (y + z - 1);

			derivatives[0, 1] = x * y * (x - 1) + (x * (2 * y - 1) * (x - 1)) / 2;
			derivatives[1, 1] = 0.0;
			derivatives[2, 1] = x * (x - 1) * (y + z - 1) + (x * (x - 1) * (2 * y + 2 * z - 1)) / 2;
			derivatives[3, 1] = x * y * (x + 1) + (x * (2 * y - 1) * (x + 1)) / 2;
			derivatives[4, 1] = 0.0;
			derivatives[5, 1] = x * (x + 1) * (y + z - 1) + (x * (x + 1) * (2 * y + 2 * z - 1)) / 2;
			derivatives[6, 1] = 2 * x * z * (x - 1);
			derivatives[7, 1] = -2 * x * (x - 1) * (y + z - 1) - 2 * x * y * (x - 1);
			derivatives[8, 1] = -2 * y * (x * x - 1) - (x * x - 1) * (2 * y - 1);
			derivatives[9, 1] = -2 * x * z * (x - 1);
			derivatives[10, 1] = 0.0;
			derivatives[11, 1] = -(x * x - 1) * (2 * y + 2 * z - 1) - 2 * (x * x - 1) * (y + z - 1);
			derivatives[12, 1] = 2 * x * z * (x + 1);
			derivatives[13, 1] = -2 * x * (x + 1) * (y + z - 1) - 2 * x * y * (x + 1);
			derivatives[14, 1] = -2 * x * z * (x + 1);
			derivatives[15, 1] = -4 * z * (x * x - 1);
			derivatives[16, 1] = 4 * y * (x * x - 1) + 4 * (x * x - 1) * (y + z - 1);
			derivatives[17, 1] = 4 * z * (x * x - 1);

			derivatives[0, 2] = 0.0;
			derivatives[1, 2] = x * z * (x - 1) + (x * (2 * z - 1) * (x - 1)) / 2;
			derivatives[2, 2] = x * (x - 1) * (y + z - 1) + (x * (x - 1) * (2 * y + 2 * z - 1)) / 2;
			derivatives[3, 2] = 0.0;
			derivatives[4, 2] = x * z * (x + 1) + (x * (2 * z - 1) * (x + 1)) / 2;
			derivatives[5, 2] = x * (x + 1) * (y + z - 1) + (x * (x + 1) * (2 * y + 2 * z - 1)) / 2;
			derivatives[6, 2] = 2 * x * y * (x - 1);
			derivatives[7, 2] = -2 * x * y * (x - 1);
			derivatives[8, 2] = 0.0;
			derivatives[9, 2] = -2 * x * (x - 1) * (y + z - 1) - 2 * x * z * (x - 1);
			derivatives[10, 2] = -2 * z * (x * x - 1) - (x * x - 1) * (2 * z - 1);
			derivatives[11, 2] = -(x * x - 1) * (2 * y + 2 * z - 1) - 2 * (x * x - 1) * (y + z - 1);
			derivatives[12, 2] = 2 * x * y * (x + 1);
			derivatives[13, 2] = -2 * x * y * (x + 1);
			derivatives[14, 2] = -2 * x * (x + 1) * (y + z - 1) - 2 * x * z * (x + 1);
			derivatives[15, 2] = -4 * y * (x * x - 1);
			derivatives[16, 2] = 4 * y * (x * x - 1);
			derivatives[17, 2] = 4 * z * (x * x - 1) + 4 * (x * x - 1) * (y + z - 1);

			return derivatives;
		}
	}
}
