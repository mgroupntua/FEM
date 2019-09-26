using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a hexahedral finite element with 27 nodes. Quadratic shape functions.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InterpolationHexa27 : IsoparametricInterpolation3DBase
	{
		private static readonly InterpolationHexa27 uniqueInstance = new InterpolationHexa27();

		private InterpolationHexa27() : base(27)
		{
			NodalNaturalCoordinates = new NaturalPoint[]
			{
				new NaturalPoint(-1, -1, -1),
				new NaturalPoint(1, -1, -1),
				new NaturalPoint(1, 1, -1),
				new NaturalPoint(-1, 1, -1),

				new NaturalPoint(-1, -1, 1),
				new NaturalPoint(1, -1, 1),
				new NaturalPoint(1, 1, 1),
				new NaturalPoint(-1, 1, 1),

				new NaturalPoint(0, -1, -1),
				new NaturalPoint(-1, 0, -1),
				new NaturalPoint(-1, -1, 0),
				new NaturalPoint(1, 0, -1),

				new NaturalPoint(1, -1, 0),
				new NaturalPoint(0, 1, -1),
				new NaturalPoint(1, 1, 0),
				new NaturalPoint(-1, 1, 0),

				new NaturalPoint(0, -1, 1),
				new NaturalPoint(-1, 0, 1),
				new NaturalPoint(1, 0, 1),
				new NaturalPoint(0, 1, 1),

				new NaturalPoint(0, 0, -1),
				new NaturalPoint(0, -1, 0),
				new NaturalPoint(-1, 0, 0),
				new NaturalPoint(1, 0, 0),

				new NaturalPoint(0, 1, 0),
				new NaturalPoint(0, 0, 1),
				new NaturalPoint(0, 0, 0),
			};
		}

		/// <summary>
		/// The coordinate of the finite element's nodes in the natural (element local) coordinate system. The order of these
		/// nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique <see cref="InterpolationHexa27"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationHexa27 UniqueInstance => uniqueInstance;

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation2D.CheckElementNodes(IReadOnlyList{Node})"/>
		/// </summary>
		public override void CheckElementNodes(IReadOnlyList<Node> nodes)
		{
			if (nodes.Count != 27) throw new ArgumentException(
				$"A Hexa27 finite element has 27 nodes, but {nodes.Count} nodes were provided.");
			// TODO: Also check the order of the nodes too and perhaps even the shape
		}

		/// <summary>
		/// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> node)
			=> throw new NotImplementedException("Iterative procedure needed");

		/// <summary>
		/// Evaluates Hexa27 shape functions according to <see cref="https://www.code-aster.org/V2/doc/v13/en/man_r/r3/r3.01.01.pdf">this</see>
		/// </summary>
		/// <param name="xi"></param>
		/// <param name="eta"></param>
		/// <param name="zeta"></param>
		/// <returns></returns>
		protected sealed override double[] EvaluateAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;
			var values = new double[27];
			values[0] = 1 / 8.0 * x * (x - 1) * y * (y - 1) * z * (z - 1);
			values[1] = 1 / 8.0 * x * (x + 1) * y * (y - 1) * z * (z - 1);
			values[2] = 1 / 8.0 * x * (x + 1) * y * (y + 1) * z * (z - 1);
			values[3] = 1 / 8.0 * x * (x - 1) * y * (y + 1) * z * (z - 1);

			values[4] = 1 / 8.0 * x * (x - 1) * y * (y - 1) * z * (z + 1);
			values[5] = 1 / 8.0 * x * (x + 1) * y * (y - 1) * z * (z + 1);
			values[6] = 1 / 8.0 * x * (x + 1) * y * (y + 1) * z * (z + 1);
			values[7] = 1 / 8.0 * x * (x - 1) * y * (y + 1) * z * (z + 1);

			values[8] = 1 / 4.0 * (1 - x * x) * y * (y - 1) * z * (z - 1);
			values[9] = 1 / 4.0 * x * (x - 1) * (1 - y * y) * z * (z - 1);
			values[10] = 1 / 4.0 * x * (x - 1) * y * (y - 1) * (1 - z * z);
			values[11] = 1 / 4.0 * x * (x + 1) * (1 - y * y) * z * (z - 1);

			values[12] = 1 / 4.0 * x * (x + 1) * y * (y - 1) * (1 - z * z);
			values[13] = 1 / 4.0 * (1 - x * x) * y * (y + 1) * z * (z - 1);
			values[14] = 1 / 4.0 * x * (x + 1) * y * (y + 1) * (1 - z * z);
			values[15] = 1 / 4.0 * x * (x - 1) * y * (y + 1) * (1 - z * z);

			values[16] = 1 / 4.0 * (1 - x * x) * y * (y - 1) * z * (z + 1);
			values[17] = 1 / 4.0 * x * (x - 1) * (1 - y * y) * z * (z + 1);
			values[18] = 1 / 4.0 * x * (x + 1) * (1 - y * y) * z * (z + 1);
			values[19] = 1 / 4.0 * (1 - x * x) * y * (y + 1) * z * (z + 1);

			values[20] = 1 / 2.0 * (1 - x * x) * (1 - y * y) * z * (z - 1);
			values[21] = 1 / 2.0 * (1 - x * x) * y * (y - 1) * (1 - z * z);
			values[22] = 1 / 2.0 * x * (x - 1) * (1 - y * y) * (1 - z * z);
			values[23] = 1 / 2.0 * x * (x + 1) * (1 - y * y) * (1 - z * z);

			values[24] = 1 / 2.0 * (1 - x * x) * y * (y + 1) * (1 - z * z);
			values[25] = 1 / 2.0 * (1 - x * x) * (1 - y * y) * z * (z + 1);
			values[26] = (1 - x * x) * (1 - y * y) * (1 - z * z);
			return values;
		}

		protected override Matrix EvaluateGradientsAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;

			var derivatives = Matrix.CreateZero(27, 3);

			derivatives[0, 0] = (x * y * z * (y - 1) * (z - 1)) / 8 + (y * z * (x - 1) * (y - 1) * (z - 1)) / 8;
			derivatives[1, 0] = (x * y * z * (y - 1) * (z - 1)) / 8 + (y * z * (x + 1) * (y - 1) * (z - 1)) / 8;
			derivatives[2, 0] = (x * y * z * (y + 1) * (z - 1)) / 8 + (y * z * (x + 1) * (y + 1) * (z - 1)) / 8;
			derivatives[3, 0] = (x * y * z * (y + 1) * (z - 1)) / 8 + (y * z * (x - 1) * (y + 1) * (z - 1)) / 8;
			derivatives[4, 0] = (x * y * z * (y - 1) * (z + 1)) / 8 + (y * z * (x - 1) * (y - 1) * (z + 1)) / 8;
			derivatives[5, 0] = (x * y * z * (y - 1) * (z + 1)) / 8 + (y * z * (x + 1) * (y - 1) * (z + 1)) / 8;
			derivatives[6, 0] = (x * y * z * (y + 1) * (z + 1)) / 8 + (y * z * (x + 1) * (y + 1) * (z + 1)) / 8;
			derivatives[7, 0] = (x * y * z * (y + 1) * (z + 1)) / 8 + (y * z * (x - 1) * (y + 1) * (z + 1)) / 8;
			derivatives[8, 0] = -(x * y * z * (y - 1) * (z - 1)) / 2;
			derivatives[9, 0] = -(x * z * (y * y - 1) * (z - 1)) / 4 - (z * (y * y - 1) * (x - 1) * (z - 1)) / 4;
			derivatives[10, 0] = -(x * y * (z * z - 1) * (y - 1)) / 4 - (y * (z * z - 1) * (x - 1) * (y - 1)) / 4;
			derivatives[11, 0] = -(x * z * (y * y - 1) * (z - 1)) / 4 - (z * (y * y - 1) * (x + 1) * (z - 1)) / 4;
			derivatives[12, 0] = -(x * y * (z * z - 1) * (y - 1)) / 4 - (y * (z * z - 1) * (x + 1) * (y - 1)) / 4;
			derivatives[13, 0] = -(x * y * z * (y + 1) * (z - 1)) / 2;
			derivatives[14, 0] = -(x * y * (z * z - 1) * (y + 1)) / 4 - (y * (z * z - 1) * (x + 1) * (y + 1)) / 4;
			derivatives[15, 0] = -(x * y * (z * z - 1) * (y + 1)) / 4 - (y * (z * z - 1) * (x - 1) * (y + 1)) / 4;
			derivatives[16, 0] = -(x * y * z * (y - 1) * (z + 1)) / 2;
			derivatives[17, 0] = -(x * z * (y * y - 1) * (z + 1)) / 4 - (z * (y * y - 1) * (x - 1) * (z + 1)) / 4;
			derivatives[18, 0] = -(x * z * (y * y - 1) * (z + 1)) / 4 - (z * (y * y - 1) * (x + 1) * (z + 1)) / 4;
			derivatives[19, 0] = -(x * y * z * (y + 1) * (z + 1)) / 2;
			derivatives[20, 0] = x * z * (y * y - 1) * (z - 1);
			derivatives[21, 0] = x * y * (z * z - 1) * (y - 1);
			derivatives[22, 0] = (x * (y * y - 1) * (z * z - 1)) / 2 + ((y * y - 1) * (z * z - 1) * (x - 1)) / 2;
			derivatives[23, 0] = (x * (y * y - 1) * (z * z - 1)) / 2 + ((y * y - 1) * (z * z - 1) * (x + 1)) / 2;
			derivatives[24, 0] = x * y * (z * z - 1) * (y + 1);
			derivatives[25, 0] = x * z * (y * y - 1) * (z + 1);
			derivatives[26, 0] = -2 * x * (y * y - 1) * (z * z - 1);

			derivatives[0, 1] = (x * y * z * (x - 1) * (z - 1)) / 8 + (x * z * (x - 1) * (y - 1) * (z - 1)) / 8;
			derivatives[1, 1] = (x * y * z * (x + 1) * (z - 1)) / 8 + (x * z * (x + 1) * (y - 1) * (z - 1)) / 8;
			derivatives[2, 1] = (x * y * z * (x + 1) * (z - 1)) / 8 + (x * z * (x + 1) * (y + 1) * (z - 1)) / 8;
			derivatives[3, 1] = (x * y * z * (x - 1) * (z - 1)) / 8 + (x * z * (x - 1) * (y + 1) * (z - 1)) / 8;
			derivatives[4, 1] = (x * y * z * (x - 1) * (z + 1)) / 8 + (x * z * (x - 1) * (y - 1) * (z + 1)) / 8;
			derivatives[5, 1] = (x * y * z * (x + 1) * (z + 1)) / 8 + (x * z * (x + 1) * (y - 1) * (z + 1)) / 8;
			derivatives[6, 1] = (x * y * z * (x + 1) * (z + 1)) / 8 + (x * z * (x + 1) * (y + 1) * (z + 1)) / 8;
			derivatives[7, 1] = (x * y * z * (x - 1) * (z + 1)) / 8 + (x * z * (x - 1) * (y + 1) * (z + 1)) / 8;
			derivatives[8, 1] = -z * (x * x - 1) / 4 * (y - 1) * (z - 1) - y * z * (x * x - 1) / 4 * (z - 1);
			derivatives[9, 1] = -(x * y * z * (x - 1) * (z - 1)) / 2;
			derivatives[10, 1] = -(x * y * (z * z - 1) * (x - 1)) / 4 - (x * (z * z - 1) * (x - 1) * (y - 1)) / 4;
			derivatives[11, 1] = -(x * y * z * (x + 1) * (z - 1)) / 2;
			derivatives[12, 1] = -(x * y * (z * z - 1) * (x + 1)) / 4 - (x * (z * z - 1) * (x + 1) * (y - 1)) / 4;
			derivatives[13, 1] = -z * (x * x - 1) / 4 * (y + 1) * (z - 1) - y * z * (x * x - 1) / 4 * (z - 1);
			derivatives[14, 1] = -(x * y * (z * z - 1) * (x + 1)) / 4 - (x * (z * z - 1) * (x + 1) * (y + 1)) / 4;
			derivatives[15, 1] = -(x * y * (z * z - 1) * (x - 1)) / 4 - (x * (z * z - 1) * (x - 1) * (y + 1)) / 4;
			derivatives[16, 1] = -z * (x * x - 1) / 4 * (y - 1) * (z + 1) - y * z * (x * x - 1) / 4 * (z + 1);
			derivatives[17, 1] = -(x * y * z * (x - 1) * (z + 1)) / 2;
			derivatives[18, 1] = -(x * y * z * (x + 1) * (z + 1)) / 2;
			derivatives[19, 1] = -z * (x * x - 1) / 4 * (y + 1) * (z + 1) - y * z * (x * x - 1) / 4 * (z + 1);
			derivatives[20, 1] = 2 * y * z * (x * x - 1) / 2 * (z - 1);
			derivatives[21, 1] = (z * z - 1) * (x * x - 1) / 2 * (y - 1) + y * (z * z - 1) * (x * x - 1) / 2;
			derivatives[22, 1] = x * y * (z * z - 1) * (x - 1);
			derivatives[23, 1] = x * y * (z * z - 1) * (x + 1);
			derivatives[24, 1] = (z * z - 1) * (x * x - 1) / 2 * (y + 1) + y * (z * z - 1) * (x * x - 1) / 2;
			derivatives[25, 1] = 2 * y * z * (x * x - 1) / 2 * (z + 1);
			derivatives[26, 1] = -2 * y * (x * x - 1) * (z * z - 1);

			derivatives[0, 2] = (x * y * z * (x - 1) * (y - 1)) / 8 + (x * y * (x - 1) * (y - 1) * (z - 1)) / 8;
			derivatives[1, 2] = (x * y * z * (x + 1) * (y - 1)) / 8 + (x * y * (x + 1) * (y - 1) * (z - 1)) / 8;
			derivatives[2, 2] = (x * y * z * (x + 1) * (y + 1)) / 8 + (x * y * (x + 1) * (y + 1) * (z - 1)) / 8;
			derivatives[3, 2] = (x * y * z * (x - 1) * (y + 1)) / 8 + (x * y * (x - 1) * (y + 1) * (z - 1)) / 8;
			derivatives[4, 2] = (x * y * z * (x - 1) * (y - 1)) / 8 + (x * y * (x - 1) * (y - 1) * (z + 1)) / 8;
			derivatives[5, 2] = (x * y * z * (x + 1) * (y - 1)) / 8 + (x * y * (x + 1) * (y - 1) * (z + 1)) / 8;
			derivatives[6, 2] = (x * y * z * (x + 1) * (y + 1)) / 8 + (x * y * (x + 1) * (y + 1) * (z + 1)) / 8;
			derivatives[7, 2] = (x * y * z * (x - 1) * (y + 1)) / 8 + (x * y * (x - 1) * (y + 1) * (z + 1)) / 8;
			derivatives[8, 2] = -y * (x * x - 1) / 4 * (y - 1) * (z - 1) - y * z * (x * x - 1) / 4 * (y - 1);
			derivatives[9, 2] = -(x * z * (y * y - 1) * (x - 1)) / 4 - (x * (y * y - 1) * (x - 1) * (z - 1)) / 4;
			derivatives[10, 2] = -(x * y * z * (x - 1) * (y - 1)) / 2;
			derivatives[11, 2] = -(x * z * (y * y - 1) * (x + 1)) / 4 - (x * (y * y - 1) * (x + 1) * (z - 1)) / 4;
			derivatives[12, 2] = -(x * y * z * (x + 1) * (y - 1)) / 2;
			derivatives[13, 2] = -y * (x * x - 1) / 4 * (y + 1) * (z - 1) - y * z * (x * x - 1) / 4 * (y + 1);
			derivatives[14, 2] = -(x * y * z * (x + 1) * (y + 1)) / 2;
			derivatives[15, 2] = -(x * y * z * (x - 1) * (y + 1)) / 2;
			derivatives[16, 2] = -y * (x * x - 1) / 4 * (y - 1) * (z + 1) - y * z * (x * x - 1) / 4 * (y - 1);
			derivatives[17, 2] = -(x * z * (y * y - 1) * (x - 1)) / 4 - (x * (y * y - 1) * (x - 1) * (z + 1)) / 4;
			derivatives[18, 2] = -(x * z * (y * y - 1) * (x + 1)) / 4 - (x * (y * y - 1) * (x + 1) * (z + 1)) / 4;
			derivatives[19, 2] = -y * (x * x - 1) / 4 * (y + 1) * (z + 1) - y * z * (x * x - 1) / 4 * (y + 1);
			derivatives[20, 2] = (y * y - 1) * (x * x - 1) / 2 * (z - 1) + z * (y * y - 1) * (x * x - 1) / 2;
			derivatives[21, 2] = 2 * y * z * (x * x - 1) / 2 * (y - 1);
			derivatives[22, 2] = x * z * (y * y - 1) * (x - 1);
			derivatives[23, 2] = x * z * (y * y - 1) * (x + 1);
			derivatives[24, 2] = (z * z - 1) * (x * x - 1) / 2 * (y + 1) + y * (z * z - 1) * (x * x - 1) / 2;
			derivatives[25, 2] = (y * y - 1) * (x * x - 1) / 2 * (z + 1) + z * (y * y - 1) * (x * x - 1) / 2;
			derivatives[26, 2] = -2 * z * (x * x - 1) * (y * y - 1);

			return derivatives;
		}
	}
}
