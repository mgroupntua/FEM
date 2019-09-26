using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a hexahedral finite element with 8 nodes. Quadratic shape functions.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InterpolationHexa20 : IsoparametricInterpolation3DBase
	{
		private static readonly InterpolationHexa20 uniqueInstance = new InterpolationHexa20();

		private InterpolationHexa20() : base(20)
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
			};
		}

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
		/// nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique <see cref="InterpolationHexa20"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationHexa20 UniqueInstance => uniqueInstance;

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation2D.CheckElementNodes(IReadOnlyList{Node})"/>
		/// </summary>
		public override void CheckElementNodes(IReadOnlyList<Node> nodes)
		{
			if (nodes.Count != 20) throw new ArgumentException(
				$"A Hexa20 finite element has 20 nodes, but {nodes.Count} nodes were provided.");
			// TODO: Also check the order of the nodes too and perhaps even the shape
		}

		/// <summary>
		/// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> node)
			=> throw new NotImplementedException("Iterative procedure needed");

		protected sealed override double[] EvaluateAt(double xi, double eta, double zeta)
		{
			var values = new double[20];
			values[0] = 1 / 8.0 * (1 - xi) * (1 - eta) * (1 - zeta) * (-2 - xi - eta - zeta);
			values[1] = 1 / 8.0 * (1 + xi) * (1 - eta) * (1 - zeta) * (-2 + xi - eta - zeta);
			values[2] = 1 / 8.0 * (1 + xi) * (1 + eta) * (1 - zeta) * (-2 + xi + eta - zeta);
			values[3] = 1 / 8.0 * (1 - xi) * (1 + eta) * (1 - zeta) * (-2 - xi + eta - zeta);

			values[4] = 1 / 8.0 * (1 - xi) * (1 - eta) * (1 + zeta) * (-2 - xi - eta + zeta);
			values[5] = 1 / 8.0 * (1 + xi) * (1 - eta) * (1 + zeta) * (-2 + xi - eta + zeta);
			values[6] = 1 / 8.0 * (1 + xi) * (1 + eta) * (1 + zeta) * (-2 + xi + eta + zeta);
			values[7] = 1 / 8.0 * (1 - xi) * (1 + eta) * (1 + zeta) * (-2 - xi + eta + zeta);

			values[8] = 1 / 4.0 * (1 - xi * xi) * (1 - eta) * (1 - zeta);
			values[9] = 1 / 4.0 * (1 - eta * eta) * (1 - xi) * (1 - zeta);
			values[10] = 1 / 4.0 * (1 - zeta * zeta) * (1 - xi) * (1 - eta);
			values[11] = 1 / 4.0 * (1 - eta * eta) * (1 + xi) * (1 - zeta);

			values[12] = 1 / 4.0 * (1 - zeta * zeta) * (1 + xi) * (1 - eta);
			values[13] = 1 / 4.0 * (1 - xi * xi) * (1 + eta) * (1 - zeta);
			values[14] = 1 / 4.0 * (1 - zeta * zeta) * (1 + xi) * (1 + eta);
			values[15] = 1 / 4.0 * (1 - zeta * zeta) * (1 - xi) * (1 + eta);

			values[16] = 1 / 4.0 * (1 - xi * xi) * (1 - eta) * (1 + zeta);
			values[17] = 1 / 4.0 * (1 - eta * eta) * (1 - xi) * (1 + zeta);
			values[18] = 1 / 4.0 * (1 - eta * eta) * (1 + xi) * (1 + zeta);
			values[19] = 1 / 4.0 * (1 - xi * xi) * (1 + eta) * (1 + zeta);
			return values;
		}

		protected sealed override Matrix EvaluateGradientsAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;

			var derivatives = Matrix.CreateZero(20, 3);

			derivatives[0, 0] = (x - 1) / 8 * (y - 1) * (z - 1) + ((y - 1) * (z - 1) * (x + y + z + 2)) / 8;
			derivatives[1, 0] = (x + 1) / 8 * (y - 1) * (z - 1) - ((y - 1) * (z - 1) * (y - x + z + 2)) / 8;
			derivatives[2, 0] = -((y + 1) * (z - 1) * (x + y - z - 2)) / 8 - (x + 1) / 8 * (y + 1) * (z - 1);
			derivatives[3, 0] = -((y + 1) * (z - 1) * (x - y + z + 2)) / 8 - (x - 1) / 8 * (y + 1) * (z - 1);
			derivatives[4, 0] = -((y - 1) * (z + 1) * (x + y - z + 2)) / 8 - (x - 1) / 8 * (y - 1) * (z + 1);
			derivatives[5, 0] = -((y - 1) * (z + 1) * (x - y + z - 2)) / 8 - (x + 1) / 8 * (y - 1) * (z + 1);
			derivatives[6, 0] = (x + 1) / 8 * (y + 1) * (z + 1) + ((y + 1) * (z + 1) * (x + y + z - 2)) / 8;
			derivatives[7, 0] = (x - 1) / 8 * (y + 1) * (z + 1) + ((y + 1) * (z + 1) * (x - y - z + 2)) / 8;
			derivatives[8, 0] = -(x * (y - 1) * (z - 1)) / 2;
			derivatives[9, 0] = -(y * y - 1) / 4 * (z - 1);
			derivatives[10, 0] = -(z * z - 1) / 4 * (y - 1);
			derivatives[11, 0] = (y * y - 1) / 4 * (z - 1);
			derivatives[12, 0] = (z * z - 1) / 4 * (y - 1);
			derivatives[13, 0] = (x * (y + 1) * (z - 1)) / 2;
			derivatives[14, 0] = -(z * z - 1) / 4 * (y + 1);
			derivatives[15, 0] = (z * z - 1) / 4 * (y + 1);
			derivatives[16, 0] = (x * (y - 1) * (z + 1)) / 2;
			derivatives[17, 0] = (y * y - 1) / 4 * (z + 1);
			derivatives[18, 0] = -(y * y - 1) / 4 * (z + 1);
			derivatives[19, 0] = -(x * (y + 1) * (z + 1)) / 2;

			derivatives[0, 1] = (x - 1) / 8 * (z - 1) * (x + y + z + 2) + (x - 1) / 8 * (y - 1) * (z - 1);
			derivatives[1, 1] = -(x + 1) / 8 * (y - 1) * (z - 1) - (x + 1) / 8 * (z - 1) * (y - x + z + 2);
			derivatives[2, 1] = -(x + 1) / 8 * (y + 1) * (z - 1) - (x + 1) / 8 * (z - 1) * (x + y - z - 2);
			derivatives[3, 1] = (x - 1) / 8 * (y + 1) * (z - 1) - (x - 1) / 8 * (z - 1) * (x - y + z + 2);
			derivatives[4, 1] = -(x - 1) / 8 * (y - 1) * (z + 1) - (x - 1) / 8 * (z + 1) * (x + y - z + 2);
			derivatives[5, 1] = (x + 1) / 8 * (y - 1) * (z + 1) - (x + 1) / 8 * (z + 1) * (x - y + z - 2);
			derivatives[6, 1] = (x + 1) / 8 * (z + 1) * (x + y + z - 2) + (x + 1) / 8 * (y + 1) * (z + 1);
			derivatives[7, 1] = (x - 1) / 8 * (z + 1) * (x - y - z + 2) - (x - 1) / 8 * (y + 1) * (z + 1);
			derivatives[8, 1] = -(x * x - 1) / 4 * (z - 1);
			derivatives[9, 1] = -(y * (x - 1) * (z - 1)) / 2;
			derivatives[10, 1] = -(z * z - 1) / 4 * (x - 1);
			derivatives[11, 1] = (y * (x + 1) * (z - 1)) / 2;
			derivatives[12, 1] = (z * z - 1) / 4 * (x + 1);
			derivatives[13, 1] = (x * x - 1) / 4 * (z - 1);
			derivatives[14, 1] = -(z * z - 1) / 4 * (x + 1);
			derivatives[15, 1] = (z * z - 1) / 4 * (x - 1);
			derivatives[16, 1] = (x * x - 1) / 4 * (z + 1);
			derivatives[17, 1] = (y * (x - 1) * (z + 1)) / 2;
			derivatives[18, 1] = -(y * (x + 1) * (z + 1)) / 2;
			derivatives[19, 1] = -(x * x - 1) / 4 * (z + 1);

			derivatives[0, 2] = (x - 1) / 8 * (y - 1) * (x + y + z + 2) + (x - 1) / 8 * (y - 1) * (z - 1);
			derivatives[1, 2] = -(x + 1) / 8 * (y - 1) * (z - 1) - (x + 1) / 8 * (y - 1) * (y - x + z + 2);
			derivatives[2, 2] = (x + 1) / 8 * (y + 1) * (z - 1) - (x + 1) / 8 * (y + 1) * (x + y - z - 2);
			derivatives[3, 2] = -(x - 1) / 8 * (y + 1) * (z - 1) - (x - 1) / 8 * (y + 1) * (x - y + z + 2);
			derivatives[4, 2] = (x - 1) / 8 * (y - 1) * (z + 1) - (x - 1) / 8 * (y - 1) * (x + y - z + 2);
			derivatives[5, 2] = -(x + 1) / 8 * (y - 1) * (z + 1) - (x + 1) / 8 * (y - 1) * (x - y + z - 2);
			derivatives[6, 2] = (x + 1) / 8 * (y + 1) * (x + y + z - 2) + (x + 1) / 8 * (y + 1) * (z + 1);
			derivatives[7, 2] = (x - 1) / 8 * (y + 1) * (x - y - z + 2) - (x - 1) / 8 * (y + 1) * (z + 1);
			derivatives[8, 2] = -(x * x - 1) / 4 * (y - 1);
			derivatives[9, 2] = -(y * y - 1) / 4 * (x - 1);
			derivatives[10, 2] = -(z * (x - 1) * (y - 1)) / 2;
			derivatives[11, 2] = (y * y - 1) / 4 * (x + 1);
			derivatives[12, 2] = (z * (x + 1) * (y - 1)) / 2;
			derivatives[13, 2] = (x * x - 1) / 4 * (y + 1);
			derivatives[14, 2] = -(z * (x + 1) * (y + 1)) / 2;
			derivatives[15, 2] = (z * (x - 1) * (y + 1)) / 2;
			derivatives[16, 2] = (x * x - 1) / 4 * (y - 1);
			derivatives[17, 2] = (y * y - 1) / 4 * (x - 1);
			derivatives[18, 2] = -(y * y - 1) / 4 * (x + 1);
			derivatives[19, 2] = -(x * x - 1) / 4 * (y + 1);
			return derivatives;
		}
	}
}
