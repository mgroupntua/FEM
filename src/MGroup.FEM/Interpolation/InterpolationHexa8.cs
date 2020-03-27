using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Mesh;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a hexahedral finite element with 8 nodes. Linear shape functions.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InterpolationHexa8 : IsoparametricInterpolation3DBase
	{
		public override CellType CellType { get; } = CellType.Hexa8;

		private const double oneOverEight = 0.125;

		private static readonly InterpolationHexa8 uniqueInstance = new InterpolationHexa8();

		private InterpolationHexa8() : base(8)
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
			};
		}

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
		/// nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique <see cref="InterpolationHexa8"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationHexa8 UniqueInstance => uniqueInstance;

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation2D.CheckElementNodes(IReadOnlyList{Node})"/>
		/// </summary>
		public override void CheckElementNodes(IReadOnlyList<Node> nodes)
		{
			if (nodes.Count != 8) throw new ArgumentException(
				$"A Hexa8 finite element has 8 nodes, but {nodes.Count} nodes were provided.");
			// TODO: Also check the order of the nodes too and perhaps even the shape
		}

		/// <summary>
		/// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> nodes) =>
			new InverseInterpolationHexa8(nodes);

		protected sealed override double[] EvaluateAt(double xi, double eta, double zeta)
		{
			var values = new double[8];
			values[0] = oneOverEight * (1 - xi) * (1 - eta) * (1 - zeta);
			values[1] = oneOverEight * (1 + xi) * (1 - eta) * (1 - zeta);
			values[2] = oneOverEight * (1 + xi) * (1 + eta) * (1 - zeta);
			values[3] = oneOverEight * (1 - xi) * (1 + eta) * (1 - zeta);
			values[4] = oneOverEight * (1 - xi) * (1 - eta) * (1 + zeta);
			values[5] = oneOverEight * (1 + xi) * (1 - eta) * (1 + zeta);
			values[6] = oneOverEight * (1 + xi) * (1 + eta) * (1 + zeta);
			values[7] = oneOverEight * (1 - xi) * (1 + eta) * (1 + zeta);
			return values;
		}

		protected sealed override Matrix EvaluateGradientsAt(double xi, double eta, double zeta)
		{
			var derivatives = Matrix.CreateZero(8, 3);

			derivatives[0, 0] = -oneOverEight * (1 - eta) * (1 - zeta);
			derivatives[1, 0] = +oneOverEight * (1 - eta) * (1 - zeta);
			derivatives[2, 0] = +oneOverEight * (1 + eta) * (1 - zeta);
			derivatives[3, 0] = -oneOverEight * (1 + eta) * (1 - zeta);
			derivatives[4, 0] = -oneOverEight * (1 - eta) * (1 + zeta);
			derivatives[5, 0] = +oneOverEight * (1 - eta) * (1 + zeta);
			derivatives[6, 0] = +oneOverEight * (1 + eta) * (1 + zeta);
			derivatives[7, 0] = -oneOverEight * (1 + eta) * (1 + zeta);

			derivatives[0, 1] = -oneOverEight * (1 - xi) * (1 - zeta);
			derivatives[1, 1] = -oneOverEight * (1 + xi) * (1 - zeta);
			derivatives[2, 1] = +oneOverEight * (1 + xi) * (1 - zeta);
			derivatives[3, 1] = +oneOverEight * (1 - xi) * (1 - zeta);
			derivatives[4, 1] = -oneOverEight * (1 - xi) * (1 + zeta);
			derivatives[5, 1] = -oneOverEight * (1 + xi) * (1 + zeta);
			derivatives[6, 1] = +oneOverEight * (1 + xi) * (1 + zeta);
			derivatives[7, 1] = +oneOverEight * (1 - xi) * (1 + zeta);

			derivatives[0, 2] = -oneOverEight * (1 - xi) * (1 - eta);
			derivatives[1, 2] = -oneOverEight * (1 + xi) * (1 - eta);
			derivatives[2, 2] = -oneOverEight * (1 + xi) * (1 + eta);
			derivatives[3, 2] = -oneOverEight * (1 - xi) * (1 + eta);
			derivatives[4, 2] = +oneOverEight * (1 - xi) * (1 - eta);
			derivatives[5, 2] = +oneOverEight * (1 + xi) * (1 - eta);
			derivatives[6, 2] = +oneOverEight * (1 + xi) * (1 + eta);
			derivatives[7, 2] = +oneOverEight * (1 - xi) * (1 + eta);

			#region untested
			//var x = xi;
			//var y = eta;
			//var z = zeta;

			//derivatives[0, 0] = -(y - 1) * (z - 1) / 8;
			//derivatives[1, 0] = (y - 1) * (z - 1) / 8;
			//derivatives[2, 0] = -(y + 1) * (z - 1) / 8;
			//derivatives[3, 0] = (y + 1) * (z - 1) / 8;
			//derivatives[4, 0] = (y - 1) * (z + 1) / 8;
			//derivatives[5, 0] = -(y - 1) * (z + 1) / 8;
			//derivatives[6, 0] = (y + 1) * (z + 1) / 8;
			//derivatives[7, 0] = -(y + 1) * (z + 1) / 8;

			//derivatives[0, 1] = -(x - 1) * (z - 1) / 8;
			//derivatives[1, 1] = (x + 1) * (z - 1) / 8;
			//derivatives[2, 1] = -(x + 1) * (z - 1) / 8;
			//derivatives[3, 1] = (x - 1) * (z - 1) / 8;
			//derivatives[4, 1] = (x - 1) * (z + 1) / 8;
			//derivatives[5, 1] = -(x + 1) * (z + 1) / 8;
			//derivatives[6, 1] = (x + 1) * (z + 1) / 8;
			//derivatives[7, 1] = -(x - 1) * (z + 1) / 8;

			//derivatives[0, 2] = -(x - 1) * (y - 1) / 8;
			//derivatives[1, 2] = (x + 1) * (y - 1) / 8;
			//derivatives[2, 2] = -(x + 1) * (y + 1) / 8;
			//derivatives[3, 2] = (x - 1) * (y + 1) / 8;
			//derivatives[4, 2] = (x - 1) * (y - 1) / 8;
			//derivatives[5, 2] = -(x + 1) * (y - 1) / 8;
			//derivatives[6, 2] = (x + 1) * (y + 1) / 8;
			//derivatives[7, 2] = -(x - 1) * (y + 1) / 8;
			#endregion

			return derivatives;
		}
	}
}
