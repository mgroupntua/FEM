using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a pyramid finite element with 13 nodes. Quadratic shape functions.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InterpolationPyra13 : IsoparametricInterpolation3DBase
	{
		private static readonly InterpolationPyra13 uniqueInstance = new InterpolationPyra13();

		private InterpolationPyra13() : base(13)
		{
			NodalNaturalCoordinates = new NaturalPoint[]
			{
				new NaturalPoint(1,0,0),
				new NaturalPoint(0,1,0),
				new NaturalPoint(-1,0,0),
				new NaturalPoint(0,-1,0),
				new NaturalPoint(0,0,1),

				new NaturalPoint(0.5,0.5,0),
				new NaturalPoint(0.5,-0.5,0),
				new NaturalPoint(0.5,0,0.5),
				new NaturalPoint(-0.5,0.5,0),
				new NaturalPoint(0,0.5,0.5),

				new NaturalPoint(-0.5,-0.5,0),
				new NaturalPoint(-0.5,0,0.5),
				new NaturalPoint(0,-0.5,0.5),
			};
		}

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order
		/// of these nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique instance <see cref="InterpolationPyra13"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationPyra13 UniqueInstance => uniqueInstance;

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation2D.CheckElementNodes(IReadOnlyList{Node})"/>
		/// </summary>
		public override void CheckElementNodes(IReadOnlyList<Node> nodes)
		{
			if (nodes.Count != 13) throw new ArgumentException(
				$"A Pyra13 finite element has 13 nodes, but {nodes.Count} nodes were provided.");
			// TODO: Also check the order of the nodes too and perhaps even the shape
		}

		/// <summary>
		/// The reverse mapping for this interpolation, namely from global cartesian coordinates to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> node)
			=> throw new NotImplementedException("Iterative procedure required");

		protected override double[] EvaluateAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;

			var values = new double[13];
			values[0] = (Math.Abs(z - 1) < 10e-10) ? 0 : (-x + y + z - 1) * (-x - y + z - 1) * (x - 0.5) / (2 * (1 - z));
			values[1] = (Math.Abs(z - 1) < 10e-10) ? 0 : (-x - y + z - 1) * (x - y + z - 1) * (y - 0.5) / (2 * (1 - z));
			values[2] = (Math.Abs(z - 1) < 10e-10) ? 0 : (x - y + z - 1) * (x + y + z - 1) * (-x - 0.5) / (2 * (1 - z));
			values[3] = (Math.Abs(z - 1) < 10e-10) ? 0 : (x + y + z - 1) * (-x + y + z - 1) * (-y - 0.5) / (2 * (1 - z));
			values[4] = 2 * z * (z - 0.5);
			values[5] = (Math.Abs(z - 1) < 10e-10) ? 0 : -(-x + y + z - 1) * (-x - y + z - 1) * (x - y + z - 1) / (2 * (1 - z));
			values[6] = (Math.Abs(z - 1) < 10e-10) ? 0 : -(x + y + z - 1) * (-x + y + z - 1) * (-x - y + z - 1) / (2 * (1 - z));
			values[7] = (Math.Abs(z - 1) < 10e-10) ? 0 : z * (-x + y + z - 1) * (-x - y + z - 1) / (1 - z);
			values[8] = (Math.Abs(z - 1) < 10e-10) ? 0 : -(-x - y + z - 1) * (x - y + z - 1) * (x + y + z - 1) / (2 * (1 - z));
			values[9] = (Math.Abs(z - 1) < 10e-10) ? 0 : z * (-x - y + z - 1) * (x - y + z - 1) / (1 - z);
			values[10] = (Math.Abs(z - 1) < 10e-10) ? 0 : -(x - y + z - 1) * (x + y + z - 1) * (-x + y + z - 1) / (2 * (1 - z));
			values[11] = (Math.Abs(z - 1) < 10e-10) ? 0 : z * (x - y + z - 1) * (x + y + z - 1) / (1 - z);
			values[12] = (Math.Abs(z - 1) < 10e-10) ? 0 : z * (x + y + z - 1) * (-x + y + z - 1) / (1 - z);
			return values;
		}


		// TODO: verify derivatives of Pyra13
		protected override Matrix EvaluateGradientsAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;

			var derivatives = Matrix.CreateZero(13, 3);

			derivatives[0, 0] = -((x - 1 / 2.0) * (x - y - z + 1)) / (2 * z - 2) -
								((x + y - z + 1) * (x - y - z + 1)) / (2 * z - 2) -
								((x - 1 / 2.0) * (x + y - z + 1)) / (2 * z - 2);
			derivatives[1, 0] = ((y - 1 / 2.0) * (x + y - z + 1)) / (2 * z - 2) +
								((y - 1 / 2.0) * (x - y + z - 1)) / (2 * z - 2);
			derivatives[2, 0] = ((x + 1 / 2.0) * (x + y + z - 1)) / (2 * z - 2) +
								((x - y + z - 1) * (x + y + z - 1)) / (2 * z - 2) +
								((x + 1 / 2.0) * (x - y + z - 1)) / (2 * z - 2);
			derivatives[3, 0] = -((y + 1 / 2.0) * (x - y - z + 1)) / (2 * z - 2) -
								((y + 1 / 2.0) * (x + y + z - 1)) / (2 * z - 2);
			derivatives[4, 0] = 0.0;
			derivatives[5, 0] = ((x + y - z + 1) * (x - y + z - 1)) / (2 * z - 2) +
								((x + y - z + 1) * (x - y - z + 1)) / (2 * z - 2) +
								((x - y + z - 1) * (x - y - z + 1)) / (2 * z - 2);
			derivatives[6, 0] = ((x - y - z + 1) * (x + y + z - 1)) / (2 * z - 2) +
								((x + y - z + 1) * (x - y - z + 1)) / (2 * z - 2) +
								((x + y - z + 1) * (x + y + z - 1)) / (2 * z - 2);
			derivatives[7, 0] = -(z * (x + y - z + 1)) / (z - 1) - (z * (x - y - z + 1)) / (z - 1);
			derivatives[8, 0] = -((x + y - z + 1) * (x - y + z - 1)) / (2 * z - 2) -
								((x + y - z + 1) * (x + y + z - 1)) / (2 * z - 2) -
								((x - y + z - 1) * (x + y + z - 1)) / (2 * z - 2);
			derivatives[9, 0] = (z * (x + y - z + 1)) / (z - 1) + (z * (x - y + z - 1)) / (z - 1);
			derivatives[10, 0] = -((x - y - z + 1) * (x + y + z - 1)) / (2 * z - 2) -
								 ((x - y + z - 1) * (x - y - z + 1)) / (2 * z - 2) -
								 ((x - y + z - 1) * (x + y + z - 1)) / (2 * z - 2);
			derivatives[11, 0] = -(z * (x - y + z - 1)) / (z - 1) - (z * (x + y + z - 1)) / (z - 1);
			derivatives[12, 0] = (z * (x - y - z + 1)) / (z - 1) + (z * (x + y + z - 1)) / (z - 1);

			derivatives[0, 1] = ((x - 1 / 2.0) * (x + y - z + 1)) / (2 * z - 2) -
								((x - 1 / 2.0) * (x - y - z + 1)) / (2 * z - 2);
			derivatives[1, 1] = ((x + y - z + 1) * (x - y + z - 1)) / (2 * z - 2) -
								((y - 1 / 2.0) * (x + y - z + 1)) / (2 * z - 2) +
								((y - 1 / 2.0) * (x - y + z - 1)) / (2 * z - 2);
			derivatives[2, 1] = ((x + 1 / 2.0) * (x - y + z - 1)) / (2 * z - 2) -
								((x + 1 / 2.0) * (x + y + z - 1)) / (2 * z - 2);
			derivatives[3, 1] = ((y + 1 / 2.0) * (x + y + z - 1)) / (2 * z - 2) -
								((y + 1 / 2.0) * (x - y - z + 1)) / (2 * z - 2) -
								((x - y - z + 1) * (x + y + z - 1)) / (2 * z - 2);
			derivatives[4, 1] = 0.0;
			derivatives[5, 1] = ((x - y + z - 1) * (x - y - z + 1)) / (2 * z - 2) -
								((x + y - z + 1) * (x - y - z + 1)) / (2 * z - 2) -
								((x + y - z + 1) * (x - y + z - 1)) / (2 * z - 2);
			derivatives[6, 1] = ((x - y - z + 1) * (x + y + z - 1)) / (2 * z - 2) +
								((x + y - z + 1) * (x - y - z + 1)) / (2 * z - 2) -
								((x + y - z + 1) * (x + y + z - 1)) / (2 * z - 2);
			derivatives[7, 1] = (z * (x + y - z + 1)) / (z - 1) - (z * (x - y - z + 1)) / (z - 1);
			derivatives[8, 1] = ((x + y - z + 1) * (x + y + z - 1)) / (2 * z - 2) -
								((x + y - z + 1) * (x - y + z - 1)) / (2 * z - 2) -
								((x - y + z - 1) * (x + y + z - 1)) / (2 * z - 2);
			derivatives[9, 1] = (z * (x - y + z - 1)) / (z - 1) - (z * (x + y - z + 1)) / (z - 1);
			derivatives[10, 1] = ((x - y - z + 1) * (x + y + z - 1)) / (2 * z - 2) -
								 ((x - y + z - 1) * (x - y - z + 1)) / (2 * z - 2) +
								 ((x - y + z - 1) * (x + y + z - 1)) / (2 * z - 2);
			derivatives[11, 1] = (z * (x + y + z - 1)) / (z - 1) - (z * (x - y + z - 1)) / (z - 1);
			derivatives[12, 1] = (z * (x - y - z + 1)) / (z - 1) - (z * (x + y + z - 1)) / (z - 1);

			derivatives[0, 2] = ((x - 1 / 2.0) * (x - y - z + 1)) / (2 * z - 2) +
								((x - 1 / 2.0) * (x + y - z + 1)) / (2 * z - 2) +
								(2 * (x - 1 / 2.0) * (x + y - z + 1) * (x - y - z + 1)) / Math.Pow(2 * z - 2, 2);
			derivatives[1, 2] = ((y - 1 / 2.0) * (x + y - z + 1)) / (2 * z - 2) -
								((y - 1 / 2.0) * (x - y + z - 1)) / (2 * z - 2) -
								(2 * (y - 1 / 2.0) * (x + y - z + 1) * (x - y + z - 1)) / Math.Pow(2 * z - 2, 2);
			derivatives[2, 2] = ((x + 1 / 2.0) * (x + y + z - 1)) / (2 * z - 2) +
								((x + 1 / 2.0) * (x - y + z - 1)) / (2 * z - 2) -
								(2 * (x + 1 / 2.0) * (x - y + z - 1) * (x + y + z - 1)) / Math.Pow(2 * z - 2, 2);
			derivatives[3, 2] = ((y + 1 / 2.0) * (x + y + z - 1)) / (2 * z - 2) -
								((y + 1 / 2.0) * (x - y - z + 1)) / (2 * z - 2) +
								(2 * (y + 1 / 2.0) * (x - y - z + 1) * (x + y + z - 1)) / Math.Pow(2 * z - 2, 2);
			derivatives[4, 2] = 4 * z - 1;
			derivatives[5, 2] = ((x + y - z + 1) * (x - y - z + 1)) / (2 * z - 2) -
								((x + y - z + 1) * (x - y + z - 1)) / (2 * z - 2) -
								((x - y + z - 1) * (x - y - z + 1)) / (2 * z - 2) -
								(2 * (x + y - z + 1) * (x - y + z - 1) * (x - y - z + 1)) / Math.Pow(2 * z - 2, 2);
			derivatives[6, 2] = ((x + y - z + 1) * (x - y - z + 1)) / (2 * z - 2) -
								((x - y - z + 1) * (x + y + z - 1)) / (2 * z - 2) -
								((x + y - z + 1) * (x + y + z - 1)) / (2 * z - 2) -
								(2 * (x + y - z + 1) * (x - y - z + 1) * (x + y + z - 1)) / Math.Pow(2 * z - 2, 2);
			derivatives[7, 2] = (z * (x + y - z + 1)) / (z - 1) - ((x + y - z + 1) * (x - y - z + 1)) / (z - 1) +
								(z * (x - y - z + 1)) / (z - 1) +
								(z * (x + y - z + 1) * (x - y - z + 1)) / Math.Pow(z - 1, 2);
			derivatives[8, 2] = ((x - y + z - 1) * (x + y + z - 1)) / (2 * z - 2) -
								((x + y - z + 1) * (x + y + z - 1)) / (2 * z - 2) -
								((x + y - z + 1) * (x - y + z - 1)) / (2 * z - 2) +
								(2 * (x + y - z + 1) * (x - y + z - 1) * (x + y + z - 1)) / Math.Pow(2 * z - 2, 2);
			derivatives[9, 2] = (z * (x + y - z + 1)) / (z - 1) - (z * (x - y + z - 1)) / (z - 1) +
								((x + y - z + 1) * (x - y + z - 1)) / (z - 1) -
								(z * (x + y - z + 1) * (x - y + z - 1)) / Math.Pow(z - 1, 2);
			derivatives[10, 2] = ((x - y + z - 1) * (x + y + z - 1)) / (2 * z - 2) -
								 ((x - y + z - 1) * (x - y - z + 1)) / (2 * z - 2) -
								 ((x - y - z + 1) * (x + y + z - 1)) / (2 * z - 2) +
								 (2 * (x - y + z - 1) * (x - y - z + 1) * (x + y + z - 1)) / Math.Pow(2 * z - 2, 2);
			derivatives[11, 2] = (z * (x - y + z - 1) * (x + y + z - 1)) / Math.Pow(z - 1, 2) -
								 ((x - y + z - 1) * (x + y + z - 1)) / (z - 1) - (z * (x + y + z - 1)) / (z - 1) -
								 (z * (x - y + z - 1)) / (z - 1);
			derivatives[12, 2] = (z * (x - y - z + 1)) / (z - 1) - (z * (x + y + z - 1)) / (z - 1) +
								 ((x - y - z + 1) * (x + y + z - 1)) / (z - 1) -
								 (z * (x - y - z + 1) * (x + y + z - 1)) / Math.Pow(z - 1, 2);

			return derivatives;
		}
	}
}
