using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of pyramid finite element with 5 nodes. Linear shape functions.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InterpolationPyra5 : IsoparametricInterpolation3DBase
	{
		private static readonly InterpolationPyra5 uniqueInstance = new InterpolationPyra5();

		private InterpolationPyra5() : base(5)
		{
			NodalNaturalCoordinates = new NaturalPoint[]
			{
				new NaturalPoint(1,0,0),
				new NaturalPoint(0,1,0),
				new NaturalPoint(-1,0,0),
				new NaturalPoint(0,-1,0),
				new NaturalPoint(0,0,1)
			};
		}

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order
		/// of these nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique instance <see cref="InterpolationPyra5"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationPyra5 UniqueInstance => uniqueInstance;

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation2D.CheckElementNodes(IReadOnlyList{Node})"/>
		/// </summary>
		public override void CheckElementNodes(IReadOnlyList<Node> nodes)
		{
			if (nodes.Count != 5) throw new ArgumentException(
				$"A Pyra5 finite element has 5 nodes, but {nodes.Count} nodes were provided.");
			// TODO: Also check the order of the nodes too and perhaps even the shape
		}

		/// <summary>
		/// The reverse mapping for this interpolation, namely from global cartesian coordinates to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> node)
			=> throw new NotImplementedException();

		protected sealed override double[] EvaluateAt(double xi, double eta, double zeta)
		{
			var values = new double[5];
			var x = xi;
			var y = eta;
			var z = zeta;

			values[0] = (Math.Abs(z - 1) < 10e-10) ? 0 : (-x + y + z - 1) * (-x - y + z - 1) / (4 * (1 - z));
			values[1] = (Math.Abs(z - 1) < 10e-10) ? 0 : (-x - y + z - 1) * (x - y + z - 1) / (4 * (1 - z));
			values[2] = (Math.Abs(z - 1) < 10e-10) ? 0 : (x + y + z - 1) * (x - y + z - 1) / (4 * (1 - z));
			values[3] = (Math.Abs(z - 1) < 10e-10) ? 0 : (x + y + z - 1) * (-x + y + z - 1) / (4 * (1 - z));
			values[4] = z;

			return values;
		}

		protected override Matrix EvaluateGradientsAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;

			var derivatives = Matrix.CreateZero(5, 3);
			derivatives[0, 0] = -(x + y - z + 1) / (4 * z - 4) - (x - y - z + 1) / (4 * z - 4);
			derivatives[1, 0] = (x + y - z + 1) / (4 * z - 4) + (x - y + z - 1) / (4 * z - 4);
			derivatives[2, 0] = -(x - y + z - 1) / (4 * z - 4) - (x + y + z - 1) / (4 * z - 4);
			derivatives[3, 0] = (x - y - z + 1) / (4 * z - 4) + (x + y + z - 1) / (4 * z - 4);
			derivatives[4, 0] = 0.0;

			derivatives[0, 1] = (x + y - z + 1) / (4 * z - 4) - (x - y - z + 1) / (4 * z - 4);
			derivatives[1, 1] = (x - y + z - 1) / (4 * z - 4) - (x + y - z + 1) / (4 * z - 4);
			derivatives[2, 1] = (x + y + z - 1) / (4 * z - 4) - (x - y + z - 1) / (4 * z - 4);
			derivatives[3, 1] = (x - y - z + 1) / (4 * z - 4) - (x + y + z - 1) / (4 * z - 4);
			derivatives[4, 1] = 0.0;

			derivatives[0, 2] = (x + y - z + 1) / (4 * z - 4) + (x - y - z + 1) / (4 * z - 4) +
								(4 * (x + y - z + 1) * (x - y - z + 1)) / Math.Pow(4 * z - 4, 2);
			derivatives[1, 2] = (x + y - z + 1) / (4 * z - 4) - (x - y + z - 1) / (4 * z - 4) -
								(4 * (x + y - z + 1) * (x - y + z - 1)) / Math.Pow(4 * z - 4, 2);
			derivatives[2, 2] = (4 * (x - y + z - 1) * (x + y + z - 1)) / Math.Pow(4 * z - 4, 2) -
								(x + y + z - 1) / (4 * z - 4) - (x - y + z - 1) / (4 * z - 4);
			derivatives[3, 2] = (x - y - z + 1) / (4 * z - 4) - (x + y + z - 1) / (4 * z - 4) -
								(4 * (x - y - z + 1) * (x + y + z - 1)) / Math.Pow(4 * z - 4, 2);
			derivatives[4, 2] = 1.0;

			return derivatives;
		}
	}
}
