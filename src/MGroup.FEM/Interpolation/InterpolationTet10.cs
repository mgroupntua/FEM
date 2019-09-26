using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a tetrahedral finite element with 10 nodes. Quadratic shape function.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis.
	/// </summary>
	public class InterpolationTet10 : IsoparametricInterpolation3DBase
	{
		private static readonly InterpolationTet10 uniqueInstance = new InterpolationTet10();

		private InterpolationTet10() : base(10)
		{
			NodalNaturalCoordinates = new NaturalPoint[]
			{
				new NaturalPoint(0,0,0),
				new NaturalPoint(1,0,0),
				new NaturalPoint(0,1,0),
				new NaturalPoint(0,0,1),

				new NaturalPoint(0.5,0,0),
				new NaturalPoint(0.5,0.5,0),
				new NaturalPoint(0,0.5,0),

				new NaturalPoint(0,0,0.5),
				new NaturalPoint(0,0.5,0.5),
				new NaturalPoint(0.5,0,0.5),
			};
		}

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order
		/// of these nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique instance <see cref="InterpolationTet10"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationTet10 UniqueInstance => uniqueInstance;

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation2D.CheckElementNodes(IReadOnlyList{Node})"/>
		/// </summary>
		public override void CheckElementNodes(IReadOnlyList<Node> nodes)
		{
			if (nodes.Count != 10) throw new ArgumentException(
				$"A Tetra10 finite element has 10 nodes, but {nodes.Count} nodes were provided.");
			// TODO: Also check the order of the nodes too and perhaps even the shape
		}

		/// <summary>
		/// The reverse mapping for this interpolation, namely from global cartesian coordinates to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> node)
			=> throw new NotImplementedException("Requires iterative procedure.");


		/// <summary>
		/// Returns the shape functions a tetrahedral quadratic element evaluated on a single point.
		/// Implementation is based on <see cref="https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch10.d/AFEM.Ch10.pdf">Carlos Felippa - Introduction to Finite Element Methods</see>
		/// </summary>
		/// <param name="xi"></param>
		/// <param name="eta"></param>
		/// <param name="zeta"></param>
		/// <returns></returns>
		protected sealed override double[] EvaluateAt(double xi, double eta, double zeta)
		{
			var z1 = xi;
			var z2 = eta;
			var z3 = zeta;
			var z4 = 1 - z1 - z2 - z3;

			var values = new double[10];
			values[0] = z4 * (2 * z4 - 1);
			values[1] = z1 * (2 * z1 - 1);
			values[2] = z2 * (2 * z2 - 1);
			values[3] = z3 * (2 * z3 - 1);
			values[4] = 4 * z1 * z4;
			values[5] = 4 * z1 * z2;
			values[6] = 4 * z2 * z4;
			values[7] = 4 * z3 * z4;
			values[8] = 4 * z2 * z3;
			values[9] = 4 * z3 * z1;

			return values;
		}


		protected sealed override Matrix EvaluateGradientsAt(double xi, double eta, double zeta)
		{
			var x = xi;
			var y = eta;
			var z = zeta;

			var derivatives = Matrix.CreateZero(10, 3);

			derivatives[0, 0] = 4 * x + 4 * y + 4 * z - 3;
			derivatives[1, 0] = 4 * x - 1;
			derivatives[2, 0] = 0.0;
			derivatives[3, 0] = 0.0;
			derivatives[4, 0] = 4 - 4 * y - 4 * z - 8 * x;
			derivatives[5, 0] = 4 * y;
			derivatives[6, 0] = -4 * y;
			derivatives[7, 0] = -4 * z;
			derivatives[8, 0] = 0.0;
			derivatives[9, 0] = 4 * z;

			derivatives[0, 1] = 4 * x + 4 * y + 4 * z - 3;
			derivatives[1, 1] = 0.0;
			derivatives[2, 1] = 4 * y - 1;
			derivatives[3, 1] = 0.0;
			derivatives[4, 1] = -4 * x;
			derivatives[5, 1] = 4 * x;
			derivatives[6, 1] = 4 - 8 * y - 4 * z - 4 * x;
			derivatives[7, 1] = -4 * z;
			derivatives[8, 1] = 4 * z;
			derivatives[9, 1] = 0.0;

			derivatives[0, 2] = 4 * x + 4 * y + 4 * z - 3;
			derivatives[1, 2] = 0.0;
			derivatives[2, 2] = 0.0;
			derivatives[3, 2] = 4 * z - 1;
			derivatives[4, 2] = -4 * x;
			derivatives[5, 2] = 0.0;
			derivatives[6, 2] = -4 * y;
			derivatives[7, 2] = 4 - 4 * y - 8 * z - 4 * x;
			derivatives[8, 2] = 4 * y;
			derivatives[9, 2] = 4 * x;

			return derivatives;
		}
	}
}
