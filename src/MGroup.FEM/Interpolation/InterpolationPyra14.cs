using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// Isoparamateric interpolation of a pyramid finite element with 14 nodes. Quadratic shape functions.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InterpolationPyra14 : IsoparametricInterpolation3DBase
	{
		private static readonly InterpolationPyra14 uniqueInstance = new InterpolationPyra14();

		private InterpolationPyra14() : base(14)
		{
			NodalNaturalCoordinates = new NaturalPoint[]
			{
				new NaturalPoint(-1,-1,-1),
				new NaturalPoint(1,-1,-1),
				new NaturalPoint(1,1,-1),
				new NaturalPoint(-1,1,-1),
				new NaturalPoint(0,0,1),

				new NaturalPoint(0,-1,-1),
				new NaturalPoint(-1,0,-1),
				new NaturalPoint(-0.5,-0.5,0),
				new NaturalPoint(1,0,-1),
				new NaturalPoint(0.5,-0.5,0),

				new NaturalPoint(0,1,-1),
				new NaturalPoint(0.5,0.5,0),
				new NaturalPoint(-0.5,0.5,0),
				new NaturalPoint(0, 0, -1),
			};
		}

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order
		/// of these nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique instance <see cref="InterpolationPyra14"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationPyra14 UniqueInstance => uniqueInstance;

		/// <summary>
		/// See <see cref="IIsoparametricInterpolation2D.CheckElementNodes(IReadOnlyList{Node})"/>
		/// </summary>
		public override void CheckElementNodes(IReadOnlyList<Node> nodes)
		{
			if (nodes.Count != 14) throw new ArgumentException(
				$"A Pyra14 finite element has 14 nodes, but {nodes.Count} nodes were provided.");
			// TODO: Also check the order of the nodes too and perhaps even the shape
		}

		/// <summary>
		/// The reverse mapping for this interpolation, namely from global cartesian coordinates to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> node)
			=> throw new NotImplementedException("Iterative procedure required");

		// based on https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch12.d/AFEM.Ch12.pdf
		protected override double[] EvaluateAt(double xi, double eta, double zeta)
		{
			var Nc = 0.5 * (1 - Math.Pow(xi, 2)) * (1 - Math.Pow(eta, 2)) * (1 - zeta);

			var values = new double[14];
			values[0] = -1 / 16.0 * (1 - xi) * (1 - eta) * (1 - zeta) * (4 + 3 * xi + 3 * eta + 2 * xi * eta + 2 * zeta + xi * zeta + eta * zeta + 2 * xi * eta * zeta) + 1 / 4.0 * Nc;
			values[1] = -1 / 16.0 * (1 + xi) * (1 - eta) * (1 - zeta) * (4 - 3 * xi + 3 * eta - 2 * xi * eta + 2 * zeta - xi * zeta + eta * zeta - 2 * xi * eta * zeta) + 1 / 4.0 * Nc;
			values[2] = -1 / 16.0 * (1 + xi) * (1 + eta) * (1 - zeta) * (4 - 3 * xi - 3 * eta + 2 * xi * eta + 2 * zeta - xi * zeta - eta * zeta + 2 * xi * eta * zeta) + 1 / 4.0 * Nc;
			values[3] = -1 / 16.0 * (1 - xi) * (1 + eta) * (1 - zeta) * (4 + 3 * xi - 3 * eta - 2 * xi * eta + 2 * zeta + xi * zeta - eta * zeta - 2 * xi * eta * zeta) + 1 / 4.0 * Nc;
			values[4] = 0.5 * zeta * (1 + zeta);
			values[5] = 1 / 8.0 * (1 - Math.Pow(xi, 2)) * (1 - eta) * (1 - zeta) * (2 + eta + eta * zeta) - 0.5 * Nc;
			values[6] = 1 / 8.0 * (1 - xi) * (1 - Math.Pow(eta, 2)) * (1 - zeta) * (2 + xi + xi * zeta) - 0.5 * Nc;
			values[7] = 1 / 4.0 * (1 - xi) * (1 - eta) * (1 - Math.Pow(zeta, 2));
			values[8] = 1 / 8.0 * (1 + xi) * (1 - Math.Pow(eta, 2)) * (1 - zeta) * (2 - xi - xi * zeta) - 0.5 * Nc;
			values[9] = 1 / 4.0 * (1 + xi) * (1 - eta) * (1 - Math.Pow(zeta, 2));
			values[10] = 1 / 8.0 * (1 - Math.Pow(xi, 2)) * (1 + eta) * (1 - zeta) * (2 - eta - eta * zeta) - 0.5 * Nc;
			values[11] = 1 / 4.0 * (1 + xi) * (1 + eta) * (1 - Math.Pow(zeta, 2));
			values[12] = 1 / 4.0 * (1 - xi) * (1 + eta) * (1 - Math.Pow(zeta, 2));
			values[13] = Nc;

			return values;
		}


		protected override Matrix EvaluateGradientsAt(double xi, double eta, double zeta)
		{
			var Ncx = -xi * (1 - eta * eta) * (1 - zeta);
			var Ncy = -eta * (1 - xi * xi) * (1 - zeta);
			var Ncz = -0.5 * (1 - xi * xi) * (1 - eta * eta);

			var derivatives = Matrix.CreateZero(14, 3);

			derivatives[0, 0] = 1 / 16.0 * (1 - eta) * (1 - zeta) * (1 + 6 * xi + eta + 4 * xi * eta + zeta + 2 * xi * zeta - eta * zeta + 4 * xi * eta * zeta) + 1 / 4.0 * Ncx;
			derivatives[1, 0] = -1 / 16.0 * (1 - eta) * (1 - zeta) * (1 - 6 * xi + eta - 4 * xi * eta + zeta - 2 * xi * zeta - eta * zeta - 4 * xi * eta * zeta) + 1 / 4.0 * Ncx;
			derivatives[2, 0] = -1 / 16.0 * (1 + eta) * (1 - zeta) * (1 - 6 * xi - eta + 4 * xi * eta + zeta - 2 * xi * zeta + eta * zeta + 4 * xi * eta * zeta) + 1 / 4.0 * Ncx;
			derivatives[3, 0] = 1 / 16.0 * (1 + eta) * (1 - zeta) * (1 + 6 * xi - eta - 4 * xi * eta + zeta + 2 * xi * zeta + eta * zeta - 4 * xi * eta * zeta) + 1 / 4.0 * Ncx;
			derivatives[4, 0] = 0.0;
			derivatives[5, 0] = -1 / 4.0 * xi * (1 - eta) * (1 - zeta) * (2 + eta + eta * zeta) - 0.5 * Ncx;
			derivatives[6, 0] = -1 / 8.0 * (1 - eta * eta) * (1 - zeta) * (1 + 2 * xi - zeta + 2 * xi * zeta) - 0.5 * Ncx;
			derivatives[7, 0] = -1 / 4.0 * (1 - eta) * (1 - zeta * zeta);
			derivatives[8, 0] = 1 / 8.0 * (1 - eta * eta) * (1 - zeta) * (1 - 2 * xi - zeta - 2 * xi * zeta) - 0.5 * Ncx;
			derivatives[9, 0] = 1 / 4.0 * (1 - eta) * (1 - zeta * zeta);
			derivatives[10, 0] = -1 / 4.0 * xi * (1 + eta) * (1 - zeta) * (2 - eta - eta * zeta) - 0.5 * Ncx;
			derivatives[11, 0] = 1 / 4.0 * (1 + eta) * (1 - zeta * zeta);
			derivatives[12, 0] = -1 / 4.0 * (1 + eta) * (1 - zeta * zeta);
			derivatives[13, 0] = Ncx;

			derivatives[0, 1] = 1 / 16.0 * (1 - xi) * (1 + zeta) * (1 + xi + 6 * eta + 4 * xi * eta + zeta - xi * zeta + 2 * eta * zeta + 4 * xi * eta * zeta) + 1 / 4.0 * Ncy;
			derivatives[1, 1] = 1 / 16.0 * (1 + xi) * (1 - zeta) * (1 - xi + 6 * eta - 4 * xi * eta + zeta + xi * zeta + 2 * eta * zeta - 4 * xi * eta * zeta) + 1 / 4.0 * Ncy;
			derivatives[2, 1] = -1 / 16.0 * (1 + xi) * (1 - zeta) * (1 - xi - 6 * eta + 4 * xi * eta + zeta + xi * zeta - 2 * eta * zeta + 4 * xi * eta * zeta) + 1 / 4.0 * Ncy;
			derivatives[3, 1] = -1 / 16.0 * (1 - xi) * (1 + zeta) * (1 + xi - 6 * eta - 4 * xi * eta + zeta - xi * zeta - 2 * eta * zeta - 4 * xi * eta * zeta) + 1 / 4.0 * Ncy;
			derivatives[4, 1] = 0.0;
			derivatives[5, 1] = -1 / 8.0 * (1 - xi * xi) * (1 - zeta) * (1 + 2 * eta - zeta + 2 * eta * zeta) - 0.5 * Ncy;
			derivatives[6, 1] = -1 / 4.0 * (1 - xi) * eta * (1 - zeta) * (2 + xi + xi * zeta) - 0.5 * Ncy;
			derivatives[7, 1] = -1 / 4.0 * (1 - xi) * (1 - zeta * zeta);
			derivatives[8, 1] = -1 / 4.0 * (1 + xi) * eta * (1 - zeta) * (2 - xi - xi * zeta) - 0.5 * Ncy;
			derivatives[9, 1] = -1 / 4.0 * (1 + xi) * (1 - zeta * zeta); ;
			derivatives[10, 1] = 1 / 8.0 * (1 - xi * xi) * (1 - zeta) * (1 - 2 * eta - zeta - 2 * eta * zeta) - 0.5 * Ncy;
			derivatives[11, 1] = 1 / 4.0 * (1 + xi) * (1 - zeta * zeta); ;
			derivatives[12, 1] = 1 / 4.0 * (1 - xi) * (1 - zeta * zeta); ;
			derivatives[13, 1] = Ncy;

			derivatives[0, 2] = 1 / 8.0 * (1 - xi) * (1 - eta) * (1 + xi + eta + 2 * zeta + xi * zeta + eta * zeta + 2 * xi * eta * zeta) + 1 / 4.0 * Ncz;
			derivatives[1, 2] = 1 / 8.0 * (1 + xi) * (1 - eta) * (1 - xi + eta + 2 * zeta - xi * zeta + eta * zeta - 2 * xi * eta * zeta) + 1 / 4.0 * Ncz; ;
			derivatives[2, 2] = 1 / 8.0 * (1 + xi) * (1 + eta) * (1 - xi - eta + 2 * zeta - xi * zeta - eta * zeta + 2 * xi * eta * zeta) + 1 / 4.0 * Ncz; ;
			derivatives[3, 2] = 1 / 8.0 * (1 - xi) * (1 + eta) * (1 + xi - eta + 2 * zeta + xi * zeta - eta * zeta - 2 * xi * eta * zeta) + 1 / 4.0 * Ncz; ;
			derivatives[4, 2] = 0.5 + zeta;
			derivatives[5, 2] = -1 / 4.0 * (1 - xi * xi) * (1 - eta) * (1 + eta * zeta) - 0.5 * Ncz;
			derivatives[6, 2] = -1 / 4.0 * (1 - xi) * (1 - eta * eta) * (1 + xi * zeta) - 0.5 * Ncz; ;
			derivatives[7, 2] = -0.5 * (1 - xi) * (1 - eta) * zeta;
			derivatives[8, 2] = -1 / 4.0 * (1 + xi) * (1 - eta * eta) * (1 - xi * zeta) - 0.5 * Ncz;
			derivatives[9, 2] = -0.5 * (1 + xi) * (1 - eta) * zeta;
			derivatives[10, 2] = -1 / 4.0 * (1 - xi * xi) * (1 + eta) * (1 - eta * zeta) - 0.5 * Ncz; ;
			derivatives[11, 2] = -0.5 * (1 + xi) * (1 + eta) * zeta; ;
			derivatives[12, 2] = -0.5 * (1 - xi) * (1 + eta) * zeta; ;
			derivatives[13, 2] = Ncz;

			return derivatives;
		}
	}
}
