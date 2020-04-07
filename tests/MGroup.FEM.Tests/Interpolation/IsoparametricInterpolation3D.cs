using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Interpolation;
using MGroup.FEM.Interpolation;
using MGroup.MSolve.Geometry.Coordinates;
using Xunit;

namespace MGroup.FEM.Tests.Interpolation
{
	/// <summary>
	/// Unit testing implementations of <see cref="IIsoparametricInterpolation3D_OLD"/>
	/// </summary>
	public class IsoparametricInterpolation3D
	{
		private const int numRandomPoints = 10;
		private delegate NaturalPoint[] GenerateRandomPoints();

		public static readonly IEnumerable<object[]> interpolations = new List<object[]>()
		{
			new object[]{InterpolationTet4.UniqueInstance},
			new object[]{InterpolationTet10.UniqueInstance},
			new object[]{InterpolationHexa8.UniqueInstance},
			new object[]{InterpolationHexa20.UniqueInstance},
			new object[]{InterpolationHexa27.UniqueInstance},
			new object[]{InterpolationWedge6.UniqueInstance},
			new object[]{InterpolationWedge15.UniqueInstance},
			new object[]{InterpolationWedge18.UniqueInstance},
			new object[]{InterpolationPyra5.UniqueInstance},
			new object[]{InterpolationPyra13.UniqueInstance},
			//TODO: new object[]{InterpolationPyra14.UniqueInstance}
		};

		private static readonly Dictionary<IIsoparametricInterpolation3D, GenerateRandomPoints> pointGenerators =
			new Dictionary<IIsoparametricInterpolation3D, GenerateRandomPoints>()
			{
				{InterpolationTet4.UniqueInstance, GenerateRandomPointsInTetrahedron},
				{InterpolationTet10.UniqueInstance, GenerateRandomPointsInTetrahedron},
				{InterpolationHexa8.UniqueInstance, GenerateRandomPointsInCube},
				{InterpolationHexa20.UniqueInstance, GenerateRandomPointsInCube},
				{InterpolationHexa27.UniqueInstance, GenerateRandomPointsInCube},
				{InterpolationWedge6.UniqueInstance, GenerarandomPointsInWedge},
				{InterpolationWedge15.UniqueInstance, GenerarandomPointsInWedge},
				{InterpolationWedge18.UniqueInstance, GenerarandomPointsInWedge},
				{InterpolationPyra5.UniqueInstance, GenerateRandomPointsInPyramid},
				{InterpolationPyra13.UniqueInstance, GenerateRandomPointsInPyramid},
				//TODO: {InterpolationPyra14.UniqueInstance, GenerateRandomPointsInPyramid}
			};


		[Theory]
		[MemberData(nameof(interpolations))]
		private static void TestPartitionOfUnity(IIsoparametricInterpolation3D interpolation)
		{
			double tolerance = 1e-10;
			NaturalPoint[] points = pointGenerators[interpolation]();
			for (int p = 0; p < points.Length; p++)
			{
				double[] shapeFunctions = interpolation.EvaluateFunctionsAt(points[p]);
				double sum = 0.0;
				for (int f = 0; f < interpolation.NumFunctions; f++) sum += shapeFunctions[f];
				Assert.True(Utilities.AreValuesEqual(1.0, sum, tolerance));
			}
		}

		[Theory]
		[MemberData(nameof(interpolations))]
		private static void TestValuesAtNodes(IIsoparametricInterpolation3D interpolation)
		{
			double tolerance = 1e-10;
			for (int n = 0; n < interpolation.NodalNaturalCoordinates.Count; n++)
			{
				double[] shapeFunctions = interpolation.EvaluateFunctionsAt(interpolation.NodalNaturalCoordinates[n]);
				for (int f = 0; f < interpolation.NumFunctions; f++)
				{
					if (f == n) Assert.True(Utilities.AreValuesEqual(1.0, shapeFunctions[f], tolerance));
					else Assert.True(Utilities.AreValuesEqual(0.0, shapeFunctions[f], tolerance));
				}
			}
		}


		/// <summary>
		/// Generates random points in the tetrahedron.
		/// </summary>
		/// <returns></returns>
		private static NaturalPoint[] GenerateRandomPointsInTetrahedron()
		{
			var rand = new Random();
			var randomPoints = new NaturalPoint[numRandomPoints];
			for (int i = 0; i < numRandomPoints; i++)
			{
				double xi = rand.NextDouble();
				double eta = rand.NextDouble() * xi;
				double zeta = rand.NextDouble() * eta;
				randomPoints[i] = new NaturalPoint(xi, eta, zeta);
			}

			return randomPoints;
		}

		/// <summary>
		/// Generates random points in a wedge
		/// </summary>
		/// <returns></returns>
		private static NaturalPoint[] GenerarandomPointsInWedge()
		{
			var rand = new Random();
			var randomPoints = new NaturalPoint[numRandomPoints];
			for (int i = 0; i < numRandomPoints; i++)
			{
				double xi = -1 + rand.NextDouble() * 2.0;
				double eta = rand.NextDouble();
				double zeta = rand.NextDouble() * eta;
				randomPoints[i] = new NaturalPoint(xi, eta, zeta);
			}

			return randomPoints;
		}

		/// <summary>
		/// Generates random points in the parent cube.
		/// </summary>
		/// <returns></returns>
		private static NaturalPoint[] GenerateRandomPointsInCube()
		{
			var rand = new Random();
			var randomPoints = new NaturalPoint[numRandomPoints];
			for (int i = 0; i < numRandomPoints; i++)
			{
				double xi = -1 + rand.NextDouble() * 2.0;
				double eta = -1 + rand.NextDouble() * 2.0;
				double zeta = -1 + rand.NextDouble() * 2.0;
				randomPoints[i] = new NaturalPoint(xi, eta, zeta);
			}

			return randomPoints;
		}

		/// <summary>
		/// Generates random points inside a pyramid <see cref="https://people.sc.fsu.edu/~jburkardt/m_src/pyramid_grid/pyramid_grid.html"/>
		/// </summary>
		/// <returns></returns>
		private static NaturalPoint[] GenerateRandomPointsInPyramid()
		{
			var rand = new Random();
			var randomPoints = new NaturalPoint[numRandomPoints];
			for (int i = 0; i < numRandomPoints; i++)
			{
				double zeta = rand.NextDouble();
				double xi = (-1 + rand.NextDouble() * 2.0) * (1 - zeta);
				double eta = (-1 + rand.NextDouble() * 2.0) * xi;
				randomPoints[i] = new NaturalPoint(xi, eta, zeta);
			}

			return randomPoints;
		}
	}
}
