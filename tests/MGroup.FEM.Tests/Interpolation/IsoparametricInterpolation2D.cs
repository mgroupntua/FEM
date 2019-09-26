using System;
using System.Collections.Generic;
using MGroup.FEM.Interpolation;
using MGroup.MSolve.Geometry.Coordinates;
using Xunit;

namespace MGroup.FEM.Tests.Interpolation
{
	/// <summary>
	/// Unit testing for implentations of <see cref="IIsoparametricInterpolation2D"/>.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class IsoparametricInterpolation2D
	{
		private const int numRandomPoints = 10;
		private delegate NaturalPoint[] GenerateRandomPoints();

		public static readonly IEnumerable<object[]> interpolations = new List<object[]>
		{
			new object[] { InterpolationQuad4.UniqueInstance },
			new object[] { InterpolationQuad8.UniqueInstance },
			new object[] { InterpolationQuad9.UniqueInstance },
			new object[] { InterpolationTri3.UniqueInstance },
			new object[] { InterpolationTri6.UniqueInstance }
		};

		private static readonly Dictionary<IIsoparametricInterpolation2D, GenerateRandomPoints> pointGenerators =
			new Dictionary<IIsoparametricInterpolation2D, GenerateRandomPoints>
			{
				{ InterpolationQuad4.UniqueInstance, GenerateRandomPointsInSquare },
				{ InterpolationQuad8.UniqueInstance, GenerateRandomPointsInSquare },
				{ InterpolationQuad9.UniqueInstance, GenerateRandomPointsInSquare },
				{ InterpolationTri3.UniqueInstance, GenerateRandomPointsInTriangle },
				{ InterpolationTri6.UniqueInstance, GenerateRandomPointsInTriangle }
			};

		[Theory]
		[MemberData(nameof(interpolations))]
		private static void TestPartitionOfUnity(IIsoparametricInterpolation2D interpolation)
		{
			double tolerance = 1e-10;
			NaturalPoint[] points = pointGenerators[interpolation]();
			for (int p = 0; p < points.Length; ++p)
			{
				double[] shapeFuncs = interpolation.EvaluateFunctionsAt(points[p]);
				double sum = 0.0;
				for (int f = 0; f < interpolation.NumFunctions; ++f) sum += shapeFuncs[f];
				Assert.True(Utilities.AreValuesEqual(1.0, sum, tolerance));
			}
		}

		[Theory]
		[MemberData(nameof(interpolations))]
		private static void TestValuesAtNodes(IIsoparametricInterpolation2D interpolation)
		{
			double tolerance = 1e-10;
			for (int n = 0; n < interpolation.NodalNaturalCoordinates.Count; ++n)
			{
				double[] shapeFuncs = interpolation.EvaluateFunctionsAt(interpolation.NodalNaturalCoordinates[n]);
				for (int f = 0; f < interpolation.NumFunctions; ++f)
				{
					if (f == n) Assert.True(Utilities.AreValuesEqual(1.0, shapeFuncs[f], tolerance));
					else Assert.True(Utilities.AreValuesEqual(0.0, shapeFuncs[f], tolerance));
				}
			}
		}

		/// <summary>
		/// Generates random point in the square: -1 &lt;= xi &lt; 1 , -1 &lt;= eta &lt; 1
		/// </summary>
		/// <returns></returns>
		private static NaturalPoint[] GenerateRandomPointsInSquare()
		{
			var rand = new Random();
			var randomPoints = new NaturalPoint[numRandomPoints];
			for (int i = 0; i < numRandomPoints; ++i)
			{
				double xi = -1.0 + rand.NextDouble() * 2.0;
				double eta = -1.0 + rand.NextDouble() * 2.0;
				randomPoints[i] = new NaturalPoint(xi, eta);
			}
			return randomPoints;
		}

		/// <summary>
		/// Generates random point in the triangle: 0 &lt;= xi &lt; 1 , 0 &lt;= eta &lt; xi
		/// </summary>
		/// <returns></returns>
		private static NaturalPoint[] GenerateRandomPointsInTriangle()
		{
			var rand = new Random();
			var randomPoints = new NaturalPoint[numRandomPoints];
			for (int i = 0; i < numRandomPoints; ++i)
			{
				double xi = rand.NextDouble();
				double eta = rand.NextDouble() * xi;
				randomPoints[i] = new NaturalPoint(xi, eta);
			}
			return randomPoints;
		}
	}
}
