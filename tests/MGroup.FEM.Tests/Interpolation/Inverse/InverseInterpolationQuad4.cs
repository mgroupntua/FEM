using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.MSolve.Discretization.Commons;
using MGroup.MSolve.Geometry.Coordinates;
using Xunit;

namespace MGroup.FEM.Tests.Interpolation.Inverse
{
	public static class InverseInterpolationQuad4
	{
		private static readonly ValueComparer comparer = new ValueComparer(1E-10);

		/// <summary>
		/// Random shape, not too distorted.
		/// </summary>
		private static readonly IReadOnlyList<Node> nodeSet = new Node[]
		{
			new Node( id: 0, x: 0.7, y:  2.0 ),
			new Node( id: 1, x: 0.2, y:  0.3 ),
			new Node( id: 2, x: 2.0, y:  0.9 ),
			new Node( id: 3, x: 3.0, y:  2.7 )
		};

		private static bool Coincide(NaturalPoint point1, NaturalPoint point2)
			=> comparer.AreEqual(point1.Xi, point2.Xi) && comparer.AreEqual(point1.Eta, point2.Eta);

		/// <summary>
		/// Reorders the nodes such that the 1st one becomes the 2nd, the 2nd one becomes the 3rd, etc.
		/// </summary>
		/// <param name="originalOrder"></param>
		private static Node[] CycleNodes(IReadOnlyList<Node> originalOrder)
		{
			var cycled = new Node[originalOrder.Count];
			cycled[0] = originalOrder[originalOrder.Count - 1];
			for (int i = 0; i < originalOrder.Count - 1; ++i) cycled[i + 1] = originalOrder[i];
			return cycled;
		}

		/// <summary>
		/// Generates random points in the square: -1 &lt;= xi &lt; 1 , -1 &lt;= eta &lt; 1
		/// </summary>
		/// <returns></returns>
		private static NaturalPoint[] GenerateRandomPointsInSquare(int numRandomPoints)
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

		[Fact]
		private static void TestInverseMapping()
		{
			var directMapping = InterpolationQuad4.UniqueInstance;
			int numRandomPoints = 10;
			NaturalPoint[] naturalPoints = GenerateRandomPointsInSquare(numRandomPoints);
			IReadOnlyList<Node> elementNodes = nodeSet;

			for (int i = 0; i < 4; ++i)
			{
				IInverseInterpolation2D inverseMapping = directMapping.CreateInverseMappingFor(elementNodes);
				foreach (NaturalPoint originalPoint in naturalPoints)
				{
					CartesianPoint cartesianPoint = directMapping.TransformNaturalToCartesian(elementNodes, originalPoint);
					NaturalPoint remappedPoint = inverseMapping.TransformPointCartesianToNatural(cartesianPoint);
					Assert.True(Coincide(originalPoint, remappedPoint));
				}

				elementNodes = CycleNodes(elementNodes); // The next iteration will use a different node order
			}
		}
	}
}
