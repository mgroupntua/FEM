using System.Collections.Generic;
using MGroup.MSolve.Discretization.Entities;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.FEM.Thermal.Isoparametric;
using MGroup.Constitutive.Thermal;
using Xunit;

namespace MGroup.FEM.Thermal.Tests.Elements
{
	public static class ThermalTri3
	{
		private static double thickness = 1.0;
		private static double thermalConductivity = 1.0;
		private static double density = 12.0;
		private static double specialHeatCoeff = 2.0;

		/// <summary>
		/// Random shape, not too distorted.
		/// </summary>
		private static readonly IReadOnlyList<Node> nodeSet0 = new Node[]
		{
			new Node( id: 0, x: 0.0, y:  0.0 ),
			new Node( id: 1, x: 1.0, y:  0.0 ),
			new Node( id: 2, x: 0.0, y:  1.0 )
		};

		[Fact]
		private static void TestCapacity()
		{
			var factory = new ThermalElement2DFactory(thickness, new ThermalProperties(density, specialHeatCoeff, thermalConductivity));
			ThermalElement2D element = factory.CreateElement(CellType.Tri3, nodeSet0);
			IMatrix M = element.BuildCapacityMatrix();

			var expectedM = Matrix.CreateFromArray(new double[,]
			{
				{2, 1, 1 },
				{1, 2, 1 },
				{1, 1, 2 }

			});

			Assert.True(expectedM.Equals(M, 1e-10));
		}

		[Fact]
		private static void TestConductivity()
		{
			var factory = new ThermalElement2DFactory(thickness, new ThermalProperties(density, specialHeatCoeff, thermalConductivity));
			ThermalElement2D element = factory.CreateElement(CellType.Tri3, nodeSet0);
			IMatrix K = element.BuildConductivityMatrix();

			var expectedK = Matrix.CreateFromArray(new double[,]
			{
				{ 1.0, -0.5, -0.5 },
				{ -0.5, 0.5,  0.0 },
				{ -0.5, 0.0,  0.5 }
			});

			Assert.True(expectedK.Equals(K, 1e-10));
		}
	}
}
