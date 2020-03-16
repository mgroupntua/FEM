using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Mesh;
using Xunit;

namespace MGroup.FEM.Tests.Elements
{
	using Constitutive.Thermal;
	using Thermal.Elements;

	public static class ThermalQuad4
	{
		private static double thickness = 1.0;
		private static double thermalConductivity = 1.0;
		private static double density = 1.0;
		private static double specialHeatCoeff = 1.0;

		/// <summary>
		/// Random shape, not too distorted.
		/// </summary>
		private static readonly IReadOnlyList<Node> nodeSet0 = new Node[]
		{
			new Node( id: 0, x: 0.0, y:  0.0 ),
			new Node( id: 1, x: 1.0, y:  0.0 ),
			new Node( id: 2, x: 1.0, y:  1.0 ),
			new Node( id: 3, x: 0.0, y:  1.0 )
		};

		[Fact]
		private static void TestCapacity()
		{
			var factory = new ThermalElement2DFactory(thickness, new ThermalMaterial(density, specialHeatCoeff, thermalConductivity));
			ThermalElement2D element = factory.CreateElement(CellType.Quad4, nodeSet0);
			IMatrix M = element.BuildCapacityMatrix();

			var expectedM = Matrix.CreateFromArray(new double[,]
			{
				{4, 2, 1, 2 },
				{2, 4, 2, 1 },
				{1, 2, 4, 2 },
				{2, 1, 2, 4 }
			});
			expectedM.ScaleIntoThis(density * specialHeatCoeff / 36);

			Assert.True(expectedM.Equals(M, 1e-10));
		}

		[Fact]
		private static void TestConductivity()
		{
			var factory = new ThermalElement2DFactory(thickness, new ThermalMaterial(density, specialHeatCoeff, thermalConductivity));
			ThermalElement2D element = factory.CreateElement(CellType.Quad4, nodeSet0);
			IMatrix K = element.BuildConductivityMatrix();

			var expectedK = Matrix.CreateFromArray(new double[,]
			{
				{ 2.0/3, -1.0/6, -1.0/3, -1/6.0},
				{-1/6.0,  2.0/3, -1.0/6, -1.0/3},
				{-1.0/3, -1.0/6,  2.0/3, -1.0/6},
				{-1/6.0, -1.0/3, -1.0/6,  2.0/3}

			});

			Assert.True(expectedK.Equals(K, 1e-10));
		}
	}
}
