namespace MGroup.FEM.Structural.Tests.Elements
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;
	using DotNumerics.Optimization.TN;

	using MGroup.Constitutive.Structural.Cohesive;
	using MGroup.Constitutive.Structural.Continuum;
	using MGroup.Constitutive.Structural.Line;
	using MGroup.FEM.Structural.Embedding;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Numerics.Integration.Quadratures;

	using Xunit;

	public class CohesiveBeam3DToBeam3DTest
	{
		[Fact]
		private static void RunTest()
		{
			var nodeSet = new List<INode>()
			{
				new Node(0, -11.261562, 21.900388, 1.716324),
				new Node(1, -42.483085, -17.153164, 1.525223),
				new Node(2, -11.261562, 21.900388, 1.716324),
				new Node(3, -42.483085, -17.153164, 1.525223),
			};
			var elementNodesBeam = new List<INode>()
			{
				new Node(0, -11.261562, 21.900388, 1.716324),
				new Node(1, -42.483085, -17.153164, 1.525223),
			};

			IIsotropicContinuumMaterial3D matrixMaterial = new ElasticMaterial3D(youngModulus: 3.5, poissonRatio: 0.4);
			var K_el = 10;
			var K_pl = 1.0;
			var T_max = 0.10;
			var cohesiveMaterial = new BondSlipMaterial(K_el, K_pl, 100.0, T_max, new double[2], new double[2], 1e-3);
			
			double area = 694.77;
			double inertiaY = 100.18;
			double inertiaZ = 100.18;
			double torsionalInertia = 68.77;
			var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, 0, 0);
			
			double mi = 8.0;
			double ni = 8.0;
			double thickness_CNT = 0.34;
			double a = 0.241;
			double diameter_CNT = (a / Math.PI) * Math.Sqrt(Math.Pow(ni, 2) + ni * mi + Math.Pow(mi, 2));
			double radius_CNT = diameter_CNT / 2.0;
			double radius_CNT_outer = radius_CNT + (thickness_CNT / 2);
			double CntPerimeter = 2.0 * Math.PI * radius_CNT_outer;
			var element = new CohesiveBeam3DToBeam3D(nodeSet, cohesiveMaterial, GaussLegendre1D.GetQuadratureWithOrder(2), elementNodesBeam, elementNodesBeam, matrixMaterial, 1, beamSection, CntPerimeter);

			var elementStiffnessMatrix = element.StiffnessMatrix();

			double[] localDisplacements =
			{
				-0.22817359334790202, -0.093157552807095936, -0.00032526368182175715, 1.3688506028227088E-05, 0.000259653566909844,
				4.9007290948375906E-05, -0.8453951256523183, 0.067674668134558316, -0.0095639398099483467, -3.1091507972934296E-05,
				0.00035233238898999089, -0.00012733831428290391, -0.23387978688156119, -0.10029516508465312, -0.00036019020640142024,
				1.3688021032242603E-05, 0.00025965431209030315, 4.9005687219101045E-05, -0.83968893211865792, 0.074812280412115614,
				-0.0095290132853686309, -3.1091022976958875E-05, 0.00035233164380952571, -0.00012733671055360516,
			};

			Tuple<double[], double[]> results = element.CalculateResponse(localDisplacements);

			var forces = element.CalculateResponseIntegral();


			var expectedStrains = new double[] { 0.0052756089927531976, 6.11479059092248E-05, 1.3659736733557528E-06 };
			var expectedStresses = new double[] { 0.052756089927531974, 0.000611479059092248, 0.00013659736733557528 };
			var expectedForces = new double[]
			{
				2.1064736285205248, 2.6349929239668777, 0.0050331374077453254, 2.4373019262327917E-05, -3.7448346514700273E-05,
				8.0593913678409413E-05, -2.1064736285244603, -2.6349929239668777, -0.0050331374078813225, -2.4373018806862645E-05,
				3.7448346817095162E-05, -8.0593914878794527E-05, -2.1064736285205248, -2.6349929239668777, -0.0050331374077453254,
				-2.4373019262327917E-05, 3.7448346514700273E-05, -8.0593913678409413E-05, 2.1064736285244603, 2.6349929239668777,
				0.0050331374078813225, 2.4373018806862645E-05, -3.7448346817095162E-05, 8.0593914878794527E-05,
			};

			Assert.Equal(results.Item1, expectedStrains);
			Assert.Equal(results.Item2, expectedStresses);
			Assert.Equal(forces, expectedForces);
		}
	}
}
