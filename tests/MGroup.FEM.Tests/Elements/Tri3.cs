using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation;
using MGroup.FEM.Interpolation.GaussPointExtrapolation;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Integration;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Mesh;
using Xunit;

//TODO: Add tests for wrong node orders, too distorted shapes, etc.
//TODO: Add tests presented in https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf
namespace MGroup.FEM.Tests.Elements
{
	using Constitutive.Structural;
	using Constitutive.Structural.PlanarElements;
	using MSolve.Constitutive;
	using Structural.Elements;

	/// <summary>
	/// Tests 3-noded triangular instances of <see cref="ContinuumElement2D"/> against Abaqus and the notes of the excellent
	/// University of Colorado at Boulder FEM course.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class Tri3
	{
		private static double thickness = 1.0;

		private static readonly ElasticMaterial2D material0 = new ElasticMaterial2D(StressState2D.PlaneStress)
		{
			YoungModulus = 2.1e5,
			PoissonRatio = 0.3
		};

		private static readonly DynamicMaterial dynamicMaterial = new DynamicMaterial(78.5, 0, 0);

		/// <summary>
		/// Random shape, not too distorted.
		/// </summary>
		private static readonly IReadOnlyList<Node> nodeSet0 = new Node[]
		{
			new Node( id: 0, x: 1.5, y:  3.8 ),
			new Node( id: 1, x: 1.0, y:  1.0 ),
			new Node( id: 2, x: 4.0, y:  0.8 )
		};

		/// <summary>
		/// Right triangle.
		/// </summary>
		private static readonly IReadOnlyList<Node> nodeSet1 = new Node[]
		{
			new Node( id: 0, x: 0.0, y:  0.0 ),
			new Node( id: 1, x: 1.0, y:  0.0 ),
			new Node( id: 2, x: 0.0, y:  1.0 )
		};

		/// <summary>
		/// Consistent mass with as many integration points as there are nodes.
		/// The reference solution can be found from: http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf, (eq 31.27) 
		/// </summary>
		[Fact]
		private static void TestConsistentMass0()
		{
			TestConsistentMassParametric(nodeSet0, TriangleQuadratureSymmetricGaussian.Order2Points3, false);
		}

		/// <summary>
		/// Consistent mass with less integration points than nodes (just for testing, do not do that).
		/// The reference solution can be found from: http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf, (eq 31.27) 
		/// </summary>
		[Fact]
		private static void TestConsistentMass1()
		{
			TestConsistentMassParametric(nodeSet1, TriangleQuadratureSymmetricGaussian.Order1Point1, true);
		}

		private static void TestConsistentMassParametric(IReadOnlyList<Node> nodeSet, IQuadrature2D quadratureForMass,
			bool reducedQuadrature)
		{
			var materialsAtGaussPoints = new List<IContinuumMaterial2D>();
			foreach (GaussPoint gaussPoint in quadratureForMass.IntegrationPoints)
			{
				materialsAtGaussPoints.Add(material0.Clone());
			}
			var tri3 = new ContinuumElement2D(thickness, nodeSet, InterpolationTri3.UniqueInstance,
				TriangleQuadratureSymmetricGaussian.Order1Point1, quadratureForMass,
				ExtrapolationGaussTriangular1Point.UniqueInstance,
				materialsAtGaussPoints, dynamicMaterial);
			IMatrix M = tri3.BuildConsistentMassMatrix();

			// Reference: http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf, (eq 31.27) 
			Matrix expectedM = Matrix.CreateFromArray(new double[,]
			{
				{ 2, 0, 1, 0, 1, 0 },
				{ 0, 2, 0, 1, 0, 1 },
				{ 1, 0, 2, 0, 1, 0 },
				{ 0, 1, 0, 2, 0, 1 },
				{ 1, 0, 1, 0, 2, 0 },
				{ 0, 1, 0, 1, 0, 2 }
			});
			double area = CalcTriangleArea(nodeSet);
			double scalar = dynamicMaterial.Density * thickness * area / 12.0;
			expectedM.ScaleIntoThis(scalar);

			if (reducedQuadrature) Assert.False(M.Equals(expectedM, 1e-10));
			else Assert.True(M.Equals(expectedM, 1e-10));

		}

		/// <summary>
		/// To reproduce the reference solution, run the Abaqus job tri3_test0.inp and look at 
		/// TRI3_TEST0_STIFFNESS_MATRICES.mtx. 
		/// </summary>
		[Fact]
		private static void TestStiffness0()
		{
			var factory = new ContinuumElement2DFactory(thickness, material0, dynamicMaterial);
			ContinuumElement2D tri3 = factory.CreateElement(CellType.Tri3, nodeSet0);
			IMatrix K = tri3.BuildStiffnessMatrix();
			double[,] expectedK = new double[,]
			{
				{ 43303.16742081400, 5294.11764705880, -43778.28054298600, -44796.38009049800, 475.11312217194, 39502.26244343900 },
				{ 5294.11764705880, 122361.99095023000, -39027.14932126700, -104660.63348416000, 33733.03167420800, -17701.35746606300 },
				{ -43778.28054298600, -39027.14932126700, 151866.51583710000, 66176.47058823500, -108088.23529412000, -27149.32126696800 },
				{ -44796.38009049800, -104660.63348416000, 66176.47058823500, 127601.80995475000, -21380.09049773800, -22941.17647058800 },
				{ 475.11312217194, 33733.03167420800, -108088.23529412000, -21380.09049773800, 107613.12217195000, -12352.94117647100 },
				{ 39502.26244343900, -17701.35746606300, -27149.32126696800, -22941.17647058800, -12352.94117647100, 40642.53393665200 }
			}; // from Abaqus
			Assert.True(K.Equals(Matrix.CreateFromArray(expectedK), 1e-10));
		}

		/// <summary>
		/// To reproduce the reference solution, run the Abaqus job tri3_test0.inp and look at tri3_tes0.dat. Tri3 is a constant 
		/// strain element(CST), therefore the strains / stresses are constant throughout the element (integration points and 
		/// nodes).
		/// </summary>
		[Fact]
		public static void TestStrainsStresses0()
		{
			var factory = new ContinuumElement2DFactory(thickness, material0, null);
			ContinuumElement2D tri3 = factory.CreateElement(CellType.Tri3, nodeSet0);

			// Abaqus results
			double[] displacements =
			{
				0.0, 0.0,                    // Node 1
                0.0, 0.0,                    // Node 2
                -1.4632E-03, -1.2747E-02     // Node 3  
            };

			// There is only 1 Gauss point.
			double[][] expectedStrainsAtGPs =
			{
				new double[] { -4.8201E-04,  7.4983E-04, -4.1130E-03 }  // Gauss point 1
            };
			double[][] expectedStressesAtGPs =
			{
				new double[] { -59.32, 139.7, -332.2 }  // Gauss point 1
            };

			// The order of the nodes is the same (at least in this job)
			double[][] expectedStrainsAtNodes =
			{
				new double[] { -4.8201E-04,  7.4983E-04, -4.1130E-03 },  // Node 1
                new double[] { -4.8201E-04,  7.4983E-04, -4.1130E-03 },  // Node 2
                new double[] { -4.8201E-04,  7.4983E-04, -4.1130E-03 }   // Node 3
            };
			double[][] expectedStressesAtNodes =
			{
				new double[] { -59.32, 139.7, -332.2 },  // Node 1
                new double[] { -59.32, 139.7, -332.2 },  // Node 2
                new double[] { -59.32, 139.7, -332.2 }   // Node 3
            };

			(IReadOnlyList<double[]> strainsAtGPs, IReadOnlyList<double[]> stressesAtGPs) =
				tri3.UpdateStrainsStressesAtGaussPoints(displacements);
			IReadOnlyList<double[]> strainsAtNodes =
				tri3.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(strainsAtGPs, tri3.Interpolation);
			IReadOnlyList<double[]> stressesAtNodes =
				tri3.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(stressesAtGPs, tri3.Interpolation);

			Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtGPs, strainsAtGPs, 1e-4));
			Assert.True(Utilities.AreTensorsEqual(expectedStressesAtGPs, stressesAtGPs, 1e-3));
			Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtNodes, strainsAtNodes, 1e-4));
			Assert.True(Utilities.AreTensorsEqual(expectedStressesAtNodes, stressesAtNodes, 1e-3));
		}

		private static double CalcTriangleArea(IReadOnlyList<Node> nodes)
		{
			return (nodes[0].X * (nodes[1].Y - nodes[2].Y)
				+ nodes[1].X * (nodes[2].Y - nodes[0].Y)
				+ nodes[2].X * (nodes[0].Y - nodes[1].Y)) / 2.0;
		}
	}
}
