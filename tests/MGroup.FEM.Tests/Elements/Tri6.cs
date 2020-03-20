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
	/// Tests 6-noded triangular instances of <see cref="ContinuumElement2D"/> against Abaqus.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class Tri6
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
			new Node( id: 2, x: 4.0, y:  0.8 ),
			new Node( id: 3, x: 1.4, y:  2.3 ),
			new Node( id: 4, x: 2.6, y:  1.2 ),
			new Node( id: 5, x: 2.9, y:  2.6 )
		};

		/// <summary>
		/// To reproduce the reference solution, run the Abaqus job tri6_test0.inp and look at 
		/// TRI6_TEST0_MASS_MATRICES.mtx.
		/// </summary>
		[Fact]
		private static void TestConsistentMass0()
		{
			IQuadrature2D quadratureForMass = TriangleQuadratureSymmetricGaussian.Order4Points6;

			var materialsAtGaussPoints = new List<IContinuumMaterial2D>();
			foreach (GaussPoint gaussPoint in quadratureForMass.IntegrationPoints)
			{
				materialsAtGaussPoints.Add(material0.Clone());
			}
			var tri6 = new ContinuumElement2D(thickness, nodeSet0, InterpolationTri6.UniqueInstance,
				TriangleQuadratureSymmetricGaussian.Order2Points3, quadratureForMass,
				ExtrapolationGaussTriangular3Points.UniqueInstance,
				materialsAtGaussPoints, dynamicMaterial);

			IMatrix M = tri6.BuildConsistentMassMatrix();
			Matrix expectedM = Matrix.CreateFromArray(new double[,]
			{
				{ 11.83599396, 0, -1.829621484, 0, -1.968668415, 0, 0.702983761, 0, -6.922261694, 0, 0.623796095, 0,  },
				{ 0, 11.83599396, 0, -1.829621484, 0, -1.968668415, 0, 0.702983761, 0, -6.922261694, 0, 0.623796095,  },
				{ -1.829621484, 0, 9.926994738, 0, -1.638252825, 0, -0.618678599, 0, -0.46271906, 0, -7.392556104, 0, },
				{ 0, -1.829621484, 0, 9.926994738, 0, -1.638252825, 0, -0.618678599, 0, -0.46271906, 0, -7.392556104, },
				{ -1.968668415, 0, -1.638252825, 0, 10.74103411, 0, -7.234180772, 0, 9.35E-02, 0, -0.141678539, 0,    },
				{ 0, -1.968668415, 0, -1.638252825, 0, 10.74103411, 0, -7.234180772, 0, 9.35E-02, 0, -0.141678539,    },
				{ 0.702983761, 0, -0.618678599, 0, -7.234180772, 0, 57.86118359, 0, 28.31288493, 0, 29.25347375, 0,   },
				{ 0, 0.702983761, 0, -0.618678599, 0, -7.234180772, 0, 57.86118359, 0, 28.31288493, 0, 29.25347375,   },
				{ -6.922261694, 0, -0.46271906, 0, 9.35E-02, 0, 28.31288493, 0, 56.00999156, 0, 28.6296356, 0,        },
				{ 0, -6.922261694, 0, -0.46271906, 0, 9.35E-02, 0, 28.31288493, 0, 56.00999156, 0, 28.6296356,        },
				{ 0.623796095, 0, -7.392556104, 0, -0.141678539, 0, 29.25347375, 0, 28.6296356, 0, 58.49121809, 0,    },
				{ 0, 0.623796095, 0, -7.392556104, 0, -0.141678539, 0, 29.25347375, 0, 28.6296356, 0, 58.49121809,    }
			}); // from Abaqus
			Assert.True(M.Equals(expectedM, 1e-3));
		}

		/// <summary>
		/// To reproduce the reference solution, run the Abaqus job tri6_test0.inp and look at 
		/// TRI6_TEST0_STIFFNESS_MATRICES.mtx.
		/// </summary>
		[Fact]
		private static void TestStiffness0()
		{
			var factory = new ContinuumElement2DFactory(thickness, material0, dynamicMaterial);
			ContinuumElement2D tri6 = factory.CreateElement(CellType.Tri6, nodeSet0);
			IMatrix K = tri6.BuildStiffnessMatrix();
			double[,] expectedK = new double[,]
			{
				{ 51173.75966, -20778.82606, 4273.968122, 6461.400457, 3161.87986, -13008.35802, -16897.34468, -27025.3615, 9688.11919, 8859.762175, -51400.38215, 45491.38295      },
				{ -20778.82606, 122356.4459, 4538.323534, 34691.5801, -11085.2811, 5376.486541, -19333.05381, -142941.0483, 8859.762175, 6529.450147, 37799.07526, -26012.91439     },
				{ 4273.968122, 4538.323534, 88582.98224, 45240.07515, 32066.51649, 7151.333545, -34941.29312, -27750.00033, -109199.8665, -39522.63295, 19217.69278, 10342.90105    },
				{ 6461.400457, 34691.5801, 45240.07515, 104048.7171, 5228.256622, 3589.586543, -35442.30802, -149252.8235, -31830.32526, 1377.198943, 10342.90105, 5545.740823      },
				{ 3161.87986, -11085.2811, 32066.51649, 5228.256622, 112494.273, -16055.03055, 1567.844178, 1895.723659, -172177.9559, -18317.78738, 22887.44239, 38334.11875       },
				{ -13008.35802, 5376.486541, 7151.333545, 3589.586543, -16055.03055, 44329.93416, 1895.723659, 5503.165063, -26010.09507, -38438.07765, 46026.42644, -20361.09466   },
				{ -16897.34468, -19333.05381, -34941.29312, -35442.30802, 1567.844178, 1895.723659, 363079.0153, 24675.86222, -73599.74689, 108516.9545, -239208.4748, -80313.17851 },
				{ -27025.3615, -142941.0483, -27750.00033, -149252.8235, 1895.723659, 5503.165063, 24675.86222, 449052.1476, 108516.9545, -110192.2562, -80313.17851, -52169.18472  },
				{ 9688.11919, 8859.762175, -109199.8665, -31830.32526, -172177.9559, -26010.09507, -73599.74689, 108516.9545, 485631.597, 60752.92074, -140342.1468, -120289.217    },
				{ 8859.762175, 6529.450147, -39522.63295, 1377.198943, -18317.78738, -38438.07765, 108516.9545, -110192.2562, 60752.92074, 406726.0715, -120289.217, -266002.3867   },
				{ -51400.38215, 37799.07526, 19217.69278, 10342.90105, 22887.44239, 46026.42644, -239208.4748, -80313.17851, -140342.1468, -120289.217, 388845.8685, 106433.9928    },
				{ 45491.38295, -26012.91439, 10342.90105, 5545.740823, 38334.11875, -20361.09466, -80313.17851, -52169.18472, -120289.217, -266002.3867, 106433.9928, 358999.8396   }
			}; // from Abaqus
			Assert.True(K.Equals(Matrix.CreateFromArray(expectedK), 1e-9));
		}

		/// <summary>
		/// To reproduce the reference solution, run the Abaqus job tri6_test0.inp and look at tri6_tes0.dat.
		/// nodes).
		/// </summary>
		[Fact]
		public static void TestStrainsStresses0()
		{
			var factory = new ContinuumElement2DFactory(thickness, material0, null);
			ContinuumElement2D tri6 = factory.CreateElement(CellType.Tri6, nodeSet0);

			// Abaqus results
			double[] displacements =
			{
				0.0, 0.0,                   // Node 1
                0.0, 0.0,                   // Node 2
                -1.7903E-02, -4.2697E-02,   // Node 3
                1.5326E-03, -1.0672E-03,    // Node 4
                -8.8920E-03, -1.1654E-02,   // Node 5
                3.5728E-03, -1.2996E-02     // Node 6
            };

			double[][] expectedStrainsAtGPs =
			{
				new double[] {  3.3735E-03,  8.8893E-04, -4.6302E-03 },  // Gauss point 1
                new double[] { -4.4781E-03,  5.7927E-04,  3.6832E-04 },  // Gauss point 2
                new double[] { -1.4936E-03,  3.0967E-03, -7.1318E-03 }   // Gauss point 3
            };
			double[][] expectedStressesAtGPs =
			{
				new double[] {  840.0,  438.7, -374.00 },   // Gauss point 1
                new double[] { -993.3, -176.3,   29.75 },   // Gauss point 2
                new double[] { -130.3,  611.2, -576.00 }    // Gauss point 3
            };

			// The order of the nodes is the same (at least in this job)
			double[][] expectedStrainsAtNodes =
			{
				new double[] {  7.6130E-03,  2.5621E-04, -5.4625E-03 }, // Node 1
                new double[] { -8.0901E-03, -3.6310E-04,  4.5346E-03 }, // Node 2
                new double[] { -2.1212E-03,  4.6718E-03, -1.0466E-02 }, // Node 3
                new double[] { -2.3853E-04, -5.3446E-05, -4.6399E-04 }, // Node 4
                new double[] { -5.1056E-03,  2.1544E-03, -2.9656E-03 }, // Node 5
                new double[] {  2.7459E-03,  2.4640E-03, -7.9641E-03 }  // Node 6
            };
			double[][] expectedStressesAtNodes =
			{
				new double[] {  1775.00,  586.20, -441.20 }, // Node 1
                new double[] { -1892.00, -643.90,  366.30 }, // Node 2
                new double[] {  -166.10,  931.30, -845.30 }, // Node 3
                new double[] {   -58.74,  -28.85,  -37.48 }, // Node 4
                new double[] { -1029.00,  143.70, -239.50 }, // Node 5
                new double[] {   804.30,  758.70, -643.30 }  // Node 6
            };

			(IReadOnlyList<double[]> strainsAtGPs, IReadOnlyList<double[]> stressesAtGPs) =
				tri6.UpdateStrainsStressesAtGaussPoints(displacements);
			IReadOnlyList<double[]> strainsAtNodes =
				tri6.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(strainsAtGPs, tri6.Interpolation);
			IReadOnlyList<double[]> stressesAtNodes =
				tri6.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(stressesAtGPs, tri6.Interpolation);

			Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtGPs, strainsAtGPs, 1e-3));
			Assert.True(Utilities.AreTensorsEqual(expectedStressesAtGPs, stressesAtGPs, 1e-3));
			Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtNodes, strainsAtNodes, 1e-3));
			Assert.True(Utilities.AreTensorsEqual(expectedStressesAtNodes, stressesAtNodes, 1e-3));
		}
	}
}
