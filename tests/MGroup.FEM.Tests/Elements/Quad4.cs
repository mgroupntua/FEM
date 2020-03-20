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
//TODO: Also add the tests presented in https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch23.d/IFEM.Ch23.pdf
namespace MGroup.FEM.Tests.Elements
{
	using Constitutive.Structural;
	using Constitutive.Structural.PlanarElements;
	using MSolve.Constitutive;
	using Structural.Elements;

	/// <summary>
	/// Tests 4-noded quadrilateral instances of <see cref="ContinuumElement2D"/> against Abaqus and the notes of the excellent 
	/// University of Colorado at Boulder FEM course.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class Quad4
	{
		#region reproducible tests
		private static double thickness = 1.0;

		private static readonly ElasticMaterial2D material0 = new ElasticMaterial2D(StressState2D.PlaneStress)
		{
			YoungModulus = 2.1e5,
			PoissonRatio = 0.3,
		};

		private static readonly DynamicMaterial dynamicMaterial = new DynamicMaterial(78.5, 0, 0);

		/// <summary>
		/// Random shape, not too distorted.
		/// </summary>
		private static readonly IReadOnlyList<Node> nodeSet0 = new Node[]
		{
			new Node( id: 0, x: 0.7, y:  2.0 ),
			new Node( id: 1, x: 0.2, y:  0.3 ),
			new Node( id: 2, x: 2.0, y:  0.9 ),
			new Node( id: 3, x: 3.0, y:  2.7 )
		};

		/// <summary>
		/// Rectangle.
		/// </summary>
		private static readonly IReadOnlyList<Node> nodeSet1 = new Node[]
		{
			new Node( id: 0, x:  0.0, y:   0.0 ),
			new Node( id: 1, x: 20.0, y:   0.0 ),
			new Node( id: 2, x: 20.0, y:  10.0 ),
			new Node( id: 3, x:  0.0, y:  10.0 )
		};

		/// <summary>
		/// Consistent mass matrix with less integration points than nodes (just for testing, do not do that).
		/// The reference solution can be found from http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf, (eq 31.29a) 
		/// </summary>
		[Fact]
		private static void TestConsistentMass0()
		{
			// reduced integration rule - bad idea
			IQuadrature2D quadratureForMass = GaussLegendre2D.GetQuadratureWithOrder(1, 1);
			var materialsAtGaussPoints = new List<IContinuumMaterial2D>();
			foreach (GaussPoint gaussPoint in quadratureForMass.IntegrationPoints)
			{
				materialsAtGaussPoints.Add(material0.Clone());
			}
			var quad4 = new ContinuumElement2D(thickness, nodeSet1, InterpolationQuad4.UniqueInstance,
				GaussLegendre2D.GetQuadratureWithOrder(2, 2), quadratureForMass,
				ExtrapolationGaussLegendre2x2.UniqueInstance, materialsAtGaussPoints, dynamicMaterial);
			IMatrix M = quad4.BuildConsistentMassMatrix();

			Matrix expectedM = Matrix.CreateFromArray(new double[,]
			{
				{ 1, 0, 1, 0, 1, 0, 1, 0 },
				{ 0, 1, 0, 1, 0, 1, 0, 1 },
				{ 1, 0, 1, 0, 1, 0, 1, 0 },
				{ 0, 1, 0, 1, 0, 1, 0, 1 },
				{ 1, 0, 1, 0, 1, 0, 1, 0 },
				{ 0, 1, 0, 1, 0, 1, 0, 1 },
				{ 1, 0, 1, 0, 1, 0, 1, 0 },
				{ 0, 1, 0, 1, 0, 1, 0, 1 }
			});
			double lengthX = nodeSet1[1].X - nodeSet1[0].X;
			double lengthY = nodeSet1[2].Y - nodeSet1[1].Y;

			// For some reason , only half the thickness is used for Quad4 elements, as shown in Fig. 31.9. Therefore the 
			// coefficient 1/16 (full thickness) became 1/32 (half thickness). Here 1/16 is used.
			double scalar = dynamicMaterial.Density * thickness * lengthX * lengthY / 16.0;
			expectedM.ScaleIntoThis(scalar);

			Assert.True(M.Equals(expectedM, 1e-10));
		}

		/// <summary>
		/// Consistent mass matrix with as many integration points as there are nodes. 
		/// The reference solution can be found from http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf, (eq 31.29b).
		/// </summary>
		[Fact]
		private static void TestConsistentMass1()
		{
			// full integration rule
			IQuadrature2D quadratureForMass = GaussLegendre2D.GetQuadratureWithOrder(2, 2);
			var materialsAtGaussPoints = new List<IContinuumMaterial2D>();
			foreach (GaussPoint gaussPoint in quadratureForMass.IntegrationPoints)
			{
				materialsAtGaussPoints.Add(material0.Clone());
			}
			var quad4 = new ContinuumElement2D(thickness, nodeSet1, InterpolationQuad4.UniqueInstance,
				GaussLegendre2D.GetQuadratureWithOrder(2, 2), quadratureForMass,
				ExtrapolationGaussLegendre2x2.UniqueInstance, materialsAtGaussPoints, dynamicMaterial);
			IMatrix M = quad4.BuildConsistentMassMatrix();

			Matrix expectedM = Matrix.CreateFromArray(new double[,]
			{
				{ 4, 0, 2, 0, 1, 0, 2, 0 },
				{ 0, 4, 0, 2, 0, 1, 0, 2 },
				{ 2, 0, 4, 0, 2, 0, 1, 0 },
				{ 0, 2, 0, 4, 0, 2, 0, 1 },
				{ 1, 0, 2, 0, 4, 0, 2, 0 },
				{ 0, 1, 0, 2, 0, 4, 0, 2 },
				{ 2, 0, 1, 0, 2, 0, 4, 0 },
				{ 0, 2, 0, 1, 0, 2, 0, 4 }
			});
			double lengthX = nodeSet1[1].X - nodeSet1[0].X;
			double lengthY = nodeSet1[2].Y - nodeSet1[1].Y;
			// For some reason , only half the thickness is used for Quad4 elements, as shown in Fig. 31.9. Therefore the 
			// coefficient 1/36 (full thickness) became 1/72 (half thickness). Here 1/36 is used.
			double scalar = dynamicMaterial.Density * thickness * lengthX * lengthY / 36.0;
			expectedM.ScaleIntoThis(scalar);

			Assert.True(M.Equals(expectedM, 1e-10));
		}

		/// <summary>
		/// To reproduce the reference solution, run the Abaqus job quad4_test0.inp and look at 
		/// QUAD4_TEST0_STIFFNESS_MATRICES.mtx
		/// </summary>
		[Fact]
		private static void TestStiffness0()
		{
			var factory = new ContinuumElement2DFactory(thickness, material0, dynamicMaterial);
			ContinuumElement2D quad4 = factory.CreateElement(CellType.Quad4, nodeSet0);
			IMatrix K = quad4.BuildStiffnessMatrix();
			double[,] expectedK = new double[,]
			{
				{ 181603.19122884000, -89089.52288016200, -4991.10465483030, 7519.91442892909, -126789.16986973000, 70773.21914644800, -49822.91670427800, 10796.38930478500 },
				{ -89089.52288016200, 210532.99196037000, 13289.14519816000, -85869.32985046300, 70773.21914644800, -146868.68100448000, 5027.15853555400, 22205.01889457400 },
				{ -4991.10465483030, 13289.14519816000, 73155.65851293700, 4056.21872780570, -66433.22620703800, 10577.55360637600, -1731.32765106900, -27922.91753234200 },
				{ 7519.91442892909, -85869.32985046300, 4056.21872780570, 85360.61861633800, 16346.78437560700, 2912.80345339600, -27922.91753234200, -2404.09221927070 },
				{ -126789.16986973000, 70773.21914644800, -66433.22620703800, 16346.78437560700, 200705.04715702000, -95472.47721160800, -7482.65108024430, 8352.47368955240 },
				{ 70773.21914644800, -146868.68100448000, 10577.55360637600, 2912.80345339600, -95472.47721160800, 232719.03971773000, 14121.70445878300, -88763.16216664000 },
				{ -49822.91670427800, 5027.15853555400, -1731.32765106900, -27922.91753234200, -7482.65108024430, 14121.70445878300, 59036.89543559100, 8774.05453800470 },
				{ 10796.38930478500, 22205.01889457400, -27922.91753234200, -2404.09221927070, 8352.47368955240, -88763.16216664000, 8774.05453800470, 68962.23549133600 }
			}; // from Abaqus
			Assert.True(K.Equals(Matrix.CreateFromArray(expectedK), 1e-10));
		}

		/// <summary>
		/// To reproduce the reference solution, run the Abaqus job quad4_test0.inp and look at quad4_tes0.dat
		/// </summary>
		[Fact]
		public static void TestStrainsStresses0()
		{
			var factory = new ContinuumElement2DFactory(thickness, material0, null);
			ContinuumElement2D quad4 = factory.CreateElement(CellType.Quad4, nodeSet0);

			// Abaqus results
			double[] displacements =
			{
				0.0, 0.0,                   // Node 1
                0.0, 0.0,                   // Node 2
                -3.1797E-03, -8.6560E-03,   // Node 3  
                4.4274E-03, -1.8570E-02     // Node 4
            };

			// The Gauss points in Abaqus have the same order as in GaussLegendre2D.Order2x2: Xi major, Eta minor
			double[][] expectedStrainsAtGPs =
			{
				new double[] {  1.1178E-03,  1.5988E-03, -7.4619E-03 },  // Gauss point 1
                new double[] { -1.2757E-03,  8.6753E-04, -4.5415E-03 },  // Gauss point 2
                new double[] {  2.8614E-04, -7.3510E-04, -4.0517E-03 },  // Gauss point 3
                new double[] { -2.3014E-03, -1.8646E-03, -5.0422E-04 }   // Gauss point 4
            };
			double[][] expectedStressesAtGPs =
			{
				new double[] {  368.6,  446.3, -602.7 },  // Gauss point 1
                new double[] { -234.3,  111.9, -366.8 },  // Gauss point 2
                new double[] {  15.14, -149.8, -327.3 },  // Gauss point 3
                new double[] { -660.2, -589.6,  -40.73 }  // Gauss point 4
            };

			// However the order of nodes is different, since in that Abaqus job the preprocessor numbered them differently from   
			// the counter-clockwise order of Quad4 nodes.
			double[][] expectedStrainsAtNodes =
			{
				new double[] {  2.2723E-03,  2.6674E-03, -9.6950E-03 },  // Node 1
                new double[] { -1.7504E-03,  1.6532E-03, -5.0343E-03 },  // Node 2
                new double[] { -3.6499E-03, -3.3314E-03,  2.3560E-03 },  // Node 3
                new double[] {  9.5484E-04, -1.1226E-03, -4.1860E-03 }   // Node 4
            };
			double[][] expectedStressesAtNodes =
			{
				new double[] {   709.0,   772.9, -783.1 },  // Node 1
                new double[] {  -289.5,   260.3, -406.6 },  // Node 2
                new double[] { -1073.0, -1021.0,  190.3 },  // Node 3
                new double[] {   142.6,  -193.0, -338.1 }   // Node 4
            };

			(IReadOnlyList<double[]> strainsAtGPs, IReadOnlyList<double[]> stressesAtGPs) =
				quad4.UpdateStrainsStressesAtGaussPoints(displacements);
			IReadOnlyList<double[]> strainsAtNodes =
				quad4.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(strainsAtGPs, quad4.Interpolation);
			IReadOnlyList<double[]> stressesAtNodes =
				quad4.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(stressesAtGPs, quad4.Interpolation);

			Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtGPs, strainsAtGPs, 1e-4));
			Assert.True(Utilities.AreTensorsEqual(expectedStressesAtGPs, stressesAtGPs, 1e-3));
			Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtNodes, strainsAtNodes, 1e-4));
			Assert.True(Utilities.AreTensorsEqual(expectedStressesAtNodes, stressesAtNodes, 1e-3));
		}
		#endregion

		#region older tests (source has been long lost)
		private static readonly ElasticMaterial2D material1 = new ElasticMaterial2D(StressState2D.PlaneStress)
		{
			YoungModulus = 2e6,
			PoissonRatio = 0.3
		};

		private static readonly ElasticMaterial2D material2 = new ElasticMaterial2D(StressState2D.PlaneStress)
		{
			YoungModulus = 1.0,
			PoissonRatio = 0.25
		};

		private static readonly IReadOnlyList<Node> nodeSet2 = new Node[]
		{
			new Node( id: 0, x: -1.0, y:  -1.0 ),
			new Node( id: 1, x: +1.0, y:  -1.0 ),
			new Node( id: 2, x: +1.0, y:  +1.0 ),
			new Node( id: 3, x: -1.0, y:  +1.0 )
		};

		private static readonly IReadOnlyList<Node> nodeSet3 = new Node[]
		{
			new Node( id: 0, x: 0.2, y:  0.3 ),
			new Node( id: 1, x: 2.2, y:  1.5 ),
			new Node( id: 2, x: 3.0, y:  2.7 ),
			new Node( id: 3, x: 0.7, y:  2.0 )
		};

		[Fact]
		private static void TestStiffness1()
		{
			var factory = new ContinuumElement2DFactory(thickness, material1, dynamicMaterial);
			ContinuumElement2D quad4 = factory.CreateElement(CellType.Quad4, nodeSet2);
			IMatrix K = quad4.BuildStiffnessMatrix();
			double[,] expectedK = new double[,]
			{
				{  989010.98901099,  357142.85714286, -604395.60439561, -27472.527472528, -494505.49450549, -357142.85714286,  109890.10989011,  27472.527472528 },
				{  357142.85714286,  989010.98901099,  27472.527472528,  109890.10989011, -357142.85714286, -494505.49450549, -27472.527472528, -604395.60439561 },
				{ -604395.60439561,  27472.527472528,  989010.98901099, -357142.85714286,  109890.10989011, -27472.527472528, -494505.49450549,  357142.85714286 },
				{ -27472.527472528,  109890.10989011, -357142.85714286,  989010.98901099,  27472.527472528, -604395.60439561,  357142.85714286, -494505.49450549 },
				{ -494505.49450549, -357142.85714286,  109890.10989011,  27472.527472528,  989010.98901099,  357142.85714286, -604395.60439561, -27472.527472528 },
				{ -357142.85714286, -494505.49450549, -27472.527472528, -604395.60439561,  357142.85714286,  989010.98901099,  27472.527472528,  109890.10989011 },
				{  109890.10989011, -27472.527472528, -494505.49450549,  357142.85714286, -604395.60439561,  27472.527472528,  989010.98901099, -357142.85714286 },
				{  27472.527472528, -604395.60439561,  357142.85714286, -494505.49450549, -27472.527472528,  109890.10989011, -357142.85714286,  989010.98901099 }
			}; // from Abaqus
			Assert.True(K.Equals(Matrix.CreateFromArray(expectedK), 1e-10));
		}

		[Fact]
		private static void TestStiffness2()
		{
			var factory = new ContinuumElement2DFactory(thickness, material1, dynamicMaterial);
			ContinuumElement2D quad4 = factory.CreateElement(CellType.Quad4, nodeSet1);
			IMatrix K = quad4.BuildStiffnessMatrix();
			double[,] expectedK = new double[,]
			{
				{  879120.87912088,  357142.85714286, -109890.10989011, -27472.527472527, -439560.43956044, -357142.85714286, -329670.32967033,  27472.527472527 },
				{  357142.85714286,  1593406.5934066,  27472.527472528,  604395.60439560, -357142.85714286, -796703.29670329, -27472.527472528, -1401098.9010989 },
				{ -109890.10989011,  27472.527472528,  879120.87912088, -357142.85714286, -329670.32967033, -27472.527472527, -439560.43956044,  357142.85714286 },
				{ -27472.527472527,  604395.60439560, -357142.85714286,  1593406.5934066,  27472.527472528, -1401098.9010989,  357142.85714286, -796703.29670329 },
				{ -439560.43956044, -357142.85714286, -329670.32967033,  27472.527472528,  879120.87912088,  357142.85714286, -109890.10989011, -27472.527472527 },
				{ -357142.85714286, -796703.29670329, -27472.527472527, -1401098.9010989,  357142.85714286,  1593406.5934066,  27472.527472528,  604395.60439560 },
				{ -329670.32967033, -27472.527472528, -439560.43956044,  357142.85714286, -109890.10989011,  27472.527472528,  879120.87912088, -357142.85714286 },
				{  27472.527472527, -1401098.9010989,  357142.85714286, -796703.29670329, -27472.527472527,  604395.60439560, -357142.85714286,  1593406.5934066 }
			}; // from Abaqus
			Assert.True(K.Equals(Matrix.CreateFromArray(expectedK), 1e-10));
		}

		[Fact]
		private static void TestStiffness3()
		{
			var factory = new ContinuumElement2DFactory(thickness, material1, dynamicMaterial);
			ContinuumElement2D quad4 = factory.CreateElement(CellType.Quad4, nodeSet3);
			IMatrix K = quad4.BuildStiffnessMatrix();
			double[,] expectedK = new double[,]
			{
				{  514198.06499808, -10764.170693892, -403744.28140248,  6179.6240659003,  136202.13267487, -257206.34711690, -246655.91627047,  261790.89374489 },
				{ -10764.170693892,  877762.34956189,  61124.679010955,  241708.83734231, -257206.34711690, -50430.336321830,  206845.83879984, -1069040.8505824 },
				{ -403744.28140248,  61124.679010955,  2618366.6203953, -1268741.0140479, -648721.49301375,  372500.34071924, -1565900.8459791,  835115.99431770 },
				{  6179.6240659003,  241708.83734231, -1268741.0140479,  3119360.6646083,  427445.39566429, -1580482.4587671,  835115.99431770, -1780587.0431835 },
				{  136202.13267487, -257206.34711690, -648721.49301375,  427445.39566429,  691579.93709070, -83847.039187745, -179060.57675182, -86392.009359646 },
				{ -257206.34711690, -50430.336321830,  372500.34071924, -1580482.4587671, -83847.039187745,  1103398.3531728, -31446.954414592,  527514.44191615 },
				{ -246655.91627047,  206845.83879984, -1565900.8459791,  835115.99431770, -179060.57675182, -31446.954414592,  1991617.3390014, -1010514.8787030 },
				{  261790.89374489, -1069040.8505824,  835115.99431770, -1780587.0431835, -86392.009359646,  527514.44191615, -1010514.8787030,  2322113.4518497 }
			}; // from Abaqus
			Assert.True(K.Equals(Matrix.CreateFromArray(expectedK), 1e-10));
		}

		private static void TestStiffness4()
		{
			var factory = new ContinuumElement2DFactory(thickness, material1, dynamicMaterial);
			ContinuumElement2D quad4 = factory.CreateElement(CellType.Quad4, nodeSet3);
			IMatrix K = quad4.BuildStiffnessMatrix();
			double[,] expectedK = new double[,]
			{
				{ 0.2592059570522744,       -0.005023279657148771,  -0.1906544880785286,    -0.017628995948735127,  0.06474697564228854,    -0.12002962865455298,   -0.13329844461603435,   0.1426819042604369      },
				{ -0.005023279657148771,    0.42886928984872097,    0.049037670717931546,   0.11055696733570626,    -0.12002962865455298,   -0.022348176556173396,  0.0760152375937702,     -0.5170780806282538     },
				{ -0.1906544880785286,      0.049037670717931546,   1.3012408988907096,     -0.5920791398890163,    -0.33356025755449087,   0.15332067182282197,    -0.77702615325769,      0.3897207973482628      },
				{ -0.017628995948735127,    0.11055696733570626,    -0.5920791398890163,    1.5350381195234322,     0.21998733848948865,    -0.7683820415727372,    0.3897207973482628,     -0.8772130452864011     },
				{ 0.06474697564228854,      -0.12002962865455298,   -0.33356025755449087,   0.21998733848948865,    0.3475567568780642,     -0.03912861828761283,   -0.07874347496586184,   -0.060829091547322856   },
				{ -0.12002962865455298,     -0.022348176556173396,  0.15332067182282197,    -0.7683820415727372,    -0.03912861828761283,   0.5397386843830521,     0.005837575119343828,   0.25099153374585864     },
				{ -0.13329844461603435,     0.0760152375937702,     -0.77702615325769,      0.3897207973482628,     -0.07874347496586184,   0.005837575119343828,   0.9890680728395862,     -0.4715736100613768     },
				{ 0.1426819042604369,       -0.5170780806282538,    0.3897207973482628,     -0.8772130452864011,    -0.060829091547322856,  0.25099153374585864,    -0.4715736100613768,    1.1432995921687963      }
			}; // from Solverize
			   //Assert.True(Utilities.AreMatricesEqual(K, new Matrix(expectedK), 1e-10));
		}

		[Fact]
		public static void TestStrainsStresses1()
		{
			var factory = new ContinuumElement2DFactory(thickness, material0, null);
			ContinuumElement2D quad4 = factory.CreateElement(CellType.Quad4, nodeSet2);

			// Abaqus results
			double[] displacements =
			{
				  0.0,          0.0,          // Node 1
                  0.0,          0.0,          // Node 2
                  1.0,       -499.614E-03,    // Node 4  
                986.100E-03,    1.0           // Node 3
            };

			// The Gauss points in Abaqus have the same order as in GaussLegendre2D.Order2x2: Xi major, Eta minor
			double[][] expectedStrainsAtGPs =
			{
				new double[] { 1.46867E-03,  341.547E-03,  336.066E-03 },  // Gauss point 1
                new double[] { 1.46867E-03,  -91.3541E-03, 340.078E-03 },  // Gauss point 2
                new double[] { 5.48114E-03,  341.547E-03,  -96.8352E-03 }, // Gauss point 3
                new double[] { 5.48114E-03,  -91.3541E-03, -92.8228E-03 }  // Gauss point 4
            };
			double[][] expectedStressesAtGPs =
			{
				new double[] { 23.9845E+03,   78.9202E+03,  27.1438E+03 },  // Gauss point 1
                new double[] { -5.98559E+03, -20.9800E+03,  27.4679E+03 },  // Gauss point 2
                new double[] { 24.9104E+03,   79.1980E+03,  -7.82131E+03 }, // Gauss point 3
                new double[] { -5.05964E+03, -20.7023E+03,  -7.49722E+03 }  // Gauss point 4
            };

			// The order of the nodes is also the same (at least in this job)
			double[][] expectedStrainsAtNodes =
			{
				new double[] { 58.2077E-12,  500.000E-03,  493.050E-03 },  // Node 1
                new double[] {  0.0,        -249.807E-03,  500.0E-03 },    // Node 2
                new double[] { 6.94981E-03, -249.807E-03, -249.807E-03 },  // Node 4
                new double[] { 6.94981E-03,  500.000E-03, -256.757E-03 }   // Node 3
            };
			double[][] expectedStressesAtNodes =
			{
				new double[] {  34.6154E+03, 115.385E+03,   39.8233E+03 },  // Node 1
                new double[] { -17.2943E+03, -57.6478E+03,  40.3846E+03 },  // Node 2
                new double[] { -15.6905E+03, -57.1666E+03, -20.1767E+03 },  // Node 4
                new double[] {  36.2192E+03, 115.866E+03,  -20.7380E+03 }   // Node 3
            };

			(IReadOnlyList<double[]> strainsAtGPs, IReadOnlyList<double[]> stressesAtGPs) =
				quad4.UpdateStrainsStressesAtGaussPoints(displacements);
			IReadOnlyList<double[]> strainsAtNodes =
				quad4.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(strainsAtGPs, quad4.Interpolation);
			IReadOnlyList<double[]> stressesAtNodes =
				quad4.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(stressesAtGPs, quad4.Interpolation);

			Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtGPs, strainsAtGPs, 1e-4));
			Assert.True(Utilities.AreTensorsEqual(expectedStressesAtGPs, stressesAtGPs, 1e-4));
			Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtNodes, strainsAtNodes, 1e-4));
			Assert.True(Utilities.AreTensorsEqual(expectedStressesAtNodes, stressesAtNodes, 1e-4));
		}
		#endregion
	}
}
