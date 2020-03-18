using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Commons;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Loads;
using MGroup.Solvers.Direct;
using Xunit;

namespace MGroup.FEM.Tests.IntegrationTests
{
	using Constitutive.Structural;
	using Constitutive.Structural.ShellElements;
	using NumericalAnalyzers;
	using NumericalAnalyzers.Logging;
	using NumericalAnalyzers.NonLinear;
	using Structural.Elements;

	public static class Shell8andCohesiveNonLinear
	{
		private const int subdomainID = 0;

		[Fact]
		public static void RunTest()
		{
			IReadOnlyList<Dictionary<int, double>> expectedDisplacements = GetExpectedDisplacements();
			//IncrementalDisplacementsLog computedDisplacements = SolveModel();
			TotalDisplacementsPerIterationLog computedDisplacements = SolveModel();
			Assert.True(AreDisplacementsSame(expectedDisplacements, computedDisplacements));
		}

		private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements,
			TotalDisplacementsPerIterationLog computedDisplacements)
		{
			var comparer = new ValueComparer(1E-10); // for node major dof order and skyline solver
													 //var comparer = new ValueComparer(1E-3); // for other solvers. It may require adjusting after visual inspection
			for (int iter = 0; iter < expectedDisplacements.Count; ++iter)
			{
				foreach (int dof in expectedDisplacements[iter].Keys)
				{
					if (!comparer.AreEqual(expectedDisplacements[iter][dof], computedDisplacements.GetTotalDisplacement(iter, subdomainID, dof)))
					{
						return false;
					}
				}
			}
			return true;
		}

		private static IReadOnlyList<Dictionary<int, double>> GetExpectedDisplacements()
		{
			var expectedDisplacements = new Dictionary<int, double>[5]; //TODO: this should be 11 EINAI ARRAY APO DICTIONARIES

			expectedDisplacements[0] = new Dictionary<int, double> {
	{ 0,-1.501306714739351400e-05 }, {11,4.963733738129490800e-06 }, {23,-1.780945407868029400e-05 }, {35,-1.499214801866540600e-05 }, {39,-5.822833969672272200e-05 }};
			expectedDisplacements[1] = new Dictionary<int, double> {
	{ 0,-1.500991892603005000e-05 }, {11,4.962619842302796000e-06 }, {23,-1.780557361553905700e-05 }, {35,-1.498958552758854400e-05 }, {39,-5.821676140520536400e-05 }};
			expectedDisplacements[2] = new Dictionary<int, double> {
	{ 0,-3.001954880280401800e-05 }, {11,9.925100656477526600e-06 }, {23,-3.561116405104391700e-05 }, {35,-2.997946837566090700e-05 }, {39,-1.164336113147322500e-04 }};
			expectedDisplacements[3] = new Dictionary<int, double> {
	{ 0,-3.074327250558424700e-05 }, {11,1.064972618932890100e-05 }, {23,-3.846410374898863100e-05 }, {35,-3.069783728664514200e-05 }, {39,-1.191612724600880000e-04 }};
			expectedDisplacements[4] = new Dictionary<int, double> {
	{ 0,-3.074281618479765600e-05 }, {11,1.064926767853693300e-05 }, {23,-3.846254167901110600e-05 }, {35,-3.069737876082750600e-05 }, {39,-1.191596225034872200e-04 }};


			return expectedDisplacements;
		}

		private static TotalDisplacementsPerIterationLog SolveModel()
		{
			var model = new Model();
			//model.dofOrderer = (subdomain) => (new SimpleDofOrderer()).OrderDofs(model);
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));
			ShellAndCohesiveRAM_11tlkShellPaktwsh(model);

			// Solver
			var solverBuilder = new SkylineSolver.Builder();
			//var solverBuilder = new Solvers.Dense.DenseMatrixSolver.Builder();
			//solverBuilder.DofOrderer = new DofOrderer(new SimpleDofOrderingStrategy(), new NullReordering());
			var solver = solverBuilder.BuildSolver(model);

			// Problem type
			var provider = new ProblemStructural(model, solver);

			// Analyzers
			int increments = 2;
			var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
			childAnalyzerBuilder.ResidualTolerance = 1E-8;
			childAnalyzerBuilder.MaxIterationsPerIncrement = 100;
			childAnalyzerBuilder.NumIterationsForMatrixRebuild = 1;
			//childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
			LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();

			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			// Output
			var watchDofs = new Dictionary<int, int[]>();
			watchDofs.Add(subdomainID, new int[5] { 0, 11, 23, 35, 39 });
			var log1 = new TotalDisplacementsPerIterationLog(watchDofs);
			childAnalyzer.TotalDisplacementsPerIterationLog = log1;

			// Run the anlaysis 
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			return log1;
		}

		private static void ShellAndCohesiveRAM_11tlkShellPaktwsh(Model model)
		{
			//Origin: dhmiourgithike kata to ParadeigmataElegxwnBuilder.ShellAndCohesiveRAM_11ShellPaktwsh(model);
			// allaxame to cohesive element
			// gewmetria
			double Tk = 0.5;

			int nodeID = 1;

			double startX = 0;
			double startY = 0;
			double startZ = 0;
			for (int l = 0; l < 3; l++)
			{
				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: startX, y: startY + l * 0.25, z: startZ));
				nodeID++;
			}

			startX = 0.25;
			for (int l = 0; l < 2; l++)
			{
				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: startX, y: startY + l * 0.5, z: startZ));
				nodeID++;
			}

			startX = 0.5;
			for (int l = 0; l < 3; l++)
			{
				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: startX, y: startY + l * 0.25, z: startZ));
				nodeID++;
			}

			// katw strwsh pou tha paktwthei

			startX = 0;
			for (int l = 0; l < 3; l++)
			{
				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: startX, y: startY + l * 0.25, z: startZ - 0.5 * Tk));
				nodeID++;
			}

			startX = 0.25;
			for (int l = 0; l < 2; l++)
			{
				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: startX, y: startY + l * 0.5, z: startZ - 0.5 * Tk));
				nodeID++;
			}

			startX = 0.5;
			for (int l = 0; l < 3; l++)
			{
				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: startX, y: startY + l * 0.25, z: startZ - 0.5 * Tk));
				nodeID++;
			}

			double[][] VH = new double[8][];

			for (int j = 0; j < 8; j++)
			{
				VH[j] = new double[3];
				VH[j][0] = 0;
				VH[j][1] = 0;
				VH[j][2] = 1;
			}
			// perioxh gewmetrias ews edw

			// constraints

			nodeID = 9;
			for (int j = 0; j < 8; j++)
			{
				model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
				model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
				nodeID++;
			}
			//perioxh constraints ews edw

			// perioxh materials 
			var material1 = new BenzeggaghKenaneCohesiveMaterial()
			{
				T_o_3 = 57, // New load case argurhs NR_shell_coh.m
				D_o_3 = 5.7e-5,
				D_f_3 = 0.0098245610,
				T_o_1 = 57,
				D_o_1 = 5.7e-5,
				D_f_1 = 0.0098245610,
				n_curve = 1.4,
			};

			//IContinuumMaterial3D material2 = new IContinuumMaterial3D()
			//{
			//    YoungModulus = 1353000,
			//    PoissonRatio = 0.3,
			//};
			var material2 = new ShellElasticMaterial3D()
			{
				YoungModulus = 1353000,
				PoissonRatio = 0.3,
				ShearCorrectionCoefficientK = 5 / 6,
			};
			// perioxh materials ews edw


			//eisagwgh tou shell element
			double[] Tk_vec = new double[8];
			for (int j = 0; j < 8; j++)
			{
				Tk_vec[j] = Tk;
			}

			Element e1;
			e1 = new Element()
			{
				ID = 1,
				ElementType = new Shell8NonLinear(material2, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))
				{
					oVn_i = VH,
					tk = Tk_vec,
				}
			};
			e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
			e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
			e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
			e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
			e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
			e1.NodesDictionary.Add(2, model.NodesDictionary[2]);
			e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
			e1.NodesDictionary.Add(7, model.NodesDictionary[7]);

			model.ElementsDictionary.Add(e1.ID, e1);
			model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
			//eisagwgh shell ews edw

			// eisagwgh tou cohesive element
			int[] coh_global_nodes;
			coh_global_nodes = new int[] { 8, 3, 1, 6, 5, 2, 4, 7, 16, 11, 9, 14, 13, 10, 12, 15 };

			Element e2;
			e2 = new Element()
			{
				ID = 2,
				ElementType = new CohesiveShell8ToHexa20(material1, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
				{
					oVn_i = VH,
					tk = Tk_vec,
					ShellElementSide = 0,
				}
			};

			for (int j = 0; j < 16; j++)
			{
				e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
			}

			model.ElementsDictionary.Add(e2.ID, e2);
			model.SubdomainsDictionary[subdomainID].Elements.Add(e2);
			// eisagwgh cohesive ews edw

			// perioxh loads
			double value_ext;
			value_ext = 2 * 2.5 * 0.5;

			int[] points_with_negative_load;
			points_with_negative_load = new int[] { 1, 3, 6, 8 };
			int[] points_with_positive_load;
			points_with_positive_load = new int[] { 2, 4, 5, 7 };

			Load load1;
			Load load2;

			// LOADCASE '' orthi ''
			//for (int j = 0; j < 4; j++)
			//{
			//    load1 = new Load()
			//    {
			//        Node = model.NodesDictionary[points_with_negative_load[j]],
			//        DOF = DOFType.Z,
			//        Amount = -0.3333333 * value_ext,
			//    };
			//    model.Loads.Add(load1);

			//    load2 = new Load()
			//    {
			//        Node = model.NodesDictionary[points_with_positive_load[j]],
			//        DOF = DOFType.Z,
			//        Amount = 1.3333333 * value_ext,
			//    };
			//    model.Loads.Add(load2);
			//}

			// LOADCASE '' orthi '' dixws ta duo prwta fortia  (-0.3333) kai (1.3333)
			for (int j = 0; j < 3; j++)
			{
				load1 = new Load()
				{
					Node = model.NodesDictionary[points_with_negative_load[j + 1]],
					DOF = StructuralDof.TranslationZ,
					Amount = -0.3333333 * value_ext,
				};
				model.Loads.Add(load1);

				load2 = new Load()
				{
					Node = model.NodesDictionary[points_with_positive_load[j + 1]],
					DOF = StructuralDof.TranslationZ,
					Amount = 1.3333333 * value_ext,
				};
				model.Loads.Add(load2);
			}


			// perioxh loads ews edw
		}
	}

}
