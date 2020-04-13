using System.Collections.Generic;
using ISAAR.MSolve.Logging;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.ContinuumElements;
using MGroup.FEM.Entities;
using MGroup.FEM.Structural.Elements;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Commons;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;
using MGroup.MSolve.Solution;
using MGroup.NumericalAnalyzers;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.NonLinear;
using MGroup.Solvers.Direct;
using Xunit;
using System.Collections.Generic;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using Xunit;
using System.Linq;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Discretization.Loads;
using ISAAR.MSolve.FEM.Elements;
using MGroup.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation;

namespace ISAAR.MSolve.FEM.Tests.Elements
{
	

	public class CohesiveElementHexa8Interface
	{
		[Fact]
		public static void SolveCantilever()
		{
			Model model1 = new Model();
			model1.SubdomainsDictionary.Add(1, new Subdomain(1));
			Example_cohesive_hexa_orthi_constr_anw_benc1(model1, 8);
			var log1 = RunAnalysis(model1);

			double[] displacements1 = new double[4] { 0.00019727473744350047, 0.00019727473744365822, 0.00019727473744356587, 0.0001972747374435206 };

			double[] expected1 = new double[] { log1.GetTotalDisplacement(4, 1, 0), log1.GetTotalDisplacement(4, 1, 1), log1.GetTotalDisplacement(4, 1, 2), log1.GetTotalDisplacement(4, 1, 3) };
			Model model2 = new Model();
			model2.SubdomainsDictionary.Add(1, new Subdomain(1));
			Example_cohesive_hexa_orthi_constr_anw_benc1(model2, 13.75);
			var log2 = RunAnalysis(model2);


			double[] displacements2 = new double[4] { 0.00045445796734459886, 0.00045445796734467687, 0.00045445796734467741, 0.00045445796734453625 };
			double[] expected2 = new double[] { log2.GetTotalDisplacement(4, 1, 0), log2.GetTotalDisplacement(4, 1, 1), log2.GetTotalDisplacement(4, 1, 2), log2.GetTotalDisplacement(4, 1, 3) };


			Assert.True(AreDisplacementsSame(expected1, displacements1));
			Assert.True(AreDisplacementsSame(expected2, displacements2));

		}

		private static bool AreDisplacementsSame(double[] expectedDisplacements, double[] computedDisplacements)//.
		{
			var comparer = new ValueComparer(1E-13);
			for (int i1 = 0; i1 < expectedDisplacements.Length; i1++)
			{
				if (!comparer.AreEqual(expectedDisplacements[i1], computedDisplacements[i1]))
				{
					return false;
				}
			}//.
			return true;
		}

		public static TotalDisplacementsPerIterationLog RunAnalysis(Model model)
		{
			// Solver
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);

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
			watchDofs.Add(1, new int[] { 0, 1, 2, 3 });
			var log1 = new TotalDisplacementsPerIterationLog(watchDofs);
			childAnalyzer.TotalDisplacementsPerIterationLog = log1;

			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			return log1;
		}

		//original:C:\Users\acivi\Documents\notes_elegxoi_2\develop_cohesive_tet\MSolve-0f36677a620c9801e578ad6e2ff6fc5ff80b2992_apo_links_internet
		//plhrofories: sto idio link apo panw
		public static void Example_cohesive_hexa_mixed(Model model)
		{
			CohesiveElementHexa8Interface.Example2Hexa8NL1Cohesive8node(model);
			CohesiveElementHexa8Interface.Example2Hexa8NL1Cohesive8nodeConstraintsMixed(model);
			CohesiveElementHexa8Interface.Example2Hexa8NL1Cohesive8nodeLoadsMIxed(model, 1);
		}

		public static void Example_cohesive_hexa_orthi_elastic(Model model)
		{
			CohesiveElementHexa8Interface.Example2Hexa8NL1Cohesive8node(model); // me 1353000
			CohesiveElementHexa8Interface.Example2Hexa8NL1Cohesive8nodeConstraintsOrthiElastic(model);
			CohesiveElementHexa8Interface.Example2Hexa8NL1Cohesive8nodeLoadsElasticOrthi(model, 3.47783);
		}

		public static void Example_cohesive_hexa_orthi_constr_anw_benc1(Model model, double load_value)
		{
			CohesiveElementHexa8Interface.Example2Hexa8NL1Cohesive8node(model); // me 135300
			CohesiveElementHexa8Interface.Example2Hexa8NL1Cohesive8nodeConstraintsBenc1(model);
			//ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8nodeLoadsBenc1(model, 8);// gia elastiko klado 
			//ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8nodeLoadsBenc1(model, 13.75); // gia metelastiko
			CohesiveElementHexa8Interface.Example2Hexa8NL1Cohesive8nodeLoadsBenc1(model, load_value);// gia elastiko klado 
		}



		public static void Example2Hexa8NL1Cohesive8node(Model model)
		{
			ElasticMaterial3D material1 = new ElasticMaterial3D()
			{
				YoungModulus = 135300, // 1353000 gia to allo paradeigma
				PoissonRatio = 0.3,
			};
			BenzeggaghKenaneCohesiveMaterial material2 = new BenzeggaghKenaneCohesiveMaterial()
			{
				T_o_3 = 57,// N / mm2
				D_o_3 = 0.000057, // mm
				D_f_3 = 0.0098245610,  // mm

				T_o_1 = 57,// N / mm2
				D_o_1 = 0.000057, // mm
				D_f_1 = 0.0098245610,  // mm

				n_curve = 1.4
			};


			double[,] nodeData = new double[,] {
			{0.500000,0.000000,1.000000},
			{0.500000,0.500000,1.000000},
			{0.000000,0.000000,1.000000},
			{0.000000,0.500000,1.000000},
			{0.500000,0.000000,0.500000},
			{0.500000,0.500000,0.500000},
			{0.000000,0.000000,0.500000},
			{0.000000,0.500000,0.500000},
			{0.500000,0.000000,0.000000},
			{0.500000,0.500000,0.000000},
			{0.000000,0.000000,0.000000},
			{0.000000,0.500000,0.000000},
			{0.500000,0.000000,0.500000},
			{0.500000,0.500000,0.500000},
			{0.000000,0.000000,0.500000},
			{0.000000,0.500000,0.500000} };

			int[,] elementData = new int[,] { { 1, 1, 2, 4, 3, 5, 6, 8, 7 }, { 2, 13, 14, 16, 15, 9, 10, 12, 11 }, };

			int[] cohesive8Nodes = new int[] { 5, 6, 8, 7, 13, 14, 16, 15, };

			// orismos shmeiwn
			for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
			{
				model.NodesDictionary.Add(nNode + 1, new Node(id: nNode + 1, x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));

			}

			// orismos elements 
			Element e1;
			int subdomainID = 1;
			for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
			{
				List<Node> nodeSet = new List<Node>(8);
				for (int j = 0; j < 8; j++)
				{
					int nodeID = elementData[nElement, j];
					nodeSet.Add((Node)model.NodesDictionary[nodeID]);
				}
				e1 = new Element()
				{
					ID = nElement + 1,
					ElementType = new ContinuumElement3DNonLinear(nodeSet, material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3), InterpolationHexa8.UniqueInstance) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
				};
				for (int j = 0; j < 8; j++)
				{
					e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
				}
				model.ElementsDictionary.Add(e1.ID, e1);
				model.SubdomainsDictionary[subdomainID].Elements.Add( e1);
			}

			// kai to cohesive
			e1 = new Element()
			{
				ID = 3,
				ElementType = new cohesiveElement(material2, GaussLegendre2D.GetQuadratureWithOrder(3, 3), InterpolationQuad4.UniqueInstance) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
			};
			for (int j = 0; j < 8; j++)
			{
				e1.NodesDictionary.Add(cohesive8Nodes[j], model.NodesDictionary[cohesive8Nodes[j]]);
			}
			model.ElementsDictionary.Add(e1.ID, e1);
			model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
		}

		public static void Example2Hexa8NL1Cohesive8nodeConstraintsOrthiElastic(Model model)
		{
			for (int k = 13; k < 17; k++)
			{
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			}
		}

		public static void Example2Hexa8NL1Cohesive8nodeConstraintsMixed(Model model)
		{
			for (int k = 13; k < 17; k++)
			{
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			}
		}

		public static void Example2Hexa8NL1Cohesive8nodeConstraintsBenc1(Model model)
		{
			for (int k = 1; k < 5; k++)
			{
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			}
			for (int k = 5; k < 9; k++)
			{
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			}
			for (int k = 9; k < 17; k++)
			{
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			}
		}

		public static void Example2Hexa8NL1Cohesive8nodeLoadsMIxed(Model model, double load_value)
		{
			Load load1;
			for (int k = 5; k < 9; k++)
			{
				load1 = new Load()
				{
					Node = model.NodesDictionary[k],
					DOF = StructuralDof.TranslationX,
					Amount = 1 * load_value

				};
				model.Loads.Add(load1);
				load1 = new Load()
				{
					Node = model.NodesDictionary[k],
					DOF = StructuralDof.TranslationZ,
					Amount = 0.2 * load_value

				};
				model.Loads.Add(load1);
			}
		}

		public static void Example2Hexa8NL1Cohesive8nodeLoadsElasticOrthi(Model model, double load_value)
		{
			Load load1;
			for (int k = 1; k < 5; k++)
			{
				load1 = new Load()
				{
					Node = model.NodesDictionary[k],
					DOF = StructuralDof.TranslationZ,
					Amount = 1 * load_value

				};
				model.Loads.Add(load1);
			}
		}

		public static void Example2Hexa8NL1Cohesive8nodeLoadsBenc1(Model model, double load_value)
		{
			Load load1;
			for (int k = 5; k < 9; k++)
			{
				load1 = new Load()
				{
					Node = model.NodesDictionary[k],
					DOF = StructuralDof.TranslationZ,
					Amount = 1 * load_value

				};
				model.Loads.Add(load1);
			}
		}
	}
}
