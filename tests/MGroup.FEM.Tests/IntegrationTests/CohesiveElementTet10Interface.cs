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
using MGroup.LinearAlgebra.Vectors;

namespace MTwin.Engine.NumericalAnalyzers.IntegrationTests
{

	public class CohesiveElementTet10Interface
	{
		[Fact]
		public void SolveInterface()
		{
			int subdomainID = 1;
			const int increments = 2;
			const double nodalDisplacement = -5.0;
			
			// Create Model
			Model model = new Model();

			// Create Subdomain
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

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


			int[] nodeIds = new int[] { 5, 6, 7, 8, 17, 18, 19, 20, 38, 39, 40, 41,42 };

			double[,] nodeData = new double[,]
				{{-0.2500000000000000,0.2500000000000000,1.0000000000000000},
				{0.2500000000000000,0.2500000000000000,1.0000000000000000},
				{0.2500000000000000,-0.2500000000000000,1.0000000000000000},
				{-0.2500000000000000,-0.2500000000000000,1.0000000000000000},

				{0.0000000000000000,0.2500000000000000,1.0000000000000000},
				{0.2500000000000000,0.0000000000000000,1.0000000000000000},
				{0.0000000000000000,-0.2500000000000000,1.0000000000000000},
				{-0.2500000000000000,0.0000000000000000,1.0000000000000000},

				{0.0000000000000000,0.0000000000000000,1.0000000000000000},
				{-0.1250000000000000,0.1250000000000000,1.0000000000000000},
				{0.1250000000000000,0.1250000000000000,1.0000000000000000},
				{0.1250000000000000,-0.1250000000000000,1.0000000000000000},
				{-0.1250000000000000,-0.1250000000000000,1.0000000000000000},};

			int[,] elementData = new int[,]
			{
					{6,5,38,17,39,40 },
					{5,8,38,20,42,39 },
					{7,6,38,18,40,41 },
					{8,7,38,19,41,42 }
			};

			int[] upperSideSuplicateNodeIds = nodeIds.Select(x => x + 100).ToArray();

			// orismos shmeiwn
			for (int nNode = 0; nNode < nodeIds.GetLength(0); nNode++)
			{
				model.NodesDictionary.Add(nodeIds[nNode], new Node(id: nodeIds[nNode], x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));

			}
			for (int nNode = 0; nNode < upperSideSuplicateNodeIds.GetLength(0); nNode++)
			{
				model.NodesDictionary.Add(upperSideSuplicateNodeIds[nNode], new Node(id: upperSideSuplicateNodeIds[nNode], x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));

			}

			// orismos elements 
			Element e1;
			
			for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
			{
				e1 = new Element()
				{
					ID = nElement+1,
					ElementType = new cohesiveElement(material2, /*GaussLegendre2D.GetQuadratureWithOrder(3, 3)*/ TriangleQuadratureSymmetricGaussian.Order2Points3, InterpolationTri6.UniqueInstance) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
				};
				for (int j = 0; j < 6; j++)
				{
					e1.NodesDictionary.Add(elementData[nElement, j] + 100, model.NodesDictionary[elementData[nElement, j] + 100]);
				}
				for (int j = 0; j < 6; j++)
				{
					e1.NodesDictionary.Add(elementData[nElement, j], model.NodesDictionary[elementData[nElement, j]]);
				}
				model.ElementsDictionary.Add(e1.ID, e1);
				model.SubdomainsDictionary[subdomainID].Elements.Add( e1);
			}

			// Boundary Condtitions
			for (int i1 = 0; i1< nodeIds.GetLength(0); i1++)
			{
				model.NodesDictionary[nodeIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[nodeIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
				model.NodesDictionary[nodeIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			}

			//for (int i1 = 0; i1 < upperSideSuplicateNodeIds.GetLength(0)-1; i1++)
			//{
			//	model.NodesDictionary[upperSideSuplicateNodeIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			//	model.NodesDictionary[upperSideSuplicateNodeIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			//	model.NodesDictionary[upperSideSuplicateNodeIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ, Amount = 0.0001 });
			//}

			//model.NodesDictionary[upperSideSuplicateNodeIds[upperSideSuplicateNodeIds.GetLength(0) - 1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			//model.NodesDictionary[upperSideSuplicateNodeIds[upperSideSuplicateNodeIds.GetLength(0) - 1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			//model.NodesDictionary[upperSideSuplicateNodeIds[upperSideSuplicateNodeIds.GetLength(0) - 1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ, Amount = 0.0001 });

			var solution = Vector.CreateFromArray(new double[upperSideSuplicateNodeIds.Count() * 3]);
			var dSolution = Vector.CreateFromArray(new double[upperSideSuplicateNodeIds.Count() * 3]);

			//model.SubdomainsDictionary.ElementAt(0).Value.GetRhsFromSolution(solution, dSolution);

			for (int i1 = 0; i1 < upperSideSuplicateNodeIds.GetLength(0); i1++)
			{
				solution[3 * i1 + 2] = 0.00019727473744350047;
			}


			//model.ConnectDataStructures();



			// Choose linear equation system solver
			var solverBuilder = new SkylineSolver.Builder();
			SkylineSolver solver = solverBuilder.BuildSolver(model);

			// Choose the provider of the problem -> here a structural problem
			var provider = new ProblemStructural(model, solver);

			// Choose child analyzer -> Child: DisplacementControlAnalyzer
			var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) };
			var childAnalyzerBuilder = new DisplacementControlAnalyzer.Builder(model, solver, provider, increments)
			{
				MaxIterationsPerIncrement = 100,
				NumIterationsForMatrixRebuild = 1,
				ResidualTolerance = 1E-03
			};
			var childAnalyzer = childAnalyzerBuilder.Build();

			// Choose parent analyzer -> Parent: Static
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			// Request output
			childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { 0 });

			var watchDofs = new Dictionary<int, int[]>();
			watchDofs.Add(subdomainID, new int[] { 0, 1 });
			var log1 = new TotalLoadsDisplacementsPerIncrementLog(model.SubdomainsDictionary[subdomainID], 2, model.NodesDictionary[142], StructuralDof.TranslationZ, " ");
			childAnalyzer.IncrementalLogs.Add(1, log1);

			// Run the analysis
			parentAnalyzer.Initialize();
			var forces = model.SubdomainsDictionary.ElementAt(0).Value.GetRhsFromSolution(solution, dSolution);


			double totalForce = forces.CopyToArray().Sum();
			double expectedTotalForce = 14.04535167391627;
			Assert.True(((totalForce - expectedTotalForce) / expectedTotalForce) < 1e-8);

		}
	}
}
