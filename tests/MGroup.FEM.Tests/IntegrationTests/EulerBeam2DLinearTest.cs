using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Loads;
using Xunit;

namespace MGroup.FEM.Tests.IntegrationTests
{
	using Constitutive.Structural;
	using MSolve.Solution;
	using NumericalAnalyzers;
	using Solvers.Direct;
	using Structural.Elements;

	public class EulerBeam2DLinearTest
	{
		[Fact]
		public void TestEulerBeam2DLinearBendingExample()
		{
			double youngModulus = 21000.0;
			double poissonRatio = 0.3;
			double nodalLoad = 2000.0;
			int nElems = 2;
			int monitorNode = 3;

			// Create new material
			var material = new ElasticMaterial()
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio,
			};

			// Node creation
			IList<Node> nodes = new List<Node>();
			Node node1 = new Node(id: 1, x: 0.0, y: 0.0, z: 0.0);
			Node node2 = new Node(id: 2, x: 100.0, y: 0.0, z: 0.0);
			Node node3 = new Node(id: 3, x: 200.0, y: 0.0, z: 0.0);
			nodes.Add(node1);
			nodes.Add(node2);
			nodes.Add(node3);

			// Model creation
			Model model = new Model();

			// Add a single subdomain to the model
			model.SubdomainsDictionary.Add(1, new Subdomain(1));

			// Add nodes to the nodes dictonary of the model
			for (int i = 0; i < nodes.Count; ++i)
			{
				model.NodesDictionary.Add(i + 1, nodes[i]);
			}

			// Constrain bottom nodes of the model
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });
			// Generate elements of the structure
			int iNode = 1;
			for (int iElem = 0; iElem < nElems; iElem++)
			{
				// Create new Beam2D section and element
				var beam = new EulerBeam2D(youngModulus)
				{
					Density = 7.85,
					SectionArea = 91.04,
					MomentOfInertia = 8091.00,
				};

				// Create elements
				var element = new Element()
				{
					ID = iElem + 1,
					ElementType = beam
				};
				// Add nodes to the created element
				element.AddNode(model.NodesDictionary[iNode]);
				element.AddNode(model.NodesDictionary[iNode + 1]);

				var a = beam.StiffnessMatrix(element);

				// Add Hexa element to the element and subdomains dictionary of the model
				model.ElementsDictionary.Add(element.ID, element);
				model.SubdomainsDictionary[1].Elements.Add(element);
				iNode++;
			}

			// Add nodal load values at the top nodes of the model
			model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = StructuralDof.TranslationY });

			// Solver
			//var solverBuilder = new MSolve.Solvers.Iterative.PcgSolver.Builder(
			//    (new LinearAlgebra.Iterative.ConjugateGradient.PcgAlgorithm.Builder()).Build());
			//var solver = solverBuilder.BuildSolver(model);
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);

			// Problem type
			var provider = new ProblemStructural(model, solver);

			// Analyzers
			var childAnalyzer = new LinearAnalyzer(model, solver, provider);
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			// Run the anlaysis 
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			double solutionNorm = solver.LinearSystems[1].Solution.Norm2();
			double rhsNorm = solver.LinearSystems[1].RhsVector.Norm2();

			Assert.Equal(31.388982074929341, solver.LinearSystems[1].Solution[4], 12);
		}
	}
}
