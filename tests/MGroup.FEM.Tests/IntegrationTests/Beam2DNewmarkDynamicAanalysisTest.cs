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
	using NumericalAnalyzers.Dynamic;
	using Solvers.Direct;
	using Structural.Elements;

	public class Beam2DNewmarkDynamicAanalysisTest
	{
		public void LinearElasticBeam2DNewmarkDynamicAnalysisTest()
		{
			double youngModulus = 21000;
			double poissonRatio = 0.3;
			double area = 91.04;
			double inertiaY = 2843.0;
			double inertiaZ = 8091.0;
			double density = 7.85;
			double nodalLoad = 1000.0;
			int totalNodes = 2;

			var material = new ElasticMaterial()
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio,
			};

			// Node creation
			IList<Node> nodes = new List<Node>();
			Node node1 = new Node(id: 1, x: 0.0, y: 0.0, z: 0.0);
			Node node2 = new Node(id: 2, x: 300.0, y: 0.0, z: 0.0);
			nodes.Add(node1);
			nodes.Add(node2);

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
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });

			// Create a new Beam2D element
			var beam = new EulerBeam2D(youngModulus)
			{
				Density = density,
				SectionArea = area,
				MomentOfInertia = inertiaZ
			};

			var element = new Element()
			{
				ID = 1,
				ElementType = beam
			};

			// Add nodes to the created element
			element.AddNode(model.NodesDictionary[1]);
			element.AddNode(model.NodesDictionary[2]);

			// Element Stiffness Matrix
			var a = beam.StiffnessMatrix(element);
			var b = beam.MassMatrix(element);

			// Add Hexa element to the element and subdomains dictionary of the model
			model.ElementsDictionary.Add(element.ID, element);
			model.SubdomainsDictionary[1].Elements.Add(element);

			// define loads
			model.Loads.Add(new Load { Amount = nodalLoad, Node = model.NodesDictionary[totalNodes], DOF = StructuralDof.TranslationY });

			// Solver
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);

			// Problem type
			var provider = new ProblemStructural(model, solver);

			// Analyzers
			var childAnalyzer = new LinearAnalyzer(model, solver, provider);
			var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 0.28, 3.36);
			parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
			NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			Assert.Equal(2.2840249264795207, solver.LinearSystems[1].Solution[1], 8);
		}
	}
}
