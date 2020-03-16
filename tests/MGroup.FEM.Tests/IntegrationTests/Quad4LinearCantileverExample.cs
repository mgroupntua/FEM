using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using MGroup.LinearAlgebra.Iterative.Termination;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;
using Xunit;

namespace MGroup.FEM.Tests.IntegrationTests
{
	using Constitutive.Structural;
	using Constitutive.Structural.PlanarElements;
	using MSolve.Constitutive;
	using NumericalAnalyzers;
	using NumericalAnalyzers.Logging;
	using Solvers.Iterative;
	using Structural.Elements;

	public static class Quad4LinearCantileverExample
	{
		[Fact]
		private static void TestQuad4LinearCantileverExample()
		{
			// Model & subdomains
			var model = new Model();
			int subdomainID = 0;
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

			// Materials
			double youngModulus = 3.76;
			double poissonRatio = 0.3779;
			double thickness = 1.0;
			double nodalLoad = 500.0;
			var material = new ElasticMaterial2D(StressState2D.PlaneStress)
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio
			};

			// Nodes
			var nodes = new Node[]
			{
				new Node( id: 1, x:  0.0, y:   0.0, z: 0.0 ),
				new Node( id: 2, x: 10.0, y:   0.0, z: 0.0 ),
				new Node( id: 3, x: 10.0, y:  10.0, z: 0.0 ),
				new Node( id: 4, x:  0.0, y:  10.0, z: 0.0 )
			};
			for (int i = 0; i < nodes.Length; ++i) model.NodesDictionary.Add(i, nodes[i]);


			// Elements
			var factory = new ContinuumElement2DFactory(thickness, material, null);

			var elementWrapper = new Element()
			{
				ID = 0,
				ElementType = factory.CreateElement(CellType.Quad4, nodes)
			};
			elementWrapper.AddNodes(nodes);
			model.ElementsDictionary.Add(elementWrapper.ID, elementWrapper);
			model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);

			//var a = quad.StiffnessMatrix(element);

			// Prescribed displacements
			model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
			model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
			model.NodesDictionary[3].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
			model.NodesDictionary[3].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });

			// Nodal loads
			model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[1], DOF = StructuralDof.TranslationX });
			model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[2], DOF = StructuralDof.TranslationX });

			// Solver
			var pcgBuilder = new PcgAlgorithm.Builder();
			pcgBuilder.ResidualTolerance = 1E-6;
			pcgBuilder.MaxIterationsProvider = new PercentageMaxIterationsProvider(0.5);
			var solverBuilder = new PcgSolver.Builder();
			solverBuilder.PcgAlgorithm = pcgBuilder.Build();
			PcgSolver solver = solverBuilder.BuildSolver(model);

			// Problem type
			var provider = new ProblemStructural(model, solver);

			// Analyzers
			var childAnalyzer = new LinearAnalyzer(model, solver, provider);
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);
			//NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);

			// Request output
			childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { 0 });

			// Run the anlaysis 
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check output
			DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
			Assert.Equal(253.132375961535, log.DOFValues[0], 8);
		}
	}
}
