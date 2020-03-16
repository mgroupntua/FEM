using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Commons;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;
using Xunit;

namespace MGroup.FEM.Tests.IntegrationTests
{
	using Constitutive.Thermal;
	using NumericalAnalyzers;
	using Solvers.Direct;
	using Thermal.Elements;

	public class ThermalBenchmarkProvatidis_11_2
	{
		private const int subdomainID = 0;

		[Fact]
		private static void RunTest()
		{
			Model Model = CreateModel();
			IVectorView solution = SolveModel(Model);
			Assert.True(CompareResults(solution));
		}

		private static bool CompareResults(IVectorView solution)
		{
			var comparer = new ValueComparer(1E-8);

			//                                                   dofs:   1,   2,   4,   5,   7,   8
			var expectedSolution = Vector.CreateFromArray(new double[] { 150, 200, 150, 200, 150, 200 });
			int numFreeDofs = 6;
			if (solution.Length != 6) return false;
			for (int i = 0; i < numFreeDofs; ++i)
			{
				if (!comparer.AreEqual(expectedSolution[i], solution[i])) return false;
			}
			return true;
		}

		private static Model CreateModel()
		{
			var model = new Model();

			// Subdomains
			model.SubdomainsDictionary.Add(0, new Subdomain(subdomainID));

			// Material
			double density = 1.0;
			double k = 1.0;
			double c = 1.0;

			// Nodes
			int numNodes = 9;
			var nodes = new Node[numNodes];
			nodes[0] = new Node(id: 0, x: 0.0, y: 0.0);
			nodes[1] = new Node(id: 1, x: 1.0, y: 0.0);
			nodes[2] = new Node(id: 2, x: 2.0, y: 0.0);
			nodes[3] = new Node(id: 3, x: 0.0, y: 1.0);
			nodes[4] = new Node(id: 4, x: 1.0, y: 1.0);
			nodes[5] = new Node(id: 5, x: 2.0, y: 1.0);
			nodes[6] = new Node(id: 6, x: 0.0, y: 2.0);
			nodes[7] = new Node(id: 7, x: 1.0, y: 2.0);
			nodes[8] = new Node(id: 8, x: 2.0, y: 2.0);

			for (int i = 0; i < numNodes; ++i) model.NodesDictionary[i] = nodes[i];

			// Elements
			int numElements = 4;
			var elementFactory = new ThermalElement2DFactory(1.0, new ThermalMaterial(density, c, k));
			var elements = new ThermalElement2D[4];
			elements[0] = elementFactory.CreateElement(CellType.Quad4, new Node[] { nodes[0], nodes[1], nodes[4], nodes[3] });
			elements[1] = elementFactory.CreateElement(CellType.Quad4, new Node[] { nodes[1], nodes[2], nodes[5], nodes[4] });
			elements[2] = elementFactory.CreateElement(CellType.Quad4, new Node[] { nodes[3], nodes[4], nodes[7], nodes[6] });
			elements[3] = elementFactory.CreateElement(CellType.Quad4, new Node[] { nodes[4], nodes[5], nodes[8], nodes[7] });

			for (int i = 0; i < numElements; ++i)
			{
				var elementWrapper = new Element() { ID = i, ElementType = elements[i] };
				foreach (var node in elements[i].Nodes) elementWrapper.AddNode(node);
				model.ElementsDictionary[i] = elementWrapper;
				model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
			}

			// Dirichlet BC
			model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });
			model.NodesDictionary[3].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });
			model.NodesDictionary[6].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });

			// Neumann BC
			double q = 50.0;
			model.Loads.Add(new Load() { Amount = q / 2.0, Node = model.NodesDictionary[2], DOF = ThermalDof.Temperature });
			model.Loads.Add(new Load() { Amount = q, Node = model.NodesDictionary[5], DOF = ThermalDof.Temperature });
			model.Loads.Add(new Load() { Amount = q / 2.0, Node = model.NodesDictionary[8], DOF = ThermalDof.Temperature });

			return model;
		}

		private static IVectorView SolveModel(Model model)
		{
			SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);
			var provider = new ProblemThermal(model, solver);

			var childAnalyzer = new LinearAnalyzer(model, solver, provider);
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			return solver.LinearSystems[subdomainID].Solution;
		}
	}
}
