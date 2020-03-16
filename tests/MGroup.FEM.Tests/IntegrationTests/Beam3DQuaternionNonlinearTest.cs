using System.Collections.Generic;
using MGroup.FEM.Elements.SupportiveClasses;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Loads;
using Xunit;

namespace MGroup.FEM.Tests.IntegrationTests
{
	using Constitutive.Structural;
	using Constitutive.Structural.ContinuumElements;
	using MGroup.MSolve.Solution;
	using MSolve.Constitutive;
	using NumericalAnalyzers;
	using NumericalAnalyzers.Logging;
	using NumericalAnalyzers.NonLinear;
	using Solvers.Direct;
	using Structural.Elements;

	public class Beam3DQuaternionNonlinearTest
	{
		private const string output = @"C:\Users\Serafeim\Desktop\output.txt";

		public void SolveNLBeam()
		{
			var m = new Model();
			m.NodesDictionary.Add(1, new Node(id: 1, x: 0, y: 0, z: 0));
			m.NodesDictionary.Add(2, new Node(id: 2, x: 5, y: 0, z: 0));
			m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
			m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
			m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });
			m.ElementsDictionary.Add(1, new Element()
			{
				ID = 1,
				ElementType = new Beam3DCorotationalQuaternion(m.Nodes, new ElasticMaterial3D() { YoungModulus = 2.1e6, PoissonRatio = 0.2 }, 1,
				new BeamSection3D(0.06, 0.0002, 0.00045, 0.000818, 0.05, 0.05))
			});
			m.ElementsDictionary[1].AddNodes(m.Nodes);
			m.SubdomainsDictionary.Add(1, new Subdomain(1));
			m.SubdomainsDictionary[1].Elements.Add(m.ElementsDictionary[1]);
			m.Loads.Add(new Load() { Node = m.NodesDictionary[2], Amount = 100, DOF = StructuralDof.TranslationY });

			// Solver
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(m);

			// Problem type
			var provider = new ProblemStructural(m, solver);

			// Analyzers
			int increments = 10;
			var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(m, solver, provider, increments);
			//childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
			LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
			var parentAnalyzer = new StaticAnalyzer(m, solver, provider, childAnalyzer);

			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();
		}

		[Fact]
		public void CantileverYBeam3DQuaternionNonlinearTest()
		{
			double youngModulus = 21000.0;
			double poissonRatio = 0.3;
			double nodalLoad = 20000.0;
			double area = 91.04;
			double inertiaY = 2843.0;
			double inertiaZ = 8091.0;
			double torsionalInertia = 76.57;
			double effectiveAreaY = 91.04;
			double effectiveAreaZ = 91.04;
			int nNodes = 3;
			int nElems = 2;
			int monitorNode = 3;

			// Create new 3D material
			var material = new ElasticMaterial3D()
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
				// element nodes
				IList<Node> elementNodes = new List<Node>();
				elementNodes.Add(model.NodesDictionary[iNode]);
				elementNodes.Add(model.NodesDictionary[iNode + 1]);

				// Create new Beam3D section and element
				var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
				var beam = new Beam3DCorotationalQuaternion(elementNodes, material, 7.85, beamSection);

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
				//var writer = new FullMatrixWriter();
				//writer.WriteToFile(a, output);


				// Add beam element to the element and subdomains dictionary of the model
				model.ElementsDictionary.Add(element.ID, element);
				model.SubdomainsDictionary[1].Elements.Add(element);
				iNode++;
			}

			// Add nodal load values at the top nodes of the model
			model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = StructuralDof.TranslationY });

			// Solver
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);

			// Problem type
			var provider = new ProblemStructural(model, solver);

			// Analyzers
			int increments = 10;
			var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
			childAnalyzerBuilder.ResidualTolerance = 1E-3;
			//childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
			LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			// Request output
			childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 7 });

			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check output
			DOFSLog log = (DOFSLog)childAnalyzer.Logs[1][0]; //There is a list of logs for each subdomain and we want the first one
			Assert.Equal(148.936792350562, log.DOFValues[7], 2);
		}

		[Fact]
		public void PlaneFrameTest()
		{
			double youngModulus = 21000.0;
			double poissonRatio = 0.3;
			double nodalLoad = 500000.0;
			double area = 91.04;
			double inertiaY = 2843.0;
			double inertiaZ = 8091.0;
			double torsionalInertia = 76.57;
			double effectiveAreaY = 91.04;
			double effectiveAreaZ = 91.04;
			int nNodes = 4;
			int nElems = 3;
			int monitorNode = 2;

			// Create new 3D material
			var material = new ElasticMaterial3D()
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio,
			};

			// Node creation
			IList<Node> nodes = new List<Node>();
			Node node1 = new Node(id: 1, x: 0.0, y: 0.0, z: 0.0);
			Node node2 = new Node(id: 2, x: 0.0, y: 100.0, z: 0.0);
			Node node3 = new Node(id: 3, x: 100.0, y: 100.0, z: 0.0);
			Node node4 = new Node(id: 4, x: 100.0, y: 0.0, z: 0.0);

			nodes.Add(node1);
			nodes.Add(node2);
			nodes.Add(node3);
			nodes.Add(node4);

			// Model creation
			Model model = new Model();

			// Add a single subdomain to the model
			model.SubdomainsDictionary.Add(1, new Subdomain(1));

			// Add nodes to the nodes dictonary of the model
			for (int i = 0; i < nodes.Count; ++i)
			{
				model.NodesDictionary.Add(i + 1, nodes[i]);
			}

			// Constrain first and last nodes of the model
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });
			model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
			model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
			model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });

			// Generate elements of the structure
			int iNode = 1;
			for (int iElem = 0; iElem < nElems; iElem++)
			{
				// element nodes
				IList<Node> elementNodes = new List<Node>();
				elementNodes.Add(model.NodesDictionary[iNode]);
				elementNodes.Add(model.NodesDictionary[iNode + 1]);

				// Create new Beam3D section and element
				var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
				var beam = new Beam3DCorotationalQuaternion(elementNodes, material, 7.85, beamSection);

				// Create elements
				var element = new Element()
				{
					ID = iElem + 1,
					ElementType = beam
				};

				var a = beam.StiffnessMatrix(element);

				// Add nodes to the created element
				element.AddNode(model.NodesDictionary[iNode]);
				element.AddNode(model.NodesDictionary[iNode + 1]);

				// Add beam element to the element and subdomains dictionary of the model
				model.ElementsDictionary.Add(element.ID, element);
				model.SubdomainsDictionary[1].Elements.Add(element);
				iNode++;
			}

			// Add nodal load values at the top nodes of the model
			model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = StructuralDof.TranslationX });

			// Solver
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);

			// Problem type
			var provider = new ProblemStructural(model, solver);

			// Analyzers
			int increments = 10;
			var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
			childAnalyzerBuilder.ResidualTolerance = 1E-3;
			//childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
			LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			// Request output
			childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 0 });

			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();


			// Check output
			DOFSLog log = (DOFSLog)childAnalyzer.Logs[1][0]; //There is a list of logs for each subdomain and we want the first one
			Assert.Equal(120.1108698752, log.DOFValues[0], 2);
		}
	}
}
