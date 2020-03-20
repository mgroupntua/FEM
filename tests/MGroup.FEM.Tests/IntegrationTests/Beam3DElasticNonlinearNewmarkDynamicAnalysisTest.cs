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
	using NumericalAnalyzers.Dynamic;
	using NumericalAnalyzers.NonLinear;
	using Solvers.Direct;
	using Structural.Elements;

	public class Beam3DElasticNonlinearNewmarkDynamicAnalysisTest
	{
		private static void TestBeam3DElasticNonlinearNewmarkDynamicAnalysisExample()
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
			double density = 7.85;
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
			Node node2 = new Node(id: 2, x: 300.0, y: 0.0, z: 0.0);
			Node node3 = new Node(id: 3, x: 600.0, y: 0.0, z: 0.0);

			nodes.Add(node1);
			nodes.Add(node2);
			nodes.Add(node3);

			// Model creation
			Model model = new Model();

			// Add a single subdomain to the model
			int subdomainID = 0;
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

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
				var beam = new Beam3DCorotationalQuaternion(elementNodes, material, density, beamSection);

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
				var b = beam.MassMatrix(element);

				// Add beam element to the element and subdomains dictionary of the model
				model.ElementsDictionary.Add(element.ID, element);
				model.SubdomainsDictionary[subdomainID].Elements.Add(element);
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
			childAnalyzerBuilder.MaxIterationsPerIncrement = 120;
			childAnalyzerBuilder.NumIterationsForMatrixRebuild = 500;
			//childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
			LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
			var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 0.28, 3.36);
			parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
			NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			Assert.Equal(148.936792350562, solver.LinearSystems[subdomainID].Solution[7], 12);
		}
	}
}
