using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Tests.IntegrationTests;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Loads;
using Xunit;

namespace MGroup.FEM.Tests.IntegrationTests
{
	using Constitutive.Structural;
	using Constitutive.Structural.PlanarElements;
	//using MGroup.FEM.Tests.IntegrationTests.Supportive_Classes;
	using MSolve.Constitutive;
	using MSolve.Solution;
	using NumericalAnalyzers;
	using Solvers.Direct;
	using Structural.Elements;

	public class IntegrationTests
	{
		//[Fact]
		//public void TestSolveHexaCantileverBeam()
		//{
		//	var model = new Model();
		//	model.SubdomainsDictionary.Add(1, new Subdomain(1));

		//	HexaSimpleCantileverBeam.MakeCantileverBeam(model, 0, 0, 0, model.NodesDictionary.Count + 1, model.ElementsDictionary.Count + 1, 1);

		//	model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[16], DOF = StructuralDof.TranslationZ });
		//	model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[17], DOF = StructuralDof.TranslationZ });
		//	model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[18], DOF = StructuralDof.TranslationZ });
		//	model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[19], DOF = StructuralDof.TranslationZ });

		//	var solverBuilder = new SkylineSolver.Builder();
		//	ISolver solver = solverBuilder.BuildSolver(model);
		//	var provider = new ProblemStructural(model, solver);
		//	var childAnalyzer = new LinearAnalyzer(model, solver, provider);
		//	var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);
		//	//childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 47 });

		//	parentAnalyzer.Initialize();
		//	parentAnalyzer.Solve();

		//	double[] expectedDisplacements = new double[]
		//	{
		//		-0.0000025899520106, -0.0000004898560318, -0.0000031099520106, -0.0000025899520106, 0.0000004898560318,
		//		-0.0000031099520106, 0.0000025899520106, 0.0000004898560318, -0.0000031099520106, 0.0000025899520106,
		//		-0.0000004898560318, -0.0000031099520106, -0.0000045673419128, -0.0000002423136749, -0.0000107872459340,
		//		-0.0000045673419128, 0.0000002423136749, -0.0000107872459340, 0.0000045673419128, 0.0000002423136749,
		//		-0.0000107872459340, 0.0000045673419128, -0.0000002423136749, -0.0000107872459340, -0.0000057299058132,
		//		-0.0000001253780263, -0.0000216044936601, -0.0000057299058132, 0.0000001253780263, -0.0000216044936601,
		//		0.0000057299058132, 0.0000001253780263, -0.0000216044936601, 0.0000057299058132, -0.0000001253780263,
		//		-0.0000216044936601, -0.0000061325564473, -0.0000000425738760, -0.0000339869559207, -0.0000061325564473,
		//		0.0000000425738760, -0.0000339869559207, 0.0000061325564473, 0.0000000425738760, -0.0000339869559207,
		//		0.0000061325564473, -0.0000000425738760, -0.0000339869559207
		//	};

		//	for (int i = 0; i < expectedDisplacements.Length; i++)
		//		Assert.Equal(expectedDisplacements[i], solver.LinearSystems[1].Solution[i], 10);
		//}

		[Fact]
		public void SolveCantileverBeam2D()
		{
			double youngModulus = 2.0e08;
			double poissonRatio = 0.3;
			double nodalLoad = 10.0;

			var material = new ElasticMaterial()
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio,
			};

			// Node creation
			IList<Node> nodes = new List<Node>();
			Node node1 = new Node(id: 1, x: 0.0, y: 0.0, z: 0.0);
			Node node2 = new Node(id: 2, x: 5.0, y: 0.0, z: 0.0);
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
				SectionArea = 1,
				MomentOfInertia = .1
			};

			var element = new Element()
			{
				ID = 1,
				ElementType = beam
			};

			// Add nodes to the created element
			element.AddNode(model.NodesDictionary[1]);
			element.AddNode(model.NodesDictionary[2]);

			var a = beam.StiffnessMatrix(element);

			// Add Hexa element to the element and subdomains dictionary of the model
			model.ElementsDictionary.Add(element.ID, element);
			model.SubdomainsDictionary[1].Elements.Add(element);

			// Add nodal load values at the top nodes of the model
			model.Loads.Add(new Load() { Amount = -nodalLoad, Node = model.NodesDictionary[2], DOF = StructuralDof.TranslationY });

			// Solvers, providers, analyzers
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);
			var provider = new ProblemStructural(model, solver);
			var childAnalyzer = new LinearAnalyzer(model, solver, provider);
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			Assert.Equal(-2.08333333333333333e-5, solver.LinearSystems[1].Solution[1], 10);
		}

		[Fact]
		private static void SolveQuadCantileverDecompositionTest1()
		{
			#region Quad Cantilever Model
			double youngModulus = 3.0e07;
			double poissonRatio = 0.3;
			double nodalLoad = 1000;

			// Create a new elastic 2D material
			var material = new ElasticMaterial2D(StressState2D.PlaneStress)
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio
			};

			// Model creation
			var model = new Model();

			// Add a single subdomain to the model
			model.SubdomainsDictionary.Add(0, new Subdomain(0));

			// Add nodes to the nodes dictonary of the model
			int indexNode = 0;
			for (int i = 0; i < 25; i++)
			{
				for (int j = 0; j < 5; j++)
				{
					model.NodesDictionary.Add(indexNode, new Node(id: indexNode++, x: i, y: j, z: 0.0));
				}
			}

			// Constrain left nodes of the model
			for (int i = 0; i < 5; i++)
			{
				model.NodesDictionary[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			}

			int indexElement = 0;
			for (int i = 0; i < 24; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					var element = new Element()
					{
						ID = indexElement,
						ElementType = new Quad4(material)
					};
					element.AddNode(model.NodesDictionary[i * 5 + j]);
					element.AddNode(model.NodesDictionary[(i + 1) * 5 + j]);
					element.AddNode(model.NodesDictionary[(i + 1) * 5 + j + 1]);
					element.AddNode(model.NodesDictionary[i * 5 + j + 1]);
					model.ElementsDictionary.Add(indexElement, element);
					model.SubdomainsDictionary[0].Elements.Add(element);
					indexElement++;
				}
			}
			// Add nodal load values to node 3
			model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[124], DOF = StructuralDof.TranslationY });

			#endregion

			// Needed in order to make all the required data structures
			model.ConnectDataStructures();

			var domainDecomposer = new AutomaticDomainDecomposer(model, 3);
			domainDecomposer.UpdateModel();


			Dictionary<int, int[]> expectedSubdomains = new Dictionary<int, int[]>()
			{
				{ 0, new int[] {0,4,1,5,8,9,2,6,10,12,13,14,3,7,11,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}},
				{ 1, new int[] {32,36,33,37,40,41,34,38,42,44,45,46,35,39,43,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63}},
				{ 2, new int[] {64,68,65,69,72,73,66,70,74,76,77,78,67,71,75,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95}}
			};

			for (int i = 0; i < expectedSubdomains.Count; i++)
			{
				var subdomainElements = model.SubdomainsDictionary[i].Elements;
				Assert.Equal(expectedSubdomains[i].Length, model.SubdomainsDictionary[i].Elements.Count);
				for (int j = 0; j < expectedSubdomains[i].Length; j++)
				{
					Assert.Equal(expectedSubdomains[i][j], subdomainElements[j].ID);
				}
			}
		}

		[Fact]
		private static void SolveQuadCantileverDecompositionTest2()
		{
			#region Quad Cantilever Model
			double youngModulus = 3.0e07;
			double poissonRatio = 0.3;
			double nodalLoad = 1000;

			// Create a new elastic 2D material
			var material = new ElasticMaterial2D(StressState2D.PlaneStress)
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio
			};

			// Model creation
			Model model = new Model();

			// Add a single subdomain to the model
			model.SubdomainsDictionary.Add(0, new Subdomain(0));

			// Add nodes to the nodes dictonary of the model
			int indexNode = 0;
			for (int i = 0; i < 25; i++)
			{
				for (int j = 0; j < 5; j++)
				{
					model.NodesDictionary.Add(indexNode, new Node(id: indexNode++, x: i, y: j, z: 0.0));
				}
			}

			// Constrain left nodes of the model
			for (int i = 0; i < 5; i++)
			{
				model.NodesDictionary[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			}

			int indexElement = 0;
			for (int i = 0; i < 24; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					var element = new Element()
					{
						ID = indexElement,
						ElementType = new Quad4(material)
					};
					element.AddNode(model.NodesDictionary[i * 5 + j]);
					element.AddNode(model.NodesDictionary[(i + 1) * 5 + j]);
					element.AddNode(model.NodesDictionary[(i + 1) * 5 + j + 1]);
					element.AddNode(model.NodesDictionary[i * 5 + j + 1]);
					model.ElementsDictionary.Add(indexElement, element);
					model.SubdomainsDictionary[0].Elements.Add(element);
					indexElement++;
				}
			}
			// Add nodal load values to node 3
			model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[124], DOF = StructuralDof.TranslationY });

			// Needed in order to make all the required data structures
			model.ConnectDataStructures();

			#endregion

			var domainDecomposer = new AutomaticDomainDecomposer(model, 8);
			domainDecomposer.UpdateModel();


			Dictionary<int, int[]> expectedSubdomains = new Dictionary<int, int[]>()
			{
				{ 0, new int[] {0,4,1,5,8,9,2,6,10,12,13,14}},
				{ 1, new int[] {3,7,11,15,18,19,17,21,22,23,16,20}},
				{ 2, new int[] {24,28,25,29,32,33,26,30,34,36,37,38}},
				{ 3, new int[] {27,31,35,39,42,43,41,45,46,47,40,44}},
				{ 4, new int[] {48,52,49,53,56,57,50,54,58,60,61,62}},
				{ 5, new int[] {51,55,59,63,66,67,65,69,70,71,64,68}},
				{ 6, new int[] {72,76,73,77,80,81,74,78,82,84,85,86}},
				{ 7, new int[] {75,79,83,87,90,91,89,93,94,95,88,92}}
			};

			for (int i = 0; i < expectedSubdomains.Count; i++)
			{
				var subdomainElements = model.SubdomainsDictionary[i].Elements;
				Assert.Equal(expectedSubdomains[i].Length, model.SubdomainsDictionary[i].Elements.Count);
				for (int j = 0; j < expectedSubdomains[i].Length; j++)
				{
					Assert.Equal(expectedSubdomains[i][j], subdomainElements[j].ID);
				}
			}
		}

		[Fact]
		private static void SolveQuadCantileverDecompositionTest3()
		{
			#region Quad Cantilever Model
			double youngModulus = 3.0e07;
			double poissonRatio = 0.3;
			double nodalLoad = 1000;

			// Create a new elastic 2D material
			var material = new ElasticMaterial2D(StressState2D.PlaneStress)
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio
			};

			// Model creation
			Model model = new Model();

			// Add a single subdomain to the model
			model.SubdomainsDictionary.Add(0, new Subdomain(0));

			// Add nodes to the nodes dictonary of the model
			int indexNode = 0;
			for (int i = 0; i < 6; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					model.NodesDictionary.Add(indexNode, new Node(id: indexNode++, x: i, y: j, z: 0.0));
				}
			}

			// Constrain left nodes of the model
			for (int i = 0; i < 3; i++)
			{
				model.NodesDictionary[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			}

			int indexElement = 0;
			for (int i = 0; i < 5; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					var element = new Element()
					{
						ID = indexElement,
						ElementType = new Quad4(material)
					};
					element.AddNode(model.NodesDictionary[i * 3 + j]);
					element.AddNode(model.NodesDictionary[(i + 1) * 3 + j]);
					element.AddNode(model.NodesDictionary[(i + 1) * 3 + j + 1]);
					element.AddNode(model.NodesDictionary[i * 3 + j + 1]);
					model.ElementsDictionary.Add(indexElement, element);
					model.SubdomainsDictionary[0].Elements.Add(element);
					indexElement++;
				}
			}
			// Add nodal load values to node 3
			model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[17], DOF = StructuralDof.TranslationY });

			// Needed in order to make all the required data structures
			model.ConnectDataStructures();

			#endregion

			var domainDecomposer = new AutomaticDomainDecomposer(model, 3);
			domainDecomposer.UpdateModel();


			Dictionary<int, int[]> expectedSubdomains = new Dictionary<int, int[]>()
			{
				{ 0, new int[] {0,2,1,3}},
				{ 1, new int[] {4,6,5,7}},
				{ 2, new int[] {8,9}}
			};

			for (int i = 0; i < expectedSubdomains.Count; i++)
			{
				var subdomainElements = model.SubdomainsDictionary[i].Elements;
				Assert.Equal(expectedSubdomains[i].Length, model.SubdomainsDictionary[i].Elements.Count);
				for (int j = 0; j < expectedSubdomains[i].Length; j++)
				{
					Assert.Equal(expectedSubdomains[i][j], subdomainElements[j].ID);
				}
			}
		}

		[Fact]
		private static void SolveLinearTrussExample()
		{
			#region CreateGeometry

			IList<Node> nodes = new List<Node>();
			Node node1 = new Node(id: 1, x: 0, y: 0);
			Node node2 = new Node(id: 2, x: 0, y: 40);
			Node node3 = new Node(id: 3, x: 40, y: 40);

			nodes.Add(node1);
			nodes.Add(node2);
			nodes.Add(node3);

			double youngMod = 10e6;
			double poisson = 0.3;
			double loadX = 500;
			double loadY = 300;
			double sectionArea = 1.5;

			var trussModel = new Model();

			trussModel.SubdomainsDictionary.Add(0, new Subdomain(0));

			for (int i = 0; i < nodes.Count; i++)
			{
				trussModel.NodesDictionary.Add(i + 1, nodes[i]);
			}

			trussModel.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			trussModel.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			trussModel.NodesDictionary[2].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			trussModel.NodesDictionary[2].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });


			var element1 = new Element() { ID = 1, ElementType = new Rod2D(youngMod) { Density = 1, SectionArea = sectionArea } };
			var element2 = new Element() { ID = 2, ElementType = new Rod2D(youngMod) { Density = 1, SectionArea = sectionArea } };

			element1.AddNode(trussModel.NodesDictionary[1]);
			element1.AddNode(trussModel.NodesDictionary[3]);

			element2.AddNode(trussModel.NodesDictionary[2]);
			element2.AddNode(trussModel.NodesDictionary[3]);

			trussModel.ElementsDictionary.Add(element1.ID, element1);
			trussModel.ElementsDictionary.Add(element2.ID, element2);

			trussModel.SubdomainsDictionary[0].Elements.Add(element1);
			trussModel.SubdomainsDictionary[0].Elements.Add(element2);

			trussModel.Loads.Add(new Load() { Amount = loadX, Node = trussModel.NodesDictionary[3], DOF = StructuralDof.TranslationX });
			trussModel.Loads.Add(new Load() { Amount = loadY, Node = trussModel.NodesDictionary[3], DOF = StructuralDof.TranslationY });
			#endregion

			// Solvers, providers, analyzers
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(trussModel);
			var provider = new ProblemStructural(trussModel, solver);
			var childAnalyzer = new LinearAnalyzer(trussModel, solver, provider);
			var parentAnalyzer = new StaticAnalyzer(trussModel, solver, provider, childAnalyzer);

			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			Assert.Equal(0.00053333333333333336, solver.LinearSystems[0].Solution[0], 10);
			Assert.Equal(0.0017294083664636196, solver.LinearSystems[0].Solution[1], 10);
		}



	}
}
