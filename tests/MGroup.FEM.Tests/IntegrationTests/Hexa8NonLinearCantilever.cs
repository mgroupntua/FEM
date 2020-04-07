using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Commons;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Loads;
using Xunit;

namespace MGroup.FEM.Tests.IntegrationTests
{
	using Constitutive.Structural;
	using Constitutive.Structural.ContinuumElements;
	using ISAAR.MSolve.FEM.Elements;
	using ISAAR.MSolve.FEM.Interpolation;
	using MGroup.FEM.Structural.Elements.supportiveClasses;
	using MGroup.MSolve.Discretization.Mesh;
	using MSolve.Constitutive;
	using MSolve.Solution;
	using NumericalAnalyzers;
	using NumericalAnalyzers.Logging;
	using NumericalAnalyzers.NonLinear;
	using Solvers.Direct;
	using Structural.Elements;

	public static class Hexa8NonLinearCantilever
	{
		private const int subdomainID = 0;

		[Fact]
		private static void RunTest()
		{
			IReadOnlyList<Dictionary<int, double>> expectedDisplacements = GetExpectedDisplacements();
			TotalDisplacementsPerIterationLog computedDisplacements = SolveModel();
			Assert.True(AreDisplacementsSame(expectedDisplacements, computedDisplacements));
		}

		private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements,
			TotalDisplacementsPerIterationLog computedDisplacements)
		{
			var comparer = new ValueComparer(1E-13);
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
			var expectedDisplacements = new Dictionary<int, double>[11]; //TODO: this should be 11 EINAI ARRAY APO DICTIONARIES

			expectedDisplacements[0] = new Dictionary<int, double> {
				{ 0, 0.039075524153873623}, {11, -0.032541895181220408}, {23, -0.057387148941853101}, {35, -0.071994381984550326}, {47, -0.077053554770404833}
			};

			expectedDisplacements[0] = new Dictionary<int, double> {
	{ 0,3.907552415387362300e-02 }, {11,-3.254189518122040800e-02 }, {23,-5.738714894185310100e-02 }, {35,-7.199438198455032600e-02 }, {47,-7.705355477040483300e-02 }};
			expectedDisplacements[1] = new Dictionary<int, double> {
	{ 0,4.061313406968563400e-02 }, {11,-3.418876666892714500e-02 }, {23,-6.682708262609965400e-02 }, {35,-9.647418428408424700e-02 }, {47,-1.214556593711370000e-01 }};
			expectedDisplacements[2] = new Dictionary<int, double> {
	{ 0,4.036171804663909300e-02 }, {11,-3.396515033613205900e-02 }, {23,-6.665084050819490600e-02 }, {35,-9.713633946904017000e-02 }, {47,-1.236631490430697600e-01 }};
			expectedDisplacements[3] = new Dictionary<int, double> {
	{ 0,4.032905162001462800e-02 }, {11,-3.393260905426281900e-02 }, {23,-6.657423779424630200e-02 }, {35,-9.701032579889114200e-02 }, {47,-1.234941821043235900e-01 }};
			expectedDisplacements[4] = new Dictionary<int, double> {
	{ 0,4.032900093364350700e-02 }, {11,-3.393255831972321500e-02 }, {23,-6.657411965268195100e-02 }, {35,-9.701012513482368300e-02 }, {47,-1.234939001150344400e-01 }};
			expectedDisplacements[5] = new Dictionary<int, double> {
	{ 0,8.095088461395548400e-02 }, {11,-6.826589092291023000e-02 }, {23,-1.393261307096994000e-01 }, {35,-2.129883579558797000e-01 }, {47,-2.840192458274605800e-01 }};
			expectedDisplacements[6] = new Dictionary<int, double> {
	{ 0,8.179065808895391600e-02 }, {11,-6.914910025670165100e-02 }, {23,-1.449912527358244700e-01 }, {35,-2.283048858573358000e-01 }, {47,-3.126785624370127000e-01 }};
			expectedDisplacements[7] = new Dictionary<int, double> {
	{ 0,8.008398180684392400e-02 }, {11,-6.747544383562544000e-02 }, {23,-1.408463169597064000e-01 }, {35,-2.210877012127209200e-01 }, {47,-3.022981704019522300e-01 }};
			expectedDisplacements[8] = new Dictionary<int, double> {
	{ 0,7.976397887674688300e-02 }, {11,-6.715673915988762400e-02 }, {23,-1.400151566610138300e-01 }, {35,-2.195056794855129700e-01 }, {47,-2.998365539162924900e-01 }};
			expectedDisplacements[9] = new Dictionary<int, double> {
	{ 0,7.975945236918889600e-02 }, {11,-6.715223199537226400e-02 }, {23,-1.400036710136937400e-01 }, {35,-2.194845023343510200e-01 }, {47,-2.998046100841828000e-01 }};
			expectedDisplacements[10] = new Dictionary<int, double> {
	{ 0,7.975944951878896600e-02 }, {11,-6.715222916021290600e-02 }, {23,-1.400036636464831200e-01 }, {35,-2.194844883932760600e-01 }, {47,-2.998045884933974200e-01 }};


			return expectedDisplacements;
		}

		private static TotalDisplacementsPerIterationLog SolveModel()
		{
			var model = new Model();
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

			BuildCantileverModel(model, 850);

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
			watchDofs.Add(subdomainID, new int[5] { 0, 11, 23, 35, 47 });
			var log1 = new TotalDisplacementsPerIterationLog(watchDofs);
			childAnalyzer.TotalDisplacementsPerIterationLog = log1;

			// Run the anlaysis 
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			return log1;
		}

		private static void BuildCantileverModel(Model model, double load_value)
		{
			//xrhsimopoiithike to  ParadeigmataElegxwnBuilder.HexaCantileverBuilder(Model model, double load_value)
			// allagh tou element kai tou material

			//IContinuumMaterial3DTemp material1 = new IContinuumMaterial3DTemp()
			//{
			//    YoungModulus = 1353000,
			//    PoissonRatio = 0.3,
			//};


			//VonMisesMaterial3D material1 = new VonMisesMaterial3D(1353000, 0.30, 1353000, 0.15);
			var material1 = new ElasticMaterial3D() { PoissonRatio = 0.3, YoungModulus = 1353000 };

			double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
			{0.250000,-0.250000,-1.000000},
			{-0.250000,0.250000,-1.000000},
			{0.250000,0.250000,-1.000000},
			{-0.250000,-0.250000,-0.500000},
			{0.250000,-0.250000,-0.500000},
			{-0.250000,0.250000,-0.500000},
			{0.250000,0.250000,-0.500000},
			{-0.250000,-0.250000,0.000000},
			{0.250000,-0.250000,0.000000},
			{-0.250000,0.250000,0.000000},
			{0.250000,0.250000,0.000000},
			{-0.250000,-0.250000,0.500000},
			{0.250000,-0.250000,0.500000},
			{-0.250000,0.250000,0.500000},
			{0.250000,0.250000,0.500000},
			{-0.250000,-0.250000,1.000000},
			{0.250000,-0.250000,1.000000},
			{-0.250000,0.250000,1.000000},
			{0.250000,0.250000,1.000000}};

			int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
			{2,12,11,9,10,8,7,5,6},
			{3,16,15,13,14,12,11,9,10},
			{4,20,19,17,18,16,15,13,14}, };

			// orismos shmeiwn
			for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
			{
				model.NodesDictionary.Add(nNode + 1, new Node(id: nNode + 1, x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));

			}

			// orismos elements 
			Element e1;
			DynamicMaterial DynamicMaterial = new DynamicMaterial(1, 0, 0);
			var factory = new ContinuumElement3DFactory(material1, DynamicMaterial);
			
			int subdomainID = Hexa8NonLinearCantilever.subdomainID;
			for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
			{
				List<Node> nodeSet = new List<Node>(8);
				for (int j = 0; j < 8; j++)
				{
					int nodeID = elementData[nElement, j+1];
					nodeSet.Add((Node)model.NodesDictionary[nodeID]);
				}
				e1 = new Element()
				{
					ID = nElement + 1,
					ElementType  //= factory.CreateNonLinearElement(CellType.Hexa8, nodeSet, material1, DynamicMaterial)
								 //nnew Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8                    
					= new ContinuumElement3DNonLinear(nodeSet, material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3), InterpolationHexa8.UniqueInstance),
				};
				for (int j = 0; j < 8; j++)
				{
					e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
				}
				model.ElementsDictionary.Add(e1.ID, e1);
				model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
			}

			// constraint vashh opou z=-1
			for (int k = 1; k < 5; k++)
			{
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
				model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			}

			// fortish korufhs
			Load load1;
			for (int k = 17; k < 21; k++)
			{
				load1 = new Load()
				{
					Node = model.NodesDictionary[k],
					DOF = StructuralDof.TranslationX,
					Amount = 1 * load_value
				};
				model.Loads.Add(load1);
			}
		}
	}

}
