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
	using Constitutive.Structural.ShellElements;
	using MSolve.Solution;
	using NumericalAnalyzers;
	using NumericalAnalyzers.Logging;
	using NumericalAnalyzers.NonLinear;
	using Solvers.Direct;
	using Structural.Elements;

	public static class Shell8NonLinearPlate
	{
		private const int subdomainID = 1;

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
			var comparer = new ValueComparer(1E-11);
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
			var expectedDisplacements = new Dictionary<int, double>[4]; //TODO: this should be 11 EINAI ARRAY APO DICTIONARIES

			// MIKRH AKRIVEIA TWN GAUSS POINT
			//expectedDisplacements[0] = new Dictionary<int, double> {
			//{ 0,1.341690670391555900e-21 }, {11,5.426715620119622300e-22 }, {23,3.510351163876216100e-19 }, {35,6.417697126233693400e-23 }, {47,9.703818759449350300e-06 }};
			//expectedDisplacements[1] = new Dictionary<int, double> {
			//{ 0,-1.766631410417230800e-10 }, {11,-1.766623168621910600e-10 }, {23,-1.384381502841035000e-16 }, {35,6.203787968704410500e-17 }, {47,9.703818433005175000e-06 }};
			//expectedDisplacements[2] = new Dictionary<int, double> {
			//{ 0,-5.299892610611580200e-10 }, {11,-5.299884368826366100e-10 }, {23,-1.389870365102331700e-16 }, {35,6.203729982771704200e-17 }, {47,1.940763607980495300e-05 }};
			//expectedDisplacements[3] = new Dictionary<int, double> {
			//{ 0,-7.066521335188173400e-10 }, {11,-7.066513887514916600e-10 }, {23,6.902764393184828500e-16 }, {35,-6.400975456584840300e-20 }, {47,1.940763490956924100e-05 }};

			// MEGALH AKRIVEIA TWN GAUSS POINT
			expectedDisplacements[0] = new Dictionary<int, double> {
	{ 0,-2.609907246226515400e-22 }, {11,1.413528855321975800e-22 }, {23,-2.458757134937532800e-19 }, {35,2.289334051771179900e-22 }, {47,9.703818759449467200e-06 }};
			expectedDisplacements[1] = new Dictionary<int, double> {
	{ 0,-1.766635527587974300e-10 }, {11,-1.766634689918627600e-10 }, {23,9.157257237104792300e-17 }, {35,-2.559311444145733000e-16 }, {47,9.703818432907000400e-06 }};
			expectedDisplacements[2] = new Dictionary<int, double> {
	{ 0,-5.299896727797873100e-10 }, {11,-5.299895890111070100e-10 }, {23,9.124758457251682500e-17 }, {35,-2.559323845353319000e-16 }, {47,1.940763607970688300e-05 }};
			expectedDisplacements[3] = new Dictionary<int, double> {
	{ 0,-7.066535263910381200e-10 }, {11,-7.066531664241640700e-10 }, {23,4.128219398586412200e-16 }, {35,2.340064775305142000e-18 }, {47,1.940763490936303600e-05 }};

			return expectedDisplacements;
		}

		private static TotalDisplacementsPerIterationLog SolveModel()
		{
			var model = new Model();
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

			ShellPlateBuilder(model, 1);

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

		private static void ShellPlateBuilder(Model model, double load_value)
		{
			// Origin: branch master idio onoma ParadeigmataElegxwnBuilder2.ShellPlateBuilder(Model model, double load_value)
			var material1 = new ShellElasticMaterial3D()
			{
				YoungModulus = 135300,
				PoissonRatio = 0.3,
				ShearCorrectionCoefficientK = 5 / 6,
			};

			double[,] nodeData = new double[,] {{10.000000,10.000000,0.000000},
			{7.500000,10.000000,0.000000},
			{5.000000,10.000000,0.000000},
			{2.500000,10.000000,0.000000},
			{0.000000,10.000000,0.000000},
			{10.000000,7.500000,0.000000},
			{7.500000,7.500000,0.000000},
			{5.000000,7.500000,0.000000},
			{2.500000,7.500000,0.000000},
			{0.000000,7.500000,0.000000},
			{10.000000,5.000000,0.000000},
			{7.500000,5.000000,0.000000},
			{5.000000,5.000000,0.000000},
			{2.500000,5.000000,0.000000},
			{0.000000,5.000000,0.000000},
			{10.000000,2.500000,0.000000},
			{7.500000,2.500000,0.000000},
			{5.000000,2.500000,0.000000},
			{2.500000,2.500000,0.000000},
			{0.000000,2.500000,0.000000},
			{10.000000,0.000000,0.000000},
			{7.500000,0.000000,0.000000},
			{5.000000,0.000000,0.000000},
			{2.500000,0.000000,0.000000},
			{0.000000,0.000000,0.000000},
			{8.750000,10.000000,0.000000},
			{6.250000,10.000000,0.000000},
			{3.750000,10.000000,0.000000},
			{1.250000,10.000000,0.000000},
			{8.750000,7.500000,0.000000},
			{6.250000,7.500000,0.000000},
			{3.750000,7.500000,0.000000},
			{1.250000,7.500000,0.000000},
			{8.750000,5.000000,0.000000},
			{6.250000,5.000000,0.000000},
			{3.750000,5.000000,0.000000},
			{1.250000,5.000000,0.000000},
			{8.750000,2.500000,0.000000},
			{6.250000,2.500000,0.000000},
			{3.750000,2.500000,0.000000},
			{1.250000,2.500000,0.000000},
			{8.750000,0.000000,0.000000},
			{6.250000,0.000000,0.000000},
			{3.750000,0.000000,0.000000},
			{1.250000,0.000000,0.000000},
			{10.000000,8.750000,0.000000},
			{10.000000,6.250000,0.000000},
			{10.000000,3.750000,0.000000},
			{10.000000,1.250000,0.000000},
			{7.500000,8.750000,0.000000},
			{7.500000,6.250000,0.000000},
			{7.500000,3.750000,0.000000},
			{7.500000,1.250000,0.000000},
			{5.000000,8.750000,0.000000},
			{5.000000,6.250000,0.000000},
			{5.000000,3.750000,0.000000},
			{5.000000,1.250000,0.000000},
			{2.500000,8.750000,0.000000},
			{2.500000,6.250000,0.000000},
			{2.500000,3.750000,0.000000},
			{2.500000,1.250000,0.000000},
			{0.000000,8.750000,0.000000},
			{0.000000,6.250000,0.000000},
			{0.000000,3.750000,0.000000},
			{0.000000,1.250000,0.000000}, };

			int[,] elementData = new int[,] { {1,1,2,7,6,26,50,30,46},
			{2,2,3,8,7,27,54,31,50},
			{3,3,4,9,8,28,58,32,54},
			{4,4,5,10,9,29,62,33,58},
			{5,6,7,12,11,30,51,34,47},
			{6,7,8,13,12,31,55,35,51},
			{7,8,9,14,13,32,59,36,55},
			{8,9,10,15,14,33,63,37,59},
			{9,11,12,17,16,34,52,38,48},
			{10,12,13,18,17,35,56,39,52},
			{11,13,14,19,18,36,60,40,56},
			{12,14,15,20,19,37,64,41,60},
			{13,16,17,22,21,38,53,42,49},
			{14,17,18,23,22,39,57,43,53},
			{15,18,19,24,23,40,61,44,57},
			{16,19,20,25,24,41,65,45,61},};

			// orismos shmeiwn
			for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
			{
				model.NodesDictionary.Add(nNode + 1, new Node(id: nNode + 1, x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));
			}

			// orismos elements 
			Element e1;
			int subdomainID = 1;
			double tk_shell_plate = 0.5;
			for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
			{
				e1 = new Element()
				{
					ID = nElement + 1,
					ElementType = new Shell8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 2))//ElementType = new Shell8dispCopyGet(material2, 3, 3, 3)
					{
						//oVn_i= new double[][] { new double [] {ElementID, ElementID }, new double [] { ElementID, ElementID } },
						oVn_i = new double[][] { new double[] { 0,0,1 },
												 new double[] { 0,0,1 },
												 new double[] { 0,0,1 },
												 new double[] { 0,0,1 },
												 new double[] { 0,0,1 },
												 new double[] { 0,0,1 },
												 new double[] { 0,0,1 },
												 new double[] { 0,0,1 },},
						tk = new double[] { tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate },
					}
				};
				for (int j = 0; j < 8; j++)
				{
					e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
				}
				model.ElementsDictionary.Add(e1.ID, e1);
				model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
			}

			// constraint paaktwsh gurw gurw plevres
			int pointID;
			int[] cnstrnd = new int[] { 21, 22, 23, 24, 25, 26, 27, 28, 29, 1, 2, 3, 4, 5, 42, 43, 44, 45, 46, 47, 48, 49, 6, 11, 16, 10, 15, 20, 62, 63, 64, 65 };
			for (int k = 0; k < cnstrnd.GetLength(0); k++)
			{
				pointID = cnstrnd[k];
				model.NodesDictionary[pointID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[pointID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
				model.NodesDictionary[pointID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
				model.NodesDictionary[pointID].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
				model.NodesDictionary[pointID].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
			}

			// fortish korufhs
			Load load1;

			load1 = new Load()
			{
				Node = model.NodesDictionary[13],
				DOF = StructuralDof.TranslationZ,
				Amount = 1 * load_value
			};
			model.Loads.Add(load1);

		}
	}

}
