using System;
using System.Collections.Generic;
using System.Linq;
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
	using Constitutive.Structural.ShellElements;
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
	using Structural.Embedding;

	public static class EmbeddingNLRVE
	{
		private const int subdomainID = 0;

		[Fact]
		public static void RunTest()
		{
			IReadOnlyList<Dictionary<int, double>> expectedDisplacements = GetExpectedDisplacements();
			TotalDisplacementsPerIterationLog computedDisplacements = SolveModel();
			Assert.True(AreDisplacementsSame(expectedDisplacements, computedDisplacements));
		}

		private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements,
			TotalDisplacementsPerIterationLog computedDisplacements)
		{
			var comparer = new ValueComparer(1E-9); // for node major dof order and skyline solver
													 //var comparer = new ValueComparer(1E-1); // for other solvers. It may require adjusting after visual inspection
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
			var expectedDisplacements = new Dictionary<int, double>[6]; //TODO: this should be 11 EINAI ARRAY APO DICTIONARIES

			expectedDisplacements[0] = new Dictionary<int, double> {
	{ 0,2.329148123666882600e-01 }, {11,1.745892242785442000e-04 }, {23,1.984890825094956000e-01 }, {35,-1.201717100714760900e-04 }, {47,-9.593766953486826400e-02 }};
			expectedDisplacements[1] = new Dictionary<int, double> {
	{ 0,2.312129749268104200e-01 }, {11,1.832976227505586300e-04 }, {23,1.971828675655811200e-01 }, {35,-9.890823091917878000e-05 }, {47,-9.559024674611843500e-02 }};
			expectedDisplacements[2] = new Dictionary<int, double> {
	{ 0,2.312128564190369100e-01 }, {11,1.833039372311667600e-04 }, {23,1.971829073503008300e-01 }, {35,-9.889271565360632300e-05 }, {47,-9.559022884187393100e-02 }};
			expectedDisplacements[3] = new Dictionary<int, double> {
	{ 0,4.607459982538425500e-01 }, {11,3.749252050030026000e-04 }, {23,3.930823598857800000e-01 }, {35,-1.772467305888356800e-04 }, {47,-1.908369235967369300e-01 }};
			expectedDisplacements[4] = new Dictionary<int, double> {
	{ 0,4.591105757986753100e-01 }, {11,3.830623301676184900e-04 }, {23,3.918152286679007500e-01 }, {35,-1.584825922987456900e-04 }, {47,-1.905002821626486100e-01 }};
			expectedDisplacements[5] = new Dictionary<int, double> {
	{ 0,4.591104745433216600e-01 }, {11,3.830652621728760000e-04 }, {23,3.918153586144528200e-01 }, {35,-1.584690208711596200e-04 }, {47,-1.905002649158537600e-01 }};

			return expectedDisplacements;
		}

		private static TotalDisplacementsPerIterationLog SolveModel()
		{
			var model = new Model();
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));
			Reference2RVEExample1000ddm_test_for_Msolve_release_version(model);

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

		public static void Reference2RVEExample1000ddm_test_for_Msolve_release_version(Model model)
		{
			// Origin: RVEkanoninkhsGewmetriasBuilder.Reference2RVEExample1000ddm(....)
			double[,] Dq; //TODOGerasimos this will be used TOUSE to use ox rotated creation entos MSOLVE
			Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
			rveMatrixParameters mp;
			grapheneSheetParameters gp;
			int graphene_sheets_number = 1;
			//string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_50_000_renu_new_multiple_algorithms_check_develop_1000ddm\REF_new_total_numbering.txt";
			int[] renumbering_vector =
			{30,1,4,5,7,11,
			12,8,13,14,15,17,
			18,24,31,25,32,33,
			16,19,20,26,34,35,
			27,36,37,2,3,9,
			6,10,21,38,28,22,
			29,23,39,40,41,42,
			43,44,45,46,47,48,
			49,50,51,52,53,54,
			55,56,57,58,59,60,
			61,62,63,64,65,66, };
			//string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_50_000_renu_new_multiple_algorithms_check_develop_1000ddm\Fxk_p_komvoi_rve.txt";
			double[] Fxk_p_komvoi_rve = new double[]
			{-20,0,0,0,0,0,
			 20,0,0,-40,0,0,
			 0,0,0,40,0,0,
			 -20,0,0,0,0,0,
			 20,0,0,-40,0,0,
			 0,0,0,40,0,0,
			 -80,0,0,80,0,0,
			 -40,0,0,0,0,0,
			 40,0,0,-20,0,0,
			 0,0,0,20,0,0,
			 -40,0,0,0,0,0,
			 40,0,0,-20,0,0,
			 0,0,0,20,0,0, };
			//string o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_50_000_renu_new_multiple_algorithms_check_develop_1000ddm\o_xsunol_gs_{0}.txt";
			double[][] o_xsunol = new double[graphene_sheets_number][];
			o_xsunol[0] = new double[]
			{-5.00000000000000270000,-25.98076211353316000000,-9.99999999999999820000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			-15.00000000000000200000,-8.66025403784438550000,-9.99999999999999820000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			-25.00000000000000000000,8.66025403784438910000,-9.99999999999999820000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			2.49999999999999820000,-21.65063509461096600000,-4.99999999999999910000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			-17.50000000000000000000,12.99038105676658200000,-4.99999999999999910000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			9.99999999999999820000,-17.32050807568877500000,0.00000000000000000000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			0.00000000000000000000,0.00000000000000000000,0.00000000000000000000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			-9.99999999999999820000,17.32050807568877500000,0.00000000000000000000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			17.50000000000000000000,-12.99038105676658200000,4.99999999999999910000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			-2.49999999999999820000,21.65063509461096600000,4.99999999999999910000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			25.00000000000000000000,-8.66025403784438910000,9.99999999999999820000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			15.00000000000000200000,8.66025403784438550000,9.99999999999999820000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,
			5.00000000000000270000,25.98076211353316000000,9.99999999999999820000,-0.43301270189221930000,-0.24999999999999994000,0.86602540378443871000,};
			int subdiscr1 = 1;
			int discr1 = 2;
			// int discr2 dn xrhsimopoieitai
			int discr3 = 2;
			int subdiscr1_shell = 2;
			int discr1_shell = 1;
			mpgp = GetReferenceKanonikhGewmetriaRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
			mp = mpgp.Item1; //mp.hexa1 = 6; mp.hexa2 = 6; mp.hexa3 = 6;
			gp = mpgp.Item2; gp.elem1 = 2; gp.elem2 = 1;

			o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
			double[][] ekk_xyz = new double[graphene_sheets_number][];
			//double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };


			Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
			HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector);

			int hexaElementsNumber = model.ElementsDictionary.Count();

			IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
			List<int> EmbeddedElementsIDs = new List<int>();
			int element_counter_after_Adding_sheet;
			element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
			int shellElementsNumber;

			for (int j = 0; j < graphene_sheets_number; j++)
			{
				//string file_no = (j + 1).ToString();
				//string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
				AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector, o_xsunol[j]);
				shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
																												   //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
																												   //embdeddedGroup.Concat(embdeddedGroup_adittion);
				for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
				{
					EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
				}
				element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

			}

			// model: add loads
			AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve, renumbering_vector);
			//RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
			// model: add constraints
			AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector);
			//RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

			int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
			IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
			var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
		}


		public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetReferenceKanonikhGewmetriaRveExampleParameters(int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell)
		{
			rveMatrixParameters mp;
			mp = new rveMatrixParameters()
			{
				E_disp = 3.5, //Gpa
				ni_disp = 0.4, // stather Poisson
				L01 = 95, //150, // diastaseis
				L02 = 95, //150,
				L03 = 95, //40,
				hexa1 = discr1 * subdiscr1,// diakritopoihsh
				hexa2 = discr1 * subdiscr1,
				hexa3 = discr3
			};

			grapheneSheetParameters gp;
			gp = new grapheneSheetParameters()
			{
				// parametroi shell
				E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
				ni_shell = 0.0607, // stathera poisson
				elem1 = discr1_shell * subdiscr1_shell,
				elem2 = discr1_shell * subdiscr1_shell,
				L1 = 50,// nm  // DIORTHOSI 2 graphene sheets
				L2 = 50,// nm
				L3 = 112.5096153846, // nm
				a1_shell = 0, // nm
				tk = 0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

				//parametroi cohesive epifaneias
				T_o_3 = 0.05,// Gpa = 1000Mpa = 1000N / mm2
				D_o_3 = 0.5, // nm
				D_f_3 = 4, // nm
				T_o_1 = 0.05,// Gpa
				D_o_1 = 0.5, // nm
				D_f_1 = 4, // nm
				n_curve = 1.4
			};

			Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
			return gpmp;
		}

		// temporarily commented out static method
		//public static void HexaElementsOnlyRVEwithRenumbering(Model model, rveMatrixParameters mp, double[,] Dq, string renumberingVectorPath)
		//{
		//    // Perioxh renumbering initialization 
		//    renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
		//    // perioxh renumbering initialization ews edw 

		//    // Perioxh parametroi Rve Matrix
		//    double E_disp = mp.E_disp; //Gpa
		//    double ni_disp = mp.ni_disp; // stather Poisson
		//    double L01 = mp.L01; // diastaseis
		//    double L02 = mp.L02;
		//    double L03 = mp.L03;
		//    int hexa1 = mp.hexa1;// diakritopoihsh
		//    int hexa2 = mp.hexa2;
		//    int hexa3 = mp.hexa3;
		//    // Perioxh parametroi Rve Matrix ews edw


		//    // Perioxh Gewmetria shmeiwn
		//    int nodeCounter = 0;

		//    int nodeID;
		//    double nodeCoordX;
		//    double nodeCoordY;
		//    double nodeCoordZ;
		//    int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
		//    int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
		//    int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

		//    for (int h1 = 0; h1 < hexa1 + 1; h1++)
		//    {
		//        for (int h2 = 0; h2 < hexa2 + 1; h2++)
		//        {
		//            for (int h3 = 0; h3 < hexa3 + 1; h3++)
		//            {
		//                nodeID = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
		//                nodeCoordX = -0.5 * L01 + (h1 + 1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
		//                nodeCoordy: -0.5 * L02 + (h2 + 1 - 1) * (L02 / hexa2);
		//                nodeCoordZ = -0.5 * L03 + (h3 + 1 - 1) * (L03 / hexa3);

		//                model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ });
		//                nodeCounter++;
		//            }
		//        }
		//    }
		//    // Perioxh Gewmetria shmeiwn ews edw

		//    //Perioxh Eisagwgh elements
		//    int elementCounter = 0;

		//    IContinuumMaterial3D material1 = new IContinuumMaterial3D()
		//    {
		//        YoungModulus = E_disp,
		//        PoissonRatio = ni_disp,
		//    };
		//    Element e1;
		//    int ElementID;
		//    int[] globalNodeIDforlocalNode_i = new int[8];

		//    for (int h1 = 0; h1 < hexa1; h1++)
		//    {
		//        for (int h2 = 0; h2 < hexa2; h2++)
		//        {
		//            for (int h3 = 0; h3 < hexa3; h3++)
		//            {
		//                ElementID = h1 + 1 + (h2 + 1 - 1) * hexa1 + (h3 + 1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based
		//                globalNodeIDforlocalNode_i[0] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
		//                globalNodeIDforlocalNode_i[1] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
		//                globalNodeIDforlocalNode_i[2] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
		//                globalNodeIDforlocalNode_i[3] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
		//                globalNodeIDforlocalNode_i[4] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
		//                globalNodeIDforlocalNode_i[5] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
		//                globalNodeIDforlocalNode_i[6] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
		//                globalNodeIDforlocalNode_i[7] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));

		//                e1 = new Element()
		//                {
		//                    ID = ElementID,
		//                    ElementType = new Hexa8NLRAM_1(material1, 3, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
		//                };

		//                for (int j = 0; j < 8; j++)
		//                {
		//                    e1.NodesDictionary.Add(globalNodeIDforlocalNode_i[j], model.NodesDictionary[globalNodeIDforlocalNode_i[j]]);
		//                }
		//                model.ElementsDictionary.Add(e1.ID, e1);
		//                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
		//                elementCounter++;
		//            }
		//        }
		//    }
		//    //Perioxh Eisagwgh elements

		//    //Tuple<int, int> nodeElementCounters = new Tuple<int, int>(nodeCounter, elementCounter);
		//    // change one tuple value
		//    //nodeElementCounters = new Tuple<int, int>(nodeElementCounters.Item1, elementCounter);
		//    // get one tuple value
		//    //elementCounter = nodeElementCounters.Item2;            
		//    //return nodeElementCounters;

		//    int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
		//    int f_komvoi_rve = kuvos;
		//    int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
		//    int komvos;
		//    //Dq = new double[9, 3 * p_komvoi_rve];
		//    for (int j = 0; j < p_komvoi_rve; j++)
		//    {
		//        komvos = renumbering.GetNewNodeNumbering(f_komvoi_rve + j + 1);
		//        Dq[0, 3 * j + 0] = model.NodesDictionary[komvos].X;
		//        Dq[1, 3 * j + 1] = model.NodesDictionary[komvos].Y;
		//        Dq[2, 3 * j + 2] = model.NodesDictionary[komvos].Z;
		//        Dq[3, 3 * j + 0] = model.NodesDictionary[komvos].Y;
		//        Dq[4, 3 * j + 1] = model.NodesDictionary[komvos].Z;
		//        Dq[5, 3 * j + 2] = model.NodesDictionary[komvos].X;
		//        Dq[6, 3 * j + 0] = model.NodesDictionary[komvos].Z;
		//        Dq[7, 3 * j + 1] = model.NodesDictionary[komvos].X;
		//        Dq[8, 3 * j + 2] = model.NodesDictionary[komvos].Y;
		//    }


		//}

		public static void HexaElementsOnlyRVEwithRenumbering(Model model, rveMatrixParameters mp, double[,] Dq, int[] renumberingVector)
		{
			// Perioxh renumbering initialization 
			renumbering renumbering = new renumbering(renumberingVector);
			// perioxh renumbering initialization ews edw 

			// Perioxh parametroi Rve Matrix
			double E_disp = mp.E_disp; //Gpa
			double ni_disp = mp.ni_disp; // stather Poisson
			double L01 = mp.L01; // diastaseis
			double L02 = mp.L02;
			double L03 = mp.L03;
			int hexa1 = mp.hexa1;// diakritopoihsh
			int hexa2 = mp.hexa2;
			int hexa3 = mp.hexa3;
			// Perioxh parametroi Rve Matrix ews edw


			// Perioxh Gewmetria shmeiwn
			int nodeCounter = 0;

			int nodeID;
			double nodeCoordX;
			double nodeCoordY;
			double nodeCoordZ;
			int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
			int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
			int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

			for (int h1 = 0; h1 < hexa1 + 1; h1++)
			{
				for (int h2 = 0; h2 < hexa2 + 1; h2++)
				{
					for (int h3 = 0; h3 < hexa3 + 1; h3++)
					{
						nodeID = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
						nodeCoordX = -0.5 * L01 + (h1 + 1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
						nodeCoordY = -0.5 * L02 + (h2 + 1 - 1) * (L02 / hexa2);
						nodeCoordZ = -0.5 * L03 + (h3 + 1 - 1) * (L03 / hexa3);

						model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
						nodeCounter++;
					}
				}
			}
			// Perioxh Gewmetria shmeiwn ews edw

			//Perioxh Eisagwgh elements
			int elementCounter = 0;

			var material1 = new ElasticMaterial3D
			{
				YoungModulus = E_disp,
				PoissonRatio = ni_disp,
			};
			Element e1;
			int ElementID;
			int[] globalNodeIDforlocalNode_i = new int[8];

			DynamicMaterial DynamicMaterial = new DynamicMaterial(1, 0, 0);
			var factory = new ContinuumElement3DFactory(material1, DynamicMaterial);

			for (int h1 = 0; h1 < hexa1; h1++)
			{
				for (int h2 = 0; h2 < hexa2; h2++)
				{
					for (int h3 = 0; h3 < hexa3; h3++)
					{
						ElementID = h1 + 1 + (h2 + 1 - 1) * hexa1 + (h3 + 1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based
						globalNodeIDforlocalNode_i[0] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
						globalNodeIDforlocalNode_i[1] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
						globalNodeIDforlocalNode_i[2] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
						globalNodeIDforlocalNode_i[3] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
						globalNodeIDforlocalNode_i[4] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
						globalNodeIDforlocalNode_i[5] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
						globalNodeIDforlocalNode_i[6] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
						globalNodeIDforlocalNode_i[7] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));

						List<Node> nodeSet = new List<Node>(8);
						for (int j = 0; j < 8; j++)
						{
							int nodeID1 = globalNodeIDforlocalNode_i[ j ];
							nodeSet.Add((Node)model.NodesDictionary[nodeID1]);
						}

						e1 = new Element()
						{
							ID = ElementID,
							ElementType //= factory.CreateNonLinearElement(CellType.Hexa8, nodeSet, material1, DynamicMaterial)
							=new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
							//= new ContinummElement3DNonLinear(nodeSet, material1, GaussLegendre3D.GetQuadratureWithOrder(3,3,3), InterpolationHexa8.UniqueInstance)
					};

						for (int j = 0; j < 8; j++)
						{
							e1.NodesDictionary.Add(globalNodeIDforlocalNode_i[j], model.NodesDictionary[globalNodeIDforlocalNode_i[j]]);
						}
						model.ElementsDictionary.Add(e1.ID, e1);
						model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
						elementCounter++;
					}
				}
			}
			//Perioxh Eisagwgh elements

			//Tuple<int, int> nodeElementCounters = new Tuple<int, int>(nodeCounter, elementCounter);
			// change one tuple value
			//nodeElementCounters = new Tuple<int, int>(nodeElementCounters.Item1, elementCounter);
			// get one tuple value
			//elementCounter = nodeElementCounters.Item2;            
			//return nodeElementCounters;

			int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
			int f_komvoi_rve = kuvos;
			int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
			int komvos;
			//Dq = new double[9, 3 * p_komvoi_rve];
			for (int j = 0; j < p_komvoi_rve; j++)
			{
				komvos = renumbering.GetNewNodeNumbering(f_komvoi_rve + j + 1);
				Dq[0, 3 * j + 0] = model.NodesDictionary[komvos].X;
				Dq[1, 3 * j + 1] = model.NodesDictionary[komvos].Y;
				Dq[2, 3 * j + 2] = model.NodesDictionary[komvos].Z;
				Dq[3, 3 * j + 0] = model.NodesDictionary[komvos].Y;
				Dq[4, 3 * j + 1] = model.NodesDictionary[komvos].Z;
				Dq[5, 3 * j + 2] = model.NodesDictionary[komvos].X;
				Dq[6, 3 * j + 0] = model.NodesDictionary[komvos].Z;
				Dq[7, 3 * j + 1] = model.NodesDictionary[komvos].X;
				Dq[8, 3 * j + 2] = model.NodesDictionary[komvos].Y;
			}


		}

		//temporarily comented out method
		//public static void AddGrapheneSheet_with_o_x_Input_withRenumbering(Model model, grapheneSheetParameters gp, double[] ekk_xyz, o_x_parameters o_x_parameters, string renumberingVectorPath, string o_xsunol_input_path)
		//{
		//    // Perioxh renumbering initialization 
		//    renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
		//    // perioxh renumbering initialization ews edw 

		//    // Perioxh parametroi Graphene sheet
		//    // parametroi shell
		//    double E_shell = gp.E_shell; // GPa = 1000Mpa = 1000N / mm2
		//    double ni_shell = gp.ni_shell; // stathera poisson
		//    int elem1 = gp.elem1;
		//    int elem2 = gp.elem2;
		//    double L1 = gp.L1;// nm
		//    double L2 = gp.L2;// nm
		//    double L3 = gp.L3; // nm
		//    double a1_shell = gp.a1_shell; // nm
		//    double tk = gp.tk;  // 0.0125016478913782nm

		//    //parametroi cohesive epifaneias
		//    //T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
		//    double T_o_3 = gp.T_o_3;// Gpa = 1000Mpa = 1000N / mm2
		//    double D_o_3 = gp.D_o_3; // nm
		//    double D_f_3 = gp.D_f_3; // nm

		//    double T_o_1 = gp.T_o_1;// Gpa
		//    double D_o_1 = gp.D_o_1; // nm
		//    double D_f_1 = gp.D_f_1; // nm

		//    double n_curve = gp.n_curve;
		//    // Perioxh parametroi Graphene sheet ews edw


		//    int eswterikosNodeCounter = 0;
		//    int eswterikosElementCounter = 0;
		//    int PreviousElementsNumberValue = model.ElementsDictionary.Count();
		//    int PreviousNodesNumberValue = model.NodesDictionary.Count();


		//    // Perioxh gewmetrias (orismos nodes) meshs epifaneias
		//    int new_rows = 2 * elem1 + 1;
		//    int new_lines = 2 * elem2 + 1;
		//    double[] o_xsunol;
		//    int NodeID;
		//    double nodeCoordX;
		//    double nodeCoordY;
		//    double nodeCoordZ;

		//    //o_xsunol = ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, L1, L2, elem1, elem2, a1_shell, ekk_xyz, o_x_parameters);
		//    o_xsunol = PrintUtilities.ReadVector(o_xsunol_input_path);

		//    for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
		//    {
		//        NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
		//        nodeCoordX = o_xsunol[6 * nNode + 0];
		//        nodeCoordy: o_xsunol[6 * nNode + 1];
		//        nodeCoordZ = o_xsunol[6 * nNode + 2];

		//        model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ });
		//        eswterikosNodeCounter++;
		//    }
		//    int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;
		//    // perioxh gewmetrias meshs epifaneias ews edw


		//    // perioxh orismou shell elements
		//    IContinuumMaterial3D material2 = new IContinuumMaterial3D()
		//    {
		//        YoungModulus = E_shell,
		//        PoissonRatio = ni_shell,
		//    };

		//    int elements = elem1 * elem2;
		//    int fdof_8 = 5 * (elem1 * (3 * elem2 + 2) + 2 * elem2 + 1);
		//    int komvoi_8 = fdof_8 / 5;
		//    int[,] t_shell;
		//    t_shell = topologia_shell_coh(elements, elem1, elem2, komvoi_8); // ta stoixeia tou einai 1 based to idio einai 0 based

		//    double[] Tk_vec = new double[8];
		//    double[][] VH = new double[8][];
		//    int[] midsurfaceNodeIDforlocalShellNode_i = new int[8];
		//    Element e2;
		//    int ElementID;


		//    for (int j = 0; j < 8; j++) // paxos idio gia ola telements
		//    {
		//        Tk_vec[j] = tk;
		//    }

		//    for (int nElement = 0; nElement < elements; nElement++)
		//    {
		//        ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
		//        // ta dianusmata katefthunshs allazoun analoga to element 
		//        for (int j1 = 0; j1 < 8; j1++)
		//        {
		//            midsurfaceNodeIDforlocalShellNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
		//            VH[j1] = new double[3];
		//            VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 3];
		//            VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 4];
		//            VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 5];
		//        }

		//        e2 = new Element()
		//        {
		//            ID = ElementID,
		//            //
		//            ElementType = new Shell8dispCopyGetRAM_1(material2, 3, 3, 3) //ElementType = new Shell8dispCopyGetRAM_1(material2, 3, 3, 3)
		//            {
		//                //oVn_i= new double[][] { new double [] {ElementID, ElementID }, new double [] { ElementID, ElementID } },
		//                oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 5] },},
		//                tk = Tk_vec,
		//            }
		//        };
		//        for (int j1 = 0; j1 < 8; j1++)
		//        {
		//            e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue)]);
		//        }
		//        model.ElementsDictionary.Add(e2.ID, e2);
		//        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
		//        eswterikosElementCounter++;
		//    }
		//    int arithmosShellElements = eswterikosElementCounter;
		//    // perioxh orismou shell elements ews edw

		//    // orismos shmeiwn katw strwshs
		//    for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
		//    {
		//        NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
		//        nodeCoordX = o_xsunol[6 * nNode + 0] - 0.5 * tk * o_xsunol[6 * nNode + 3];
		//        nodeCoordy: o_xsunol[6 * nNode + 1] - 0.5 * tk * o_xsunol[6 * nNode + 4];
		//        nodeCoordZ = o_xsunol[6 * nNode + 2] - 0.5 * tk * o_xsunol[6 * nNode + 5];

		//        model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ });
		//        eswterikosNodeCounter++;
		//    }
		//    //

		//    //orismos elements katw strwshs
		//    BenzeggaghKenaneCohMat material3 = new Materials.BenzeggaghKenaneCohMat()
		//    {
		//        T_o_3 = T_o_3,
		//        D_o_3 = D_o_3,
		//        D_f_3 = D_f_3,
		//        T_o_1 = T_o_1,
		//        D_o_1 = D_o_1,
		//        D_f_1 = D_f_1,
		//        n_curve = n_curve,
		//    };

		//    int[] midsurfaceNodeIDforlocalCohesiveNode_i = new int[8];
		//    for (int nElement = 0; nElement < elements; nElement++)
		//    {
		//        ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
		//        // ta dianusmata katefthunshs allazoun analoga to element 
		//        for (int j1 = 0; j1 < 8; j1++)
		//        {
		//            midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
		//            VH[j1] = new double[3];
		//            VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
		//            VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
		//            VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
		//        }

		//        e2 = new Element()
		//        {
		//            ID = ElementID,
		//            ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_11_tlk(material3, 3, 3) //ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_1(material3, 3, 3)
		//            {
		//                oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
		//                tk = Tk_vec,
		//                endeixi_element_2 = 0,
		//            }
		//        };
		//        for (int j1 = 0; j1 < 8; j1++)
		//        {
		//            e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
		//        }
		//        for (int j1 = 0; j1 < 8; j1++)
		//        {
		//            e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface),
		//                model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface)]);
		//        }
		//        model.ElementsDictionary.Add(e2.ID, e2);
		//        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
		//        eswterikosElementCounter++;
		//    }
		//    // orismos elements katw strwshs ews edw

		//    // orismos shmeiwn anw strwshs
		//    for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
		//    {
		//        NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
		//        nodeCoordX = o_xsunol[6 * nNode + 0] + 0.5 * tk * o_xsunol[6 * nNode + 3];
		//        nodeCoordy: o_xsunol[6 * nNode + 1] + 0.5 * tk * o_xsunol[6 * nNode + 4];
		//        nodeCoordZ = o_xsunol[6 * nNode + 2] + 0.5 * tk * o_xsunol[6 * nNode + 5];

		//        model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ });
		//        eswterikosNodeCounter++;
		//    }
		//    //
		//    //orismos elements anw strwshs 
		//    for (int nElement = 0; nElement < elements; nElement++)
		//    {
		//        ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
		//        // ta dianusmata katefthunshs allazoun analoga to element 
		//        for (int j1 = 0; j1 < 8; j1++)
		//        {
		//            midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
		//            VH[j1] = new double[3];
		//            VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
		//            VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
		//            VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
		//        }

		//        e2 = new Element()
		//        {
		//            ID = ElementID,
		//            ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_11_tlk(material3, 3, 3) //ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_1(material3, 3, 3)
		//            {
		//                oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
		//                                         new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
		//                tk = Tk_vec,
		//                endeixi_element_2 = 1,
		//            }
		//        };
		//        for (int j1 = 0; j1 < 8; j1++)
		//        {
		//            e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
		//        }
		//        for (int j1 = 0; j1 < 8; j1++)
		//        {
		//            e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface),
		//                model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface)]);
		//        }
		//        model.ElementsDictionary.Add(e2.ID, e2);
		//        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
		//        eswterikosElementCounter++;
		//    }
		//    // orismos elements anw strwshs ews edw

		//}

		public static void AddGrapheneSheet_with_o_x_Input_withRenumbering(Model model, grapheneSheetParameters gp, double[] ekk_xyz, o_x_parameters o_x_parameters, int[] renumberingVector, double[] o_xsunol)
		{
			// Perioxh renumbering initialization 
			renumbering renumbering = new renumbering(renumberingVector);
			// perioxh renumbering initialization ews edw 

			// Perioxh parametroi Graphene sheet
			// parametroi shell
			double E_shell = gp.E_shell; // GPa = 1000Mpa = 1000N / mm2
			double ni_shell = gp.ni_shell; // stathera poisson
			int elem1 = gp.elem1;
			int elem2 = gp.elem2;
			double L1 = gp.L1;// nm
			double L2 = gp.L2;// nm
			double L3 = gp.L3; // nm
			double a1_shell = gp.a1_shell; // nm
			double tk = gp.tk;  // 0.0125016478913782nm

			//parametroi cohesive epifaneias
			//T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
			double T_o_3 = gp.T_o_3;// Gpa = 1000Mpa = 1000N / mm2
			double D_o_3 = gp.D_o_3; // nm
			double D_f_3 = gp.D_f_3; // nm

			double T_o_1 = gp.T_o_1;// Gpa
			double D_o_1 = gp.D_o_1; // nm
			double D_f_1 = gp.D_f_1; // nm

			double n_curve = gp.n_curve;
			// Perioxh parametroi Graphene sheet ews edw


			int eswterikosNodeCounter = 0;
			int eswterikosElementCounter = 0;
			int PreviousElementsNumberValue = model.ElementsDictionary.Count();
			int PreviousNodesNumberValue = model.NodesDictionary.Count();


			// Perioxh gewmetrias (orismos nodes) meshs epifaneias
			int new_rows = 2 * elem1 + 1;
			int new_lines = 2 * elem2 + 1;
			//double[] o_xsunol;
			int NodeID;
			double nodeCoordX;
			double nodeCoordY;
			double nodeCoordZ;

			//o_xsunol = ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, L1, L2, elem1, elem2, a1_shell, ekk_xyz, o_x_parameters);
			//o_xsunol = PrintUtilities.ReadVector(o_xsunol_input_path);

			for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
			{
				NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
				nodeCoordX = o_xsunol[6 * nNode + 0];
				nodeCoordY = o_xsunol[6 * nNode + 1];
				nodeCoordZ = o_xsunol[6 * nNode + 2];

				model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
				eswterikosNodeCounter++;
			}
			int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;
			// perioxh gewmetrias meshs epifaneias ews edw


			// perioxh orismou shell elements
			//IContinuumMaterial3DTemp material2 = new IContinuumMaterial3DTemp()
			//{
			//    YoungModulus = E_shell,
			//    PoissonRatio = ni_shell,
			//};
			var material2 = new ShellElasticMaterial3D()
			{
				YoungModulus = E_shell,
				PoissonRatio = ni_shell,
				ShearCorrectionCoefficientK = 5 / 6,
			};

			int elements = elem1 * elem2;
			int fdof_8 = 5 * (elem1 * (3 * elem2 + 2) + 2 * elem2 + 1);
			int komvoi_8 = fdof_8 / 5;
			int[,] t_shell;
			t_shell = topologia_shell_coh(elements, elem1, elem2, komvoi_8); // ta stoixeia tou einai 1 based to idio einai 0 based

			double[] Tk_vec = new double[8];
			double[][] VH = new double[8][];
			int[] midsurfaceNodeIDforlocalShellNode_i = new int[8];
			Element e2;
			int ElementID;

			for (int j = 0; j < 8; j++) // paxos idio gia ola telements
			{
				Tk_vec[j] = tk;
			}

			for (int nElement = 0; nElement < elements; nElement++)
			{
				ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
				// ta dianusmata katefthunshs allazoun analoga to element 
				for (int j1 = 0; j1 < 8; j1++)
				{
					midsurfaceNodeIDforlocalShellNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
					VH[j1] = new double[3];
					VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 3];
					VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 4];
					VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 5];
				}

				e2 = new Element()
				{
					ID = ElementID,
					ElementType = new Shell8NonLinear(material2, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) //ElementType = new Shell8dispCopyGetRAM_1(material2, 3, 3, 3)
					{
						//oVn_i= new double[][] { new double [] {ElementID, ElementID }, new double [] { ElementID, ElementID } },
						oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 5] },},
						tk = Tk_vec,
					}
				};
				for (int j1 = 0; j1 < 8; j1++)
				{
					e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue)]);
				}
				model.ElementsDictionary.Add(e2.ID, e2);
				model.SubdomainsDictionary[subdomainID].Elements.Add(e2);
				eswterikosElementCounter++;
			}
			int arithmosShellElements = eswterikosElementCounter;
			// perioxh orismou shell elements ews edw

			// orismos shmeiwn katw strwshs
			for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
			{
				NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
				nodeCoordX = o_xsunol[6 * nNode + 0] - 0.5 * tk * o_xsunol[6 * nNode + 3];
				nodeCoordY = o_xsunol[6 * nNode + 1] - 0.5 * tk * o_xsunol[6 * nNode + 4];
				nodeCoordZ = o_xsunol[6 * nNode + 2] - 0.5 * tk * o_xsunol[6 * nNode + 5];

				model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
				eswterikosNodeCounter++;
			}
			//

			//orismos elements katw strwshs
			var material3 = new BenzeggaghKenaneCohesiveMaterial()
			{
				T_o_3 = T_o_3,
				D_o_3 = D_o_3,
				D_f_3 = D_f_3,
				T_o_1 = T_o_1,
				D_o_1 = D_o_1,
				D_f_1 = D_f_1,
				n_curve = n_curve,
			};

			int[] midsurfaceNodeIDforlocalCohesiveNode_i = new int[8];
			for (int nElement = 0; nElement < elements; nElement++)
			{
				ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
				// ta dianusmata katefthunshs allazoun analoga to element 
				for (int j1 = 0; j1 < 8; j1++)
				{
					midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
					VH[j1] = new double[3];
					VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
					VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
					VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
				}

				e2 = new Element()
				{
					ID = ElementID,
					ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
					{
						oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
						tk = Tk_vec,
						ShellElementSide = 0,
					}
				};
				for (int j1 = 0; j1 < 8; j1++)
				{
					e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
				}
				for (int j1 = 0; j1 < 8; j1++)
				{
					e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface),
						model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface)]);
				}
				model.ElementsDictionary.Add(e2.ID, e2);
				model.SubdomainsDictionary[subdomainID].Elements.Add(e2);
				eswterikosElementCounter++;
			}
			// orismos elements katw strwshs ews edw

			// orismos shmeiwn anw strwshs
			for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
			{
				NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
				nodeCoordX = o_xsunol[6 * nNode + 0] + 0.5 * tk * o_xsunol[6 * nNode + 3];
				nodeCoordY = o_xsunol[6 * nNode + 1] + 0.5 * tk * o_xsunol[6 * nNode + 4];
				nodeCoordZ = o_xsunol[6 * nNode + 2] + 0.5 * tk * o_xsunol[6 * nNode + 5];

				model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
				eswterikosNodeCounter++;
			}
			//
			//orismos elements anw strwshs 
			for (int nElement = 0; nElement < elements; nElement++)
			{
				ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
				// ta dianusmata katefthunshs allazoun analoga to element 
				for (int j1 = 0; j1 < 8; j1++)
				{
					midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
					VH[j1] = new double[3];
					VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
					VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
					VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
				}

				e2 = new Element()
				{
					ID = ElementID,
					ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
					{
						oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
												 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
						tk = Tk_vec,
						ShellElementSide = 1,
					}
				};
				for (int j1 = 0; j1 < 8; j1++)
				{
					e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
				}
				for (int j1 = 0; j1 < 8; j1++)
				{
					e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface),
						model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface)]);
				}
				model.ElementsDictionary.Add(e2.ID, e2);
				model.SubdomainsDictionary[subdomainID].Elements.Add(e2);
				eswterikosElementCounter++;
			}
			// orismos elements anw strwshs ews edw

		}

		//temporarily comented out method
		//public static void AddLoadsOnRveFromFile_withRenumbering(Model model, int hexa1, int hexa2, int hexa3, string vectorpath, string renumberingVectorPath)
		//{
		//    // Perioxh renumbering initialization 
		//    renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
		//    // perioxh renumbering initialization ews edw 

		//    int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
		//    double[] Fxk_p_komvoi_rve;
		//    //Fxk_p_komvoi_rve = PrintUtilities.ReadVector(@"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\elegxos_alalgwn_fe2_tax_me1_arxiko_chol_dixws_me1_OneElementRVECheckExample\Fxk_p_komvoi_rve.txt");
		//    Fxk_p_komvoi_rve = PrintUtilities.ReadVector(vectorpath);
		//    int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
		//    int f_komvoi_rve = kuvos;
		//    int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
		//    int komvos;

		//    Load load_i;
		//    for (int j = 0; j < p_komvoi_rve; j++)
		//    {
		//        komvos = f_komvoi_rve + j + 1;
		//        load_i = new Load()
		//        {
		//            Node = model.NodesDictionary[renumbering.GetNewNodeNumbering(komvos)],
		//            DOF = DOFType.X,
		//            Amount = Fxk_p_komvoi_rve[3 * (j) + 0]
		//        };
		//        model.Loads.Add(load_i);

		//        load_i = new Load()
		//        {
		//            Node = model.NodesDictionary[renumbering.GetNewNodeNumbering(komvos)],
		//            DOF = DOFType.Y,
		//            Amount = Fxk_p_komvoi_rve[3 * (j) + 1]
		//        };
		//        model.Loads.Add(load_i);

		//        load_i = new Load()
		//        {
		//            Node = model.NodesDictionary[renumbering.GetNewNodeNumbering(komvos)],
		//            DOF = DOFType.Z,
		//            Amount = Fxk_p_komvoi_rve[3 * (j) + 2]
		//        };
		//        model.Loads.Add(load_i);
		//    }

		//    // Afairesh fortiwn apo tous desmevmenous vathmous eleftherias 
		//    int nodeID;
		//    int[] supportedDOFs = new int[9];
		//    int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
		//    int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
		//    nodeID = Topol_rve(1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
		//    supportedDOFs[0] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
		//    supportedDOFs[1] = 3 * (nodeID - f_komvoi_rve - 1) + 1;
		//    supportedDOFs[2] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

		//    nodeID = Topol_rve(hexa1 + 1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
		//    supportedDOFs[3] = 3 * (nodeID - f_komvoi_rve - 1) + 1;
		//    supportedDOFs[4] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

		//    nodeID = Topol_rve(1, hexa2 + 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
		//    supportedDOFs[5] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
		//    supportedDOFs[6] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

		//    nodeID = Topol_rve(1, 1, hexa3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
		//    supportedDOFs[7] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
		//    supportedDOFs[8] = 3 * (nodeID - f_komvoi_rve - 1) + 1;

		//    for (int j = 0; j < 9; j++)
		//    {
		//        model.Loads.RemoveAt(supportedDOFs[8 - j]); // afairoume apo pisw pros ta mpros gia na mh xalaei h thesh twn epomenwn pou tha afairethoun
		//    }

		//}

		public static void AddLoadsOnRveFromFile_withRenumbering(Model model, int hexa1, int hexa2, int hexa3, double[] Fxk_p_komvoi_rve, int[] renumberingVector)
		{
			// Perioxh renumbering initialization 
			renumbering renumbering = new renumbering(renumberingVector);
			// perioxh renumbering initialization ews edw 

			int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
			//double[] Fxk_p_komvoi_rve;
			//Fxk_p_komvoi_rve = PrintUtilities.ReadVector(@"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\elegxos_alalgwn_fe2_tax_me1_arxiko_chol_dixws_me1_OneElementRVECheckExample\Fxk_p_komvoi_rve.txt");
			//Fxk_p_komvoi_rve = PrintUtilities.ReadVector(vectorpath);
			int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
			int f_komvoi_rve = kuvos;
			int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
			int komvos;

			Load load_i;
			for (int j = 0; j < p_komvoi_rve; j++)
			{
				komvos = f_komvoi_rve + j + 1;
				load_i = new Load()
				{
					Node = model.NodesDictionary[renumbering.GetNewNodeNumbering(komvos)],
					DOF = StructuralDof.TranslationX,
					Amount = Fxk_p_komvoi_rve[3 * (j) + 0]
				};
				model.Loads.Add(load_i);

				load_i = new Load()
				{
					Node = model.NodesDictionary[renumbering.GetNewNodeNumbering(komvos)],
					DOF = StructuralDof.TranslationY,
					Amount = Fxk_p_komvoi_rve[3 * (j) + 1]
				};
				model.Loads.Add(load_i);

				load_i = new Load()
				{
					Node = model.NodesDictionary[renumbering.GetNewNodeNumbering(komvos)],
					DOF = StructuralDof.TranslationZ,
					Amount = Fxk_p_komvoi_rve[3 * (j) + 2]
				};
				model.Loads.Add(load_i);
			}

			// Afairesh fortiwn apo tous desmevmenous vathmous eleftherias 
			int nodeID;
			int[] supportedDOFs = new int[9];
			int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
			int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
			nodeID = Topol_rve(1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
			supportedDOFs[0] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
			supportedDOFs[1] = 3 * (nodeID - f_komvoi_rve - 1) + 1;
			supportedDOFs[2] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

			nodeID = Topol_rve(hexa1 + 1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
			supportedDOFs[3] = 3 * (nodeID - f_komvoi_rve - 1) + 1;
			supportedDOFs[4] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

			nodeID = Topol_rve(1, hexa2 + 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
			supportedDOFs[5] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
			supportedDOFs[6] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

			nodeID = Topol_rve(1, 1, hexa3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
			supportedDOFs[7] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
			supportedDOFs[8] = 3 * (nodeID - f_komvoi_rve - 1) + 1;

			for (int j = 0; j < 9; j++)
			{
				model.Loads.RemoveAt(supportedDOFs[8 - j]); // afairoume apo pisw pros ta mpros gia na mh xalaei h thesh twn epomenwn pou tha afairethoun
			}

		}

		//temporarily comented out method
		//public static void AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(Model model, int hexa1, int hexa2, int hexa3, string renumberingVectorPath)
		//{
		//    // Perioxh renumbering initialization 
		//    renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
		//    // perioxh renumbering initialization ews edw 

		//    int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
		//    int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
		//    int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
		//    int nodeID;

		//    nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
		//    model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
		//    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
		//    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);

		//    nodeID = renumbering.GetNewNodeNumbering(Topol_rve(hexa1 + 1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
		//    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
		//    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);

		//    nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, hexa2 + 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
		//    model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
		//    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);

		//    nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, 1, hexa3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
		//    model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
		//    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);

		//}

		public static void AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(Model model, int hexa1, int hexa2, int hexa3, int[] renumberingVector)
		{
			// Perioxh renumbering initialization 
			renumbering renumbering = new renumbering(renumberingVector);
			// perioxh renumbering initialization ews edw 

			int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
			int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
			int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
			int nodeID;

			nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
			model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });

			nodeID = renumbering.GetNewNodeNumbering(Topol_rve(hexa1 + 1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
			model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });

			nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, hexa2 + 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
			model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });

			nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, 1, hexa3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
			model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });

		}

		public static int Topol_rve(int h1, int h2, int h3, int hexa1, int hexa2, int hexa3, int kuvos, int endiam_plaka, int katw_plaka)
		{
			int arith;
			if (h3 == 1)
			{ arith = h1 + (h2 - 1) * (hexa1 + 1) + kuvos; }
			else
			{
				if (h3 == hexa3 + 1)
				{ arith = hexa3 * (hexa1 + 1) * (hexa2 + 1) + h1 + (h2 - 1) * (hexa1 + 1); }
				else
				{
					if (h2 == 1)
					{ arith = (h3 - 2) * endiam_plaka + kuvos + katw_plaka + h1; }
					else
					{
						if (h2 == hexa2 + 1)
						{ arith = (h3 - 2) * endiam_plaka + kuvos + katw_plaka + (hexa1 + 1) + 2 * (hexa2 - 1) + h1; }
						else
						{
							if (h1 == 1)
							{ arith = kuvos + katw_plaka + (h3 - 2) * endiam_plaka + (hexa1 + 1) + (h2 - 2) * 2 + 1; }
							else
							{
								if (h1 == hexa1 + 1)
								{ arith = kuvos + katw_plaka + (h3 - 2) * endiam_plaka + (hexa1 + 1) + (h2 - 2) * 2 + 2; }
								else
								{ arith = (h1 - 1) + (h2 - 2) * (hexa1 - 1) + (h3 - 2) * (hexa1 - 1) * (hexa2 - 1); }
							}
						}
					}

				}
			}
			return arith;
		}

		private static int[,] topologia_shell_coh(int elements, int elem1, int elem2, object komvoi_8)
		{
			int elem;
			int[,] t_shell = new int[elements, 8];
			for (int nrow = 0; nrow < elem1; nrow++)
			{
				for (int nline = 0; nline < elem2; nline++)
				{
					elem = (nrow + 1 - 1) * elem2 + nline + 1;//nrow+ 1 nline+1 einai zero based 
					t_shell[elem - 1, -1 + 1] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 3;
					t_shell[elem - 1, -1 + 8] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 2;
					t_shell[elem - 1, -1 + 4] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 1;

					t_shell[elem - 1, -1 + 5] = (nrow + 1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline + 1 - 1) * 1 + 2;
					t_shell[elem - 1, -1 + 7] = (nrow + 1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline + 1 - 1) * 1 + 1;

					t_shell[elem - 1, -1 + 2] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 3;
					t_shell[elem - 1, -1 + 6] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 2;
					t_shell[elem - 1, -1 + 3] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 1;
				}
			}
			return t_shell;

		}
	}

	public class rveMatrixParameters
	{
		public double E_disp { get; set; }
		public double ni_disp { get; set; }
		public double L01 { get; set; }
		public double L02 { get; set; }
		public double L03 { get; set; }
		public int hexa1 { get; set; }
		public int hexa2 { get; set; }
		public int hexa3 { get; set; }

		public rveMatrixParameters()
		{

		}
		public rveMatrixParameters(double E_disp, double ni_disp, double L01, double L02, double L03, int hexa1, int hexa2, int hexa3)
		{
			this.E_disp = E_disp;
			this.ni_disp = ni_disp;
			this.L01 = L01;
			this.L02 = L02;
			this.L03 = L03;
			this.hexa1 = hexa1;
			this.hexa2 = hexa2;
			this.hexa3 = hexa3;
		}
	}

	public class o_x_parameters
	{
		//public double E_disp { get; set; }
		//public double ni_disp { get; set; }
		//public double L01 { get; set; }
		//public double L02 { get; set; }
		//public double L03 { get; set; }
		//public int hexa1 { get; set; }
		//public int hexa2 { get; set; }
		//public int hexa3 { get; set; }

		public o_x_parameters()
		{

		}
		public o_x_parameters(double E_disp, double ni_disp, double L01, double L02, double L03, int hexa1, int hexa2, int hexa3)
		{
			//this.E_disp = E_disp;
			//this.ni_disp = ni_disp;
			//this.L01 = L01;
			//this.L02 = L02;
			//this.L03 = L03;
			//this.hexa1 = hexa1;
			//this.hexa2 = hexa2;
			//this.hexa3 = hexa3;
		}
	}

	public class grapheneSheetParameters
	{
		// parametroi shell
		public double E_shell; // GPa = 1000Mpa = 1000N / mm2
		public double ni_shell; // stathera poisson
		public int elem1;
		public int elem2;
		public double L1;// nm
		public double L2;// nm
		public double L3; // nm
		public double a1_shell; // nm
		public double tk;  // 0.0125016478913782nm
						   //parametroi cohesive epifaneias
		public double T_o_3;// Gpa = 1000Mpa = 1000N / mm2
		public double D_o_3; // nm
		public double D_f_3; // nm
		public double T_o_1;// Gpa
		public double D_o_1; // nm
		public double D_f_1; // nm
		public double n_curve = 1.4;

		public grapheneSheetParameters()
		{

		}
		public grapheneSheetParameters(double E_shell, double ni_shell, int elem1, int elem2, double L1, double L2, double L3, double a1_shell, double tk,
			double T_o_3, double D_o_3, double D_f_3, double T_o_1, double D_o_1, double D_f_1, double n_curve)
		{
			this.E_shell = E_shell; // GPa = 1000Mpa = 1000N / mm2
			this.ni_shell = ni_shell; // stathera poisson
			this.elem1 = elem1;
			this.elem2 = elem2;
			this.L1 = L1;// nm
			this.L2 = L2;// nm
			this.L3 = L3; // nm
			this.a1_shell = a1_shell; // nm
			this.tk = tk;  // 0.0125016478913782nm

			//parametroi cohesive epifaneias
			//T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
			this.T_o_3 = T_o_3;// Gpa = 1000Mpa = 1000N / mm2
			this.D_o_3 = D_o_3; // nm
			this.D_f_3 = D_f_3; // nm

			this.T_o_1 = T_o_1;// Gpa
			this.D_o_1 = D_o_1; // nm
			this.D_f_1 = D_f_1; // nm

			this.n_curve = n_curve;
		}
	}

	public class renumbering
	{
		public int[] sunol_nodes_numbering { get; set; }

		public renumbering()
		{

		}
		public renumbering(int[] sunol_nodes_numbering)
		{
			this.sunol_nodes_numbering = sunol_nodes_numbering;
		}

		public int GetNewNodeNumbering(int initial_node_number)
		{
			return sunol_nodes_numbering[initial_node_number - 1];
		}

	}
}
