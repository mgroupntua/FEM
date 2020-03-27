using System.Collections.Generic;
using System.Linq;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.ContinuumElements;
using MGroup.FEM.Entities;
using MGroup.FEM.Structural.Elements;
using MGroup.FEM.Structural.Elements.supportiveClasses;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Commons;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;
using MGroup.MSolve.Solution;
using MGroup.NumericalAnalyzers;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.NonLinear;
using MGroup.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    public static class Tet10ContinuumNonLinearCantilever
    {
        private const int subdomainID = 0;

        [Fact]
        private static void RunTest()
        {
            double[] expectedDisplacements = GetExpectedDisplacements();
            double[] computedDisplacements = SolveModel();
            Assert.True(AreDisplacementsSame(expectedDisplacements, computedDisplacements));
        }

        private static bool AreDisplacementsSame(double[] expectedDisplacements, 
            double[] computedDisplacements)
        {
            var comparer = new ValueComparer(1E-3);

            for (int i = 0; i < expectedDisplacements.Length; i++)
            {
                if (!comparer.AreEqual(expectedDisplacements[i], computedDisplacements[i]))
                {
                    return false;
                }
            } 
           
            return true;
        }

        private static double[]  GetExpectedDisplacements()
        {
            var expectedSolitionOfIters5And12 = new double[] { 0.56783 , 0.921567, 0.615967, 1.05899 };
            return expectedSolitionOfIters5And12;
        }

        private static double[] SolveModel()
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



            // initialize the anlaysis 
            parentAnalyzer.Initialize();

            // Output
            int[] monitoredNodes = new int[] { 6, 8 };
            int[] monitoredDofs = monitoredNodes.Select(x => model.GlobalDofOrdering.GlobalFreeDofs[model.NodesDictionary[x], StructuralDof.TranslationX]).ToArray();
            var watchDofs = new Dictionary<int, int[]>();
            watchDofs.Add(subdomainID, monitoredDofs);
            var log1 = new TotalDisplacementsPerIterationLog(watchDofs);
            childAnalyzer.TotalDisplacementsPerIterationLog = log1;
            
            //run the analysis and give output
            parentAnalyzer.Solve();



            var solutionOfIters5And12 = new double[] { log1.GetTotalDisplacement(5, 0, 3), log1.GetTotalDisplacement(12, 0, 3),
                log1.GetTotalDisplacement(5, 0, 9), log1.GetTotalDisplacement(12, 0, 9) };

            return solutionOfIters5And12;
        }

        private static void BuildCantileverModel(Model model, double load_value)
        {
            //xrhsimopoiithike to  ParadeigmataElegxwnBuilder.HexaCantileverBuilder(Model model, double load_value)
            // allagh tou element kai tou material

            //ElasticMaterial3DTemp material1 = new ElasticMaterial3DTemp()
            //{
            //    YoungModulus = 1353000,
            //    PoissonRatio = 0.3,
            //};


            //VonMisesMaterial3D material1 = new VonMisesMaterial3D(1353000, 0.30, 1353000, 0.15);
            var material1 = new ElasticMaterial3D() { PoissonRatio = 0.3, YoungModulus = 1353000 };

            double[,] nodeData = new double[,] { {-0.2500000000000000,0.2500000000000000,-1.0000000000000000},
                {0.2500000000000000,0.2500000000000000,-1.0000000000000000},
                {0.2500000000000000,-0.2500000000000000,-1.0000000000000000},
                {-0.2500000000000000,-0.2500000000000000,-1.0000000000000000},
                {-0.2500000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,-0.2500000000000000,1.0000000000000000},
                {-0.2500000000000000,-0.2500000000000000,1.0000000000000000},
                {-0.2500000000000000,0.2500000000000000,0.0000000000000000},
                {0.2500000000000000,0.2500000000000000,0.0000000000000000},
                {0.2500000000000000,-0.2500000000000000,0.0000000000000000},
                {-0.2500000000000000,-0.2500000000000000,0.0000000000000000},
                {0.0000000000000000,0.2500000000000000,-1.0000000000000000},
                {0.2500000000000000,0.0000000000000000,-1.0000000000000000},
                {0.0000000000000000,-0.2500000000000000,-1.0000000000000000},
                {-0.2500000000000000,0.0000000000000000,-1.0000000000000000},
                {0.0000000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,0.0000000000000000,1.0000000000000000},
                {0.0000000000000000,-0.2500000000000000,1.0000000000000000},
                {-0.2500000000000000,0.0000000000000000,1.0000000000000000},
                {0.0000000000000000,0.2500000000000000,0.0000000000000000},
                {0.2500000000000000,0.0000000000000000,0.0000000000000000},
                {0.0000000000000000,-0.2500000000000000,0.0000000000000000},
                {-0.2500000000000000,0.0000000000000000,0.0000000000000000},
                {-0.2500000000000000,0.2500000000000000,-0.5000000000000000},
                {0.2500000000000000,0.2500000000000000,-0.5000000000000000},
                {0.2500000000000000,-0.2500000000000000,-0.5000000000000000},
                {-0.2500000000000000,-0.2500000000000000,-0.5000000000000000},
                {-0.2500000000000000,0.2500000000000000,0.5000000000000000},
                {0.2500000000000000,0.2500000000000000,0.5000000000000000},
                {0.2500000000000000,-0.2500000000000000,0.5000000000000000},
                {-0.2500000000000000,-0.2500000000000000,0.5000000000000000},
                {0.0000000000000000,0.0000000000000000,-1.0000000000000000},
                {-0.1250000000000000,0.1250000000000000,-1.0000000000000000},
                {0.1250000000000000,0.1250000000000000,-1.0000000000000000},
                {0.1250000000000000,-0.1250000000000000,-1.0000000000000000},
                {-0.1250000000000000,-0.1250000000000000,-1.0000000000000000},
                {0.0000000000000000,0.0000000000000000,1.0000000000000000},
                {-0.1250000000000000,0.1250000000000000,1.0000000000000000},
                {0.1250000000000000,0.1250000000000000,1.0000000000000000},
                {0.1250000000000000,-0.1250000000000000,1.0000000000000000},
                {-0.1250000000000000,-0.1250000000000000,1.0000000000000000},
                {0.0000000000000000,0.0000000000000000,0.0000000000000000},
                {-0.1250000000000000,0.1250000000000000,0.0000000000000000},
                {0.1250000000000000,0.1250000000000000,0.0000000000000000},
                {0.1250000000000000,-0.1250000000000000,0.0000000000000000},
                {-0.1250000000000000,-0.1250000000000000,0.0000000000000000},
                {0.0000000000000000,0.2500000000000000,-0.5000000000000000},
                {0.1250000000000000,0.2500000000000000,-0.7500000000000000},
                {-0.1250000000000000,0.2500000000000000,-0.7500000000000000},
                {-0.1250000000000000,0.2500000000000000,-0.2500000000000000},
                {0.1250000000000000,0.2500000000000000,-0.2500000000000000},
                {0.2500000000000000,0.0000000000000000,-0.5000000000000000},
                {0.2500000000000000,-0.1250000000000000,-0.7500000000000000},
                {0.2500000000000000,0.1250000000000000,-0.7500000000000000},
                {0.2500000000000000,0.1250000000000000,-0.2500000000000000},
                {0.2500000000000000,-0.1250000000000000,-0.2500000000000000},
                {0.0000000000000000,-0.2500000000000000,-0.5000000000000000},
                {-0.1250000000000000,-0.2500000000000000,-0.7500000000000000},
                {0.1250000000000000,-0.2500000000000000,-0.7500000000000000},
                {0.1250000000000000,-0.2500000000000000,-0.2500000000000000},
                {-0.1250000000000000,-0.2500000000000000,-0.2500000000000000},
                {-0.2500000000000000,0.0000000000000000,-0.5000000000000000},
                {-0.2500000000000000,0.1250000000000000,-0.7500000000000000},
                {-0.2500000000000000,-0.1250000000000000,-0.7500000000000000},
                {-0.2500000000000000,-0.1250000000000000,-0.2500000000000000},
                {-0.2500000000000000,0.1250000000000000,-0.2500000000000000},
                {0.0000000000000000,0.2500000000000000,0.5000000000000000},
                {0.1250000000000000,0.2500000000000000,0.2500000000000000},
                {-0.1250000000000000,0.2500000000000000,0.2500000000000000},
                {-0.1250000000000000,0.2500000000000000,0.7500000000000000},
                {0.1250000000000000,0.2500000000000000,0.7500000000000000},
                {0.2500000000000000,0.0000000000000000,0.5000000000000000},
                {0.2500000000000000,-0.1250000000000000,0.2500000000000000},
                {0.2500000000000000,0.1250000000000000,0.2500000000000000},
                {0.2500000000000000,0.1250000000000000,0.7500000000000000},
                {0.2500000000000000,-0.1250000000000000,0.7500000000000000},
                {0.0000000000000000,-0.2500000000000000,0.5000000000000000},
                {-0.1250000000000000,-0.2500000000000000,0.2500000000000000},
                {0.1250000000000000,-0.2500000000000000,0.2500000000000000},
                {0.1250000000000000,-0.2500000000000000,0.7500000000000000},
                {-0.1250000000000000,-0.2500000000000000,0.7500000000000000},
                {-0.2500000000000000,0.0000000000000000,0.5000000000000000},
                {-0.2500000000000000,0.1250000000000000,0.2500000000000000},
                {-0.2500000000000000,-0.1250000000000000,0.2500000000000000},
                {-0.2500000000000000,-0.1250000000000000,0.7500000000000000},
                {-0.2500000000000000,0.1250000000000000,0.7500000000000000},
                {0.0000000000000000,0.0000000000000000,-0.5000000000000000},
                {-0.1250000000000000,0.1250000000000000,-0.7500000000000000},
                {0.1250000000000000,0.1250000000000000,-0.7500000000000000},
                {0.1250000000000000,-0.1250000000000000,-0.7500000000000000},
                {-0.1250000000000000,-0.1250000000000000,-0.7500000000000000},
                {-0.1250000000000000,0.1250000000000000,-0.2500000000000000},
                {0.1250000000000000,0.1250000000000000,-0.2500000000000000},
                {0.1250000000000000,-0.1250000000000000,-0.2500000000000000},
                {-0.1250000000000000,-0.1250000000000000,-0.2500000000000000},
                {0.0000000000000000,0.0000000000000000,-0.7500000000000000},
                {0.0000000000000000,0.1250000000000000,-0.5000000000000000},
                {0.1250000000000000,0.0000000000000000,-0.5000000000000000},
                {0.0000000000000000,-0.1250000000000000,-0.5000000000000000},
                {-0.1250000000000000,0.0000000000000000,-0.5000000000000000},
                {0.0000000000000000,0.0000000000000000,-0.2500000000000000},
                {0.0000000000000000,0.0000000000000000,0.5000000000000000},
                {-0.1250000000000000,0.1250000000000000,0.2500000000000000},
                {0.1250000000000000,0.1250000000000000,0.2500000000000000},
                {0.1250000000000000,-0.1250000000000000,0.2500000000000000},
                {-0.1250000000000000,-0.1250000000000000,0.2500000000000000},
                {-0.1250000000000000,0.1250000000000000,0.7500000000000000},
                {0.1250000000000000,0.1250000000000000,0.7500000000000000},
                {0.1250000000000000,-0.1250000000000000,0.7500000000000000},
                {-0.1250000000000000,-0.1250000000000000,0.7500000000000000},
                {0.0000000000000000,0.0000000000000000,0.2500000000000000},
                {0.0000000000000000,0.1250000000000000,0.5000000000000000},
                {0.1250000000000000,0.0000000000000000,0.5000000000000000},
                {0.0000000000000000,-0.1250000000000000,0.5000000000000000},
                {-0.1250000000000000,0.0000000000000000,0.5000000000000000},
                {0.0000000000000000,0.0000000000000000,0.7500000000000000},};

            //int[] take_msolve_nodes_from_adina_tade_nodes = new int[] { 2, 1, 3, 0, 5, 7, 8, 6, 9, 4 };

            int[,] elementData = new int[,] {{33,2,1,88,35,13,34,90,89,97},
                {33,3,2,88,36,14,35,91,90,97},
                {33,4,3,88,37,15,36,92,91,97},
                {33,1,4,88,34,16,37,89,92,97},
                {48,1,2,88,50,13,49,89,90,98},
                {48,2,10,88,49,26,52,90,94,98},
                {48,10,9,88,52,21,51,94,93,98},
                {48,9,1,88,51,25,50,93,89,98},
                {53,2,3,88,55,14,54,90,91,99},
                {53,3,11,88,54,27,57,91,95,99},
                {53,11,10,88,57,22,56,95,94,99},
                {53,10,2,88,56,26,55,94,90,99},
                {58,3,4,88,60,15,59,91,92,100},
                {58,11,3,88,61,27,60,95,91,100},
                {58,12,11,88,62,23,61,96,95,100},
                {58,4,12,88,59,28,62,92,96,100},
                {63,4,1,88,65,16,64,92,89,101},
                {63,12,4,88,66,28,65,96,92,101},
                {63,9,12,88,67,24,66,93,96,101},
                {63,1,9,88,64,25,67,89,93,101},
                {43,9,10,88,44,21,45,93,94,102},
                {43,10,11,88,45,22,46,94,95,102},
                {43,11,12,88,46,23,47,95,96,102},
                {43,12,9,88,47,24,44,96,93,102},
                {43,10,9,103,45,21,44,105,104,112},
                {43,11,10,103,46,22,45,106,105,112},
                {43,12,11,103,47,23,46,107,106,112},
                {43,9,12,103,44,24,47,104,107,112},
                {68,9,10,103,70,21,69,104,105,113},
                {68,10,6,103,69,30,72,105,109,113},
                {68,6,5,103,72,17,71,109,108,113},
                {68,5,9,103,71,29,70,108,104,113},
                {73,10,11,103,75,22,74,105,106,114},
                {73,11,7,103,74,31,77,106,110,114},
                {73,7,6,103,77,18,76,110,109,114},
                {73,6,10,103,76,30,75,109,105,114},
                {78,11,12,103,80,23,79,106,107,115},
                {78,7,11,103,81,31,80,110,106,115},
                {78,8,7,103,82,19,81,111,110,115},
                {78,12,8,103,79,32,82,107,111,115},
                {83,12,9,103,85,24,84,107,104,116},
                {83,8,12,103,86,32,85,111,107,116},
                {83,5,8,103,87,20,86,108,111,116},
                {83,9,5,103,84,29,87,104,108,116},
                {38,5,6,103,39,17,40,108,109,117},
                {38,6,7,103,40,18,41,109,110,117},
                {38,7,8,103,41,19,42,110,111,117},
                {38,8,5,103,42,20,39,111,108,117},
                 };

            // orismos shmeiwn
            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode + 1, new Node(id: nNode + 1, x: nodeData[nNode, 0], y:  nodeData[nNode, 1], z: nodeData[nNode, 2] ));

            }

            // orismos elements 
            Element e1;
			var adinaOrder = new AdinaElementLocalNodeOrdering();
            int subdomainID = Tet10ContinuumNonLinearCantilever.subdomainID;
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                DynamicMaterial DynamicMaterial = new DynamicMaterial(1, 0, 0);
                //Dictionary<int,Node3D >
                List<Node> nodeSet = new List<Node>(10);
                for (int j = 0; j < 10; j++)
                {
                    int nodeID = elementData[nElement, j];
                    nodeSet.Add((Node)model.NodesDictionary[nodeID]);
                }

                var factory = new ContinuumElement3DFactory(material1, DynamicMaterial);

				var nodeSet1 = adinaOrder.ReorderNodes(nodeSet, CellType.Tet10);
				e1 = new Element()
				{
					ID = nElement + 1,
					ElementType= factory.CreateNonLinearElement(CellType.Tet10,nodeSet1,material1,DynamicMaterial)
                };

				for (int j = 0; j < 10; j++)
                {
					int LocalNode = adinaOrder.GetNodeForLocalMsolveNode(j, CellType.Tet10);
                   e1.NodesDictionary.Add(elementData[nElement, LocalNode], model.NodesDictionary[elementData[nElement, LocalNode]]);
                }

                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].Elements.Add( e1);
            }

            int[] constrainedIds = new int[] { 1,2,3,4,13, 14, 15, 16, 33, 34, 35, 36, 37 };
            
            for (int k1 = 0; k1 < constrainedIds.Length; k1++)
            {
                int k = constrainedIds[k1];
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
            }

            // fortish korufhs
            Load load1;
            for (int k = 5; k < 9; k++)
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
