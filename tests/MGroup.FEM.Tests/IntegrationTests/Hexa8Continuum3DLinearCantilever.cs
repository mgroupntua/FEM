using System.Collections.Generic;
using ISAAR.MSolve.Logging;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.ContinuumElements;
using MGroup.FEM.Entities;
using MGroup.FEM.Structural.Elements;
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

namespace ISAAR.MSolve.FEM.Tests.Elements
{
    public static class Hexa8Continuum3DLinearCantilever
    {
        private const int subdomainID = 1;

        [Fact]
        private static void RunTest()
        {
            IReadOnlyList<Dictionary<int, double>> expectedDisplacements = GetExpectedDisplacements();
            IncrementalDisplacementsLog computedDisplacements = SolveModel();
            Assert.True(AreDisplacementsSame(expectedDisplacements, computedDisplacements,1e-9));
        }

        private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements, TotalDisplacementsPerIterationLog computedDisplacements)
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

        private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements, IncrementalDisplacementsLog computedDisplacements,double tolerance)
        {
            var comparer = new ValueComparer(tolerance);
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
            var expectedDisplacements = new Dictionary<int, double>[2]; //TODO: this should be 11 EINAI ARRAY APO DICTIONARIES

            expectedDisplacements[0] = new Dictionary<int, double> {
                { 0, 0.039075524153873623}, {11, -0.032541895181220408}, {23, -0.057387148941853101}, {35, -0.071994381984550326}, {47, -0.077053554770404833}
            };

            
            expectedDisplacements[1] = new Dictionary<int, double> {
    { 0,2* 0.039075524153873623}, {11,2*( -0.032541895181220408)}, {23,2*( -0.057387148941853101)}, {35,2*( -0.071994381984550326)}, {47,2*( -0.077053554770404833)}};
            

            return expectedDisplacements;
        }

        private static IncrementalDisplacementsLog SolveModel()
        {
            //VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
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
            var log1 = new IncrementalDisplacementsLog(watchDofs);
            childAnalyzer.IncrementalDisplacementsLog = log1;

            // Run the anlaysis 
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();



            return log1;
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
            double correction = 10;// +20;
            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode + 1, new Node( nNode + 1, nodeData[nNode, 0]+correction, nodeData[nNode, 1] + correction, nodeData[nNode, 2] + correction));

            }

            // orismos elements 
            Element e1;
            int subdomainID = Hexa8Continuum3DLinearCantilever.subdomainID;
            int[] ContinuumHexa8NodesNumbering = new int[8] { 7, 8, 5, 6, 3,4,1,2 };// { 2, 3, 7, 6, 1, 4, 8, 5 };
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                DynamicMaterial DynamicMaterial = new DynamicMaterial(1, 0, 0);
                //Dictionary<int,Node3D >
                List<Node> nodeSet = new List<Node>(8);
                for (int j = 0; j < 8; j++)
                {
                    int nodeID = elementData[nElement, j+1];
                    nodeSet.Add((Node)model.NodesDictionary[nodeID]);
                }

				
                var factory = new ContinuumElement3DFactory(material1, DynamicMaterial);

                //e1 = factory.CreateElement(CellType3D.Hexa8, nodeSet); 
                e1 = new Element()
                {
                    ID = nElement + 1,
                    ElementType = factory.CreateElement(CellType.Hexa8, nodeSet)  // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8                    
                };
                for (int j = 0; j < 8; j++)
                {
                    
                    int nodeID = elementData[nElement, j+1];
                    e1.NodesDictionary.Add(nodeID, model.NodesDictionary[nodeID]);
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
