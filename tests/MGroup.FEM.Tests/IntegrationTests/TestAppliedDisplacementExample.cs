//using System;
//using System.Collections.Generic;
//using System.Text;
//using Xunit;

//namespace ISAAR.MSolve.Tests
//{
//    public class TestAppliedDisplacementExample
//    {
//        [Fact]
//        public void AppliedDisplacementTestExample()
//        {
//            IList<Node> nodes = new List<Node>();
//            Node node1 = new Node( id: 1, x: 0.0, y:  0.0, z: 0.0 };
//            Node node2 = new Node( id: 2, x: 300.0, y:  0.0, z: 0.0 };
//            nodes.Add(node1);
//            nodes.Add(node2);

//            VectorExtensions.AssignTotalAffinityCount();
//            double youngModulus = 21000;
//            double poissonRatio = 0.3;
//            double nodalLoad = 100.0;

//            // Create a new elastic 3D material
//            ElasticMaterial material = new ElasticMaterial()
//            {
//                YoungModulus = youngModulus,
//                PoissonRatio = poissonRatio,
//            };

//            // Model creation
//            Model model = new Model();

//            // Add a single subdomain to the model
//            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

//            // Add nodes to the nodes dictonary of the model
//            for (int i = 0; i < nodes.Count; ++i)
//            {
//                model.NodesDictionary.Add(i + 1, nodes[i]);
//            }

//            // Constrain bottom nodes of the model
//            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
//            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
//            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });

//            model.NodesDictionary[2].Constraints.Add(new Constraint { DOF = DOFType.Y, Amount = -4.16666666666667E-07 });

//            //Create a new Beam2D element
//            var beam = new EulerBeam2D(youngModulus)
//            {
//                SectionArea = 1,
//                MomentOfInertia = .1
//            };

//            var element = new Element()
//            {
//                ID = 1,
//                ElementType = beam
//            };

//            //// Add nodes to the created element
//            element.AddNode(model.NodesDictionary[1]);
//            element.AddNode(model.NodesDictionary[2]);

//            var a = beam.StiffnessMatrix(element);

//            // Add Hexa element to the element and subdomains dictionary of the model
//            model.ElementsDictionary.Add(element.ID, element);
//            model.SubdomainsDictionary[0].ElementsDictionary.Add(element.ID, element);

//            // Add nodal load values at the top nodes of the model
//            //model.Loads.Add(new Load() { Amount = -nodalLoad, Node = model.NodesDictionary[2], DOF = DOFType.Y });

//            // Needed in order to make all the required data structures
//            model.ConnectDataStructures();

//            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
//            linearSystems[0] = new SkylineLinearSystem(0, model.SubdomainsDictionary[0].Forces);
//            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

//            ProblemStructural provider = new ProblemStructural(model, linearSystems);

//            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
//            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

//            // Choose dof types X, Y, Z to log for node 5
//            childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] {
//            model.NodalDOFsDictionary[2][DOFType.X],
//            model.NodalDOFsDictionary[2][DOFType.RotZ]});

//            // Analyze the problem
//            parentAnalyzer.BuildMatrices();
//            parentAnalyzer.Initialize();
//            parentAnalyzer.Solve();

//            Dictionary<int, double> results = (childAnalyzer.Logs[0][0] as DOFSLog).DOFValues;
//            double[] expected = new double[] { 0, -4.16666666666667E-07, -6.25E-07 };

//            for (int i = 0; i < expected.Length - 1; i++)
//            {
//                //if (Math.Abs(expected[i] - results[i]) > 1e-14)
//                //{
//                //    throw new SystemException("Failed beam2D test, results don't coincide for dof no: " + i + ", expected displacement: " + expected[i] + ", calculated displacement: " + results[i]);
//                //}
//                Console.WriteLine(results[i]);
//            }
//            Console.WriteLine("ran beam2d2 test");
//        }
//    }
//}
