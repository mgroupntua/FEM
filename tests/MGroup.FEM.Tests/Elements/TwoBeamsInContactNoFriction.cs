using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Xunit;

using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Mesh;
using MGroup.FEM.Structural.Contact.Elements;
using MGroup.FEM.Structural.Elements;
using MGroup.FEM.Interfaces;
using MGroup.Constitutive.Structural;
using MGroup.MSolve.Solution;
using MGroup.Solvers.Iterative;
using MGroup.MSolve.Discretization.Loads;
using MGroup.Solvers.Direct;
using MGroup.NumericalAnalyzers;
using MGroup.NumericalAnalyzers.NonLinear;
using MGroup.NumericalAnalyzers.Logging;

namespace MGroup.FEM.Tests.Elements
{
	public class TwoBeamsInContactNoFriction
	{
        static Dictionary<int, Node> nodes { get; set; }
        static Dictionary<int, Element> elements { get; set; }
        static Model model { get; set; }

        private static void CreateModel()
        {
            model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain(1));
            CreateNodes();
            for (int i = 1; i <= nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i, nodes[i]);
            }
            CreateConnectivity();

            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });
            model.NodesDictionary[6].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
            model.NodesDictionary[6].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            model.NodesDictionary[6].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });

			model.Loads.Add(new Load() { Amount = -2200000, Node = model.NodesDictionary[3], DOF = StructuralDof.TranslationY });
        }

        private static void CreateNodes()
        {
            nodes = new Dictionary<int, Node>();
            nodes[1] = new Node(1, 0.0, 0.05);
            nodes[2] = new Node(2, 0.3, 0.05);
            nodes[3] = new Node(3, 0.6, 0.05);
            nodes[4] = new Node(4, 0.6, 0.0);
            nodes[5] = new Node(5, 0.9, 0.0);
            nodes[6] = new Node(6, 1.1, 0.0);
		}

		private static void CreateConnectivity()
        {
            double E = 200.0e9;
            double A = 0.01;
            double I = 8.333e-6;
            elements = new Dictionary<int, Element>();

            for (int i=1; i<=2; i++)
			{
                elements[i] = new Element()
                {
                    ID = i,
                    ElementType = new EulerBeam2D(E) { SectionArea = A, MomentOfInertia = I }
                };
                elements[i].AddNode(nodes[i]);
                elements[i].AddNode(nodes[i + 1]);

                model.ElementsDictionary.Add(elements[i].ID, elements[i]);
                model.SubdomainsDictionary[1].Elements.Add(elements[i]);
            }
            for (int i = 3; i <= 4; i++)
            {
                elements[i] = new Element()
                {
                    ID = i,
                    ElementType = new EulerBeam2D(E) { SectionArea = A, MomentOfInertia = I }
                };
                elements[i].AddNode(nodes[i + 1]);
                elements[i].AddNode(nodes[i + 2]);

                model.ElementsDictionary.Add(elements[i].ID, elements[i]);
                model.SubdomainsDictionary[1].Elements.Add(elements[i]);
            }

            elements[5] = new Element()
            {
                ID = 5,
                ElementType = new NodeToNodeContact2D(E, A, I)
            };
            elements[5].AddNode(nodes[4]);
            elements[5].AddNode(nodes[3]);

            model.ElementsDictionary.Add(elements[5].ID, elements[5]);
            model.SubdomainsDictionary[1].Elements.Add(elements[5]);

			elements[6] = new Element()
			{
				ID = 6,
				ElementType = new NodeToNodeContact2D(E, A, I)
			};
			elements[6].AddNode(nodes[4]);
			elements[6].AddNode(nodes[5]);

			model.ElementsDictionary.Add(elements[6].ID, elements[6]);
			model.SubdomainsDictionary[1].Elements.Add(elements[6]);
		}

		[Fact]
        public static void RunNodeToNodeContact()
        {
            CreateModel();

			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);
			var provider = new ProblemStructural(model, solver);
            int increments = 15;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Request output
            //childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 5 });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }
    }
}
