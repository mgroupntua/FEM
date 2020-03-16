using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Loads;
using Moq;
using Xunit;

namespace MGroup.FEM.Tests.IntegrationTests
{
	using Constitutive.Structural;
	using MSolve.Solution;
	using NumericalAnalyzers;
	using NumericalAnalyzers.Dynamic;
	using NumericalAnalyzers.Logging;
	using Solvers.Direct;

	public static class NewmarkDynamicAnalysisTests
	{
		[Fact]
		private static void TestBatheImplicitAnalysisExample()
		{
			var model = new Model();
			int subdomainID = 0;
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

			var n = new Node(id: 0, x: double.NaN);
			var e = new Element() { ID = 0 };
			e.NodesDictionary.Add(0, n);
			var m = new Mock<IFiniteElement>();
			m.Setup(x => x.StiffnessMatrix(e)).Returns(Matrix.CreateFromArray(new double[,] { { 6, -2, }, { -2, 4 } }));
			m.Setup(x => x.MassMatrix(e)).Returns(Matrix.CreateFromArray(new double[,] { { 2, 0, }, { 0, 1 } }));
			m.Setup(x => x.DampingMatrix(e)).Returns(Matrix.CreateFromArray(new double[,] { { 0, 0, }, { 0, 0 } }));
			//m.Setup(x => x.StiffnessMatrix(e)).Returns(new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 6, -2, 4 }));
			//m.Setup(x => x.MassMatrix(e)).Returns(new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 2, 0, 1 }));
			//m.Setup(x => x.DampingMatrix(e)).Returns(new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 0, 0, 0 }));
			m.Setup(x => x.GetElementDofTypes(e)).Returns(new[] { new[] { StructuralDof.TranslationX, StructuralDof.TranslationY } });
			m.SetupGet(x => x.DofEnumerator).Returns(new GenericDofEnumerator());
			e.ElementType = m.Object;
			model.NodesDictionary.Add(0, n);
			model.ElementsDictionary.Add(0, e);
			model.SubdomainsDictionary[subdomainID].Elements.Add(e);
			model.Loads.Add(new Load() { Amount = 10, Node = n, DOF = StructuralDof.TranslationY });
			var lX = new Mock<IMassAccelerationHistoryLoad>();
			lX.SetupGet(x => x.DOF).Returns(StructuralDof.TranslationX);
			lX.SetupGet(x => x[It.IsAny<int>()]).Returns(0);
			var lY = new Mock<IMassAccelerationHistoryLoad>();
			lY.SetupGet(x => x.DOF).Returns(StructuralDof.TranslationY);
			lY.SetupGet(x => x[0]).Returns(10);
			lY.SetupGet(x => x[It.IsInRange(1, 100, Range.Inclusive)]).Returns(0);
			model.MassAccelerationHistoryLoads.Add(lX.Object);
			model.MassAccelerationHistoryLoads.Add(lY.Object);
			m.Setup(x => x.CalculateAccelerationForces(It.IsAny<Element>(), It.IsAny<IList<MassAccelerationLoad>>()))
				.Returns<Element, IList<MassAccelerationLoad>>((element, loads) =>
				{
					double[] accelerations = { loads[0].Amount, loads[1].Amount };
					var massMatrix = Matrix.CreateFromArray(new double[,] { { 2, 0, }, { 0, 1 } });
					//var massMatrix = new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 2, 0, 1 });
					return massMatrix.Multiply(accelerations);
				}
			);

			// Solver
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);

			//TODO: These overwrite the corresponding data extracted by the Model. Either set up these or the Model.
			//solver.LinearSystems[subdomainID].SetMatrix(
			//    SkylineMatrix.CreateFromArrays(2, new double[] { 6, 4, -2 }, new int[] { 0, 1, 3 }, true)); // K = [6 -2; -2 4]
			//solver.LinearSystems[subdomainID].RhsVector = Vector.CreateFromArray(new double[] { 0, 10 });

			// Problem type
			var provider = new ProblemStructural(model, solver);

			// Analyzers
			var childAnalyzer = new LinearAnalyzer(model, solver, provider);
			var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 0.28, 3.36);
			parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
																				 //parentAnalyzerBuilder.SetNewmarkParameters(0.25, 0.5); // Not necessary. This is the default
			NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

			// Request output
			childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { 0, 1 });

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			//Check output
			DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
			Assert.Equal(2.2840249264795207, log.DOFValues[0], 8);
			Assert.Equal(2.4351921891904156, log.DOFValues[1], 8);
		}
	}
}

