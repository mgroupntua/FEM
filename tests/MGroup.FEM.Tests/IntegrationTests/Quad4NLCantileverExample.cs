namespace MGroup.FEM.Tests.IntegrationTests
{
	public class Quad4NLCantileverExample
	{
		////[Fact]
		//public void TestQuad4NLCantileverExample()
		//{
		//    VectorExtensions.AssignTotalAffinityCount();
		//    double youngModulus = 3.76;
		//    double poissonRatio = 0.3779;
		//    double thickness = 1.0;
		//    double nodalLoad = 500.0;

		//    IContinuumMaterial2D material = new IContinuumMaterial2D(StressState2D.PlaneStress)
		//    {
		//        YoungModulus = youngModulus,
		//        PoissonRatio = poissonRatio
		//    };

		//    // Node creation
		//    IList<Node> nodes = new List<Node>();
		//    Node node1 = new Node( id: 1, x: 0.0, y:  0.0, z: 0.0 };
		//    Node node2 = new Node( id: 2, x: 10.0, y:  0.0, z: 0.0 };
		//    Node node3 = new Node( id: 3, x: 10.0, y:  10.0, z: 0.0 };
		//    Node node4 = new Node( id: 4, x: 0.0, y:  10.0, z: 0.0 };
		//    nodes.Add(node1);
		//    nodes.Add(node2);
		//    nodes.Add(node3);
		//    nodes.Add(node4);

		//    // Model creation
		//    Model model = new Model();

		//    // Add a single subdomain to the model
		//    model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

		//    // Add nodes to the nodes dictonary of the model
		//    for (int i = 0; i < nodes.Count; ++i)
		//    {
		//        model.NodesDictionary.Add(i + 1, nodes[i]);
		//    }

		//    // Constrain bottom nodes of the model
		//    model.NodesDictionary[1].Constraints.Add(DOFType.X);
		//    model.NodesDictionary[1].Constraints.Add(DOFType.Y);
		//    model.NodesDictionary[4].Constraints.Add(DOFType.X);
		//    model.NodesDictionary[4].Constraints.Add(DOFType.Y);

		//    // Create Quad4 element
		//    var quad = new Quad4NL(material) { Thickness = thickness };
		//    var element = new Element()
		//    {
		//        ID = 1,
		//        ElementType = quad
		//    };

		//    // Add nodes to the created element
		//    element.AddNode(model.NodesDictionary[1]);
		//    element.AddNode(model.NodesDictionary[2]);
		//    element.AddNode(model.NodesDictionary[3]);
		//    element.AddNode(model.NodesDictionary[4]);

		//    // Element Stiffness Matrix
		//    var a = quad.StiffnessMatrix(element);

		//    // Add quad element to the element and subdomains dictionary of the model
		//    model.ElementsDictionary.Add(element.ID, element);
		//    model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);

		//    // define loads
		//    model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[2], DOF = DOFType.X });
		//    model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[3], DOF = DOFType.X });
		//    model.ConnectDataStructures();
		//    var linearSystems = new Dictionary<int, ILinearSystem>();
		//    linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
		//    SolverSkyline solver = new SolverSkyline(linearSystems[1]);
		//    ProblemStructural provider = new ProblemStructural(model, linearSystems);
		//    // Choose child analyzer -> Child: Linear or NewtonRaphsonNonLinearAnalyzer
		//    //LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
		//    var linearSystemsArray = new[] { linearSystems[1] };
		//    var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
		//    var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
		//    int increments = 10;
		//    int totalDOFs = model.TotalDOFs;
		//    NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
		//    provider, increments, totalDOFs);
		//    // Choose parent analyzer -> Parent: Static or Dynamic
		//    StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
		//    //NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);
		//    parentAnalyzer.BuildMatrices();
		//    parentAnalyzer.Initialize();
		//    parentAnalyzer.Solve();
		//    Assert.Equal(8.8937311427, linearSystems[1].Solution[1], 8);
		//}
	}
}
