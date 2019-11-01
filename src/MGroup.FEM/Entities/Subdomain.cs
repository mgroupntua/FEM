using System.Collections.Generic;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Commons;
using MGroup.MSolve.Discretization.DofOrdering;

//TODO: remove code that calculates rhs vector components (nodal loads, constraints, etc). It should be moved to dedicated 
//      classes like EquivalentLoadAssembler, so that it can be reused between subdomains of different projects (FEM, IGA, XFEM).
//TODO: same for multiscale
namespace MGroup.FEM.Entities
{
	public class Subdomain : ISubdomain
	{
		private readonly List<Node> nodes = new List<Node>();

		public Subdomain(int id)
		{
			this.ID = id;
		}

		public Table<INode, IDofType, double> Constraints { get; } = new Table<INode, IDofType, double>();

		IReadOnlyList<IElement> ISubdomain.Elements => Elements;
		public List<Element> Elements { get; } = new List<Element>();

		//public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

		public int ID { get; }

		IReadOnlyList<INode> ISubdomain.Nodes => nodes;
		public IReadOnlyList<Node> Nodes => nodes;

		public ISubdomainConstrainedDofOrdering ConstrainedDofOrdering { get; set; }
		public ISubdomainFreeDofOrdering FreeDofOrdering { get; set; }

		public Vector Forces { get; set; } //TODO: this doesn't belong here

		//public bool MaterialsModified
		//{
		//    get
		//    {
		//        bool modified = false;
		//        foreach (Element element in elementsDictionary.Values)
		//            if (element.ElementType.MaterialModified)
		//            {
		//                modified = true;
		//                break;
		//            }
		//        return modified;
		//    }
		//}
		public bool StiffnessModified { get; set; } = true; // At first it is modified
		public bool ConnectivityModified { get; set; } = true; // At first it is modified

		//TODO: This belongs in EquivalentLoadsAssembler
		//TODO: the constraintScalingFactor parameter is not used.
		public double[] CalculateElementIncrementalConstraintDisplacements(IElement element, double constraintScalingFactor)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
		{
			var elementNodalDisplacements = new double[FreeDofOrdering.CountElementDofs(element)];
			SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
			return elementNodalDisplacements;
		}

		public double[] CalculateElementDisplacements(IElement element, IVectorView globalDisplacementVector)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
		{
			double[] elementNodalDisplacements = FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
			SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
			return elementNodalDisplacements;
		}

		public void ClearMaterialStresses()
		{
			foreach (Element element in Elements) element.ElementType.ClearMaterialStresses();
		}

		public void DefineNodesFromElements()
		{
			nodes.Clear();
			var nodeComparer = Comparer<Node>.Create((Node node1, Node node2) => node1.ID - node2.ID);
			var nodeSet = new SortedSet<Node>(nodeComparer);
			foreach (Element element in Elements)
			{
				foreach (Node node in element.Nodes) nodeSet.Add(node);
			}
			nodes.AddRange(nodeSet);

			//foreach (var e in modelEmbeddedNodes.Where(x => nodeIDs.IndexOf(x.Node.ID) >= 0))
			//    EmbeddedNodes.Add(e);
		}

		//TODO: constraints should not be saved inside the nodes. As it is right now (22/11/2018) the same constraint 
		//      is saved in the node, the model constraints table and the subdomain constraints table. Furthermore,
		//      displacement control analyzer updates the subdomain constraints table only (another bad design decision).  
		//      It is too easy to access the wrong instance of the constraint. 
		public void ExtractConstraintsFromGlobal(Table<INode, IDofType, double> globalConstraints)
		{
			//TODO: perhaps it is more efficient to traverse the global constraints instead of the subdomain's nodes, provided
			//      the latter are stored as a set. 
			//TODO: the next could be a Table method: Table.KeepDataOfRows(IEnumerable<TRow> rows)
			foreach (Node node in Nodes)
			{
				bool isNodeConstrained = globalConstraints.TryGetDataOfRow(node,
					out IReadOnlyDictionary<IDofType, double> constraintsOfNode);
				if (isNodeConstrained)
				{
					foreach (var dofDisplacementPair in constraintsOfNode)
					{
						Constraints[node, dofDisplacementPair.Key] = dofDisplacementPair.Value;
					}
				}
			}

			// This is probably faster but assumes that nodes store their prescribed displacements, which I hate.
			//foreach (Node node in Nodes)
			//{
			//    if (node.Constraints == null) continue;
			//    foreach (Constraint constraint in node.Constraints) Constraints[node, constraint.DOF] = constraint.Amount;
			//}
		}

		public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
		{
			var forces = Vector.CreateZero(FreeDofOrdering.NumFreeDofs); //TODO: use Vector
			foreach (Element element in Elements)
			{
				//var localSolution = GetLocalVectorFromGlobal(element, solution);//TODOMaria: This is where the element displacements are calculated //removeMaria
				//var localdSolution = GetLocalVectorFromGlobal(element, dSolution);//removeMaria

				//TODO: ElementType should operate with Vector instead of double[]. Then the ToRawArray() calls can be removed
				double[] localSolution = CalculateElementDisplacements(element, solution);
				double[] localdSolution = CalculateElementDisplacements(element, dSolution);
				element.ElementType.CalculateStresses(element, localSolution, localdSolution);
				if (element.ElementType.MaterialModified)
					element.Subdomain.StiffnessModified = true;
				var f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
				FreeDofOrdering.AddVectorElementToSubdomain(element, f, forces);
			}
			return forces;
		}

		public void ResetMaterialsModifiedProperty()
		{
			this.StiffnessModified = false;
			foreach (Element element in Elements) element.ElementType.ResetMaterialModified();
		}

		public void SaveMaterialState()
		{
			foreach (Element element in Elements) element.ElementType.SaveMaterialState();
		}

		//TODO: I am against modifying the constraints table of the subdomain. Instead the analyzer should keep a constraint
		//      displacements vector at global/subdomain scale and modify that.
		public void ScaleConstraints(double scalingFactor) => Constraints.ModifyValues((u) => scalingFactor * u);

		//ADDED1
		public IVector GetRHSFromSolutionWithInitialDisplacementsEffect(IVectorView solution, IVectorView dSolution, Dictionary<int, INode> boundaryNodes,
			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements)
		{
			var forces = Vector.CreateZero(FreeDofOrdering.NumFreeDofs); //TODO: use Vector
			foreach (Element element in Elements)
			{
				var localSolution = GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, solution);
				ImposePrescribedDisplacementsWithInitialConditionSEffect(element, localSolution, boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
				var localdSolution = GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, dSolution);
				ImposePrescribed_d_DisplacementsWithInitialConditionSEffect(element, localdSolution, boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
				element.ElementType.CalculateStresses(element, localSolution, localdSolution);
				if (element.ElementType.MaterialModified)
					element.Subdomain.StiffnessModified = true;
				var f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
				FreeDofOrdering.AddVectorElementToSubdomain(element, f, forces);
			}
			return forces;
		}

		public double[] GetLocalVectorFromGlobalWithoutPrescribedDisplacements(Element element, IVectorView globalDisplacementVector)
		{
			double[] elementNodalDisplacements = FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
			return elementNodalDisplacements;
		}


		public void ImposePrescribedDisplacementsWithInitialConditionSEffect(Element element, double[] localSolution, Dictionary<int, INode> boundaryNodes,
			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements)
		{

			var elementDOFTypes = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
			var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
			int iElementMatrixColumn = 0;
			for (int j = 0; j < elementDOFTypes.Count; j++)
			{
				INode nodeColumn = matrixAssemblyNodes[j];
				int nodalDofsNumber = elementDOFTypes[j].Count;
				if (boundaryNodes.ContainsKey(nodeColumn.ID))
				{
					Dictionary<IDofType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
					Dictionary<IDofType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
					int positionOfDofInNode = 0;
					foreach (IDofType doftype1 in elementDOFTypes[j])
					{
						if (nodalConvergedDisplacements.ContainsKey(doftype1))
						{
							localSolution[iElementMatrixColumn + positionOfDofInNode] = nodalConvergedDisplacements[doftype1] + (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
							// TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
						}
						positionOfDofInNode += 1;
					}
				}
				iElementMatrixColumn += nodalDofsNumber;
			}

		}

		public void ImposePrescribed_d_DisplacementsWithInitialConditionSEffect(Element element, double[] localSolution, Dictionary<int, INode> boundaryNodes,
			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements)
		{

			var elementDOFTypes = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
			var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
			int iElementMatrixColumn = 0;
			for (int j = 0; j < elementDOFTypes.Count; j++)
			{
				INode nodeColumn = matrixAssemblyNodes[j];
				int nodalDofsNumber = elementDOFTypes[j].Count;
				if (boundaryNodes.ContainsKey(nodeColumn.ID))
				{
					Dictionary<IDofType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
					Dictionary<IDofType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
					int positionOfDofInNode = 0;
					foreach (IDofType doftype1 in elementDOFTypes[j])
					{
						if (nodalConvergedDisplacements.ContainsKey(doftype1))
						{
							localSolution[iElementMatrixColumn + positionOfDofInNode] = (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
							// 1) den vazoume mono (1/increments) alla (nIncrement/increments) dioti metaxu aftwn twn nIncrements den exei mesolavhsei save sta material ths mikroklimakas
							// TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
						}
						positionOfDofInNode += 1;
					}
				}
				iElementMatrixColumn += nodalDofsNumber;
			}

		}
	}
}
