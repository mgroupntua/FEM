using System.Collections.Generic;
using System.Linq;
using MGroup.FEM.Interfaces;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Commons;
using MGroup.MSolve.Discretization.DofOrdering;
using MGroup.MSolve.Discretization.Loads;

//TODO: find what is going on with the dynamic loads and refactor them. That 564000000 in AssignMassAccelerationHistoryLoads()
//      cannot be correct.
//TODO: ConnectDataStructures() should not be called twice. There should be a flag that determines if it has been called. If it
//      has, the method should just return without doing anything.
//TODO: Replace all IList with IReadOnlyList. Even better, have a different class to create the model than the one used to 
//      store the entities, so that they can be accessed by solvers, analyzers & loggers. Could the latter be the same for FEM, 
//      IGA, XFEM?
namespace MGroup.FEM.Entities
{
	public class Model : IModel
	{
		//public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

		public IList<Cluster> Clusters => ClustersDictionary.Values.ToList();
		public Dictionary<int, Cluster> ClustersDictionary { get; } = new Dictionary<int, Cluster>();

		public Table<INode, IDofType, double> Constraints { get; private set; } = new Table<INode, IDofType, double>();//TODOMaria: maybe it's useless in model class

		IReadOnlyList<IElement> IModel.Elements => ElementsDictionary.Values.ToList();
		public IList<Element> Elements => ElementsDictionary.Values.ToList();
		public Dictionary<int, Element> ElementsDictionary { get; } = new Dictionary<int, Element>();

		public IList<ElementMassAccelerationHistoryLoad> ElementMassAccelerationHistoryLoads { get; }
			= new List<ElementMassAccelerationHistoryLoad>();
		public IList<ElementMassAccelerationLoad> ElementMassAccelerationLoads { get; }
			= new List<ElementMassAccelerationLoad>();
		public IList<Load> Loads { get; private set; } = new List<Load>();
		public IList<MassAccelerationLoad> MassAccelerationLoads { get; } = new List<MassAccelerationLoad>();
		public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; } = new List<IMassAccelerationHistoryLoad>();

		public IList<Node> Nodes => NodesDictionary.Values.ToList();
		IReadOnlyList<INode> IModel.Nodes => NodesDictionary.Values.ToList();
		public Dictionary<int, Node> NodesDictionary { get; } = new Dictionary<int, Node>();

		IReadOnlyList<ISubdomain> IModel.Subdomains => SubdomainsDictionary.Values.ToList();
		public IReadOnlyList<Subdomain> Subdomains => SubdomainsDictionary.Values.ToList();
		public Dictionary<int, Subdomain> SubdomainsDictionary { get; } = new Dictionary<int, Subdomain>();

		public IList<ITimeDependentNodalLoad> TimeDependentNodalLoads { get; private set; } = new List<ITimeDependentNodalLoad>();

		public IGlobalFreeDofOrdering GlobalDofOrdering { get; set; }

		public void AssignLoads(NodalLoadsToSubdomainsDistributor distributeNodalLoads)
		{
			foreach (Subdomain subdomain in SubdomainsDictionary.Values) subdomain.Forces.Clear();
			AssignNodalLoads(distributeNodalLoads);
			AssignElementMassLoads();
			AssignMassAccelerationLoads();
		}

		public void AssignMassAccelerationHistoryLoads(int timeStep)
		{
			if (MassAccelerationHistoryLoads.Count > 0)
			{
				List<MassAccelerationLoad> m = new List<MassAccelerationLoad>(MassAccelerationHistoryLoads.Count);
				foreach (IMassAccelerationHistoryLoad l in MassAccelerationHistoryLoads)
				{
					m.Add(new MassAccelerationLoad() { Amount = l[timeStep], DOF = l.DOF });
				}

				foreach (Subdomain subdomain in SubdomainsDictionary.Values)
				{
					foreach (Element element in subdomain.Elements)
					{
						double[] accelerationForces = element.ElementType.CalculateAccelerationForces(element, m);
						subdomain.FreeDofOrdering.AddVectorElementToSubdomain(element, accelerationForces, subdomain.Forces);
					}
				}
			}

			foreach (ElementMassAccelerationHistoryLoad load in ElementMassAccelerationHistoryLoads)
			{
				MassAccelerationLoad hl = new MassAccelerationLoad()
				{
					Amount = load.HistoryLoad[timeStep] * 564000000,
					DOF = load.HistoryLoad.DOF
				};
				Element element = load.Element;
				ISubdomain subdomain = element.Subdomain;
				var accelerationForces = element.ElementType.CalculateAccelerationForces(
					load.Element, (new MassAccelerationLoad[] { hl }).ToList());
				GlobalDofOrdering.SubdomainDofOrderings[subdomain].AddVectorElementToSubdomain(element, accelerationForces,
					subdomain.Forces);
			}
		}

		public void AssignNodalLoads(NodalLoadsToSubdomainsDistributor distributeNodalLoads)
		{
			var globalNodalLoads = new Table<INode, IDofType, double>();
			foreach (Load load in Loads) globalNodalLoads.TryAdd(load.Node, load.DOF, load.Amount);

			Dictionary<int, SparseVector> subdomainNodalLoads = distributeNodalLoads(globalNodalLoads);
			foreach (var idSubdomainLoads in subdomainNodalLoads)
			{
				SubdomainsDictionary[idSubdomainLoads.Key].Forces.AddIntoThis(idSubdomainLoads.Value);
			}
		}

		public void AssignTimeDependentNodalLoads(int timeStep, NodalLoadsToSubdomainsDistributor distributeNodalLoads)
		{
			var globalNodalLoads = new Table<INode, IDofType, double>();
			foreach (ITimeDependentNodalLoad load in TimeDependentNodalLoads)
			{
				globalNodalLoads.TryAdd(load.Node, load.DOF, load.GetLoadAmount(timeStep));
			}

			Dictionary<int, SparseVector> subdomainNodalLoads = distributeNodalLoads(globalNodalLoads);
			foreach (var idSubdomainLoads in subdomainNodalLoads)
			{
				SubdomainsDictionary[idSubdomainLoads.Key].Forces.AddIntoThis(idSubdomainLoads.Value);
			}
		}

		//What is the purpose of this method? If someone wanted to clear the Model, they could just create a new one.
		public void Clear()
		{
			Loads.Clear();
			ClustersDictionary.Clear();
			SubdomainsDictionary.Clear();
			ElementsDictionary.Clear();
			NodesDictionary.Clear();
			GlobalDofOrdering = null;
			Constraints.Clear();
			ElementMassAccelerationHistoryLoads.Clear();
			ElementMassAccelerationLoads.Clear();
			MassAccelerationHistoryLoads.Clear();
			MassAccelerationLoads.Clear();
		}

		// Warning: This is called by the analyzer, so that the user does not have to call it explicitly. However, it is must be 
		// called explicitly before the AutomaticDomainDecompositioner is used.
		public void ConnectDataStructures()
		{
			BuildInterconnectionData();
			AssignConstraints();
			RemoveInactiveNodalLoads();

			//TODOSerafeim: This should be called by the analyzer, which defines when the dofs are ordered and when the global vectors/matrices are built.
			//AssignLoads();
		}

		//TODO: constraints should not be saved inside the nodes. As it is right now (22/11/2018) the same constraint 
		//      is saved in the node, the model constraints table and the subdomain constraints table. Furthermore,
		//      displacement control analyzer updates the subdomain constraints table only (another bad design decision).  
		//      It is too easy to access the wrong instance of the constraint. 
		private void AssignConstraints()
		{
			foreach (Node node in NodesDictionary.Values)
			{
				if (node.Constraints == null) continue;
				foreach (Constraint constraint in node.Constraints) Constraints[node, constraint.DOF] = constraint.Amount;
			}

			foreach (Subdomain subdomain in SubdomainsDictionary.Values) subdomain.ExtractConstraintsFromGlobal(Constraints);
		}

		private void AssignElementMassLoads()
		{
			foreach (ElementMassAccelerationLoad load in ElementMassAccelerationLoads)
			{
				ISubdomain subdomain = load.Element.Subdomain;
				var accelerationForces = load.Element.ElementType.CalculateAccelerationForces(
					load.Element, MassAccelerationLoads);
				GlobalDofOrdering.SubdomainDofOrderings[subdomain].AddVectorElementToSubdomain(load.Element,
					accelerationForces, subdomain.Forces);
			}
		}

		private void AssignMassAccelerationLoads()
		{
			if (MassAccelerationLoads.Count < 1) return;

			foreach (Subdomain subdomain in SubdomainsDictionary.Values)
			{
				foreach (Element element in subdomain.Elements)
				{
					subdomain.FreeDofOrdering.AddVectorElementToSubdomain(element,
						element.ElementType.CalculateAccelerationForces(element, MassAccelerationLoads),
						subdomain.Forces);
				}
			}
		}

		private void BuildElementDictionaryOfEachNode()
		{
			foreach (Element element in ElementsDictionary.Values)
			{
				foreach (Node node in element.Nodes) node.ElementsDictionary[element.ID] = element;
			}
		}

		private void BuildInterconnectionData()//TODOMaria: maybe I have to generate the constraints dictionary for each subdomain here
		{
			BuildSubdomainOfEachElement();
			DuplicateInterSubdomainEmbeddedElements();
			BuildElementDictionaryOfEachNode();
			foreach (Node node in NodesDictionary.Values) node.BuildSubdomainDictionary();

			//BuildNonConformingNodes();

			foreach (Subdomain subdomain in SubdomainsDictionary.Values) subdomain.DefineNodesFromElements();
		}

		private void BuildSubdomainOfEachElement()
		{
			foreach (Subdomain subdomain in SubdomainsDictionary.Values)
			{
				foreach (Element element in subdomain.Elements) element.Subdomain = subdomain;
			}
		}

		private void BuildNonConformingNodes()
		{
			List<int> subIDs = new List<int>();
			foreach (Element element in ElementsDictionary.Values)
			{
				subIDs.Clear();

				foreach (Node node in element.Nodes)
				{
					foreach (int subID in node.SubdomainsDictionary.Keys)
					{
						if (!subIDs.Contains(subID)) subIDs.Add(subID);

					}
				}

				foreach (Node node in element.Nodes)
				{
					foreach (int subID in subIDs)
					{
						if (!node.SubdomainsDictionary.ContainsKey(subID))
						{
							node.NonMatchingSubdomainsDictionary.Add(subID, SubdomainsDictionary[subID]);
						}
					}
				}

			}
		}

		private void DuplicateInterSubdomainEmbeddedElements()
		{
			foreach (var e in ElementsDictionary.Values.Where(x => x.ElementType is IEmbeddedElement))
			{
				var subs = ((IEmbeddedElement)e.ElementType).EmbeddedNodes.Select(x => x.EmbeddedInElement.Subdomain).Distinct();
				foreach (var s in subs.Where(x => x.ID != e.Subdomain.ID))
					s.Elements.Add(e);
			}
		}

		private void RemoveInactiveNodalLoads()
		{
			// Static loads
			var activeLoadsStatic = new List<Load>(Loads.Count);
			foreach (Load load in Loads)
			{
				bool isConstrained = Constraints.Contains(load.Node, load.DOF);
				if (!isConstrained) activeLoadsStatic.Add(load);
			}
			Loads = activeLoadsStatic;

			// Dynamic loads
			var activeLoadsDynamic = new List<ITimeDependentNodalLoad>(TimeDependentNodalLoads.Count);
			foreach (ITimeDependentNodalLoad load in TimeDependentNodalLoads)
			{
				bool isConstrained = Constraints.Contains(load.Node, load.DOF);
				if (!isConstrained) activeLoadsDynamic.Add(load);
			}
			TimeDependentNodalLoads = activeLoadsDynamic;
		}
	}
}
