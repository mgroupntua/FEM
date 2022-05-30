using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.FEM
{
	public class AutomaticDomainDecomposer
	{
		private readonly Model Model;
		private readonly int NumberOfProcessors;
		private int numberOfElementsPerSubdomain;
		Dictionary<IElementType, List<IElementType>> ElementAdjacency;
		List<IElementType> ElementsRenumbered = new List<IElementType>();
		Dictionary<int, List<INode>> SubdomainInterfaceNodes = new Dictionary<int, List<INode>>();

		public AutomaticDomainDecomposer(Model model, int numberOfProcessors)
		{
			this.NumberOfProcessors = numberOfProcessors;
			this.Model = model;
		}


		public void UpdateModel(bool isColoringEnabled = false)
		{
			Adjacency();

			CreateSubdomains();
			AssignElementsRenumberedToSubdomains();

			if (isColoringEnabled)
			{
				var purgedElements = Purge();
				ColorDisconnectedElements(purgedElements);
			}
		}

		private void AssignElementsRenumberedToSubdomains()
		{
			Model.SubdomainsDictionary.Clear();
			var indexElement = 0;
			for (int i = 0; i < NumberOfProcessors; i++)
			{
				if (indexElement >= ElementsRenumbered.Count) break;
				Model.SubdomainsDictionary.Add(i, new Subdomain(i));
				for (int j = 0; j < numberOfElementsPerSubdomain; j++)
				{
					if (indexElement >= ElementsRenumbered.Count) break;
					Model.SubdomainsDictionary[i].Elements.Add(ElementsRenumbered[indexElement++]);
				}
			}

			UpdateModelDataStructures();
		}

		private void UpdateModelDataStructures()
		{
			foreach (Subdomain subdomain in Model.SubdomainsDictionary.Values)
			{
				foreach (IElementType element in subdomain.Elements)
					element.SubdomainID = subdomain.ID;
			}

			foreach (Node node in Model.NodesDictionary.Values)
			{
				node.Subdomains.Clear();
				node.FindAssociatedSubdomains();
			}

			foreach (Subdomain subdomain in Model.SubdomainsDictionary.Values)
			{
				//subdomain.Nodes.Clear(); // This is done automatically by the next method
				subdomain.DefineNodesFromElements();
			}
		}

		private void Adjacency()
		{
			ElementAdjacency = new Dictionary<IElementType, List<IElementType>>();
			// mask is an integer that shows if the element is used


			foreach (var element in Model.ElementsDictionary.Values)
			{
				bool[] usedElement = new bool[Model.ElementsDictionary.Count];//mask
				ElementAdjacency.Add(element, new List<IElementType>());
				usedElement[element.ID] = true;
				foreach (var node in element.Nodes)
				{
					foreach (IElementType nodeElement in node.ElementsDictionary.Values)
					{
						if (usedElement[nodeElement.ID]) continue;
						ElementAdjacency[element].Add(nodeElement);
						usedElement[nodeElement.ID] = true;
					}
				}
			}
		}

		private void CreateSubdomains()
		{
			//TODO: return IndexOutOfRangeException if nodes,elements or subdomains numbering does not start with 0.
			bool[] isInteriorBoundaryElement = new bool[Model.ElementsDictionary.Count];
			bool[] isInteriorBoundaryNode = new bool[Model.NodesDictionary.Count];

			// Number of Elements per subdomain
			numberOfElementsPerSubdomain = (Model.ElementsDictionary.Count % NumberOfProcessors == 0) ?
				Model.ElementsDictionary.Count / NumberOfProcessors : Model.ElementsDictionary.Count / NumberOfProcessors + 1;

			Dictionary<INode, int> nodeWeight = new Dictionary<INode, int>();
			foreach (Node node in Model.NodesDictionary.Values)
				nodeWeight.Add(node, node.ElementsDictionary.Count);

			var usedElementsCounter = 0;
			var mlabel = 0;
			int counterSubdomain = 0;

			do
			{
				var flag = true;
				var flagStop = true;
				#region Find Node with next minimum weight
				var finalSubdomainElement = usedElementsCounter;
				var minimumNodeWeight = int.MaxValue;
				int nodeID = 0; ;
				for (int i = 0; i < Model.NodesDictionary.Count; i++)
				{
					if (nodeWeight[Model.NodesDictionary[i]] == 0) continue;
					if (nodeWeight[Model.NodesDictionary[i]] < minimumNodeWeight)
					{
						minimumNodeWeight = nodeWeight[Model.NodesDictionary[i]];
						nodeID = i;
					}
				}
				#endregion

				// Start fill list with elements connected to node with minimum weight
				var counterSubdomainElements = 0;
				foreach (IElementType element in Model.NodesDictionary[nodeID].ElementsDictionary.Values)
				{
					var elementID = element.ID;
					if (isInteriorBoundaryElement[elementID]) continue;
					counterSubdomainElements++;
					isInteriorBoundaryElement[elementID] = true;
					ElementsRenumbered.Add(element);

					#region nomask
					//Reduce nodeWeight for all nodes connected to this element
					foreach (Node node in Model.ElementsDictionary[elementID].Nodes)
						nodeWeight[node]--;
					#endregion

					if (counterSubdomainElements == numberOfElementsPerSubdomain)
					{
						flag = false;
						break;
					}
				}

				if (flag)
				{
					// Recursively add adjacent elements to list
					do
					{
						var initialSubdomainElement = finalSubdomainElement;
						finalSubdomainElement = usedElementsCounter + counterSubdomainElements;
						var nnstart = initialSubdomainElement + 1;

						for (int i = initialSubdomainElement; i <= finalSubdomainElement - 1; i++)
						{
							int lc = 0;
							for (int j = 0; j < ElementAdjacency.First(x => x.Key.ID == ElementsRenumbered[i].ID).Value.Count; j++)
							{
								int elementID = ElementAdjacency[ElementsRenumbered[i]][j].ID;
								if (isInteriorBoundaryElement[elementID]) continue;
								lc++;
								counterSubdomainElements++;
								isInteriorBoundaryElement[elementID] = true;
								ElementsRenumbered.Add(Model.ElementsDictionary.First(x => x.Value.ID == elementID).Value);

								#region nomask
								foreach (Node node in Model.ElementsDictionary[elementID].Nodes)
									nodeWeight[node]--;
								#endregion

								if (counterSubdomainElements == numberOfElementsPerSubdomain)
								{
									flag = false;
									break;
								}
							} // 800

							if (flag)
							{
								if (lc == 0 && (usedElementsCounter + counterSubdomainElements) == finalSubdomainElement && i == finalSubdomainElement - 1)
								{
									usedElementsCounter = usedElementsCounter + counterSubdomainElements;
									flagStop = false;
									flag = false;
								}
							}

							if (!flag) break;
						}
						if (!flag) break;
					} while (counterSubdomainElements < numberOfElementsPerSubdomain);
					if (!flagStop) continue;
				}
				SubdomainInterfaceNodes.Add(counterSubdomain++,
					CalculateInterface(nodeWeight, isInteriorBoundaryNode, ElementsRenumbered, usedElementsCounter,
					counterSubdomainElements, counterSubdomain));
				usedElementsCounter = usedElementsCounter + counterSubdomainElements;
				mlabel = usedElementsCounter;
			} while (usedElementsCounter < Model.ElementsDictionary.Count);

		}

		//isInteriorBoundaryNode-> if node is on the interior interface
		private List<INode> CalculateInterface(Dictionary<INode, int> nodeWeight, bool[] isInteriorBoundaryNode,
			List<IElementType> ElementsRenumbered, int usedElementsCounter, int counterSubdomainElements, int counterSubdomain)
		{
			var locmask = new bool[Model.NodesDictionary.Count];
			List<INode> SubdomainInterfaceNodes = new List<INode>();

			for (int i = usedElementsCounter; i < usedElementsCounter + counterSubdomainElements; i++)
			{
				int elementID = ElementsRenumbered[i].ID;
				foreach (Node node in Model.ElementsDictionary[elementID].Nodes)
				{
					if ((nodeWeight[node] != 0 || isInteriorBoundaryNode[node.ID]) && !locmask[node.ID])
					{
						isInteriorBoundaryNode[node.ID] = true;
						locmask[node.ID] = true;
						SubdomainInterfaceNodes.Add(node);
					}
				}
			}
			return SubdomainInterfaceNodes;
		}

		private void ColorDisconnectedElements(List<IElementType> purgedElements)
		{
			int numberOfColors = 1;
			Dictionary<int, List<IElementType>> elementsPerColor = new Dictionary<int, List<IElementType>>();
			int indexColor = 0;
			int counterElement = 0;

			bool[] isElementUsed = new bool[Model.ElementsDictionary.Count];

			// Form the list of distinct colors
			do
			{
				elementsPerColor.Add(indexColor, new List<IElementType>());
				bool[] usedNodes = new bool[Model.NodesDictionary.Count];
				foreach (var element in purgedElements)
				{
					if (isElementUsed[element.ID]) continue;

					bool disjoint = CheckIfElementDisjoint(usedNodes, element.ID);

					if (disjoint)
					{
						counterElement++;
						elementsPerColor[indexColor].Add(Model.ElementsDictionary[element.ID]);

						#region marker
						foreach (Node node in Model.ElementsDictionary[element.ID].Nodes)
							usedNodes[node.ID] = true;
						#endregion

					}

				}
				indexColor++;
			} while (counterElement < purgedElements.Count);
			numberOfColors = indexColor;
		}

		private bool CheckIfElementDisjoint(bool[] usedNodes, int elementID)
		{
			foreach (Node node in Model.ElementsDictionary[elementID].Nodes)
				if (usedNodes[node.ID])
					return false;
			return true;
		}


		private List<IElementType> Purge()
		{
			bool[] isElementUsed = new bool[Model.ElementsDictionary.Count];

			List<INode> interfaceNodes = new List<INode>();
			foreach (var nodeList in SubdomainInterfaceNodes.Values)
				interfaceNodes.AddRange(nodeList);

			interfaceNodes = interfaceNodes.Distinct().ToList();
			int numberOfInterfaceNodes = interfaceNodes.Count;

			List<IElementType> purgedElements = new List<IElementType>();
			int numberOfPurgedElements = 0;
			for (int indexInterfaceNode = 0; indexInterfaceNode < numberOfInterfaceNodes; indexInterfaceNode++)
			{
				foreach (IElementType element in interfaceNodes[indexInterfaceNode].ElementsDictionary.Values)
				{
					if (!isElementUsed[element.ID])
					{
						isElementUsed[element.ID] = true;
						purgedElements.Add(element);
					}
				}
			}

			return purgedElements;
		}

	}
}
