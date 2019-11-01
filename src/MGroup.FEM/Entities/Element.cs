using System.Collections.Generic;
using System.Linq;
using MGroup.FEM.Interfaces;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.DofOrdering;

//TODO: Delete this ASAP. 1) Its purpose is element-node connectivity, which should be done through interfaces and inheritance,
//      2) The order of the nodes should be defined by what is now called ElementType
namespace MGroup.FEM.Entities
{
	public enum AbsorptionType
	{
		Unknown = 0,
		Compressional = 1,
		Shear = 2
	}

	public class Element : IElement
	{
		private readonly Dictionary<int, Node> nodesDictionary = new Dictionary<int, Node>();
		private readonly Dictionary<IDofType, AbsorptionType> absorptions = new Dictionary<IDofType, AbsorptionType>();
		private readonly IList<Node> embeddedNodes = new List<Node>();

		public Dictionary<int, Node> NodesDictionary => nodesDictionary;

		public Dictionary<IDofType, AbsorptionType> Absorptions => absorptions;

		IReadOnlyList<INode> IElement.Nodes => nodesDictionary.Values.ToList<INode>();
		public IList<Node> Nodes => nodesDictionary.Values.ToList();


		public IList<Node> EmbeddedNodes => embeddedNodes;

		public int ID { get; set; }

		IElementType IElement.ElementType => ElementType;
		public IFiniteElement ElementType { get; set; }

		ISubdomain IElement.Subdomain => this.Subdomain;
		public Subdomain Subdomain { get; set; }

		public void AddNode(Node node) => nodesDictionary.Add(node.ID, node);

		public void AddNodes(IList<Node> nodes)
		{
			foreach (Node node in nodes) AddNode(node);
		}
	}
}
