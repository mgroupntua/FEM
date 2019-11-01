using System;
using System.Collections.Generic;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Geometry.Coordinates;

//TODO: Elements and Subdomains should not be stored inside the node by default. E.g. for elements there should be 
//      a lightweiht unidirectional mesh class, where only elements store nodes, and a birectional mesh class, where a 
//      Dictionary of nodes -> elements is stored. Thus extra memory for these associations is not always required.
namespace MGroup.FEM.Entities
{
	/// <summary>
	/// Vertex of a finite element in a 3-dimensional space. It can also represent points in 1-dimensional or 2-dimension 
	/// spaces. Immutable.
	/// </summary>
	public class Node : CartesianPoint, INode
	{
		private readonly List<Constraint> constraints = new List<Constraint>();
		private readonly Dictionary<int, IElement> elementsDictionary = new Dictionary<int, IElement>();
		private readonly Dictionary<int, Subdomain> nonMatchingSubdomainsDictionary = new Dictionary<int, Subdomain>();
		private readonly Dictionary<int, ISubdomain> subdomainsDictionary = new Dictionary<int, ISubdomain>();

		/// <summary>
		/// Instantiates a <see cref="Node"/>.
		/// </summary>
		/// <param name="id">
		/// A unique identifier <see cref="ID"/> to differentiate this instance of <see cref="Node1D"/> from the rest. 
		/// Constraints: <see cref="ID"/> &gt;= 0.
		/// </param>
		/// <param name="x">The coordinate of the point along the single axis X.</param>
		/// <param name="y">The coordinate of the point along the single axis Y.</param>
		/// <param name="z">The coordinate of the point along the single axis Z.</param>
		public Node(int id, double x, double y = 0.0, double z = 0.0) : base(x, y, z)
		{
			if (id < 0) throw new ArgumentException("The parameter id must be non negative, but was: " + id);
			this.ID = id;
		}

		//public Element EmbeddedInElement { get; set; }

		/// <summary>
		/// A unique identifier <see cref="ID"/> to differentiate this instance of <see cref="Node"/> from the rest. 
		/// Constraints: <see cref="ID"/> &gt;= 0.
		/// </summary>
		public int ID { get; }

		public override string ToString()
		{
			var header = String.Format("{0}: ({1}, {2}, {3})", ID, X, Y, Z);
			string constraintsDescripton = string.Empty;
			constraintsDescripton = constraintsDescripton.Length > 1
				? constraintsDescripton.Substring(0, constraintsDescripton.Length - 2)
				: constraintsDescripton;

			return String.Format("{0} - Con ({1})", header, constraintsDescripton);
		}

		public List<Constraint> Constraints => constraints;
		public Dictionary<int, IElement> ElementsDictionary => elementsDictionary; //TODO: This should be IElement
		public Dictionary<int, Subdomain> NonMatchingSubdomainsDictionary => nonMatchingSubdomainsDictionary;
		public Dictionary<int, ISubdomain> SubdomainsDictionary => subdomainsDictionary;


		public void BuildSubdomainDictionary()
		{
			foreach (Element element in elementsDictionary.Values)
				if (!subdomainsDictionary.ContainsKey(element.Subdomain.ID))
					subdomainsDictionary.Add(element.Subdomain.ID, element.Subdomain);
		}

		public int CompareTo(INode other) => this.ID - other.ID;
	}
}
