using System.Collections.Generic;

namespace MGroup.FEM.Entities
{
	public class Cluster
	{
		private readonly IList<Subdomain> subdomains = new List<Subdomain>();

		public IList<Subdomain> Subdomains
		{
			get { return subdomains; }
		}

		public int ID { get; set; }
	}
}
