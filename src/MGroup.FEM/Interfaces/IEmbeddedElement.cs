using System.Collections.Generic;
using MGroup.FEM.Embedding;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization;

namespace MGroup.FEM.Interfaces
{
	public interface IEmbeddedElement
	{
		IList<EmbeddedNode> EmbeddedNodes { get; }
		Dictionary<IDofType, int> GetInternalNodalDOFs(Element element, Node node);
		double[] GetLocalDOFValues(Element hostElement, double[] hostDOFValues);
	}
}
