using MGroup.FEM.Embedding;
using MGroup.FEM.Entities;

namespace MGroup.FEM.Interfaces
{
	public interface IEmbeddedHostElement
	{
		EmbeddedNode BuildHostElementEmbeddedNode(Element element, Node node, IEmbeddedDOFInHostTransformationVector transformationVector);
		double[] GetShapeFunctionsForNode(Element element, EmbeddedNode node);
	}
}
