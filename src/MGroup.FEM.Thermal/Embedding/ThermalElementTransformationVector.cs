using System.Collections.Generic;
using MGroup.Constitutive.Thermal;
using MGroup.FEM.Embedding;
using MGroup.FEM.Interfaces;
using MGroup.MSolve.Discretization;

namespace MGroup.FEM.Thermal.Embedding
{
	public class ThermalElementTransformationVector : IEmbeddedDOFInHostTransformationVector
	{
		private readonly IDofType[] thermalDOFTypes = new IDofType[] { ThermalDof.Temperature };

		public IList<IDofType> GetDependentDOFTypes { get { return thermalDOFTypes; } }

		public IReadOnlyList<IReadOnlyList<IDofType>> GetDOFTypesOfHost(EmbeddedNode node)
		{
			return node.EmbeddedInElement.ElementType.GetElementDofTypes(node.EmbeddedInElement);
		}

		public double[][] GetTransformationVector(EmbeddedNode node)
		{
			//CheckElementType(node.EmbeddedInElement.ElementType);

			const int commonDofsPerNode = 1;
			const int hostDofsPerNode = 1;
			const int hostShapeFunctionLength = 4; //TODO: Use the interpolation for this. Probably for the next line too.
			double[] hostShapeFunctions = ((IEmbeddedHostElement)node.EmbeddedInElement.ElementType).GetShapeFunctionsForNode(node.EmbeddedInElement, node);

			var transformation = new double[commonDofsPerNode][];
			for (int j = 0; j < commonDofsPerNode; j++)
			{
				transformation[j] = new double[hostShapeFunctionLength * hostDofsPerNode];
				for (int k = 0; k < hostShapeFunctionLength; k++)
					transformation[j][hostDofsPerNode * k + j] = hostShapeFunctions[k];
			}

			return transformation;
		}
	}
}

