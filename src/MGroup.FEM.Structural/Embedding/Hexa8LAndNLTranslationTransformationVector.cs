using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Embedding;

namespace MGroup.FEM.Structural.Embedding
{
	public class Hexa8LAndNLTranslationTransformationVector : IEmbeddedDOFInHostTransformationVector
	{
		private readonly IDofType[] translationOnlyDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

		public IList<IDofType> GetDependentDOFTypes { get { return translationOnlyDOFTypes; } }

		public IReadOnlyList<IReadOnlyList<IDofType>> GetDOFTypesOfHost(EmbeddedNode node)
		{
			return node.EmbeddedInElement.GetElementDofTypes();
		}

		public double[][] GetTransformationVector(EmbeddedNode node)
		{
			CheckElementType(node.EmbeddedInElement);

			const int commonDofsPerNode = 3;
			const int hostDofsPerNode = 3;
			const int hostShapeFunctionLength = 8;
			double[] hostShapeFunctions = ((IEmbeddedHostElement)node.EmbeddedInElement).GetShapeFunctionsForNode(node.EmbeddedInElement, node);

			var transformation = new double[commonDofsPerNode][];
			for (int j = 0; j < commonDofsPerNode; j++)
			{
				transformation[j] = new double[hostShapeFunctionLength * hostDofsPerNode];
				for (int k = 0; k < hostShapeFunctionLength; k++)
					transformation[j][hostDofsPerNode * k + j] = hostShapeFunctions[k];
			}

			return transformation;
		}

		private void CheckElementType(IElementType element)
		{
			bool isHexa8= element.CellType == MSolve.Discretization.CellType.Hexa8;
			//bool validElement = element is Hexa8;
			//validElement |= element is Hexa8NonLinear;
			if (!(isHexa8)) throw new ArgumentException("Host element is not Hexa8 or Hexa8NL.");
		}
	}
}
