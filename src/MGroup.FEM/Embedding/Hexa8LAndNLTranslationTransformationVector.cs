using System;
using System.Collections.Generic;
using MGroup.FEM.Elements;
using MGroup.FEM.Interfaces;
using MGroup.MSolve.Discretization.FreedomDegrees;

namespace MGroup.FEM.Embedding
{
    public class Hexa8LAndNLTranslationTransformationVector : IEmbeddedDOFInHostTransformationVector
    {
        private readonly IDofType[] translationOnlyDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

        public IList<IDofType> GetDependentDOFTypes { get { return translationOnlyDOFTypes; } }

        public IReadOnlyList<IReadOnlyList<IDofType>> GetDOFTypesOfHost(EmbeddedNode node)
        {
            return node.EmbeddedInElement.ElementType.GetElementDofTypes(node.EmbeddedInElement);
        }

        public double[][] GetTransformationVector(EmbeddedNode node)
        {
            CheckElementType(node.EmbeddedInElement.ElementType);

            const int commonDofsPerNode = 3;
            const int hostDofsPerNode = 3;
            const int hostShapeFunctionLength = 8;
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

        private void CheckElementType(IFiniteElement element)
        {
            bool validElement = element is Hexa8;
            validElement |= element is Hexa8NonLinear;
            if (!(validElement)) throw new ArgumentException("Host element is not Hexa8 or Hexa8NL.");
        }
    }
}
