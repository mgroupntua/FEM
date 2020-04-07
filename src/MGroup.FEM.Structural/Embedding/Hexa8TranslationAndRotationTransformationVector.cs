using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.FEM.Elements;
using MGroup.FEM.Embedding;
using MGroup.FEM.Interfaces;
using MGroup.FEM.Structural.Elements;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Structural.Embedding
{
	public class Hexa8TranslationAndRotationTransformationVector : IEmbeddedDOFInHostTransformationVector
	{
		private const int commonDofsPerNode = 3;
		private const int rotationalDofsPerNode = 3;
		private const int hostDofsPerNode = 3;
		private const int hostShapeFunctionLength = 8;
		private const int shapeFunctionOffset = hostShapeFunctionLength * (commonDofsPerNode + 1);
		private readonly IDofType[] translationAndRotationDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY, StructuralDof.RotationZ };

		public IList<IDofType> GetDependentDOFTypes { get { return translationAndRotationDOFTypes; } }

		public IReadOnlyList<IReadOnlyList<IDofType>> GetDOFTypesOfHost(EmbeddedNode node)
		{
			return node.EmbeddedInElement.ElementType.GetElementDofTypes(node.EmbeddedInElement);
		}

		private Tuple<double[,], double[,]> GetJacobiansFromShapeFunctionsVector(double[] shapeFunctionsVector)
		{
			var jacobian = new double[3, 3];
			var jacobianInverse = new double[3, 3];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					jacobian[i, j] = shapeFunctionsVector[shapeFunctionOffset + i * 3 + j];
					jacobianInverse[i, j] = shapeFunctionsVector[shapeFunctionOffset + 9 + i * 3 + j];
				}

			return new Tuple<double[,], double[,]>(jacobian, jacobianInverse);
		}

		private double[][] GetTransformationVectorForTranslationsOnly(EmbeddedNode node)
		{
			if (!(node.EmbeddedInElement.ElementType.CellType == CellType.Hexa8)
				)throw new ArgumentException("Host element is not Hexa8.");

			double[] hostShapeFunctions = ((IEmbeddedHostElement)node.EmbeddedInElement.ElementType).GetShapeFunctionsForNode(node.EmbeddedInElement, node);
			var transformation = new double[commonDofsPerNode + rotationalDofsPerNode][];
			for (int j = 0; j < commonDofsPerNode; j++)
			{
				transformation[j] = new double[hostShapeFunctionLength * hostDofsPerNode];
				for (int k = 0; k < hostShapeFunctionLength; k++)
					transformation[j][hostDofsPerNode * k + j] = hostShapeFunctions[k];
			}
			for (int j = 0; j < rotationalDofsPerNode; j++)
				transformation[commonDofsPerNode + j] = new double[hostShapeFunctionLength * hostDofsPerNode];

			var jacobianAndInverse = GetJacobiansFromShapeFunctionsVector(hostShapeFunctions);
			transformation[commonDofsPerNode - 1 + 1][0] = 0;
			transformation[commonDofsPerNode - 1 + 1][1] = -jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 9] - jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 1] - jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 17];
			transformation[commonDofsPerNode - 1 + 1][2] = jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 9] + jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 1] + jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 17];
			transformation[commonDofsPerNode - 1 + 1][3] = 0;
			transformation[commonDofsPerNode - 1 + 1][4] = -jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 10] - jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 2] - jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 18];
			transformation[commonDofsPerNode - 1 + 1][5] = jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 10] + jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 2] + jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 18];
			transformation[commonDofsPerNode - 1 + 1][6] = 0;
			transformation[commonDofsPerNode - 1 + 1][7] = -jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 11] - jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 3] - jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 19];
			transformation[commonDofsPerNode - 1 + 1][8] = jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 11] + jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 3] + jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 19];
			transformation[commonDofsPerNode - 1 + 1][9] = 0;
			transformation[commonDofsPerNode - 1 + 1][10] = -jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 12] - jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 4] - jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 20];
			transformation[commonDofsPerNode - 1 + 1][11] = jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 12] + jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 4] + jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 20];
			transformation[commonDofsPerNode - 1 + 1][12] = 0;
			transformation[commonDofsPerNode - 1 + 1][13] = -jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 13] - jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 5] - jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 21];
			transformation[commonDofsPerNode - 1 + 1][14] = jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 13] + jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 5] + jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 21];
			transformation[commonDofsPerNode - 1 + 1][15] = 0;
			transformation[commonDofsPerNode - 1 + 1][16] = -jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 14] - jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 6] - jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 22];
			transformation[commonDofsPerNode - 1 + 1][17] = jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 14] + jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 6] + jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 22];
			transformation[commonDofsPerNode - 1 + 1][18] = 0;
			transformation[commonDofsPerNode - 1 + 1][19] = -jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 15] - jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 7] - jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 23];
			transformation[commonDofsPerNode - 1 + 1][20] = jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 15] + jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 7] + jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 23];
			transformation[commonDofsPerNode - 1 + 1][21] = 0;
			transformation[commonDofsPerNode - 1 + 1][22] = -jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 16] - jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 8] - jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 24];
			transformation[commonDofsPerNode - 1 + 1][23] = jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 16] + jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 8] + jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 24];
			transformation[commonDofsPerNode - 1 + 2][0] = jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 9] + jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 1] + jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 17];
			transformation[commonDofsPerNode - 1 + 2][1] = 0;
			transformation[commonDofsPerNode - 1 + 2][2] = -jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 9] - jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 1] - jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 17];
			transformation[commonDofsPerNode - 1 + 2][3] = jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 10] + jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 2] + jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 18];
			transformation[commonDofsPerNode - 1 + 2][4] = 0;
			transformation[commonDofsPerNode - 1 + 2][5] = -jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 10] - jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 2] - jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 18];
			transformation[commonDofsPerNode - 1 + 2][6] = jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 11] + jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 3] + jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 19];
			transformation[commonDofsPerNode - 1 + 2][7] = 0;
			transformation[commonDofsPerNode - 1 + 2][8] = -jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 11] - jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 3] - jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 19];
			transformation[commonDofsPerNode - 1 + 2][9] = jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 12] + jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 4] + jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 20];
			transformation[commonDofsPerNode - 1 + 2][10] = 0;
			transformation[commonDofsPerNode - 1 + 2][11] = -jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 12] - jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 4] - jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 20];
			transformation[commonDofsPerNode - 1 + 2][12] = jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 13] + jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 5] + jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 21];
			transformation[commonDofsPerNode - 1 + 2][13] = 0;
			transformation[commonDofsPerNode - 1 + 2][14] = -jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 13] - jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 5] - jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 21];
			transformation[commonDofsPerNode - 1 + 2][15] = jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 14] + jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 6] + jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 22];
			transformation[commonDofsPerNode - 1 + 2][16] = 0;
			transformation[commonDofsPerNode - 1 + 2][17] = -jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 14] - jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 6] - jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 22];
			transformation[commonDofsPerNode - 1 + 2][18] = jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 15] + jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 7] + jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 23];
			transformation[commonDofsPerNode - 1 + 2][19] = 0;
			transformation[commonDofsPerNode - 1 + 2][20] = -jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 15] - jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 7] - jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 23];
			transformation[commonDofsPerNode - 1 + 2][21] = jacobianAndInverse.Item2[2, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 16] + jacobianAndInverse.Item2[2, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 8] + jacobianAndInverse.Item2[2, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 24];
			transformation[commonDofsPerNode - 1 + 2][22] = 0;
			transformation[commonDofsPerNode - 1 + 2][23] = -jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 16] - jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 8] - jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 24];
			transformation[commonDofsPerNode - 1 + 3][0] = -jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 9] - jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 1] - jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 17];
			transformation[commonDofsPerNode - 1 + 3][1] = jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 9] + jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 1] + jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 17];
			transformation[commonDofsPerNode - 1 + 3][2] = 0;
			transformation[commonDofsPerNode - 1 + 3][3] = -jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 10] - jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 2] - jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 18];
			transformation[commonDofsPerNode - 1 + 3][4] = jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 10] + jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 2] + jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 18];
			transformation[commonDofsPerNode - 1 + 3][5] = 0;
			transformation[commonDofsPerNode - 1 + 3][6] = -jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 11] - jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 3] - jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 19];
			transformation[commonDofsPerNode - 1 + 3][7] = jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 11] + jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 3] + jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 19];
			transformation[commonDofsPerNode - 1 + 3][8] = 0;
			transformation[commonDofsPerNode - 1 + 3][9] = -jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 12] - jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 4] - jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 20];
			transformation[commonDofsPerNode - 1 + 3][10] = jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 12] + jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 4] + jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 20];
			transformation[commonDofsPerNode - 1 + 3][11] = 0;
			transformation[commonDofsPerNode - 1 + 3][12] = -jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 13] - jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 5] - jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 21];
			transformation[commonDofsPerNode - 1 + 3][13] = jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 13] + jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 5] + jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 21];
			transformation[commonDofsPerNode - 1 + 3][14] = 0;
			transformation[commonDofsPerNode - 1 + 3][15] = -jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 14] - jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 6] - jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 22];
			transformation[commonDofsPerNode - 1 + 3][16] = jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 14] + jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 6] + jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 22];
			transformation[commonDofsPerNode - 1 + 3][17] = 0;
			transformation[commonDofsPerNode - 1 + 3][18] = -jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 15] - jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 7] - jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 23];
			transformation[commonDofsPerNode - 1 + 3][19] = jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 15] + jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 7] + jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 23];
			transformation[commonDofsPerNode - 1 + 3][20] = 0;
			transformation[commonDofsPerNode - 1 + 3][21] = -jacobianAndInverse.Item2[1, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 16] - jacobianAndInverse.Item2[1, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 8] - jacobianAndInverse.Item2[1, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 24];
			transformation[commonDofsPerNode - 1 + 3][22] = jacobianAndInverse.Item2[0, 1] * hostShapeFunctions[hostShapeFunctionLength - 1 + 16] + jacobianAndInverse.Item2[0, 0] * hostShapeFunctions[hostShapeFunctionLength - 1 + 8] + jacobianAndInverse.Item2[0, 2] * hostShapeFunctions[hostShapeFunctionLength - 1 + 24];
			transformation[commonDofsPerNode - 1 + 3][23] = 0;

			for (int j = commonDofsPerNode; j < commonDofsPerNode + 3; j++)
				for (int k = 0; k < 24; k++)
					transformation[j][k] *= 0.5;

			return transformation;
		}

		public double[][] GetTransformationVector(EmbeddedNode node)
		{
			return GetTransformationVectorForTranslationsOnly(node);
		}
	}
}
