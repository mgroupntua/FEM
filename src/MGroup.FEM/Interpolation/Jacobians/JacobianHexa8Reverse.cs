using System;
using System.Collections.Generic;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;

//jacobianInverse and jacobianDeterminants can only be calculated during initialization (at the first configuration) and then cached
namespace MGroup.FEM.Interpolation.Jacobians
{
	public class JacobianHexa8Reverse
	{
		public static (Matrix[] jacobianInverse, double[] jacobianDeterminants)
			GetJ_0invHexaAndjacobianDeterminants(IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives, IReadOnlyList<INode> elementNodes, int nGaussPoints)
		{
			double[][,] initialCoordinatesReshaped; // dimension [] number of gpoints
			double[][,] jacobians;
			Matrix[] jacobianInverse;
			double[] jacobianDeterminants; // dimension [] number of gpoints

			double[][] initialCoordinates;
			initialCoordinates = new double[8][];
			for (int j = 0; j < 8; j++)
			{
				initialCoordinates[j] = new double[] { elementNodes[j].X, elementNodes[j].Y, elementNodes[j].Z, };
			}
			initialCoordinatesReshaped = new double[nGaussPoints][,];
			jacobians = new double[nGaussPoints][,];
			jacobianInverse = new Matrix[nGaussPoints];
			jacobianDeterminants = new double[nGaussPoints];

			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				// initialize matrices and calculate those known in the undeformed configuration
				initialCoordinatesReshaped[gpoint] = new double[8, 3];
				jacobians[gpoint] = new double[3, 3];
				jacobianInverse[gpoint] = Matrix.CreateZero(3, 3);

				//
				for (int m = 0; m < 8; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						initialCoordinatesReshaped[gpoint][m, n] = initialCoordinates[m][n];
					}
				}

				//
				for (int m = 0; m < 3; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						jacobians[gpoint][m, n] = 0;
						for (int p = 0; p < 8; p++)
						{
							jacobians[gpoint][m, n] += shapeFunctionNaturalDerivatives[gpoint][p,m] * initialCoordinatesReshaped[gpoint][p, n];
						}
					}
				}

				//
				double det1 = jacobians[gpoint][0, 0] *
						 ((jacobians[gpoint][1, 1] * jacobians[gpoint][2, 2]) - (jacobians[gpoint][2, 1] * jacobians[gpoint][1, 2]));
				double det2 = jacobians[gpoint][0, 1] *
							  ((jacobians[gpoint][1, 0] * jacobians[gpoint][2, 2]) - (jacobians[gpoint][2, 0] * jacobians[gpoint][1, 2]));
				double det3 = jacobians[gpoint][0, 2] *
							  ((jacobians[gpoint][1, 0] * jacobians[gpoint][2, 1]) - (jacobians[gpoint][2, 0] * jacobians[gpoint][1, 1]));
				double jacobianDeterminant = det1 - det2 + det3;
				if (jacobianDeterminant < 0)
				{
					throw new InvalidOperationException("The Jacobian Determinant is negative.");
				}
				jacobianDeterminants[gpoint] = jacobianDeterminant;

				//
				jacobianInverse[gpoint][0, 0] = ((jacobians[gpoint][1, 1] * jacobians[gpoint][2, 2]) - (jacobians[gpoint][2, 1] * jacobians[gpoint][1, 2])) *
									(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][0, 1] = ((jacobians[gpoint][2, 1] * jacobians[gpoint][0, 2]) - (jacobians[gpoint][0, 1] * jacobians[gpoint][2, 2])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][0, 2] = ((jacobians[gpoint][0, 1] * jacobians[gpoint][1, 2]) - (jacobians[gpoint][1, 1] * jacobians[gpoint][0, 2])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][1, 0] = ((jacobians[gpoint][2, 0] * jacobians[gpoint][1, 2]) - (jacobians[gpoint][1, 0] * jacobians[gpoint][2, 2])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][1, 1] = ((jacobians[gpoint][0, 0] * jacobians[gpoint][2, 2]) - (jacobians[gpoint][2, 0] * jacobians[gpoint][0, 2])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][1, 2] = ((jacobians[gpoint][1, 0] * jacobians[gpoint][0, 2]) - (jacobians[gpoint][0, 0] * jacobians[gpoint][1, 2])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][2, 0] = ((jacobians[gpoint][1, 0] * jacobians[gpoint][2, 1]) - (jacobians[gpoint][2, 0] * jacobians[gpoint][1, 1])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][2, 1] = ((jacobians[gpoint][2, 0] * jacobians[gpoint][0, 1]) - (jacobians[gpoint][2, 1] * jacobians[gpoint][0, 0])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][2, 2] = ((jacobians[gpoint][0, 0] * jacobians[gpoint][1, 1]) - (jacobians[gpoint][1, 0] * jacobians[gpoint][0, 1])) *
										(1 / jacobianDeterminants[gpoint]);
			}


			return (jacobianInverse, jacobianDeterminants);


		}
		public static (double[][,] jacobianInverse, double[] jacobianDeterminants)
			GetJacobiansAndJacobianDeterminants(double[][,] shapeFunctionNaturalDerivatives, IList<INode> elementNodes, int nGaussPoints)
		{
			double[][,] nodeCoordinatesReshaped; // dimension [] number of gpoints
			double[][,] jacobians;
			double[][,] jacobianInverse;
			double[] jacobianDeterminants; // dimension [] number of gpoints

			double[][] initialCoordinates;
			initialCoordinates = new double[8][];
			for (int j = 0; j < 8; j++)
			{
				initialCoordinates[j] = new double[] { elementNodes[j].X, elementNodes[j].Y, elementNodes[j].Z, };
			}
			nodeCoordinatesReshaped = new double[nGaussPoints][,];
			jacobians = new double[nGaussPoints][,];
			jacobianInverse = new double[nGaussPoints][,];
			jacobianDeterminants = new double[nGaussPoints];

			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				// initialize matrices and calculate those known in the undeformed configuration
				nodeCoordinatesReshaped[gpoint] = new double[8, 3];
				jacobians[gpoint] = new double[3, 3];
				jacobianInverse[gpoint] = new double[3, 3];

				//
				for (int m = 0; m < 8; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						nodeCoordinatesReshaped[gpoint][m, n] = initialCoordinates[m][n];
					}
				}

				//
				for (int m = 0; m < 3; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						jacobians[gpoint][m, n] = 0;
						for (int p = 0; p < 8; p++)
						{
							jacobians[gpoint][m, n] += shapeFunctionNaturalDerivatives[gpoint][m, p] * nodeCoordinatesReshaped[gpoint][p, n];
						}
					}
				}

				//
				double det1 = jacobians[gpoint][0, 0] *
						 ((jacobians[gpoint][1, 1] * jacobians[gpoint][2, 2]) - (jacobians[gpoint][2, 1] * jacobians[gpoint][1, 2]));
				double det2 = jacobians[gpoint][0, 1] *
							  ((jacobians[gpoint][1, 0] * jacobians[gpoint][2, 2]) - (jacobians[gpoint][2, 0] * jacobians[gpoint][1, 2]));
				double det3 = jacobians[gpoint][0, 2] *
							  ((jacobians[gpoint][1, 0] * jacobians[gpoint][2, 1]) - (jacobians[gpoint][2, 0] * jacobians[gpoint][1, 1]));
				double jacobianDeterminant = det1 - det2 + det3;
				if (jacobianDeterminant < 0)
				{
					throw new InvalidOperationException("The Jacobian Determinant is negative.");
				}
				jacobianDeterminants[gpoint] = jacobianDeterminant;

				//
				jacobianInverse[gpoint][0, 0] = ((jacobians[gpoint][1, 1] * jacobians[gpoint][2, 2]) - (jacobians[gpoint][2, 1] * jacobians[gpoint][1, 2])) *
									(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][0, 1] = ((jacobians[gpoint][2, 1] * jacobians[gpoint][0, 2]) - (jacobians[gpoint][0, 1] * jacobians[gpoint][2, 2])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][0, 2] = ((jacobians[gpoint][0, 1] * jacobians[gpoint][1, 2]) - (jacobians[gpoint][1, 1] * jacobians[gpoint][0, 2])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][1, 0] = ((jacobians[gpoint][2, 0] * jacobians[gpoint][1, 2]) - (jacobians[gpoint][1, 0] * jacobians[gpoint][2, 2])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][1, 1] = ((jacobians[gpoint][0, 0] * jacobians[gpoint][2, 2]) - (jacobians[gpoint][2, 0] * jacobians[gpoint][0, 2])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][1, 2] = ((jacobians[gpoint][1, 0] * jacobians[gpoint][0, 2]) - (jacobians[gpoint][0, 0] * jacobians[gpoint][1, 2])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][2, 0] = ((jacobians[gpoint][1, 0] * jacobians[gpoint][2, 1]) - (jacobians[gpoint][2, 0] * jacobians[gpoint][1, 1])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][2, 1] = ((jacobians[gpoint][2, 0] * jacobians[gpoint][0, 1]) - (jacobians[gpoint][2, 1] * jacobians[gpoint][0, 0])) *
										(1 / jacobianDeterminants[gpoint]);
				jacobianInverse[gpoint][2, 2] = ((jacobians[gpoint][0, 0] * jacobians[gpoint][1, 1]) - (jacobians[gpoint][1, 0] * jacobians[gpoint][0, 1])) *
										(1 / jacobianDeterminants[gpoint]);
			}


			return (jacobianInverse, jacobianDeterminants);


		}

		public static Matrix[] Get_jacobiansDeformedMatrices(int nGaussPoints, double[][] deformedCoordinates, IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives)
		{
			double[,] jacobiansDeformedMatricesb = new double[8, 3];
			Matrix[] jacobiansDeformedMatrices = new Matrix[nGaussPoints];

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				jacobiansDeformedMatrices[npoint] = Matrix.CreateZero(3, 3);
			}

			for (int m = 0; m < 8; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					//ll2[m, n] = totalDisplacements[m][n];
					jacobiansDeformedMatricesb[m, n] = deformedCoordinates[m][n];
				}
			}

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				for (int m = 0; m < 3; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						jacobiansDeformedMatrices[npoint][m, n] = 0;
						for (int p = 0; p < 8; p++)
						{
							jacobiansDeformedMatrices[npoint][m, n] += shapeFunctionNaturalDerivatives[npoint][p,m] * jacobiansDeformedMatricesb[p, n];
						}
					}
				}
			}

			return jacobiansDeformedMatrices;
		}

		public static double[][,] GetJacobiansDeformedMatrices(int nGaussPoints, double[][] deformedCoordinates, double[][,] shapeFunctionNaturalDerivatives)
		{
			double[,] jacobiansDeformedMatricesb = new double[8, 3];
			double[][,] jacobiansDeformedMatrices = new double[nGaussPoints][,];

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				jacobiansDeformedMatrices[npoint] = new double[3, 3];
			}

			for (int m = 0; m < 8; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					//ll2[m, n] = totalDisplacements[m][n];
					jacobiansDeformedMatricesb[m, n] = deformedCoordinates[m][n];
				}
			}

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				for (int m = 0; m < 3; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						jacobiansDeformedMatrices[npoint][m, n] = 0;
						for (int p = 0; p < 8; p++)
						{
							jacobiansDeformedMatrices[npoint][m, n] += shapeFunctionNaturalDerivatives[npoint][m, p] * jacobiansDeformedMatricesb[p, n];
						}
					}
				}
			}

			return jacobiansDeformedMatrices;
		}
	}
}
