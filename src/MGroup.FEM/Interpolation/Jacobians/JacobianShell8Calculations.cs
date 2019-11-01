using System;
using System.Collections.Generic;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Integration.Quadratures;

//J_0a and ll1 can only be calculated during initialization (at the first configuration) and then cached
namespace MGroup.FEM.Interpolation.Jacobians
{
	public class JacobianShell8Calculations
	{
		public static (double[][,] ll1, double[][,] J_0a)
			Getll1AndJ_0a(int nGaussPoints, double[] tk, double[][] gausscoordinates, double[][] shapeFunctions, double[][] shapeFunctionDerivatives)
		{
			double[][,] ll1;
			ll1 = new double[nGaussPoints][,];
			for (int j = 0; j < nGaussPoints; j++)
			{ ll1[j] = new double[3, 24]; }
			for (int j = 0; j < nGaussPoints; j++) //calcualte ll1 as a whole, at each gp
			{
				for (int k = 0; k < 8; k++)
				{
					ll1[j][0, 3 * k] = shapeFunctionDerivatives[k][j];
					ll1[j][0, 3 * k + 1] = 0.5 * gausscoordinates[2][j] * tk[k] * shapeFunctionDerivatives[k][j];
					ll1[j][0, 3 * k + 2] = -ll1[j][0, 3 * k + 1];
					ll1[j][1, 3 * k] = shapeFunctionDerivatives[k + 8][j];
					ll1[j][1, 3 * k + 1] = 0.5 * gausscoordinates[2][j] * tk[k] * shapeFunctionDerivatives[k + 8][j];
					ll1[j][1, 3 * k + 2] = -ll1[j][1, 3 * k + 1];
					ll1[j][2, 3 * k] = 0;
					ll1[j][2, 3 * k + 1] = 0.5 * tk[k] * shapeFunctions[k][j];
					ll1[j][2, 3 * k + 2] = -ll1[j][2, 3 * k + 1];
				}

			}

			double[][,] J_0a;//final values (not precalculated )
			J_0a = new double[nGaussPoints][,];
			for (int j = 0; j < nGaussPoints; j++)
			{ J_0a[j] = new double[3, 16]; }
			for (int j = 0; j < nGaussPoints; j++)
			{
				for (int k = 0; k < 8; k++)
				{
					J_0a[j][0, 2 * k] = ll1[j][0, 3 * k];
					J_0a[j][0, 2 * k + 1] = ll1[j][0, 3 * k + 1];
					J_0a[j][1, 2 * k] = ll1[j][1, 3 * k];
					J_0a[j][1, 2 * k + 1] = ll1[j][1, 3 * k + 1];
					J_0a[j][2, 2 * k] = ll1[j][2, 3 * k];
					J_0a[j][2, 2 * k + 1] = ll1[j][2, 3 * k + 1];
				}
			}

			return (ll1, J_0a);
		}

		public static (Matrix[] ll1, Matrix[] J_0a) Getll1AndJ_0a(IQuadrature3D quadrature, double[] tk,
			IReadOnlyList<double[]> shapeFunctions, IReadOnlyList<Matrix> shapeFunctionDerivatives)
		{
			int nGaussPoints = quadrature.IntegrationPoints.Count;
			Matrix[] ll1;
			ll1 = new Matrix[nGaussPoints];
			for (int j = 0; j < nGaussPoints; j++)
			{ ll1[j] = Matrix.CreateZero(3, 24); }
			for (int j = 0; j < nGaussPoints; j++) //calcualte ll1 as a whole, at each gp
			{
				var gaussPoint = quadrature.IntegrationPoints[j];
				for (int k = 0; k < 8; k++)
				{
					ll1[j][0, 3 * k] = shapeFunctionDerivatives[j][k, 0];
					ll1[j][0, 3 * k + 1] = 0.5 * gaussPoint.Zeta * tk[k] * shapeFunctionDerivatives[j][k, 0];
					ll1[j][0, 3 * k + 2] = -ll1[j][0, 3 * k + 1];
					ll1[j][1, 3 * k] = shapeFunctionDerivatives[j][k, 1];
					ll1[j][1, 3 * k + 1] = 0.5 * gaussPoint.Zeta * tk[k] * shapeFunctionDerivatives[j][k, 1];
					ll1[j][1, 3 * k + 2] = -ll1[j][1, 3 * k + 1];
					ll1[j][2, 3 * k] = 0;
					ll1[j][2, 3 * k + 1] = 0.5 * tk[k] * shapeFunctions[j][k];
					ll1[j][2, 3 * k + 2] = -ll1[j][2, 3 * k + 1];
				}

			}

			Matrix[] J_0a;//final values (not precalculated )
			J_0a = new Matrix[nGaussPoints];
			for (int j = 0; j < nGaussPoints; j++)
			{ J_0a[j] = Matrix.CreateZero(3, 16); }
			for (int j = 0; j < nGaussPoints; j++)
			{
				for (int k = 0; k < 8; k++)
				{
					J_0a[j][0, 2 * k] = ll1[j][0, 3 * k];
					J_0a[j][0, 2 * k + 1] = ll1[j][0, 3 * k + 1];
					J_0a[j][1, 2 * k] = ll1[j][1, 3 * k];
					J_0a[j][1, 2 * k + 1] = ll1[j][1, 3 * k + 1];
					J_0a[j][2, 2 * k] = ll1[j][2, 3 * k];
					J_0a[j][2, 2 * k + 1] = ll1[j][2, 3 * k + 1];
				}
			}

			return (ll1, J_0a);
		}

		public static (Matrix[] J_0inv, double[] detJ_0)
			GetJ_0invAndDetJ_0(Matrix[] J_0a, IReadOnlyList<INode> elementNodes, double[][] oVn_i, int nGaussPoints)
		{
			double[][] ox_i = new double[8][];
			for (int j = 0; j < 8; j++)
			{
				ox_i[j] = new double[] { elementNodes[j].X, elementNodes[j].Y, elementNodes[j].Z, };
			}

			double[,] J_0b;    //the same for all of the gauss points
			double[][,] J_0;       //different value for each gausspoint //final values (not precalculated )

			J_0b = new double[16, 3];
			for (int j = 0; j < 8; j++)
			{
				J_0b[2 * j, 0] = ox_i[j][0];
				J_0b[2 * j + 1, 0] = oVn_i[j][0];
				J_0b[2 * j, 1] = ox_i[j][1];
				J_0b[2 * j + 1, 1] = oVn_i[j][1];
				J_0b[2 * j, 2] = ox_i[j][2];
				J_0b[2 * j + 1, 2] = oVn_i[j][2];
			}

			J_0 = new double[nGaussPoints][,];
			for (int j = 0; j < nGaussPoints; j++)
			{ J_0[j] = new double[3, 3]; }
			for (int j = 0; j < nGaussPoints; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						J_0[j][k, l] = 0;
						for (int m = 0; m < 16; m++)
						{
							J_0[j][k, l] += J_0a[j][k, m] * J_0b[m, l];
						}

					}

				}
			}

			double[] detJ_0 = new double[nGaussPoints];
			for (int j = 0; j < nGaussPoints; j++)
			{
				double det1 = J_0[j][0, 0] *
					 ((J_0[j][1, 1] * J_0[j][2, 2]) - (J_0[j][2, 1] * J_0[j][1, 2]));
				double det2 = J_0[j][0, 1] *
							  ((J_0[j][1, 0] * J_0[j][2, 2]) - (J_0[j][2, 0] * J_0[j][1, 2]));
				double det3 = J_0[j][0, 2] *
							  ((J_0[j][1, 0] * J_0[j][2, 1]) - (J_0[j][2, 0] * J_0[j][1, 1]));

				double jacobianDeterminant = det1 - det2 + det3;

				if (jacobianDeterminant < 0)
				{
					throw new InvalidOperationException("The Jacobian Determinant is negative.");
				}

				detJ_0[j] = jacobianDeterminant;
			}

			Matrix[] J_0inv;
			J_0inv = new Matrix[nGaussPoints];
			for (int j = 0; j < nGaussPoints; j++) J_0inv[j] = Matrix.CreateZero(3, 3);
			for (int j = 0; j < nGaussPoints; j++)
			{
				J_0inv[j][0, 0] = ((J_0[j][1, 1] * J_0[j][2, 2]) - (J_0[j][2, 1] * J_0[j][1, 2])) *
								(1 / detJ_0[j]);
				J_0inv[j][0, 1] = ((J_0[j][2, 1] * J_0[j][0, 2]) - (J_0[j][0, 1] * J_0[j][2, 2])) *
										(1 / detJ_0[j]);
				J_0inv[j][0, 2] = ((J_0[j][0, 1] * J_0[j][1, 2]) - (J_0[j][1, 1] * J_0[j][0, 2])) *
										(1 / detJ_0[j]);
				J_0inv[j][1, 0] = ((J_0[j][2, 0] * J_0[j][1, 2]) - (J_0[j][1, 0] * J_0[j][2, 2])) *
										(1 / detJ_0[j]);
				J_0inv[j][1, 1] = ((J_0[j][0, 0] * J_0[j][2, 2]) - (J_0[j][2, 0] * J_0[j][0, 2])) *
										(1 / detJ_0[j]);
				J_0inv[j][1, 2] = ((J_0[j][1, 0] * J_0[j][0, 2]) - (J_0[j][0, 0] * J_0[j][1, 2])) *
										(1 / detJ_0[j]);
				J_0inv[j][2, 0] = ((J_0[j][1, 0] * J_0[j][2, 1]) - (J_0[j][2, 0] * J_0[j][1, 1])) *
										(1 / detJ_0[j]);
				J_0inv[j][2, 1] = ((J_0[j][2, 0] * J_0[j][0, 1]) - (J_0[j][2, 1] * J_0[j][0, 0])) *
										(1 / detJ_0[j]);
				J_0inv[j][2, 2] = ((J_0[j][0, 0] * J_0[j][1, 1]) - (J_0[j][1, 0] * J_0[j][0, 1])) *
										(1 / detJ_0[j]);
			}

			return (J_0inv, detJ_0);
		}

		public static Matrix[] Get_J_1(int nGaussPoints, double[][] tx_i, double[][] tU, Matrix[] J_0a)
		{
			double[,] J_1b;

			J_1b = new double[16, 3];
			var J_1 = new Matrix[nGaussPoints];

			for (int j = 0; j < nGaussPoints; j++)
			{
				J_1[j] = Matrix.CreateZero(3, 3);
			}

			for (int j = 0; j < 8; j++)
			{
				J_1b[2 * j, 0] = tx_i[j][0];
				J_1b[2 * j + 1, 0] = tU[j][3];
				J_1b[2 * j, 1] = tx_i[j][1];
				J_1b[2 * j + 1, 1] = tU[j][4];
				J_1b[2 * j, 2] = tx_i[j][2];
				J_1b[2 * j + 1, 2] = tU[j][5];
			}

			for (int j = 0; j < nGaussPoints; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						J_1[j][k, l] = 0;
						for (int m = 0; m < 16; m++)
						{
							J_1[j][k, l] += J_0a[j][k, m] * J_1b[m, l];
						}

					}

				}
			}
			return J_1;
		}

	}
}
