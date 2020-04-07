using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;

//TODO: Use Matrix3by3 after benchmarking it.
namespace MGroup.FEM.Interpolation.Jacobians
{
	/// <summary>
	/// This class encapsulates the determinant and inverse of the Jacobian matrix for a 3D isoparametric mapping.
	/// Let f be a mapping: x \in R^3 -> f(x) \in R^3. The Jacobian matrix of the mapping is (in numerator layout): 
	/// J = [df_1/dx_1 df_1/dx_2 df_1/dx_3; df_2/dx_1 df_2/dx_2 df_2/dx_3; df_3/dx_1 df_3/dx_2 df_3/dx_3]. 
	/// Note that some sources call the transpose of this matrix as J. In FEM we are usually interested in the determinant and 
	/// inverse of the Jacobian matrix.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class IsoparametricJacobian3D
	{
		private const double determinantTolerance = 1E-8;

		/// <summary>
		/// The caller (usually the interpolation class) assumes responsibility for matching the nodes to the shape function 
		/// derivatives.
		/// </summary>
		/// <param name="nodes">The nodes used for the interpolation.</param>
		/// <param name="naturalCoordinates">The shape function derivatives at a specific integration point.</param>
		public IsoparametricJacobian3D(IReadOnlyList<INode> nodes, Matrix naturalDerivatives)
		{
			DirectMatrix = CalculateJacobianMatrix(nodes, naturalDerivatives);
			(InverseMatrix, DirectDeterminant) = DirectMatrix.InvertAndDeterminant();
			//(InverseMatrix, DirectDeterminant) = InvertAndDeterminant(DirectMatrix);
			if (DirectDeterminant < determinantTolerance)
			{
				throw new ArgumentException("Jacobian determinant is negative or under the allowed tolerance"
					+ $" ({DirectDeterminant} < {determinantTolerance}). Check the order of nodes or the element geometry.");
			}
		}

		public IsoparametricJacobian3D(double[][] deformedStateNodeCoordinates, Matrix naturalDerivatives, bool CalculateInverseAndDetereminant)
		{
			DirectMatrix = CalculateJacobianMatrix(deformedStateNodeCoordinates, naturalDerivatives);
			if (CalculateInverseAndDetereminant)
			{
				(InverseMatrix, DirectDeterminant) = DirectMatrix.InvertAndDeterminant();
				//(InverseMatrix, DirectDeterminant) = InvertAndDeterminant(DirectMatrix);
				if (DirectDeterminant < determinantTolerance)
				{
					throw new ArgumentException("Jacobian determinant is negative or under the allowed tolerance"
						+ $" ({DirectDeterminant} < {determinantTolerance}). Check the order of nodes or the element geometry.");
				}
			}
		}


		/// <summary>
		/// The determinant of the direct Jacobian matrix <see cref="DirectMatrix"/>.
		/// </summary>
		public double DirectDeterminant { get; }

		/// <summary>
		/// The Jacobian matrix of the direct mapping. Numerator layout is used:
		/// J = [df_1/dx_1 df_1/dx_2 df_1/dx_3; df_2/dx_1 df_2/dx_2 df_2/dx_3; df_3/dx_1 df_3/dx_2 df_3/dx_3].
		/// </summary>
		public Matrix DirectMatrix { get; }

		/// <summary>
		/// The inverse of the Jacobian matrix. Numerator layout used is used:
		/// inv(J) = [dx_1/df_1 dx_1/df_2 dx_1/df_3; dx_2/df_1 dx_2/df_2 dx_2/df_3; dx_3/df_1 dx_3/df_2 dx_3/df_3]
		/// </summary>
		public Matrix InverseMatrix { get; }

		/// <summary>
		/// Transforms the gradient of a vector-valued function from the natural to the global cartesian coordinate system.
		/// </summary>
		/// <param name="naturalGradient">The gradient of a vector-valued function in the natural coordinate system. Each row 
		///     corresponds to the gradient of a single component of the vector function. Each column corresponds to the 
		///     derivatives of all components with respect to a single coordinate.</param>
		public Matrix TransformNaturalDerivativesToCartesian(Matrix naturalGradient)
			=> InverseMatrix.MultiplyLeft(naturalGradient);

		/// <summary>
		/// Transforms the gradient of a scalar-valued function from the natural to the global cartesian coordinate system.
		/// </summary>
		/// <param name="naturalGradient">The gradient of a scalar-valued function in the natural coordinate system. Each entry 
		///     corresponds to the derivative with respect to a single coordinate.</param>
		public double[] TransformNaturalDerivativesToCartesian(double[] naturalGradient) //TODO: rowVector * matrix
		{
			var result = new double[3];
			result[0] = naturalGradient[0] * InverseMatrix[0, 0] + naturalGradient[1] * InverseMatrix[1, 0] +
						naturalGradient[2] * InverseMatrix[2, 0];

			result[1] = naturalGradient[0] * InverseMatrix[0, 1] + naturalGradient[1] * InverseMatrix[1, 1] +
						naturalGradient[2] * InverseMatrix[2, 1];

			result[2] = naturalGradient[0] * InverseMatrix[0, 2] + naturalGradient[1] * InverseMatrix[1, 2] +
						naturalGradient[2] * InverseMatrix[2, 2];
			return result;
		}


		private static Matrix CalculateJacobianMatrix(IReadOnlyList<INode> nodes, Matrix naturalDerivatives)
		{
			var jacobianMatrix = Matrix.CreateZero(3, 3);
			//var jacobianMatrix = new double[3,3];

			for (int nodeIndex = 0; nodeIndex < nodes.Count; nodeIndex++)
			{
				jacobianMatrix[0, 0] += naturalDerivatives[nodeIndex, 0] * nodes[nodeIndex].X;
				jacobianMatrix[0, 1] += naturalDerivatives[nodeIndex, 1] * nodes[nodeIndex].X;
				jacobianMatrix[0, 2] += naturalDerivatives[nodeIndex, 2] * nodes[nodeIndex].X;

				jacobianMatrix[1, 0] += naturalDerivatives[nodeIndex, 0] * nodes[nodeIndex].Y;
				jacobianMatrix[1, 1] += naturalDerivatives[nodeIndex, 1] * nodes[nodeIndex].Y;
				jacobianMatrix[1, 2] += naturalDerivatives[nodeIndex, 2] * nodes[nodeIndex].Y;

				jacobianMatrix[2, 0] += naturalDerivatives[nodeIndex, 0] * nodes[nodeIndex].Z;
				jacobianMatrix[2, 1] += naturalDerivatives[nodeIndex, 1] * nodes[nodeIndex].Z;
				jacobianMatrix[2, 2] += naturalDerivatives[nodeIndex, 2] * nodes[nodeIndex].Z;
			}

			return jacobianMatrix;
			//return Matrix3by3.CreateFromArray(jacobianMatrix, false);
		}

		private Matrix CalculateJacobianMatrix(double[][] stateCoordinates, Matrix naturalDerivatives)
		{
			var stateJacobian = Matrix.CreateZero(3, 3);

			for (int m = 0; m < 3; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					for (int p = 0; p < naturalDerivatives.NumRows; p++)
					{
						stateJacobian[m, n] += naturalDerivatives[p, m] * stateCoordinates[p][n];
					}
				}
			}

			return stateJacobian;
		}

		//private static (Matrix inverse, double determinant) InvertAndDeterminant(Matrix directMatrix)
		//{
		//    double determinant = directMatrix[0, 0] *
		//                         (directMatrix[1, 1] * directMatrix[2, 2] - directMatrix[2, 1] * directMatrix[1, 2])
		//                         - directMatrix[0, 1] * (directMatrix[1, 0] * directMatrix[2, 2] -
		//                                                   directMatrix[2, 0] * directMatrix[1, 2])
		//                         + directMatrix[0, 2] * (directMatrix[1, 0] * directMatrix[2, 1] -
		//                                                   directMatrix[2, 0] * directMatrix[1, 1]);
		//    if (Math.Abs(determinant) < determinantTolerance) throw new Exception(
		//        $"|Determinant| = {Math.Abs(determinant)} < tolerance = {determinantTolerance}. The matrix is singular");

		//    var inverseJacobian = new double[3, 3];

		//    inverseJacobian[0, 0] = (directMatrix[1, 1] * directMatrix[2, 2] - directMatrix[1, 2] * directMatrix[2, 1]) / determinant;
		//    inverseJacobian[0, 1] = (directMatrix[0, 2] * directMatrix[2, 1] - directMatrix[0, 1] * directMatrix[2, 2]) / determinant;
		//    inverseJacobian[0, 2] = (directMatrix[0, 1] * directMatrix[1, 2] - directMatrix[0, 2] * directMatrix[1, 1]) / determinant;

		//    inverseJacobian[1, 0] = (directMatrix[1, 2] * directMatrix[2, 0] - directMatrix[1, 0] * directMatrix[2, 2]) / determinant;
		//    inverseJacobian[1, 1] = (directMatrix[0, 0] * directMatrix[2, 2] - directMatrix[0, 2] * directMatrix[2, 0]) / determinant;
		//    inverseJacobian[1, 2] = (directMatrix[0, 2] * directMatrix[1, 0] - directMatrix[0, 0] * directMatrix[1, 2]) / determinant;

		//    inverseJacobian[2, 0] = (directMatrix[1, 0] * directMatrix[2, 1] - directMatrix[1, 1] * directMatrix[2, 0]) / determinant;
		//    inverseJacobian[2, 1] = (directMatrix[0, 1] * directMatrix[2, 0] - directMatrix[0, 0] * directMatrix[2, 1]) / determinant;
		//    inverseJacobian[2, 2] = (directMatrix[0, 0] * directMatrix[1, 1] - directMatrix[0, 1] * directMatrix[1, 0]) / determinant;

		//    return (new Matrix2D(inverseJacobian), determinant);
		//}
	}
}
