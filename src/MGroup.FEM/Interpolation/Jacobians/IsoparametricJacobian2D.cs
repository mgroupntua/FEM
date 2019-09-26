using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Matrices;

//TODO: use Matrix2by2 after benchmarking it.
//TODO: once we know that an exception will be thrown, try to pinpoint the error: wrong node order, clockwise node order, the  
//      element's shape is too distorted, midpoints are too close to corners in quadratic elements, etc...
namespace MGroup.FEM.Interpolation.Jacobians
{
	/// <summary>
	/// This class encapsulates the determinant and inverse of the Jacobian matrix for a 2D isoparametric mapping.
	/// Let f be a mapping: x \in R^2 -> f(x) \in R^2. The Jacobian matrix of the mapping is (in numerator layout): 
	/// J = [df_1/dx_1 df_1/dx_2; df_2/dx_1 df_2/dx_2]. 
	/// Note that some sources call the transpose of this matrix as J. In FEM we are usually interested in the determinant and
	/// inverse of the Jacobian matrix. 
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class IsoparametricJacobian2D
	{
		private const double determinantTolerance = 1E-8; // This needs to be in a static settings class.

		/// <summary>
		/// The caller (usually the interpolation class) assumes responsibility for matching the nodes to the shape function 
		/// derivatives.
		/// </summary>
		/// <param name="nodes">The nodes used for the interpolation.</param>
		/// <param name="naturalDerivatives">The shape function derivatives at a specific integration point.</param>
		public IsoparametricJacobian2D(IReadOnlyList<Node> nodes, Matrix naturalDerivatives)
		{
			// The original matrix is not stored. Only the inverse and the determinant
			DirectMatrix = CalculateJacobianMatrix(nodes, naturalDerivatives);
			(InverseMatrix, DirectDeterminant) = DirectMatrix.InvertAndDeterminant();
			//(InverseMatrix, DirectDeterminant) = InvertAndDeterminant(DirectMatrix);
			if (DirectDeterminant < determinantTolerance)
			{
				throw new ArgumentException("Jacobian determinant is negative or under the allowed tolerance"
					+ $" ({DirectDeterminant} < {determinantTolerance}). Check the order of nodes or the element geometry.");
			}
		}

		/// <summary>
		/// The determinant of the direct Jacobian matrix <see cref="DirectMatrix"/>.
		/// </summary>
		public double DirectDeterminant { get; }

		/// <summary>
		/// The Jacobian matrix of the direct mapping. Numerator layout is used:
		/// J = [df_1/dx_1 df_1/dx_2; df_2/dx_1 df_2/dx_2].
		/// </summary>
		public Matrix DirectMatrix { get; }

		/// <summary>
		/// The inverse of the Jacobian matrix. Numerator layout used is used:
		/// inv(J) = [dx_1/df_1 dx_1/df_2 ; dx_2/df_1 dx_2/df_2]
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
			// naturalGradient * inverseJ = 1-by-2 * 2-by-2 = 1-by-2
			var result = new double[2];
			result[0] = naturalGradient[0] * InverseMatrix[0, 0] + naturalGradient[1] * InverseMatrix[1, 0];
			result[1] = naturalGradient[0] * InverseMatrix[0, 1] + naturalGradient[1] * InverseMatrix[1, 1];
			return result;
		}

		private static Matrix CalculateJacobianMatrix(IReadOnlyList<Node> nodes, Matrix naturalDerivatives)
		{
			//TODO: describe this as a matrix operation
			var jacobianMatrix = Matrix.CreateZero(2, 2);
			//var J = new double[2, 2];
			for (int nodeIndex = 0; nodeIndex < nodes.Count; ++nodeIndex)
			{
				double x = nodes[nodeIndex].X;
				double y = nodes[nodeIndex].Y;
				double N_xi = naturalDerivatives[nodeIndex, 0];
				double N_eta = naturalDerivatives[nodeIndex, 1];

				jacobianMatrix[0, 0] += N_xi * x;
				jacobianMatrix[0, 1] += N_eta * x;
				jacobianMatrix[1, 0] += N_xi * y;
				jacobianMatrix[1, 1] += N_eta * y;
			}
			return jacobianMatrix;
			//return Matrix2by2.CreateFromArray(J);
		}

		//private static (Matrix inverse, double det) InvertAndDeterminant(Matrix directMatrix)
		//{
		//    // Leibniz formula:
		//    double det = directMatrix[0, 0] * directMatrix[1, 1] - directMatrix[0, 1] * directMatrix[1, 0];
		//    if (Math.Abs(det) < determinantTolerance) throw new Exception(
		//        $"|Determinant| = {Math.Abs(det)} < tolerance = {determinantTolerance}. The matrix is singular");

		//    // Cramer's rule: inverse = 1/det * [a11 -a01; -a10 a00]
		//    double[,] inverse = new double[,]
		//    {
		//        { directMatrix[1, 1] / det, -directMatrix[0, 1] / det },
		//        { -directMatrix[1, 0] / det, directMatrix[0, 0] / det }
		//    };

		//    return (new Matrix2D(inverse), det);
		//}
	}
}
