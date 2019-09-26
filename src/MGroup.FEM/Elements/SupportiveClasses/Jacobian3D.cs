using System;

namespace MGroup.FEM.Elements.SupportiveClasses
{
	#region

	#endregion

	public class Jacobian3D
	{
		#region Constants and Fields

		private const int Dimensions = 3;

		private readonly double determinant;

		private readonly double[,] matrix;

		#endregion

		#region Constructors and Destructors

		public Jacobian3D(double[,] nodeCoordinates, ShapeFunctionNaturalDerivatives3D[] shapeFunctionDerivatives)
		{
			this.matrix = CalculateJacobianMatrix(nodeCoordinates, shapeFunctionDerivatives);
			this.determinant = this.CalculateJacobianDeterminant();
		}

		#endregion

		#region Properties

		public double Determinant
		{
			get
			{
				return this.determinant;
			}
		}

		public double[,] Matrix
		{
			get
			{
				return this.matrix;
			}
		}

		#endregion

		#region Public Methods

		public double[,] CalculateJacobianInverse()
		{
			double[,] jacobianInverse = new double[Dimensions, Dimensions];
			double determinantInverse = 1.0 / this.determinant;

			jacobianInverse[0, 0] = ((this.matrix[1, 1] * this.matrix[2, 2]) - (this.matrix[2, 1] * this.matrix[1, 2])) *
									determinantInverse;
			jacobianInverse[0, 1] = ((this.matrix[2, 1] * this.matrix[0, 2]) - (this.matrix[0, 1] * this.matrix[2, 2])) *
									determinantInverse;
			jacobianInverse[0, 2] = ((this.matrix[0, 1] * this.matrix[1, 2]) - (this.matrix[1, 1] * this.matrix[0, 2])) *
									determinantInverse;
			jacobianInverse[1, 0] = ((this.matrix[2, 0] * this.matrix[1, 2]) - (this.matrix[1, 0] * this.matrix[2, 2])) *
									determinantInverse;
			jacobianInverse[1, 1] = ((this.matrix[0, 0] * this.matrix[2, 2]) - (this.matrix[2, 0] * this.matrix[0, 2])) *
									determinantInverse;
			jacobianInverse[1, 2] = ((this.matrix[1, 0] * this.matrix[0, 2]) - (this.matrix[0, 0] * this.matrix[1, 2])) *
									determinantInverse;
			jacobianInverse[2, 0] = ((this.matrix[1, 0] * this.matrix[2, 1]) - (this.matrix[2, 0] * this.matrix[1, 1])) *
									determinantInverse;
			jacobianInverse[2, 1] = ((this.matrix[2, 0] * this.matrix[0, 1]) - (this.matrix[2, 1] * this.matrix[0, 0])) *
									determinantInverse;
			jacobianInverse[2, 2] = ((this.matrix[0, 0] * this.matrix[1, 1]) - (this.matrix[1, 0] * this.matrix[0, 1])) *
									determinantInverse;
			return jacobianInverse;
		}

		#endregion

		#region Methods

		private static double[,] CalculateJacobianMatrix(
			double[,] nodeCoordinates, ShapeFunctionNaturalDerivatives3D[] shapeFunctionDerivatives)
		{
			double[,] jacobianMatrix = new double[Dimensions, Dimensions];
			for (int i = 0; i < nodeCoordinates.GetLength(0); i++)
			{
				for (int j = 0; j < Dimensions; j++)
				{
					jacobianMatrix[0, j] += nodeCoordinates[i, j] * shapeFunctionDerivatives[i].Xi;
					jacobianMatrix[1, j] += nodeCoordinates[i, j] * shapeFunctionDerivatives[i].Eta;
					jacobianMatrix[2, j] += nodeCoordinates[i, j] * shapeFunctionDerivatives[i].Zeta;
				}
			}

			return jacobianMatrix;
		}

		private double CalculateJacobianDeterminant()
		{
			double det1 = this.matrix[0, 0] *
						  ((this.matrix[1, 1] * this.matrix[2, 2]) - (this.matrix[2, 1] * this.matrix[1, 2]));
			double det2 = this.matrix[0, 1] *
						  ((this.matrix[1, 0] * this.matrix[2, 2]) - (this.matrix[2, 0] * this.matrix[1, 2]));
			double det3 = this.matrix[0, 2] *
						  ((this.matrix[1, 0] * this.matrix[2, 1]) - (this.matrix[2, 0] * this.matrix[1, 1]));

			double jacobianDeterminant = det1 - det2 + det3;

			if (jacobianDeterminant < 0)
			{
				throw new InvalidOperationException("The Jacobian Determinant is negative.");
			}

			return jacobianDeterminant;
		}

		#endregion
	}
}