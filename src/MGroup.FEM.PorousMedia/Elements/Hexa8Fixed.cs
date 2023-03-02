using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using MGroup.Constitutive.Structural;
using MGroup.FEM.PorousMedia.Helpers;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Embedding;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.DataStructures;

//TODO: get rid of Hexa8.cs

namespace MGroup.FEM.PorousMedia.Elements
{
	public class Hexa8Fixed : IStructuralElementType, IEmbeddedHostElement
	{
		protected static double determinantTolerance = 0.00000001;
		protected static int iInt = 2;
		protected static int iInt2 = iInt * iInt;
		protected static int iInt3 = iInt2 * iInt;
		protected readonly static IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		protected readonly static IDofType[][] dofTypes = new IDofType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
			nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
		protected readonly IContinuumMaterial3D[] materialsAtGaussPoints;
		protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private double[][] lastStresses;

		#region Fortran imports
		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH8GAUSSMATRICES",
		//	CallingConvention = CallingConvention.Cdecl)]
		//protected static extern void CalcH8GaussMatrices(ref int iInt, [MarshalAs(UnmanagedType.LPArray)] double[,] faXYZ,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faWeight, [MarshalAs(UnmanagedType.LPArray)] double[,] faS,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,] faDS, [MarshalAs(UnmanagedType.LPArray)] double[,,] faJ,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faDetJ, [MarshalAs(UnmanagedType.LPArray)] double[,,] faB);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH8STRAINS",
		//	CallingConvention = CallingConvention.Cdecl)]
		//protected static extern void CalcH8Strains(ref int iInt,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,,] faB, [MarshalAs(UnmanagedType.LPArray)] double[] fau,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,] faStrains);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH8FORCES",
		//	CallingConvention = CallingConvention.Cdecl)]
		//protected static extern void CalcH8Forces(ref int iInt,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,,] faB, [MarshalAs(UnmanagedType.LPArray)] double[] faWeight,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,] faStresses,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faForces);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH8K",
		//	CallingConvention = CallingConvention.Cdecl)]
		//protected static extern void CalcH8K(ref int iInt, [MarshalAs(UnmanagedType.LPArray)] double[,,] faE,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,,] faB, [MarshalAs(UnmanagedType.LPArray)] double[] faWeight,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faK);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH8MLUMPED",
		//	CallingConvention = CallingConvention.Cdecl)]
		//protected static extern void CalcH8MLumped(ref int iInt, ref double fDensity,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faWeight, [MarshalAs(UnmanagedType.LPArray)] double[] faM);
		#endregion

		protected Hexa8Fixed()
		{
		}

		public Hexa8Fixed(IContinuumMaterial3D material)
		{
			lastStresses = new double[iInt3][];
			materialsAtGaussPoints = new IContinuumMaterial3D[iInt3];
			for (int i = 0; i < iInt3; i++)
				materialsAtGaussPoints[i] = (IContinuumMaterial3D)material.Clone();
		}

		public Hexa8Fixed(IContinuumMaterial3D material, IElementDofEnumerator dofEnumerator)
			: this(material)
		{
			this.dofEnumerator = dofEnumerator;
		}

		public CellType CellType { get; } = CellType.Hexa8;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		public double Density { get; set; }
		public double RayleighAlpha { get; set; }
		public double RayleighBeta { get; set; }

		protected double[,] GetCoordinates()
		{
			double[,] faXYZ = new double[dofTypes.Length, 3];
			for (int i = 0; i < dofTypes.Length; i++)
			{
				faXYZ[i, 0] = Nodes[i].X;
				faXYZ[i, 1] = Nodes[i].Y;
				faXYZ[i, 2] = Nodes[i].Z;
			}
			return faXYZ;
		}

		protected double[,] GetCoordinatesTranspose(IElementType element)
		{
			double[,] faXYZ = new double[3, dofTypes.Length];
			for (int i = 0; i < dofTypes.Length; i++)
			{
				faXYZ[0, i] = element.Nodes[i].X;
				faXYZ[1, i] = element.Nodes[i].Y;
				faXYZ[2, i] = element.Nodes[i].Z;
			}
			return faXYZ;
		}

		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public virtual IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofTypes;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		private double[] CalcH8Shape(double fXi, double fEta, double fZeta)
		{
			const double fSqC125 = 0.5;
			double fXiP = (1.0 + fXi) * fSqC125;
			double fEtaP = (1.0 + fEta) * fSqC125;
			double fZetaP = (1.0 + fZeta) * fSqC125;
			double fXiM = (1.0 - fXi) * fSqC125;
			double fEtaM = (1.0 - fEta) * fSqC125;
			double fZetaM = (1.0 - fZeta) * fSqC125;

			return new double[]
			{
				fXiM * fEtaM * fZetaM,
				fXiP * fEtaM * fZetaM,
				fXiP * fEtaP * fZetaM,
				fXiM * fEtaP * fZetaM,
				fXiM * fEtaM * fZetaP,
				fXiP * fEtaM * fZetaP,
				fXiP * fEtaP * fZetaP,
				fXiM * fEtaP * fZetaP
			};
		}

		private double[] CalcH8NablaShape(double fXi, double fEta, double fZeta)
		{
			const double fSq125 = 0.35355339059327376220042218105242;
			double fXiP = (1.0 + fXi) * fSq125;
			double fEtaP = (1.0 + fEta) * fSq125;
			double fZetaP = (1.0 + fZeta) * fSq125;
			double fXiM = (1.0 - fXi) * fSq125;
			double fEtaM = (1.0 - fEta) * fSq125;
			double fZetaM = (1.0 - fZeta) * fSq125;

			double[] faDS = new double[24];
			faDS[0] = -fEtaM * fZetaM;
			faDS[1] = -faDS[0];
			faDS[2] = fEtaP * fZetaM;
			faDS[3] = -faDS[2];
			faDS[4] = -fEtaM * fZetaP;
			faDS[5] = -faDS[4];
			faDS[6] = fEtaP * fZetaP;
			faDS[7] = -faDS[6];
			faDS[8] = -fXiM * fZetaM;
			faDS[9] = -fXiP * fZetaM;
			faDS[10] = -faDS[9];
			faDS[11] = -faDS[8];
			faDS[12] = -fXiM * fZetaP;
			faDS[13] = -fXiP * fZetaP;
			faDS[14] = -faDS[13];
			faDS[15] = -faDS[12];
			faDS[16] = -fXiM * fEtaM;
			faDS[17] = -fXiP * fEtaM;
			faDS[18] = -fXiP * fEtaP;
			faDS[19] = -fXiM * fEtaP;
			faDS[20] = -faDS[16];
			faDS[21] = -faDS[17];
			faDS[22] = -faDS[18];
			faDS[23] = -faDS[19];

			return faDS;
		}

		private Tuple<double[,], double[,], double> CalcH8JDetJ(double[,] faXYZ, double[] faDS)
		{
			double[,] faJ = new double[3, 3];
			faJ[0, 0] = faDS[0] * faXYZ[0, 0] + faDS[1] * faXYZ[0, 1] + faDS[2] * faXYZ[0, 2] + faDS[3] * faXYZ[0, 3] + faDS[4] * faXYZ[0, 4] + faDS[5] * faXYZ[0, 5] + faDS[6] * faXYZ[0, 6] + faDS[7] * faXYZ[0, 7];
			faJ[0, 1] = faDS[0] * faXYZ[1, 0] + faDS[1] * faXYZ[1, 1] + faDS[2] * faXYZ[1, 2] + faDS[3] * faXYZ[1, 3] + faDS[4] * faXYZ[1, 4] + faDS[5] * faXYZ[1, 5] + faDS[6] * faXYZ[1, 6] + faDS[7] * faXYZ[1, 7];
			faJ[0, 2] = faDS[0] * faXYZ[2, 0] + faDS[1] * faXYZ[2, 1] + faDS[2] * faXYZ[2, 2] + faDS[3] * faXYZ[2, 3] + faDS[4] * faXYZ[2, 4] + faDS[5] * faXYZ[2, 5] + faDS[6] * faXYZ[2, 6] + faDS[7] * faXYZ[2, 7];
			faJ[1, 0] = faDS[8] * faXYZ[0, 0] + faDS[9] * faXYZ[0, 1] + faDS[10] * faXYZ[0, 2] + faDS[11] * faXYZ[0, 3] + faDS[12] * faXYZ[0, 4] + faDS[13] * faXYZ[0, 5] + faDS[14] * faXYZ[0, 6] + faDS[15] * faXYZ[0, 7];
			faJ[1, 1] = faDS[8] * faXYZ[1, 0] + faDS[9] * faXYZ[1, 1] + faDS[10] * faXYZ[1, 2] + faDS[11] * faXYZ[1, 3] + faDS[12] * faXYZ[1, 4] + faDS[13] * faXYZ[1, 5] + faDS[14] * faXYZ[1, 6] + faDS[15] * faXYZ[1, 7];
			faJ[1, 2] = faDS[8] * faXYZ[2, 0] + faDS[9] * faXYZ[2, 1] + faDS[10] * faXYZ[2, 2] + faDS[11] * faXYZ[2, 3] + faDS[12] * faXYZ[2, 4] + faDS[13] * faXYZ[2, 5] + faDS[14] * faXYZ[2, 6] + faDS[15] * faXYZ[2, 7];
			faJ[2, 0] = faDS[16] * faXYZ[0, 0] + faDS[17] * faXYZ[0, 1] + faDS[18] * faXYZ[0, 2] + faDS[19] * faXYZ[0, 3] + faDS[20] * faXYZ[0, 4] + faDS[21] * faXYZ[0, 5] + faDS[22] * faXYZ[0, 6] + faDS[23] * faXYZ[0, 7];
			faJ[2, 1] = faDS[16] * faXYZ[1, 0] + faDS[17] * faXYZ[1, 1] + faDS[18] * faXYZ[1, 2] + faDS[19] * faXYZ[1, 3] + faDS[20] * faXYZ[1, 4] + faDS[21] * faXYZ[1, 5] + faDS[22] * faXYZ[1, 6] + faDS[23] * faXYZ[1, 7];
			faJ[2, 2] = faDS[16] * faXYZ[2, 0] + faDS[17] * faXYZ[2, 1] + faDS[18] * faXYZ[2, 2] + faDS[19] * faXYZ[2, 3] + faDS[20] * faXYZ[2, 4] + faDS[21] * faXYZ[2, 5] + faDS[22] * faXYZ[2, 6] + faDS[23] * faXYZ[2, 7];

			double fDet1 = faJ[0, 0] * (faJ[1, 1] * faJ[2, 2] - faJ[2, 1] * faJ[1, 2]);
			double fDet2 = -faJ[0, 1] * (faJ[1, 0] * faJ[2, 2] - faJ[2, 0] * faJ[1, 2]);
			double fDet3 = faJ[0, 2] * (faJ[1, 0] * faJ[2, 1] - faJ[2, 0] * faJ[1, 1]);
			double fDetJ = fDet1 + fDet2 + fDet3;
			if (fDetJ < determinantTolerance)
				throw new ArgumentException(String.Format("Jacobian determinant is negative or under tolerance ({0} < {1}). Check the order of nodes or the element geometry.", fDetJ, determinantTolerance));

			double fDetInv = 1.0 / fDetJ;
			double[,] faJInv = new double[3, 3];
			faJInv[0, 0] = (faJ[1, 1] * faJ[2, 2] - faJ[2, 1] * faJ[1, 2]) * fDetInv;
			faJInv[1, 0] = (faJ[2, 0] * faJ[1, 2] - faJ[1, 0] * faJ[2, 2]) * fDetInv;
			faJInv[2, 0] = (faJ[1, 0] * faJ[2, 1] - faJ[2, 0] * faJ[1, 1]) * fDetInv;
			faJInv[0, 1] = (faJ[2, 1] * faJ[0, 2] - faJ[0, 1] * faJ[2, 2]) * fDetInv;
			faJInv[1, 1] = (faJ[0, 0] * faJ[2, 2] - faJ[2, 0] * faJ[0, 2]) * fDetInv;
			faJInv[2, 1] = (faJ[2, 0] * faJ[0, 1] - faJ[2, 1] * faJ[0, 0]) * fDetInv;
			faJInv[0, 2] = (faJ[0, 1] * faJ[1, 2] - faJ[1, 1] * faJ[0, 2]) * fDetInv;
			faJInv[1, 2] = (faJ[1, 0] * faJ[0, 2] - faJ[0, 0] * faJ[1, 2]) * fDetInv;
			faJInv[2, 2] = (faJ[0, 0] * faJ[1, 1] - faJ[1, 0] * faJ[0, 1]) * fDetInv;

			return new Tuple<double[,], double[,], double>(faJ, faJInv, fDetJ);
		}

		#region Alex code
		private double[,] CalculateDeformationMatrix(
			Jacobian3D jacobian, ShapeFunctionNaturalDerivatives3D[] shapeFunctionDerivatives)
		{
			double[,] jacobianInverse = jacobian.CalculateJacobianInverse();
			double[,] b = new double[6, 24];

			for (int shapeFunction = 0; shapeFunction < 8; shapeFunction++)
			{
				b[0, (3 * shapeFunction) + 0] = (jacobianInverse[0, 0] * shapeFunctionDerivatives[shapeFunction].Xi) +
												(jacobianInverse[0, 1] * shapeFunctionDerivatives[shapeFunction].Eta) +
												(jacobianInverse[0, 2] * shapeFunctionDerivatives[shapeFunction].Zeta);
				b[1, (3 * shapeFunction) + 1] = (jacobianInverse[1, 0] * shapeFunctionDerivatives[shapeFunction].Xi) +
												(jacobianInverse[1, 1] * shapeFunctionDerivatives[shapeFunction].Eta) +
												(jacobianInverse[1, 2] * shapeFunctionDerivatives[shapeFunction].Zeta);
				b[2, (3 * shapeFunction) + 2] = (jacobianInverse[2, 0] * shapeFunctionDerivatives[shapeFunction].Xi) +
												(jacobianInverse[2, 1] * shapeFunctionDerivatives[shapeFunction].Eta) +
												(jacobianInverse[2, 2] * shapeFunctionDerivatives[shapeFunction].Zeta);
				b[3, (3 * shapeFunction) + 0] = b[1, (3 * shapeFunction) + 1];
				b[3, (3 * shapeFunction) + 1] = b[0, (3 * shapeFunction) + 0];
				b[4, (3 * shapeFunction) + 1] = b[2, (3 * shapeFunction) + 2];
				b[4, (3 * shapeFunction) + 2] = b[1, (3 * shapeFunction) + 1];
				b[5, (3 * shapeFunction) + 0] = b[2, (3 * shapeFunction) + 2];
				b[5, (3 * shapeFunction) + 2] = b[0, (3 * shapeFunction) + 0];
			}

			return b;
		}

		private ShapeFunctionNaturalDerivatives3D[] CalculateShapeDerivativeValues(
			double xi, double eta, double zeta)
		{
			ShapeFunctionNaturalDerivatives3D[] shapeFunctionDerivatives =
				new ShapeFunctionNaturalDerivatives3D[8];
			for (int shapeFunction = 0; shapeFunction < 8; shapeFunction++)
			{
				shapeFunctionDerivatives[shapeFunction] = new ShapeFunctionNaturalDerivatives3D();
			}

			const double oneOverEight = 0.125;
			double xiPlus = 1.0 + xi;
			double etaPlus = 1.0 + eta;
			double zetaPlus = 1.0 + zeta;
			double xiMinus = 1.0 - xi;
			double etaMinus = 1.0 - eta;
			double zetaMinus = 1.0 - zeta;

			shapeFunctionDerivatives[0].Xi = -oneOverEight * etaMinus * zetaMinus;
			shapeFunctionDerivatives[1].Xi = -shapeFunctionDerivatives[0].Xi;
			shapeFunctionDerivatives[2].Xi = oneOverEight * etaPlus * zetaMinus;
			shapeFunctionDerivatives[3].Xi = -shapeFunctionDerivatives[2].Xi;
			shapeFunctionDerivatives[4].Xi = -oneOverEight * etaMinus * zetaPlus;
			shapeFunctionDerivatives[5].Xi = -shapeFunctionDerivatives[4].Xi;
			shapeFunctionDerivatives[6].Xi = oneOverEight * etaPlus * zetaPlus;
			shapeFunctionDerivatives[7].Xi = -shapeFunctionDerivatives[6].Xi;

			shapeFunctionDerivatives[0].Eta = -oneOverEight * xiMinus * zetaMinus;
			shapeFunctionDerivatives[1].Eta = -oneOverEight * xiPlus * zetaMinus;
			shapeFunctionDerivatives[2].Eta = -shapeFunctionDerivatives[1].Eta;
			shapeFunctionDerivatives[3].Eta = -shapeFunctionDerivatives[0].Eta;
			shapeFunctionDerivatives[4].Eta = -oneOverEight * xiMinus * zetaPlus;
			shapeFunctionDerivatives[5].Eta = -oneOverEight * xiPlus * zetaPlus;
			shapeFunctionDerivatives[6].Eta = -shapeFunctionDerivatives[5].Eta;
			shapeFunctionDerivatives[7].Eta = -shapeFunctionDerivatives[4].Eta;

			shapeFunctionDerivatives[0].Zeta = -oneOverEight * xiMinus * etaMinus;
			shapeFunctionDerivatives[1].Zeta = -oneOverEight * xiPlus * etaMinus;
			shapeFunctionDerivatives[2].Zeta = -oneOverEight * xiPlus * etaPlus;
			shapeFunctionDerivatives[3].Zeta = -oneOverEight * xiMinus * etaPlus;
			shapeFunctionDerivatives[4].Zeta = -shapeFunctionDerivatives[0].Zeta;
			shapeFunctionDerivatives[5].Zeta = -shapeFunctionDerivatives[1].Zeta;
			shapeFunctionDerivatives[6].Zeta = -shapeFunctionDerivatives[2].Zeta;
			shapeFunctionDerivatives[7].Zeta = -shapeFunctionDerivatives[3].Zeta;

			return shapeFunctionDerivatives;
		}

		private GaussLegendrePoint3D[] CalculateGaussMatrices(double[,] nodeCoordinates)
		{
			GaussLegendrePoint1D[] integrationPointsPerAxis =
				GaussQuadrature.GetGaussLegendrePoints(iInt);
			int totalSamplingPoints = (int)Math.Pow(integrationPointsPerAxis.Length, 3);

			GaussLegendrePoint3D[] integrationPoints = new GaussLegendrePoint3D[totalSamplingPoints];

			int counter = -1;
			foreach (GaussLegendrePoint1D pointXi in integrationPointsPerAxis)
			{
				foreach (GaussLegendrePoint1D pointEta in integrationPointsPerAxis)
				{
					foreach (GaussLegendrePoint1D pointZeta in integrationPointsPerAxis)
					{
						counter += 1;
						double xi = pointXi.Coordinate;
						double eta = pointEta.Coordinate;
						double zeta = pointZeta.Coordinate;

						ShapeFunctionNaturalDerivatives3D[] shapeDerivativeValues =
							this.CalculateShapeDerivativeValues(xi, eta, zeta);
						Jacobian3D jacobian = new Jacobian3D(nodeCoordinates, shapeDerivativeValues);
						double[,] deformationMatrix = this.CalculateDeformationMatrix(jacobian, shapeDerivativeValues);
						double weightFactor = pointXi.WeightFactor * pointEta.WeightFactor * pointZeta.WeightFactor *
											  jacobian.Determinant;

						integrationPoints[counter] = new GaussLegendrePoint3D(
							xi, eta, zeta, deformationMatrix, weightFactor);
					}
				}
			}

			return integrationPoints;
		}

		public virtual IMatrix StiffnessMatrix()
		{
			double[,] coordinates = this.GetCoordinates();
			GaussLegendrePoint3D[] integrationPoints = this.CalculateGaussMatrices(coordinates);

			var stiffnessMatrix = SymmetricMatrix.CreateZero(24);

			int pointId = -1;
			foreach (GaussLegendrePoint3D intPoint in integrationPoints)
			{
				pointId++;
				IMatrixView constitutiveMatrix = materialsAtGaussPoints[pointId].ConstitutiveMatrix;
				double[,] b = intPoint.DeformationMatrix;
				for (int i = 0; i < 24; i++)
				{
					double[] eb = new double[24];
					for (int iE = 0; iE < 6; iE++)
					{
						eb[iE] = (constitutiveMatrix[iE, 0] * b[0, i]) + (constitutiveMatrix[iE, 1] * b[1, i]) +
								 (constitutiveMatrix[iE, 2] * b[2, i]) + (constitutiveMatrix[iE, 3] * b[3, i]) +
								 (constitutiveMatrix[iE, 4] * b[4, i]) + (constitutiveMatrix[iE, 5] * b[5, i]);
					}

					for (int j = i; j < 24; j++)
					{
						double stiffness = (b[0, j] * eb[0]) + (b[1, j] * eb[1]) + (b[2, j] * eb[2]) + (b[3, j] * eb[3]) +
										   (b[4, j] * eb[4]) + (b[5, j] * eb[5]);
						stiffnessMatrix[i, j] += stiffness * intPoint.WeightFactor;
					}
				}
			}

			return stiffnessMatrix;
		}
		public IMatrix PhysicsMatrix()
		{
			return StiffnessMatrix();
		}
		public IMatrix CalculateConsistentMass()
		{
			double[,] coordinates = GetCoordinates();
			GaussLegendrePoint3D[] integrationPoints = this.CalculateGaussMatrices(coordinates);

			var consistentMass = SymmetricMatrix.CreateZero(24);

			foreach (GaussLegendrePoint3D intPoint in integrationPoints)
			{
				double[] shapeFunctionValues = this.CalcH8Shape(intPoint.Xi, intPoint.Eta, intPoint.Zeta);
				double weightDensity = intPoint.WeightFactor * this.Density;
				for (int shapeFunctionI = 0; shapeFunctionI < shapeFunctionValues.Length; shapeFunctionI++)
				{
					for (int shapeFunctionJ = shapeFunctionI; shapeFunctionJ < shapeFunctionValues.Length; shapeFunctionJ++)
					{
						consistentMass[3 * shapeFunctionI, 3 * shapeFunctionJ] += shapeFunctionValues[shapeFunctionI] *
																				  shapeFunctionValues[shapeFunctionJ] *
																				  weightDensity;
					}

					for (int shapeFunctionJ = shapeFunctionI; shapeFunctionJ < shapeFunctionValues.Length; shapeFunctionJ++)
					{
						consistentMass[(3 * shapeFunctionI) + 1, (3 * shapeFunctionJ) + 1] =
							consistentMass[3 * shapeFunctionI, 3 * shapeFunctionJ];

						consistentMass[(3 * shapeFunctionI) + 2, (3 * shapeFunctionJ) + 2] =
							consistentMass[3 * shapeFunctionI, 3 * shapeFunctionJ];
					}
				}
			}

			return consistentMass;
		}

		#endregion

		public virtual IMatrix MassMatrix() => CalculateConsistentMass();


		public virtual IMatrix DampingMatrix()
		{
			IMatrix k = StiffnessMatrix();
			IMatrix m = MassMatrix();
			k.LinearCombinationIntoThis(RayleighBeta, m, RayleighAlpha);
			return k;
		}

		// TODO: Replaced dstrains with strains (goat: this might go wrong)
		public Tuple<double[], double[]> CalculateResponse(double[] localDisplacements)
		{
			double[,] nodalCoordinates = GetCoordinates();
			GaussLegendrePoint3D[] gaussMatrices = CalculateGaussMatrices(nodalCoordinates);
			int gaussPointsCount = iInt3;
			if (gaussMatrices.Count() != gaussPointsCount)
				throw new Exception("There must have been " + gaussPointsCount + " gauss points, but there are "
					+ gaussMatrices.Count());
			if (materialsAtGaussPoints.Length != gaussPointsCount)
				throw new Exception("There must have been " + gaussPointsCount + " material points, but there are "
				   + materialsAtGaussPoints.Length);

			double[] deltaStrains = new double[6];
			double[] strains = new double[6];
			for (int gp = 0; gp < gaussPointsCount; ++gp)
			{
				var deformationMatrix = Matrix.CreateZero(6, 24);
				for (int i = 0; i < 6; i++)
				{
					for (int j = 0; j < 24; j++)
					{
						deformationMatrix[i, j] = gaussMatrices[gp].DeformationMatrix[i, j];
					}
				}
				deformationMatrix.MultiplyIntoResult(localDisplacements, strains);
				//deformationMatrix.MultiplyIntoResult(localdDisplacements, deltaStrains);
				//materialsAtGaussPoints[gp].UpdateConstitutiveMatrixAndEvaluateResponse(deltaStrains);
				lastStresses[gp] = materialsAtGaussPoints[gp].UpdateConstitutiveMatrixAndEvaluateResponse(strains);
			}
			return new Tuple<double[], double[]>(strains, lastStresses[materialsAtGaussPoints.Length - 1]);
		}

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
		{
			return CalculateResponseIntegral();
		}


		public double[] CalculateResponseIntegral()
		{
			double[,] nodalCoordinates = GetCoordinates();
			GaussLegendrePoint3D[] gaussMatrices = CalculateGaussMatrices(nodalCoordinates);
			int gaussPointsCount = iInt3;
			if (gaussMatrices.Count() != gaussPointsCount)
				throw new Exception("There must have been " + gaussPointsCount + " gauss points, but there are "
					+ gaussMatrices.Count());
			if (materialsAtGaussPoints.Length != gaussPointsCount)
				throw new Exception("There must have been " + gaussPointsCount + " material points, but there are "
				   + materialsAtGaussPoints.Length);

			int dofsCount = 24; // 8 nodes * 3 dofs per node
			var internalForces = new double[dofsCount];
			for (int gp = 0; gp < gaussPointsCount; ++gp)
			{
				var stressVector = lastStresses[gp]; // Risky if stressVector is modified
																		//var stressesClone = new double[tensorComponentsCount];
																		//Array.Copy(materialsAtGaussPoints[gp].Stresses, stressesClone, tensorComponentsCount);
																		//Vector stressVector = new Vector(stressesClone);
				Matrix deformationMatrix = Matrix.CreateFromArray(gaussMatrices[gp].DeformationMatrix);
				double[] transposeBtimesStress = deformationMatrix.Multiply(stressVector, true);
				internalForces.AxpyIntoThis(transposeBtimesStress, gaussMatrices[gp].WeightFactor);
			}
			return internalForces;
		}

		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads)
		//{
		//	var accelerations = new double[24];
		//	IMatrix massMatrix = MassMatrix(element);

		//	foreach (MassAccelerationLoad load in loads)
		//	{
		//		int index = 0;
		//		foreach (IDofType[] nodalDOFTypes in dofTypes)
		//			foreach (IDofType dofType in nodalDOFTypes)
		//			{
		//				if (dofType == load.DOF) accelerations[index] += load.Amount;
		//				index++;
		//			}
		//	}

		//	return massMatrix.Multiply(accelerations);
		//}

		//public void ClearMaterialState()
		//{
		//	foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearState();
		//}

		public void SaveConstitutiveLawState(IHaveState externalState)
		{
			foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.CreateState();
		}

		//public void ClearMaterialStresses()
		//{
		//	foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
		//}

		#endregion

		#region IStructuralFiniteElement Members

		//public bool ConstitutiveLawModified
		//{
		//	get
		//	{
		//		foreach (IContinuumMaterial3D material in materialsAtGaussPoints)
		//			if (material.IsCurrentStateDifferent()) return true;
		//		return false;
		//	}
		//}

		//public void ResetConstitutiveLawModified()
		//{
		//	foreach (IContinuumMaterial3D material in materialsAtGaussPoints) material.ResetModified();
		//}

		#endregion

		#region IEmbeddedHostElement Members

		private readonly IList<INode> embeddedNodes = new List<INode>();
		public EmbeddedNode BuildHostElementEmbeddedNode(IElementType element, INode node,
			IEmbeddedDOFInHostTransformationVector transformationVector)
		{
			var points = GetNaturalCoordinates(element, node);
			if (points.Length == 0) return null;

			embeddedNodes.Add(node);
			var embeddedNode = new EmbeddedNode(node, element, transformationVector.GetDependentDOFTypes);
			for (int i = 0; i < points.Length; i++)
				embeddedNode.Coordinates.Add(points[i]);
			return embeddedNode;
		}

		public double[] GetShapeFunctionsForNode(IElementType element, EmbeddedNode node)
		{
			double[,] elementCoordinates = GetCoordinatesTranspose(element);
			var shapeFunctions = CalcH8Shape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
			var nablaShapeFunctions = CalcH8NablaShape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
			var jacobian = CalcH8JDetJ(elementCoordinates, nablaShapeFunctions);

			return new double[]
			{
				shapeFunctions[0], shapeFunctions[1], shapeFunctions[2], shapeFunctions[3], shapeFunctions[4], shapeFunctions[5], shapeFunctions[6], shapeFunctions[7],
				nablaShapeFunctions[0], nablaShapeFunctions[1], nablaShapeFunctions[2], nablaShapeFunctions[3], nablaShapeFunctions[4], nablaShapeFunctions[5], nablaShapeFunctions[6], nablaShapeFunctions[7],
				nablaShapeFunctions[8], nablaShapeFunctions[9], nablaShapeFunctions[10], nablaShapeFunctions[11], nablaShapeFunctions[12], nablaShapeFunctions[13], nablaShapeFunctions[14], nablaShapeFunctions[15],
				nablaShapeFunctions[16], nablaShapeFunctions[17], nablaShapeFunctions[18], nablaShapeFunctions[19], nablaShapeFunctions[20], nablaShapeFunctions[21], nablaShapeFunctions[22], nablaShapeFunctions[23],
				jacobian.Item1[0, 0], jacobian.Item1[0, 1], jacobian.Item1[0, 2], jacobian.Item1[1, 0], jacobian.Item1[1, 1], jacobian.Item1[1, 2], jacobian.Item1[2, 0], jacobian.Item1[2, 1], jacobian.Item1[2, 2],
				jacobian.Item2[0, 0], jacobian.Item2[0, 1], jacobian.Item2[0, 2], jacobian.Item2[1, 0], jacobian.Item2[1, 1], jacobian.Item2[1, 2], jacobian.Item2[2, 0], jacobian.Item2[2, 1], jacobian.Item2[2, 2]
			};

			//double[] coords = GetNaturalCoordinates(element, node);
			//return CalcH8Shape(coords[0], coords[1], coords[2]);

			//double fXiP = (1.0 + coords[0]) * 0.5;
			//double fEtaP = (1.0 + coords[1]) * 0.5;
			//double fZetaP = (1.0 + coords[2]) * 0.5;
			//double fXiM = (1.0 - coords[0]) * 0.5;
			//double fEtaM = (1.0 - coords[1]) * 0.5;
			//double fZetaM = (1.0 - coords[2]) * 0.5;

			//return new double[] { fXiM * fEtaM * fZetaM,
			//    fXiP * fEtaM * fZetaM,
			//    fXiP * fEtaP * fZetaM,
			//    fXiM * fEtaP * fZetaM,
			//    fXiM * fEtaM * fZetaP,
			//    fXiP * fEtaM * fZetaP,
			//    fXiP * fEtaP * fZetaP,
			//    fXiM * fEtaP * fZetaP };
		}

		private double[] GetNaturalCoordinates(IElementType element, INode node)
		{
			double[] mins = new double[] { element.Nodes[0].X, element.Nodes[0].Y, element.Nodes[0].Z };
			double[] maxes = new double[] { element.Nodes[0].X, element.Nodes[0].Y, element.Nodes[0].Z };
			for (int i = 0; i < element.Nodes.Count; i++)
			{
				mins[0] = mins[0] > element.Nodes[i].X ? element.Nodes[i].X : mins[0];
				mins[1] = mins[1] > element.Nodes[i].Y ? element.Nodes[i].Y : mins[1];
				mins[2] = mins[2] > element.Nodes[i].Z ? element.Nodes[i].Z : mins[2];
				maxes[0] = maxes[0] < element.Nodes[i].X ? element.Nodes[i].X : maxes[0];
				maxes[1] = maxes[1] < element.Nodes[i].Y ? element.Nodes[i].Y : maxes[1];
				maxes[2] = maxes[2] < element.Nodes[i].Z ? element.Nodes[i].Z : maxes[2];
			}
			//return new double[] { (node.X - mins[0]) / ((maxes[0] - mins[0]) / 2) - 1,
			//    (node.Y - mins[1]) / ((maxes[1] - mins[1]) / 2) - 1,
			//    (node.Z - mins[2]) / ((maxes[2] - mins[2]) / 2) - 1 };

			bool maybeInsideElement = node.X <= maxes[0] && node.X >= mins[0] &&
				node.Y <= maxes[1] && node.Y >= mins[1] &&
				node.Z <= maxes[2] && node.Z >= mins[2];
			if (maybeInsideElement == false) return new double[0];

			const int jacobianSize = 3;
			const int maxIterations = 1000;
			const double tolerance = 1e-10;
			int iterations = 0;
			double deltaNaturalCoordinatesNormSquare = 100;
			double[] naturalCoordinates = new double[] { 0, 0, 0 };
			const double toleranceSquare = tolerance * tolerance;

			while (deltaNaturalCoordinatesNormSquare > toleranceSquare && iterations < maxIterations)
			{
				iterations++;
				var shapeFunctions = CalcH8Shape(naturalCoordinates[0], naturalCoordinates[1], naturalCoordinates[2]);
				double[] coordinateDifferences = new double[] { 0, 0, 0 };
				for (int i = 0; i < shapeFunctions.Length; i++)
				{
					coordinateDifferences[0] += shapeFunctions[i] * element.Nodes[i].X;
					coordinateDifferences[1] += shapeFunctions[i] * element.Nodes[i].Y;
					coordinateDifferences[2] += shapeFunctions[i] * element.Nodes[i].Z;
				}
				coordinateDifferences[0] = node.X - coordinateDifferences[0];
				coordinateDifferences[1] = node.Y - coordinateDifferences[1];
				coordinateDifferences[2] = node.Z - coordinateDifferences[2];

				double[,] faXYZ = GetCoordinatesTranspose(element);
				double[] nablaShapeFunctions = CalcH8NablaShape(naturalCoordinates[0], naturalCoordinates[1], naturalCoordinates[2]);
				var inverseJacobian = CalcH8JDetJ(faXYZ, nablaShapeFunctions).Item2;

				double[] deltaNaturalCoordinates = new double[] { 0, 0, 0 };
				for (int i = 0; i < jacobianSize; i++)
					for (int j = 0; j < jacobianSize; j++)
						deltaNaturalCoordinates[i] += inverseJacobian[j, i] * coordinateDifferences[j];
				for (int i = 0; i < 3; i++)
					naturalCoordinates[i] += deltaNaturalCoordinates[i];

				deltaNaturalCoordinatesNormSquare = 0;
				for (int i = 0; i < 3; i++)
					deltaNaturalCoordinatesNormSquare += deltaNaturalCoordinates[i] * deltaNaturalCoordinates[i];
				//deltaNaturalCoordinatesNormSquare = Math.Sqrt(deltaNaturalCoordinatesNormSquare);
			}

			return naturalCoordinates.Count(x => Math.Abs(x) - 1.0 > tolerance) > 0 ? new double[0] : naturalCoordinates;
		}

		#endregion
	}
}
