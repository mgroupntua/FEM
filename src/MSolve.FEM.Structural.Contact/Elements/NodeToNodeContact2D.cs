using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;

using MSolve.FEM.Structural.Contact;
using MSolve.FEM.Structural.Contact.Helpers;

namespace MGroup.FEM.Structural.Contact.Elements
{
	public class NodeToNodeContact2D : IStructuralFiniteElement
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		private readonly double youngModulus;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		public double Density { get; set; }

		public double SectionArea { get; set; }

		public double MomentOfInertia { get; set; }

		public double RayleighAlpha { get; set; }

		public double RayleighBeta { get; set; }

		#region contactProperties
		public double[] DisplacementVector { get; set; }

		private double PenaltyFactor { get; set; }

		#endregion contactProperties

		public NodeToNodeContact2D(double youngModulus, double A, double I)
		{
			this.youngModulus = youngModulus;
			#region contact
			SectionArea = A;
			MomentOfInertia = I;
			//PenaltyFactor = youngModulus * 100.0;
			PenaltyFactor = youngModulus * 2.0 / 100;
			DisplacementVector = new double[4];
			#endregion contact

		}

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		#region IElementType Members

		public int ID => 1;

		public CellType CellType { get; } = CellType.Line;

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofs;

		public IList<Node> GetNodesForMatrixAssembly(Element element) => element.Nodes;

		public virtual IMatrix StiffnessMatrix(IElement element)
		{
			#region backup
			//double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
			//double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
			//double L = Math.Sqrt(x2 + y2);
			//double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
			//double c2 = c * c;
			//double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
			//double s2 = s * s;
			//double EL = this.youngModulus / L;
			//double EAL = EL * SectionArea;
			//double EIL = EL * MomentOfInertia;
			//double EIL2 = EIL / L;
			//double EIL3 = EIL2 / L;

			////TODO: optimize this
			//int order = 6;
			//var k = SymmetricMatrix.CreateFromPackedRowMajorArray(new double[]
			//{
			//	c2*EAL+12*s2*EIL3, c*s*EAL-12*c*s*EIL3, -6*s*EIL2, -c2*EAL-12*s2*EIL3, -c*s*EAL+12*c*s*EIL3, -6*s*EIL2,
			//	s2*EAL+12*c2*EIL3, 6*c*EIL2, -s*c*EAL+12*c*s*EIL3, -s2*EAL-12*c2*EIL3, 6*c*EIL2,
			//	4*EIL, 6*s*EIL2, -6*c*EIL2, 2*EIL,
			//	c2*EAL+12*s2*EIL3, s*c*EAL-12*c*s*EIL3, 6*s*EIL2,
			//	s2*EAL+12*c2*EIL3, -6*c*EIL2,
			//	4*EIL
			//}, order);

			//return dofEnumerator.GetTransformedMatrix(k);
			#endregion backup
			double penetration = CalculateNormalGap(element);
			if (penetration <= 0)
			{
				double[] n = CalculateNormalUnitVector(element);
				double[,] A = CalculatePositionMatrix();
				double[,] AT = MatrixOperations.Transpose(A);
				double[,] nxn = VectorOperations.VectorVectorTensorProduct(n, n);
				double[,] nxn_A = MatrixOperations.MatrixProduct(nxn, A);
				double[,] AT_nxn_A = MatrixOperations.MatrixProduct(AT, nxn_A);
				double[,] globalStiffnessMatrix = MatrixOperations.ScalarMatrixProductNew(PenaltyFactor, AT_nxn_A);
				return SymmetricMatrix.CreateFromArray(globalStiffnessMatrix);
			}
			else
			{
				double[,] globalStifnessMatrix = new double[4, 4];
				return SymmetricMatrix.CreateFromArray(globalStifnessMatrix);
			}
		}

		#region contactFunctions
		private double CalculateNormalGap(IElement element)
		{
			double[,] A = CalculatePositionMatrix();
			double[,] AT = MatrixOperations.Transpose(A);
			double[] n = CalculateNormalUnitVector(element);
			double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
			double[] xupd = new double[] {
				element.Nodes[0].X + DisplacementVector[0],
				element.Nodes[0].Y + DisplacementVector[1],
				element.Nodes[1].X + DisplacementVector[2],
				element.Nodes[1].Y + DisplacementVector[3]
			};
			double normalGap = VectorOperations.VectorDotProduct(xupd, AT_n);
			return normalGap;
		}

		private double[] CalculateNormalUnitVector(IElement element)
		{
			double X1 = element.Nodes[0].X;
			double Y1 = element.Nodes[0].Y;
			double X2 = element.Nodes[1].X;
			double Y2 = element.Nodes[1].Y;
			double[] normalVector = new double[] { X2 - X1, Y2 - Y1 };
			double normalVectorLength = VectorOperations.VectorNorm2(normalVector);
			double[] normalUnitVec = new double[] { normalVector[0] / normalVectorLength, normalVector[1] / normalVectorLength };
			return normalUnitVec;
		}

		private double[,] CalculatePositionMatrix()
		{
			double[,] aMatrix = new double[,]
				{
					{ -1,0,1,0},
					{0,-1,0,1 }
				};
			return aMatrix;
		}
		#endregion contactFunctions


		public IMatrix MassMatrix(IElement element)
		{
			#region backup
			//double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
			//double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
			//double L = Math.Sqrt(x2 + y2);
			//double L2 = L * L;
			//double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
			//double c2 = c * c;
			//double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
			//double s2 = s * s;
			//double dAL420 = Density * SectionArea * L / 420;

			//double totalMass = Density * SectionArea * L;
			//double totalMassOfDiagonalTerms = 2 * dAL420 * (140 * c2 + 156 * s2) + 2 * dAL420 * (140 * s2 + 156 * c2);
			//double scale = totalMass / totalMassOfDiagonalTerms;

			////TODO: optimize this
			//int order = 6;
			//return SymmetricMatrix.CreateFromPackedRowMajorArray(new double[]
			//{
			//	dAL420 *(140*c2+156*s2)*scale, 0, 0, 0, 0, 0,
			//	dAL420*(140*s2+156*c2)*scale, 0, 0, 0, 0,
			//	0, 0, 0, 0,
			//	dAL420*(140*c2+156*s2)*scale, 0, 0,
			//	dAL420*(140*s2+156*c2)*scale, 0,
			//	0
			//}, order);
			#endregion backup
			return SymmetricMatrix.CreateFromArray(new double[4, 4]);
		}

		public IMatrix DampingMatrix(IElement element)
		{
			#region backup
			var k = StiffnessMatrix(element);
			var m = MassMatrix(element);
			k.LinearCombinationIntoThis(RayleighBeta, m, RayleighAlpha);
			return k;
			#endregion backup
			return SymmetricMatrix.CreateFromArray(new double[4, 4]);
		}

		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements, double[] localdDisplacements)
			=> new Tuple<double[], double[]>(new double[4], new double[4]);

		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
			=> CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			DisplacementVector[0] = localDisplacements[0];
			DisplacementVector[1] = localDisplacements[1];
			DisplacementVector[2] = localDisplacements[2];
			DisplacementVector[3] = localDisplacements[3];

			double penetration = CalculateNormalGap(element);
			if (penetration <= 0)
			{
				double[,] A = CalculatePositionMatrix();
				double[,] AT = MatrixOperations.Transpose(A);
				double[] n = CalculateNormalUnitVector(element);
				double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
				double ksi = CalculateNormalGap(element);
				double[] ksi_AT_n = VectorOperations.VectorScalarProductNew(AT_n, ksi);
				double[] e_ksi_AT_n = VectorOperations.VectorScalarProductNew(ksi_AT_n, PenaltyFactor);
				return e_ksi_AT_n;
			}
			else
			{
				double[] internalGlobalForcesVector = new double[4];
				return internalGlobalForcesVector;
			}
		}

		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
			=> throw new NotImplementedException();
		#region backup
		//{
		//	var accelerations = new double[6];
		//	int index = 0;
		//	foreach (MassAccelerationLoad load in loads)
		//		foreach (IDofType[] nodalDOFTypes in dofs)
		//			foreach (IDofType dofType in nodalDOFTypes)
		//			{
		//				if (dofType == load.DOF) accelerations[index] += load.Amount;
		//				index++;
		//			}

		//	IMatrix massMatrix = MassMatrix(element);
		//	return massMatrix.Multiply(accelerations);
		//}
		#endregion backup
		public void SaveMaterialState()
		{
			// Method intentionally left empty.
		}

		#endregion

		#region IFiniteElement Members


		public bool MaterialModified => false;

		public void ResetMaterialModified()
		{
			// Method intentionally left empty.
		}

		public void ClearMaterialState()
		{
			DisplacementVector = new double[4];
		}

		public void ClearMaterialStresses() => throw new NotImplementedException();
		#endregion
	}
}
