using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.DataStructures;

namespace MGroup.FEM.Structural.Line
{
	public class EulerBeam2D : IStructuralElementType
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.RotationZ };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		private readonly double youngModulus;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		public double Density { get; set; }
		public double SectionArea { get; set; }
		public double MomentOfInertia { get; set; }
		public double RayleighAlpha { get; set; }
		public double RayleighBeta { get; set; }

		public EulerBeam2D(IReadOnlyList<INode> nodes, double youngModulus)
		{
			this.youngModulus = youngModulus;
			this.Nodes = nodes;
		}

		public EulerBeam2D(IReadOnlyList<INode> nodes, double youngModulus, IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus)
		{
			this.dofEnumerator = dofEnumerator;
		}

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public CellType CellType { get; } = CellType.Line2;

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofs;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		//[  c^2*E*A/L+12*s^2*E*I/L^3,  s*E*A/L*c-12*c*E*I/L^3*s,              -6*E*I/L^2*s, -c^2*E*A/L-12*s^2*E*I/L^3, -s*E*A/L*c+12*c*E*I/L^3*s,              -6*E*I/L^2*s]
		//[  s*E*A/L*c-12*c*E*I/L^3*s,  s^2*E*A/L+12*c^2*E*I/L^3,               6*E*I/L^2*c, -s*E*A/L*c+12*c*E*I/L^3*s, -s^2*E*A/L-12*c^2*E*I/L^3,               6*E*I/L^2*c]
		//[              -6*E*I/L^2*s,               6*E*I/L^2*c,                   4*E*I/L,               6*E*I/L^2*s,              -6*E*I/L^2*c,                   2*E*I/L]
		//[ -c^2*E*A/L-12*s^2*E*I/L^3, -s*E*A/L*c+12*c*E*I/L^3*s,               6*E*I/L^2*s,  c^2*E*A/L+12*s^2*E*I/L^3,  s*E*A/L*c-12*c*E*I/L^3*s,               6*E*I/L^2*s]
		//[ -s*E*A/L*c+12*c*E*I/L^3*s, -s^2*E*A/L-12*c^2*E*I/L^3,              -6*E*I/L^2*c,  s*E*A/L*c-12*c*E*I/L^3*s,  s^2*E*A/L+12*c^2*E*I/L^3,              -6*E*I/L^2*c]
		//[              -6*E*I/L^2*s,               6*E*I/L^2*c,                   2*E*I/L,               6*E*I/L^2*s,              -6*E*I/L^2*c,                   4*E*I/L]
		public virtual IMatrix StiffnessMatrix()
		{
			double x2 = Math.Pow(Nodes[1].X - Nodes[0].X, 2);
			double y2 = Math.Pow(Nodes[1].Y - Nodes[0].Y, 2);
			double L = Math.Sqrt(x2 + y2);
			double c = (Nodes[1].X - Nodes[0].X) / L;
			double c2 = c * c;
			double s = (Nodes[1].Y - Nodes[0].Y) / L;
			double s2 = s * s;
			double EL = this.youngModulus / L;
			double EAL = EL * SectionArea;
			double EIL = EL * MomentOfInertia;
			double EIL2 = EIL / L;
			double EIL3 = EIL2 / L;

			//TODO: optimize this
			int order = 6;
			var k = SymmetricMatrix.CreateFromPackedRowMajorArray(new double[]
			{
				c2*EAL+12*s2*EIL3, c*s*EAL-12*c*s*EIL3, -6*s*EIL2, -c2*EAL-12*s2*EIL3, -c*s*EAL+12*c*s*EIL3, -6*s*EIL2,
				s2*EAL+12*c2*EIL3, 6*c*EIL2, -s*c*EAL+12*c*s*EIL3, -s2*EAL-12*c2*EIL3, 6*c*EIL2,
				4*EIL, 6*s*EIL2, -6*c*EIL2, 2*EIL,
				c2*EAL+12*s2*EIL3, s*c*EAL-12*c*s*EIL3, 6*s*EIL2,
				s2*EAL+12*c2*EIL3, -6*c*EIL2,
				4*EIL
			}, order);

			return dofEnumerator.GetTransformedMatrix(k);
		}
		public IMatrix PhysicsMatrix()
		{
			return StiffnessMatrix();
		}
		////[ 140*c^2+156*s^2,         -16*c*s,         -22*s*L,   70*c^2+54*s^2,          16*c*s,          13*s*L]
		////[         -16*c*s, 140*s^2+156*c^2,          22*c*L,          16*c*s,   70*s^2+54*c^2,         -13*c*L]
		////[         -22*s*L,          22*c*L,           4*L^2,         -13*s*L,          13*c*L,          -3*L^2]
		////[   70*c^2+54*s^2,          16*c*s,         -13*s*L, 140*c^2+156*s^2,         -16*c*s,          22*s*L]
		////[          16*c*s,   70*s^2+54*c^2,          13*c*L,         -16*c*s, 140*s^2+156*c^2,         -22*c*L]
		////[          13*s*L,         -13*c*L,          -3*L^2,          22*s*L,         -22*c*L,           4*L^2]

		//public IMatrix2D<double> MassMatrix(Element element)
		//{
		//    double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
		//    double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
		//    double L = Math.Sqrt(x2 + y2);
		//    double L2 = L * L;
		//    double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
		//    double c2 = c * c;
		//    double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
		//    double s2 = s * s;
		//    double dAL420 = Density * SectionArea * L / 420;
		//    return new SymmetricMatrix2D<double>(new double[] { dAL420*(140*c2+156*s2), -16*dAL420*c*s, -22*dAL420*s*L, dAL420*(70*c2+54*s2), 16*dAL420*c*s, 13*dAL420*s*L,
		//        dAL420*(140*s2+156*c2), 22*dAL420*c*L, 16*dAL420*c*s, dAL420*(70*s2+54*c2), -13*dAL420*c*L,
		//        4*dAL420*L2, -13*dAL420*s*L, 13*dAL420*c*L, -3*dAL420*L2,
		//        dAL420*(140*c2+156*s2), -16*dAL420*c*s, 22*dAL420*s*L,
		//        dAL420*(140*s2+156*c2), -22*dAL420*c*L,
		//        4*dAL420*L2 });
		//}

		//[ 140*c^2+156*s^2,         -16*c*s,         -22*s*L,   70*c^2+54*s^2,          16*c*s,          13*s*L]
		//[         -16*c*s, 140*s^2+156*c^2,          22*c*L,          16*c*s,   70*s^2+54*c^2,         -13*c*L]
		//[         -22*s*L,          22*c*L,           4*L^2,         -13*s*L,          13*c*L,          -3*L^2]
		//[   70*c^2+54*s^2,          16*c*s,         -13*s*L, 140*c^2+156*s^2,         -16*c*s,          22*s*L]
		//[          16*c*s,   70*s^2+54*c^2,          13*c*L,         -16*c*s, 140*s^2+156*c^2,         -22*c*L]
		//[          13*s*L,         -13*c*L,          -3*L^2,          22*s*L,         -22*c*L,           4*L^2]
		public IMatrix MassMatrix()
		{
			double x2 = Math.Pow(Nodes[1].X - Nodes[0].X, 2);
			double y2 = Math.Pow(Nodes[1].Y - Nodes[0].Y, 2);
			double L = Math.Sqrt(x2 + y2);
			double L2 = L * L;
			double c = (Nodes[1].X - Nodes[0].X) / L;
			double c2 = c * c;
			double s = (Nodes[1].Y - Nodes[0].Y) / L;
			double s2 = s * s;
			double dAL420 = Density * SectionArea * L / 420;

			double totalMass = Density * SectionArea * L;
			double totalMassOfDiagonalTerms = 2 * dAL420 * (140 * c2 + 156 * s2) + 2 * dAL420 * (140 * s2 + 156 * c2);
			double scale = totalMass / totalMassOfDiagonalTerms;

			//TODO: optimize this
			int order = 6;
			return SymmetricMatrix.CreateFromPackedRowMajorArray(new double[]
			{
				dAL420 *(140*c2+156*s2)*scale, 0, 0, 0, 0, 0,
				dAL420*(140*s2+156*c2)*scale, 0, 0, 0, 0,
				0, 0, 0, 0,
				dAL420*(140*c2+156*s2)*scale, 0, 0,
				dAL420*(140*s2+156*c2)*scale, 0,
				0
			}, order);
		}

		public IMatrix DampingMatrix()
		{
			var k = StiffnessMatrix();
			var m = MassMatrix();
			k.LinearCombinationIntoThis(RayleighBeta, m, RayleighAlpha);
			return k;
		}

		public Tuple<double[], double[]> CalculateResponse(double[] localDisplacements)
			=> throw new NotImplementedException();

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
			=> CalculateResponseIntegral();

		public double[] CalculateResponseIntegral()
			=> throw new NotImplementedException();

		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads)
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

		public void SaveConstitutiveLawState(IHaveState externalState) => throw new NotImplementedException();

		#endregion

		#region IFiniteElement Members


		public bool ConstitutiveLawModified => false;

		public void ResetConstitutiveLawModified()
		{
			// Method intentionally left empty.
		}

		public void ClearConstitutiveLawState()
		{
			// Method intentionally left empty.
		}

		public void ClearConstitutiveLawStresses() => throw new NotImplementedException();
		#endregion
	}
}
