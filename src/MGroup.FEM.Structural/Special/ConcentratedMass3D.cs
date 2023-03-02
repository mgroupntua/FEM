using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.DataStructures;

namespace MGroup.FEM.Structural.Special
{
	public class ConcentratedMass3D : IStructuralElementType
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes };
		private readonly double massCoefficient;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public CellType CellType { get; } = CellType.Unknown;

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes()
		{
			var d = new List<List<IDofType>>();
			foreach (var node in Nodes)
			{
				var nodeDofs = new List<IDofType>();
				nodeDofs.AddRange(nodalDOFTypes);
				d.Add(nodeDofs);
			}
			return d;
		}

		public bool ConstitutiveLawModified => false;

		public ConcentratedMass3D(double massCoefficient)
		{
			this.massCoefficient = massCoefficient;
		}

		public ConcentratedMass3D(double massCoefficient, IElementDofEnumerator dofEnumerator)
			: this(massCoefficient)
		{
			this.dofEnumerator = dofEnumerator;
		}

		public IMatrix MassMatrix()
		{
			var mass = Matrix.CreateZero(3, 3);
			mass[0, 0] = massCoefficient;
			mass[1, 1] = massCoefficient;
			mass[2, 2] = massCoefficient;
			return mass;
		}

		public IMatrix StiffnessMatrix() => Matrix.CreateZero(3, 3);

		public IMatrix PhysicsMatrix()
		{
			return StiffnessMatrix();
		}

		public IMatrix DampingMatrix() => Matrix.CreateZero(3, 3);

		public void ResetConstitutiveLawModified() { }

		public Tuple<double[], double[]> CalculateResponse(double[] localDisplacements) 
			=> new Tuple<double[], double[]>(new double[6], new double[6]);

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
			=> CalculateResponseIntegral();

		public double[] CalculateResponseIntegral()
			=> new double[6];

		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads)
		//{
		//	var accelerations = new double[3];
		//	IMatrix massMatrix = MassMatrix(element);

		//	foreach (MassAccelerationLoad load in loads)
		//	{
		//		int index = 0;
		//		foreach (IDofType[] nodalDOFTypes in dofs)
		//			foreach (IDofType dofType in nodalDOFTypes)
		//			{
		//				if (dofType == load.DOF) accelerations[index] += load.Amount;
		//				index++;
		//			}
		//	}

		//	return massMatrix.Multiply(accelerations);
		//}

		public void ClearConstitutiveLawState() { }
		public void SaveConstitutiveLawState(IHaveState externalState) { }
		public void ClearConstitutiveLawStresses() { }
	}
}
