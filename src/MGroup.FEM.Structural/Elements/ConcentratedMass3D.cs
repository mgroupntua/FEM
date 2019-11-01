using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.FEM.Interfaces;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Structural.Elements
{
	public class ConcentratedMass3D : IStructuralFiniteElement
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes };
		private readonly double massCoefficient;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		public int ID => 998;
		public CellType CellType { get; } = CellType.Unknown;

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			if (element == null) return dofs;

			var d = new List<List<IDofType>>();
			foreach (var node in element.Nodes)
			{
				var nodeDofs = new List<IDofType>();
				nodeDofs.AddRange(nodalDOFTypes);
				d.Add(nodeDofs);
			}
			return d;
		}

		public bool MaterialModified => false;

		public ConcentratedMass3D(double massCoefficient)
		{
			this.massCoefficient = massCoefficient;
		}

		public ConcentratedMass3D(double massCoefficient, IElementDofEnumerator dofEnumerator)
			: this(massCoefficient)
		{
			this.dofEnumerator = dofEnumerator;
		}

		public IMatrix MassMatrix(IElement element)
		{
			var mass = Matrix.CreateZero(3, 3);
			mass[0, 0] = massCoefficient;
			mass[1, 1] = massCoefficient;
			mass[2, 2] = massCoefficient;
			return mass;
		}

		public IMatrix StiffnessMatrix(IElement element) => Matrix.CreateZero(3, 3);

		public IMatrix DampingMatrix(IElement element) => Matrix.CreateZero(3, 3);

		public void ResetMaterialModified() { }

		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements,
			double[] localdDisplacements)
			=> new Tuple<double[], double[]>(new double[6], new double[6]);

		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
			=> CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
			=> new double[6];

		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			var accelerations = new double[3];
			IMatrix massMatrix = MassMatrix(element);

			foreach (MassAccelerationLoad load in loads)
			{
				int index = 0;
				foreach (IDofType[] nodalDOFTypes in dofs)
					foreach (IDofType dofType in nodalDOFTypes)
					{
						if (dofType == load.DOF) accelerations[index] += load.Amount;
						index++;
					}
			}

			return massMatrix.Multiply(accelerations);
		}

		public void ClearMaterialState() { }
		public void SaveMaterialState() { }
		public void ClearMaterialStresses() { }
	}
}
