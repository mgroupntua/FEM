using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Structural.Elements
{
	public class Rod2D : IStructuralFiniteElement
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[2] { StructuralDof.TranslationX, StructuralDof.TranslationY };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		private readonly double youngModulus;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		public double Density { get; set; }
		public double SectionArea { get; set; }

		public Rod2D(double youngModulus)
		{
			this.youngModulus = youngModulus;
		}

		public Rod2D(double youngModulus, IElementDofEnumerator dofEnumerator)
			: this(youngModulus)
		{
			this.dofEnumerator = dofEnumerator;
		}

		public CellType CellType { get; } = CellType.Line;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		//TODO: this should be either cached, or even better the calculations should be incorporated into Stiffness()
		public IMatrix TransformationMatrix(Element element)
		{
			double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
			double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
			double L = Math.Sqrt(x2 + y2);
			double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
			double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;

			// T = [ cos sin 0 0; 0 0 cos sin]
			var transformation = Matrix.CreateZero(2, 4);
			transformation[0, 0] = c;
			transformation[0, 1] = s;
			transformation[1, 2] = c;
			transformation[1, 3] = s;
			return transformation;
		}

		/// <summary>
		/// Stress0         Stress1
		/// -> ------------ ->
		/// </summary>
		/// <param name="element"></param>
		/// <param name="localDisplacements"></param>
		/// <param name="local_d_Displacements"></param>
		/// <returns></returns>
		public double CalculateAxialStress(Element element, double[] localDisplacements, double[] local_d_Displacements)
		{
			double[] globalStresses = CalculateStresses(element, localDisplacements, local_d_Displacements).Item2; // item1 = strains
			IMatrix transformation = TransformationMatrix(element);
			double[] localStresses = transformation.Multiply(globalStresses); // In local natural system there are 2 dofs
																			  // If Stress1 = localStresses[1] > 0 => tension. Else compression
			return localStresses[1];
		}

		#region IElementType Members

		public int ID => 1;

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofs;

		public IList<Node> GetNodesForMatrixAssembly(Element element) => element.Nodes;

		public IMatrix StiffnessMatrix(IElement element)
		{
			double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
			double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
			double L = Math.Sqrt(x2 + y2);
			double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
			double c2 = c * c;
			double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
			double s2 = s * s;
			double cs = c * s;
			double E = this.youngModulus;
			double A = SectionArea;

			return dofEnumerator.GetTransformedMatrix(
				Matrix.CreateFromArray(new double[,]
				{
					{A*E*c2/L, A*E*cs/L, -A*E*c2/L, -A*E*cs/L },
					{A*E*cs/L, A*E*s2/L, -A*E*cs/L, -A*E*s2/L },
					{-A*E*c2/L, -A*E*cs/L, A*E*c2/L, A*E*cs/L },
					{-A*E*cs/L, -A*E*s2/L, A*E*cs/L, A*E*s2/L }
				}));
		}

		public IMatrix MassMatrix(IElement element)
		{
			double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
			double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
			double L = Math.Sqrt(x2 + y2);

			double totalMassOver2 = Density * SectionArea * L / 2.0;

			// Lumped mass: M = [m/2 0 0 0; 0 m/2 0 0; 0 0 m/2 0; 0 0 0 m/2]
			int order = 4;
			var lumpedMass = Matrix.CreateZero(order, order);
			for (int i = 0; i < order; ++i) lumpedMass[i, i] = totalMassOver2;
			return lumpedMass;
		}

		public IMatrix DampingMatrix(IElement element) => throw new NotImplementedException();

		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] local_Displacements,
			double[] local_d_Displacements)
		{
			// WARNING: 1) No strains are computed 2) localdDisplacements are not used.
			double[] strains = null;
			double[] forces = CalculateForces(element, local_Displacements, local_d_Displacements);
			double[] stresses = Array.ConvertAll(forces, x => x / SectionArea);
			return new Tuple<double[], double[]>(strains, stresses);
		}

		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
		{
			return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
		}

		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			IMatrix stiffness = StiffnessMatrix(element);
			return stiffness.Multiply(localdDisplacements);
		}

		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			var accelerations = new double[4];
			IMatrix massMatrix = MassMatrix(element);

			int index = 0;
			foreach (MassAccelerationLoad load in loads)
				foreach (IDofType[] nodalDOFTypes in dofs)
					foreach (IDofType dofType in nodalDOFTypes)
					{
						if (dofType == load.DOF) accelerations[index] += load.Amount;
						index++;
					}

			return massMatrix.Multiply(accelerations);
		}

		public void SaveMaterialState() { }

		#endregion

		#region IFiniteElement Members


		public bool MaterialModified => false;

		public void ResetMaterialModified() { }

		#endregion

		#region IFiniteElement Members

		public void ClearMaterialState() { }

		public void ClearMaterialStresses() => throw new NotImplementedException();

		#endregion
	}
}
