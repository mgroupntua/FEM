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
	public enum SpringDirections
	{
		X = 0, Y, Z, XY, YZ, XZ, XYZ
	}

	public class SpringDamper3D : IStructuralElementType
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		private readonly double springCoefficient, dampingCoefficient;
		private readonly SpringDirections springDirections, dampingDirections;
		private double[] currentDisplacements;
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
			var d = new List<IReadOnlyList<IDofType>>();
			foreach (var node in Nodes)
			{
				var nodeDofs = new List<IDofType>();
				nodeDofs.AddRange(nodalDOFTypes);
				d.Add(nodeDofs);
			}
			return d;
		}

		public bool ConstitutiveLawModified => false;

		public SpringDamper3D(double springCoefficient, double dampingCoefficient, SpringDirections springDirections,
			SpringDirections dampingDirections)
		{
			this.springCoefficient = springCoefficient;
			this.dampingCoefficient = dampingCoefficient;
			this.springDirections = springDirections;
			this.dampingDirections = dampingDirections;
		}

		public SpringDamper3D(double springCoefficient, double dampingCoefficient, SpringDirections springDirections,
			SpringDirections dampingDirections, IElementDofEnumerator dofEnumerator)
			: this(springCoefficient, dampingCoefficient, springDirections, dampingDirections)
		{
			this.dofEnumerator = dofEnumerator;
		}

		public IMatrix StiffnessMatrix()
		{
			double x = (springDirections == SpringDirections.X || springDirections == SpringDirections.XY || springDirections == SpringDirections.XZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
			double y = (springDirections == SpringDirections.Y || springDirections == SpringDirections.XY || springDirections == SpringDirections.YZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
			double z = (springDirections == SpringDirections.Z || springDirections == SpringDirections.XZ || springDirections == SpringDirections.YZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
			return Matrix.CreateFromArray(new double[,]
				{
				   { x, 0, 0, -x, 0, 0 },
				   { 0, y, 0, 0, -y, 0 },
				   { 0, 0, z, 0, 0, -z },
				   {-x, 0, 0, x, 0, 0 },
				   { 0,-y, 0, 0, y, 0 },
				   { 0, 0,-z, 0, 0, z }
				}
				);

			//return SymmetricMatrix.CreateFromArray(new double[] 
			//{
			//    x, 0, 0, -x, 0, 0,
			//       y, 0, 0, -y, 0, 
			//          z, 0, 0, -z,
			//             x, 0, 0,
			//                y, 0,
			//                   z
			//});
		}
		public IMatrix PhysicsMatrix()
		{
			return StiffnessMatrix();
		}
		public IMatrix MassMatrix() => Matrix.CreateZero(6, 6);

		public IMatrix DampingMatrix()
		{
			double x = (dampingDirections == SpringDirections.X || dampingDirections == SpringDirections.XY || dampingDirections == SpringDirections.XZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;
			double y = (dampingDirections == SpringDirections.Y || dampingDirections == SpringDirections.XY || dampingDirections == SpringDirections.YZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;
			double z = (dampingDirections == SpringDirections.Z || dampingDirections == SpringDirections.XZ || dampingDirections == SpringDirections.YZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;

			return Matrix.CreateFromArray(new double[,]
				{
				   { x, 0, 0, -x, 0, 0 },
				   { 0, y, 0, 0, -y, 0 },
				   { 0, 0, z, 0, 0, -z },
				   {-x, 0, 0, x, 0, 0 },
				   { 0,-y, 0, 0, y, 0 },
				   { 0, 0,-z, 0, 0, z }
				}
				);

			//return SymmetricMatrix.CreateFromArray(new double[] 
			//{
			//    x, 0, 0, -x, 0, 0,
			//    y, 0, 0, -y, 0, 
			//    z, 0, 0, -z,
			//    x, 0, 0,
			//    y, 0,
			//    z
			//});
		}

		public void ResetConstitutiveLawModified() { }

		private double[] CalculateResponseIntegral(double[] localDisplacements)
		{
			IMatrix stiffnessMatrix = StiffnessMatrix();
			return stiffnessMatrix.Multiply(localDisplacements);
		}

		public Tuple<double[], double[]> CalculateResponse(double[] localDisplacements)
		{
			if (currentDisplacements == null || currentDisplacements.Length != localDisplacements.Length)
			{
				currentDisplacements = new double[localDisplacements.Length];
			}

			Array.Copy(localDisplacements, currentDisplacements, localDisplacements.Length);
			return new Tuple<double[], double[]>(new double[6], new double[6]);
		}

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
			=> CalculateResponseIntegral();

		public double[] CalculateResponseIntegral() => CalculateResponseIntegral(currentDisplacements);

		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads) => new double[6];

		public void ClearConstitutiveLawState() { }
		public void SaveConstitutiveLawState(IHaveState externalState) { }
		public void ClearConstitutiveLawStresses() { }
	}
}
