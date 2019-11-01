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
	public enum SpringDirections
	{
		X = 0, Y, Z, XY, YZ, XZ, XYZ
	}

	public class SpringDamper3D : IStructuralFiniteElement
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		private readonly double springCoefficient, dampingCoefficient;
		private readonly SpringDirections springDirections, dampingDirections;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		public int ID => 999;
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

			var d = new List<IReadOnlyList<IDofType>>();
			foreach (var node in element.Nodes)
			{
				var nodeDofs = new List<IDofType>();
				nodeDofs.AddRange(nodalDOFTypes);
				d.Add(nodeDofs);
			}
			return d;
		}

		public bool MaterialModified => false;

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

		public IMatrix StiffnessMatrix(IElement element)
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

		public IMatrix MassMatrix(IElement element) => Matrix.CreateZero(6, 6);

		public IMatrix DampingMatrix(IElement element)
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

		public void ResetMaterialModified() { }

		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements,
			double[] localdDisplacements)
			=> new Tuple<double[], double[]>(new double[6], new double[6]);

		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
			=> CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			IMatrix stiffnessMatrix = StiffnessMatrix(element);
			return stiffnessMatrix.Multiply(localDisplacements);
		}

		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) => new double[6];

		public void ClearMaterialState() { }
		public void SaveMaterialState() { }
		public void ClearMaterialStresses() { }
	}
}
