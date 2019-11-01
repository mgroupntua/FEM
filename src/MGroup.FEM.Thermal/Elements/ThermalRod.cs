using System;
using System.Collections.Generic;
using System.Diagnostics;
using MGroup.Constitutive.Thermal;
using MGroup.FEM.Embedding;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Thermal.Elements
{
	/// <summary>
	/// Finite element for heat transfer along a single direction. Can be used for 1D, 2D and 3D problems. Does not take into 
	/// account geometric non-linearities.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class ThermalRod : IFiniteElement, IEmbeddedElement
	{
		private const int numNodes = 2;
		private const int numDofs = 2;
		private static readonly IDofType[][] dofTypes = {
			new IDofType[] { ThermalDof.Temperature }, new IDofType[] { ThermalDof.Temperature } };

		private readonly IThermalMaterial material;

		public ThermalRod(IReadOnlyList<Node> nodes, double crossSectionArea, IThermalMaterial material)
		{
			Debug.Assert(nodes.Count == 2, "Thermal rod element must have exactly 2 nodes.");
			this.material = material;
			this.Nodes = nodes;
			this.CrossSectionArea = crossSectionArea;
			this.Length = nodes[0].CalculateEuclidianDistanceFrom(nodes[1]);
		}

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public int ID => throw new NotImplementedException(
			"Element type codes should be in a settings class. Even then it's a bad design choice");
		public CellType CellType { get; } = CellType.Line;

		public double CrossSectionArea { get; }
		public double Length { get; }
		public IReadOnlyList<Node> Nodes { get; }

		public bool MaterialModified => throw new NotImplementedException();

		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

		public IMatrix MassMatrix(IElement element)
		{
			return BuildCapacityMatrix();
		}

		public Matrix BuildCapacityMatrix()
		{
			double kdAL = material.SpecialHeatCoeff * material.Density * CrossSectionArea * Length;
			double[,] capacity = { { kdAL / 3.0, kdAL / 6.0 }, { kdAL / 6.0, kdAL / 3.0 } };
			return Matrix.CreateFromArray(capacity);
		}

		public Matrix BuildConductivityMatrix()
		{

			double cAoverL = material.ThermalConductivity * CrossSectionArea / Length;
			double[,] conductivity = { { cAoverL, -cAoverL }, { -cAoverL, cAoverL } };
			return Matrix.CreateFromArray(conductivity);
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			throw new NotImplementedException();
		}

		public void SaveMaterialState()
		{
			throw new NotImplementedException();
		}

		public void ClearMaterialState()
		{
			throw new NotImplementedException();
		}

		public void ClearMaterialStresses()
		{
			throw new NotImplementedException();
		}

		public IMatrix StiffnessMatrix(IElement element)
		{
			return DofEnumerator.GetTransformedMatrix(BuildConductivityMatrix());
		}

		public IMatrix DampingMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public Dictionary<IDofType, int> GetInternalNodalDOFs(Element element, Node node)
		{
			if (node.ID == this.Nodes[0].ID) return new Dictionary<IDofType, int> { { ThermalDof.Temperature, 0 } };
			else if (node.ID == this.Nodes[1].ID) return new Dictionary<IDofType, int> { { ThermalDof.Temperature, 1 } };
			else throw new ArgumentException($"GetInternalNodalDOFs: Node {node.ID} not found in element {element.ID}.");
		}

		public double[] GetLocalDOFValues(Element hostElement, double[] hostDOFValues)
		{
			return DofEnumerator.GetTransformedDisplacementsVector(hostDOFValues);
		}
	}
}
