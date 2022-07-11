using System;
using System.Collections.Generic;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Numerics.Interpolation.GaussPointExtrapolation;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Meshes;

namespace MGroup.FEM.ConvectionDiffusion.Isoparametric
{
	public class ConvectionDiffusionElement3D : IConvectionDiffusionElementType, ICell<INode>
	{
		private readonly static IDofType[] nodalDOFTypes = new IDofType[] { ConvectionDiffusionDof.UnknownVariable};
		private readonly IDofType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
		private readonly IConvectionDiffusionProperties material;

		public ConvectionDiffusionElement3D(IReadOnlyList<INode> nodes, IIsoparametricInterpolation3D interpolation,
		IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
		IGaussPointExtrapolation3D gaussPointExtrapolation, IConvectionDiffusionProperties material)
		{
			this.material = material;
			this.GaussPointExtrapolation = gaussPointExtrapolation;
			this.Nodes = nodes;
			this.Interpolation = interpolation;
			this.QuadratureForConsistentMass = quadratureForMass;
			this.QuadratureForStiffness = quadratureForStiffness;

			dofTypes = new IDofType[nodes.Count][];
			for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new IDofType[] { ConvectionDiffusionDof.UnknownVariable };
		}

		public CellType CellType => Interpolation.CellType;

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public int ID { get; set; }

		public int SubdomainID { get; set; }

		public IReadOnlyList<INode> Nodes { get; }

		public IGaussPointExtrapolation3D GaussPointExtrapolation { get; }

		public IIsoparametricInterpolation3D Interpolation { get; }

		public IQuadrature3D QuadratureForConsistentMass { get; }

		public IQuadrature3D QuadratureForStiffness { get; }

		public bool ConstitutiveLawModified => throw new NotImplementedException();

		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		public IMatrix FirstTimeDerivativeMatrix()
		{
			return BuildFirstTimeDerivativeMatrix();
		}

		public IMatrix DiffusionMatrix()
		{
			return BuildDiffusionMatrix();
		}

		public IMatrix ConvectionMatrix()
		{
			return BuildConvectionMatrix();
		}

		public IMatrix ProductionMatrix()
		{
			return BuildProductionMatrix();
		}

		public IMatrix PhysicsMatrix()
		{
			return DiffusionMatrix().Add(ConvectionMatrix()).Add(ProductionMatrix());
		}

		public Matrix BuildFirstTimeDerivativeMatrix() //TODO: Do it
		{
			throw new NotImplementedException();
		}

		public Matrix BuildDiffusionMatrix() //TODO: Do it
		{
			throw new NotImplementedException();
		}

		public Matrix BuildConvectionMatrix() //TODO: Do it
		{
			throw new NotImplementedException();
		}

		public Matrix BuildProductionMatrix() //TODO: Do it
		{
			throw new NotImplementedException();
		}

		private Matrix BuildDeformationMatrix(Matrix shapeGradientsCartesian)
		{
			//TODO: isn't this just the transpose of [dNi/dxj]?
			var deformation = Matrix.CreateZero(3, Nodes.Count);
			for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
			{
				deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
				deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
				deformation[2, nodeIdx] = shapeGradientsCartesian[nodeIdx, 2];
			}
			return deformation;
		}

		public Matrix BuildShapeFunctionMatrix(double[] shapeFunctions)
		{
			return Matrix.CreateFromArray(shapeFunctions, 1, shapeFunctions.Length);
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofTypes;

		public void ResetConstitutiveLawModified()
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateResponse(double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateResponseIntegral()
		{
			throw new NotImplementedException();
		}

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public void SaveConstitutiveLawState()
		{
			throw new NotImplementedException();
		}

		public void ClearConstitutiveLawState()
		{
			throw new NotImplementedException();
		}

		public void ClearConstitutiveLawStresses()
		{
			throw new NotImplementedException();
		}
	}
}
