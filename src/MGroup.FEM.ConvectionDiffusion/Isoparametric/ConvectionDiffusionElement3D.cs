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

		public Matrix BuildFirstTimeDerivativeMatrix()
		{
			int numDofs = Nodes.Count;
			var capacity = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
				Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				capacity.AxpyIntoThis(partial, dA);
			}
			//TODO Ask Giannis : Should it be multiplies by surface
			capacity.ScaleIntoThis(material.FirstTimeDerivativeCoeff * Thickness);
			return capacity;
		}

		public Matrix BuildDiffusionMatrix()
		{
			int numDofs = Nodes.Count;
			var diffusion = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				Matrix shapeGradientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);

				var deformation = Matrix.CreateZero(3, Nodes.Count);
				for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
				{
					deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
					deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
					deformation[2, nodeIdx] = shapeGradientsCartesian[nodeIdx, 2];
				}

				Matrix partialK = deformation.Transpose() * deformation;

				double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
				diffusion.AxpyIntoThis(partialK, dA * material.DiffusionCoeff * Thickness);
			}

			return diffusion;
		}

		public Matrix BuildConvectionMatrix() // TODO: Check this. Cannot be the same as Capacity and production
		{
			int numDofs = Nodes.Count;
			var convection = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix_line = BuildShapeFunctionMatrix(shapeFunctions[gp]);
				double[,] shapeFunctionArray = new double[3, numDofs];

				for (var col = 0; col < numDofs; col++)
				{
					shapeFunctionArray[0, col] = shapeFunctionMatrix_line[0, col];
					shapeFunctionArray[1, col] = shapeFunctionMatrix_line[0, col];
					shapeFunctionArray[3, col] = shapeFunctionMatrix_line[0, col];
				}

				Matrix partial = Matrix.CreateFromArray(shapeFunctionArray).Transpose() * shapeGradientsNatural[gp];

				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);

				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				convection.AxpyIntoThis(partial, dA * material.ConvectionCoeff * Thickness);
			}

			return convection;
		}

		public Matrix BuildProductionMatrix()
		{
			int numDofs = Nodes.Count;
			var production = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
				Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				production.AxpyIntoThis(partial, dA);
			}

			production.ScaleIntoThis(material.IndependentSourceCoeff * Thickness);
			return production;
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
