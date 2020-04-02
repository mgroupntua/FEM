using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;
using MGroup.FEM.Interpolation;
using MGroup.FEM.Interpolation.GaussPointExtrapolation;
using MGroup.FEM.Interpolation.Jacobians;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;
using MGroup.MSolve.Geometry.Coordinates;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.FEM.Structural.Elements
{
	/// <summary>
	/// Represents a continuum finite element for 3D problems. Specific elements (e.g. Hexa8, Hexa20, ...) can be created using
	/// the appropriate <see cref="IIsoparametricInterpolation3D_OLD"/>, <see cref="IQuadrature3D"/> etc. strategies. 
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class ContinuumElement3D : IStructuralFiniteElement, ICell<Node>
	{
		private readonly static IDofType[] nodalDOFTypes = new IDofType[]
		{
			StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
		};

		private readonly IDofType[][] dofTypes;
		private IDynamicMaterial dynamicProperties;
		private readonly IReadOnlyList<IContinuumMaterial3D> materialsAtGaussPoints;

		private double[][] strainsVec;
		private double[][] strainsVecLastConverged;

		public ContinuumElement3D(IReadOnlyList<Node> nodes, IIsoparametricInterpolation3D interpolation,
			IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
			IGaussPointExtrapolation3D gaussPointExtrapolation,
			IReadOnlyList<IContinuumMaterial3D> materialsAtGaussPoints, IDynamicMaterial dynamicProperties)
		{
			this.dynamicProperties = dynamicProperties;
			this.materialsAtGaussPoints = materialsAtGaussPoints;
			this.GaussPointExtrapolation = gaussPointExtrapolation;
			this.Nodes = nodes;
			this.Interpolation = interpolation;
			this.QuadratureForConsistentMass = quadratureForMass;
			this.QuadratureForStiffness = quadratureForStiffness;

			dofTypes = new IDofType[nodes.Count][];
			for (int i = 0; i < nodes.Count; i++)
			{
				dofTypes[i] = new IDofType[]
				{
					StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
				};
			}

			strainsVec = new double[materialsAtGaussPoints.Count][];
			strainsVecLastConverged = new double[materialsAtGaussPoints.Count][];
			for (int gpoint = 0; gpoint < materialsAtGaussPoints.Count; gpoint++)
			{
				strainsVec[gpoint] = new double[6];
				strainsVecLastConverged[gpoint] = new double[6];
			}
		}

		public CellType CellType => Interpolation.CellType;
		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();
		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;
		public IGaussPointExtrapolation3D GaussPointExtrapolation { get; }
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;

		public int ID => throw new NotImplementedException(
			"Element type codes should be in a settings class. Even then it's a bad design choice");

		public IIsoparametricInterpolation3D Interpolation { get; }

		public bool MaterialModified
		{
			get
			{
				foreach (IContinuumMaterial3D material in materialsAtGaussPoints)
				{
					if (material.Modified) return true;
				}
				return false;
			}
		}

		public IReadOnlyList<Node> Nodes { get; }
		public IQuadrature3D QuadratureForConsistentMass { get; }
		public IQuadrature3D QuadratureForStiffness { get; }

		public Matrix BuildConsistentMassMatrix()
		{
			int numberOfDofs = 3 * Nodes.Count;
			var mass = Matrix.CreateZero(numberOfDofs, numberOfDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
				Matrix partial = shapeFunctionMatrix.MultiplyRight(shapeFunctionMatrix, true, false);
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				mass.AxpyIntoThis(partial, dA);
			}
			mass.ScaleIntoThis(dynamicProperties.Density);
			return mass;
		}

		public Matrix BuildLumpedMassMatrix()
		{
			int numberOfDofs = 3 * Nodes.Count;
			var lumpedMass = Matrix.CreateZero(numberOfDofs, numberOfDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			double area = 0;
			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				area += jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
			}

			double nodalMass = area * dynamicProperties.Density / Nodes.Count;
			for (int i = 0; i < numberOfDofs; i++) lumpedMass[i, i] = nodalMass;

			return lumpedMass;
		}

		public IMatrix BuildStiffnessMatrix()
		{
			int numberOfDofs = 3 * Nodes.Count;
			var stiffness = Matrix.CreateZero(numberOfDofs, numberOfDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				IMatrixView constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				Matrix shapeGradientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
				Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);

				Matrix partial = deformation.ThisTransposeTimesOtherTimesThis(constitutive);
				double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
				stiffness.AxpyIntoThis(partial, dA);
			}

			return DofEnumerator.GetTransformedMatrix(stiffness);
		}

		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			int numberOfDofs = 3 * Nodes.Count;
			var accelerations = new double[numberOfDofs];
			IMatrix massMatrix = MassMatrix(element);

			foreach (var load in loads)
			{
				int index = 0;
				foreach (var nodalDOFTypes in dofTypes)
				{
					foreach (var dofType in nodalDOFTypes)
					{
						if (dofType == load.DOF) accelerations[index] += load.Amount;
						index++;
					}
				}
			}

			return massMatrix.Multiply(accelerations);
		}

		public double[] CalculateForces(IElement element, double[] localTotalDisplacements, double[] localDisplacements)
		{
			int numberOfDofs = 3 * Nodes.Count;
			var Forces = Vector.CreateZero(numberOfDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				Vector Stresses = Vector.CreateFromArray(materialsAtGaussPoints[gp].Stresses);
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				Matrix shapeGradientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
				Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);
				Vector gpForces = deformation.Transpose() * (Stresses);
				double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
				gpForces.ScaleIntoThis(dA);
				Forces.AddIntoThis(gpForces);
			}
			return Forces.CopyToArray();
		}

		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
		{
			return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
		}

		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements,
			double[] localdDisplacements)
		{
			int numberOfDofs = 3 * Nodes.Count;
			var Forces = Vector.CreateZero(numberOfDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			//double[] strains = new double[6];
			//for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			//{
			//	strains = new double[6];
			//	var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
			//	Matrix shapeGradientsCartesian =
			//		jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
			//	Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);
			//	strains = deformation.Multiply(localDisplacements);
			//	materialsAtGaussPoints[gp].UpdateMaterial(strains);
			//}

			//double[] strains = new double[6];
			double[] strainsVecMinusLastConvergedValue = new double[6];
			for (int gpo = 0; gpo < QuadratureForStiffness.IntegrationPoints.Count; ++gpo)
			{
				//strainsVec[gpo] = new double[6];
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gpo]);
				Matrix shapeGradientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gpo]);
				Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);
				strainsVec[gpo] = deformation.Multiply(localDisplacements);
				strainsVecMinusLastConvergedValue = new double[6]
				{
					strainsVec[gpo][0]- strainsVecLastConverged[gpo][0],
					strainsVec[gpo][1]- strainsVecLastConverged[gpo][1],
					strainsVec[gpo][2]- strainsVecLastConverged[gpo][2],
					strainsVec[gpo][3]- strainsVecLastConverged[gpo][3],
					strainsVec[gpo][4]- strainsVecLastConverged[gpo][4],
					strainsVec[gpo][5]- strainsVecLastConverged[gpo][5]
				};
				materialsAtGaussPoints[gpo].UpdateMaterial(strainsVecMinusLastConvergedValue);
				//To update with total strain simplY = materialsAtGaussPoints[npoint].UpdateMaterial(strainsVec[npoint]);
			}


			return new Tuple<double[], double[]>(strainsVec[materialsAtGaussPoints.Count-1], materialsAtGaussPoints[materialsAtGaussPoints.Count - 1].Stresses);
		}

		public double CalculateVolume()
		{
			//TODO: Linear elements can use the more efficient rules for volume of polygons. Therefore this method should be 
			//      delegated to the interpolation.
			//TODO: A different integration rule should be used for integrating constant functions. For linear elements there
			//      is only 1 Gauss point (most probably), therefore the computational cost could be the same as using the 
			//      polygonal formulas.
			double volume = 0.0;
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				volume += jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight; //TODO: this is used by all methods that integrate. I should cache it.
			}
			return volume;
		}

		public void ClearMaterialState()
		{
			foreach (var material in materialsAtGaussPoints) material.ClearState();
		}

		public void ClearMaterialStresses()
		{
			foreach (var material in materialsAtGaussPoints) material.ClearStresses();
		}

		public IMatrix DampingMatrix(IElement element)
		{
			IMatrix damping = BuildStiffnessMatrix();
			damping.ScaleIntoThis(dynamicProperties.RayleighCoeffStiffness);
			damping.AxpyIntoThis(MassMatrix(element), dynamicProperties.RayleighCoeffMass);
			return damping;
		}

		/// <summary>
		/// Calculates the coordinates of the centroid of this element.
		/// </summary>
		public CartesianPoint FindCentroid()
			=> Interpolation.TransformNaturalToCartesian(Nodes, new NaturalPoint(0.0, 0.0, 0.0));


		public IMatrix MassMatrix(IElement element)
		{
			return BuildLumpedMassMatrix();
		}

		public void ResetMaterialModified()
		{
			foreach (var material in materialsAtGaussPoints) material.ResetModified();
		}

		public void SaveMaterialState()
		{
			for (int npoint = 0; npoint < materialsAtGaussPoints.Count; npoint++)
			{
				for (int i1 = 0; i1 < 6; i1++)
				{ strainsVecLastConverged[npoint][i1] = strainsVec[npoint][i1]; }
			}
			foreach (var m in materialsAtGaussPoints) m.SaveState();
		}

		public IMatrix StiffnessMatrix(IElement element) => DofEnumerator.GetTransformedMatrix(BuildStiffnessMatrix());

		public (IReadOnlyList<double[]> strains, IReadOnlyList<double[]> stresses)
			UpdateStrainStressesAtGaussPoints(double[] localDisplacements)
		{
			int numberOfGPs = QuadratureForStiffness.IntegrationPoints.Count;
			var strains = new double[numberOfGPs][];
			var stresses = new double[numberOfGPs][];
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			for (int gp = 0; gp < numberOfGPs; gp++)
			{
				IMatrixView constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				Matrix shapeGrandientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
				Matrix deformation = BuildDeformationMatrix(shapeGrandientsCartesian);

				strains[gp] = deformation.Multiply(localDisplacements);
				stresses[gp] = constitutive.Multiply(strains[gp]);
			}

			return (strains, stresses);
		}

		/// <summary>
		/// Assembles the deformation matrix of a solid element.
		/// The calculation are based on <see cref="https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch08.d/AFEM.Ch08.pdf"/>
		/// paragraph 8.4, equation 8.7
		/// </summary>
		/// <param name="shapeGradients"></param>
		/// <returns></returns>
		private Matrix BuildDeformationMatrix(Matrix shapeGradientsCartesian)
		{
			var deformation = Matrix.CreateZero(6, 3 * Nodes.Count);
			for (int nodeIdx = 0; nodeIdx < Nodes.Count; nodeIdx++)
			{
				int col0 = 3 * nodeIdx;
				int col1 = 3 * nodeIdx + 1;
				int col2 = 3 * nodeIdx + 2;

				deformation[0, col0] = shapeGradientsCartesian[nodeIdx, 0];
				deformation[1, col1] = shapeGradientsCartesian[nodeIdx, 1];
				deformation[2, col2] = shapeGradientsCartesian[nodeIdx, 2];

				deformation[3, col0] = shapeGradientsCartesian[nodeIdx, 1];
				deformation[3, col1] = shapeGradientsCartesian[nodeIdx, 0];

				deformation[4, col1] = shapeGradientsCartesian[nodeIdx, 2];
				deformation[4, col2] = shapeGradientsCartesian[nodeIdx, 1];

				deformation[5, col0] = shapeGradientsCartesian[nodeIdx, 2];
				deformation[5, col2] = shapeGradientsCartesian[nodeIdx, 0];
			}

			return deformation;
		}

		/// <summary>
		/// The shape function matrix is 2-by-2n, where n = is the number of shape functions. Row 0 corresponds to dof X, while
		/// row 1 to dof Y, etc.
		/// </summary>
		private Matrix BuildShapeFunctionMatrix(double[] shapeFunctions)
		{
			var shapeFunctionMatrix = Matrix.CreateZero(3, 3 * shapeFunctions.Length);
			for (int i = 0; i < shapeFunctions.Length; i++)
			{
				shapeFunctionMatrix[0, 3 * i] = shapeFunctions[i];
				shapeFunctionMatrix[1, 2 * i + 1] = shapeFunctions[i];
				shapeFunctionMatrix[2, 3 * i + 2] = shapeFunctions[i];
			}
			return shapeFunctionMatrix;
		}

	}
}
