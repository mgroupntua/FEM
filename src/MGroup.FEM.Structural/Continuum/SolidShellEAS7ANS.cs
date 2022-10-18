using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Constitutive.Structural.Transient;
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
using MGroup.MSolve.Geometry.Coordinates;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.FEM.Structural.Continuum
{
	/// <summary>
	/// Continuum finite Element for 3d problems with EAS7 & ANS methods
	/// for volumetric, membrane and Poisson thickness locking
	/// Authors: Kostas Margaronis.
	/// </summary>
	public class SolidShellEAS7ANS : IStructuralElementType, ICell<INode>
	{
		private readonly static IDofType[] nodalDOFTypes = new IDofType[]
		{
			StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
		};

		private readonly IDofType[][] dofTypes;
		private ITransientAnalysisProperties dynamicProperties;
		private readonly IReadOnlyList<IContinuumMaterial3D> materialsAtGaussPoints;

		private double[][] strainsVec;
		private double[][] strainsVecLastConverged;
		private double[][] lastStresses;
		private double[] internalParamVectorEAS;
		private double[] internalForceVectorEAS;
		private Matrix internalEASLMatrix;
		private Matrix internalEASDMatrix;
		private Vector displacementVectorPreviousIncrement;

		public SolidShellEAS7ANS(
			IReadOnlyList<INode> nodes,
			IContinuumMaterial3D commonMaterial,
			ITransientAnalysisProperties dynamicProperties)
		{
			IQuadrature3D integrationsForStiffness = GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2);
			IQuadrature3D integrationsForMass = GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2);
			IGaussPointExtrapolation3D extrapolation = ExtrapolationGaussLegendre2x2x2.UniqueInstance;
			var materialsGaussPoints = new IContinuumMaterial3D[integrationsForStiffness.IntegrationPoints.Count];
			for (int gp = 0; gp < integrationsForStiffness.IntegrationPoints.Count; ++gp) materialsGaussPoints[gp] = (IContinuumMaterial3D)commonMaterial.Clone();
			this.materialsAtGaussPoints = materialsGaussPoints;
			this.GaussPointExtrapolation = extrapolation;
			this.Nodes = nodes;
			this.Interpolation = InterpolationHexa8.UniqueInstance;
			this.QuadratureForConsistentMass = integrationsForMass;
			this.QuadratureForStiffness = integrationsForStiffness;
			this.dynamicProperties = dynamicProperties;
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
			lastStresses = new double[materialsAtGaussPoints.Count][];
			for (int gpoint = 0; gpoint < materialsAtGaussPoints.Count; gpoint++)
			{
				strainsVec[gpoint] = new double[6];
				strainsVecLastConverged[gpoint] = new double[6];
				lastStresses[gpoint] = new double[6];
			}

			displacementVectorPreviousIncrement = Vector.CreateFromArray(new double[24]);
			internalParamVectorEAS = new double[7];
			internalForceVectorEAS = new double[7];
			internalEASLMatrix = Matrix.CreateFromArray(new double[7, 24]);
			internalEASDMatrix = Matrix.CreateFromArray(new double[7, 7]);
		}

		public CellType CellType => Interpolation.CellType;

		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public IGaussPointExtrapolation3D GaussPointExtrapolation { get; }

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofTypes;

		public int ID { get; set; }

		public int SubdomainID { get; set; }

		public IReadOnlyList<INode> Nodes { get; }

		public IIsoparametricInterpolation3D Interpolation { get; }

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
			var LeMatrix = Matrix.CreateFromArray(new double[7, 24]);
			var deMatrix = Matrix.CreateFromArray(new double[7, 7]);
			var jacobianZero = new IsoparametricJacobian3D(Nodes,
							Interpolation.EvaluateNaturalGradientsAt(new NaturalPoint(xi: 0d, eta: 0d, zeta: 0d)));
			var detJ0 = jacobianZero.DirectDeterminant;
			var jacobianZeroInverse = jacobianZero.InverseMatrix;
			var transformationMat0 = CalculateTransformationMatrix(jacobianZeroInverse);
			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				IMatrixView constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
				var localdN = CalculateShapeFunctionsLocalDerivatives(
					QuadratureForStiffness.IntegrationPoints[gp].Coordinates);
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				var jacobianMat = jacobian.DirectMatrix;
				var detJ = jacobian.DirectDeterminant;
				var jacobianInverse = jacobian.InverseMatrix;
				var transformationMatrix = CalculateTransformationMatrix(jacobianInverse);
				var deformation = CalculateDeformationMatrix(localdN, jacobianMat, QuadratureForStiffness.IntegrationPoints[gp].Coordinates,
					transformationMatrix);
				//BuildDeformationMatrix(shapeGradientsCartesian);

				Matrix partial = deformation.ThisTransposeTimesOtherTimesThis(constitutive);
				double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
				stiffness.AxpyIntoThis(partial, dA);
				var Gamma = CalculateEnhancedStrainMatrixGamma(Matrix.CreateFromArray(transformationMat0), CreateEnhancedStrainsInterpolationMatrix(QuadratureForStiffness.IntegrationPoints[gp].Coordinates),
						detJ0, detJ);
				LeMatrix += (Gamma.Transpose().MultiplyRight(constitutive) * deformation).Scale(detJ * QuadratureForStiffness.IntegrationPoints[gp].Weight);
				deMatrix += (Gamma.ThisTransposeTimesOtherTimesThis(constitutive)).Scale(detJ * QuadratureForStiffness.IntegrationPoints[gp].Weight);
			}

			var deMatrixInv = deMatrix.Invert();
			var tangentMatrix = stiffness - LeMatrix.Transpose() * deMatrixInv * LeMatrix;
			internalEASLMatrix = Matrix.CreateFromArray(LeMatrix.CopyToArray2D());
			internalEASDMatrix = Matrix.CreateFromArray(deMatrixInv.CopyToArray2D());
			return DofEnumerator.GetTransformedMatrix(tangentMatrix);
		}

		public double[] CalculateResponseIntegral()
		{
			var vectorRe = Vector.CreateFromArray(new double[24]);
			var vectorPe = Vector.CreateFromArray(new double[7]);
			var jacobianZero = new IsoparametricJacobian3D(Nodes,
							Interpolation.EvaluateNaturalGradientsAt(new NaturalPoint(xi: 0d, eta: 0d, zeta: 0d)));
			var detJ0 = jacobianZero.DirectDeterminant;
			var jacobianZeroInverse = jacobianZero.InverseMatrix;
			var transformationMat0 = CalculateTransformationMatrix(jacobianZeroInverse);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				Vector Stresses = Vector.CreateFromArray(lastStresses[gp]);
				var localdN = CalculateShapeFunctionsLocalDerivatives(
					QuadratureForStiffness.IntegrationPoints[gp].Coordinates);
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				var jacobianMat = jacobian.DirectMatrix;
				var detJ = jacobian.DirectDeterminant;
				var jacobianInverse = jacobian.InverseMatrix;
				var transformationMatrix = CalculateTransformationMatrix(jacobianInverse);
				var deformation = CalculateDeformationMatrix(localdN, jacobianMat, QuadratureForStiffness.IntegrationPoints[gp].Coordinates,
					transformationMatrix);
				var gamma = CalculateEnhancedStrainMatrixGamma(Matrix.CreateFromArray(transformationMat0), CreateEnhancedStrainsInterpolationMatrix(QuadratureForStiffness.IntegrationPoints[gp].Coordinates),
						detJ0, detJ);
				var gpForces = deformation.Transpose() * (Stresses);
				var dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
				gpForces.ScaleIntoThis(dA);
				vectorRe.AddIntoThis(gpForces);
				vectorPe.AddIntoThis((CalculateEnhancedStrainMatrixGamma(Matrix.CreateFromArray(transformationMat0), CreateEnhancedStrainsInterpolationMatrix(QuadratureForStiffness.IntegrationPoints[gp].Coordinates),
						detJ0, detJ).Transpose() * (Stresses)).Scale(dA));
			}

			internalForceVectorEAS = vectorPe.CopyToArray();
			var forces = vectorRe - internalEASLMatrix.Transpose() * internalEASDMatrix * vectorPe;
			return forces.CopyToArray();
		}

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
		{
			return CalculateResponseIntegral();
		}

		public Tuple<double[], double[]> CalculateResponse(double[] localDisplacements)
		{
			var jacobianZero = new IsoparametricJacobian3D(Nodes,
							Interpolation.EvaluateNaturalGradientsAt(new NaturalPoint(xi: 0d, eta: 0d, zeta: 0d)));
			var detJ0 = jacobianZero.DirectDeterminant;
			var jacobianZeroInverse = jacobianZero.InverseMatrix;
			var transformationMat0 = CalculateTransformationMatrix(jacobianZeroInverse);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			var alphaVector = Vector.CreateFromArray(internalParamVectorEAS);
			for (int gpo = 0; gpo < QuadratureForStiffness.IntegrationPoints.Count; ++gpo)
			{
				var localdN = CalculateShapeFunctionsLocalDerivatives(
					QuadratureForStiffness.IntegrationPoints[gpo].Coordinates);
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gpo]);
				var jacobianMat = jacobian.DirectMatrix;
				var detJ = jacobian.DirectDeterminant;
				var jacobianInverse = jacobian.InverseMatrix;
				var transformationMatrix = CalculateTransformationMatrix(jacobianInverse);
				var deformation = CalculateDeformationMatrix(localdN, jacobianMat, QuadratureForStiffness.IntegrationPoints[gpo].Coordinates,
					transformationMatrix);
				var Gamma = CalculateEnhancedStrainMatrixGamma(Matrix.CreateFromArray(transformationMat0), CreateEnhancedStrainsInterpolationMatrix(QuadratureForStiffness.IntegrationPoints[gpo].Coordinates),
						detJ0, detJ);
				var strainsVecANS = deformation.Multiply(localDisplacements);
				var EnhStrainVector = Gamma * alphaVector;
				strainsVec[gpo] = (Vector.CreateFromArray(strainsVecANS) + EnhStrainVector).CopyToArray();
				var strainsVecMinusLastConvergedValue = new double[6]
				{
					strainsVec[gpo][0] - strainsVecLastConverged[gpo][0],
					strainsVec[gpo][1] - strainsVecLastConverged[gpo][1],
					strainsVec[gpo][2] - strainsVecLastConverged[gpo][2],
					strainsVec[gpo][3] - strainsVecLastConverged[gpo][3],
					strainsVec[gpo][4] - strainsVecLastConverged[gpo][4],
					strainsVec[gpo][5] - strainsVecLastConverged[gpo][5]
				};
				lastStresses[gpo] = materialsAtGaussPoints[gpo].UpdateConstitutiveMatrixAndEvaluateResponse(strainsVecMinusLastConvergedValue);
			}

			UpdateElementEASParameters(localDisplacements);
			return new Tuple<double[], double[]>(strainsVec[materialsAtGaussPoints.Count - 1], lastStresses[materialsAtGaussPoints.Count - 1]);
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

		//public void ClearMaterialState()
		//{
		//	foreach (var material in materialsAtGaussPoints) material.ClearState();
		//}

		//public void ClearMaterialStresses()
		//{
		//	foreach (var material in materialsAtGaussPoints) material.ClearStresses();
		//}

		public IMatrix DampingMatrix()
		{
			IMatrix damping = BuildStiffnessMatrix();
			damping.ScaleIntoThis(dynamicProperties.RayleighCoeffStiffness);
			damping.AxpyIntoThis(MassMatrix(), dynamicProperties.RayleighCoeffMass);
			return damping;
		}

		/// <summary>
		/// Calculates the coordinates of the centroid of this element.
		/// </summary>
		public CartesianPoint FindCentroid()
			=> Interpolation.TransformNaturalToCartesian(Nodes, new NaturalPoint(0.0, 0.0, 0.0));

		public IMatrix MassMatrix()
		{
			return BuildLumpedMassMatrix();
		}

		//public void ResetConstitutiveLawModified()
		//{
		//	foreach (var material in materialsAtGaussPoints) material.ResetModified();
		//}

		public void SaveConstitutiveLawState()
		{
			for (int npoint = 0; npoint < materialsAtGaussPoints.Count; npoint++)
			{
				for (int i1 = 0; i1 < 6; i1++)
				{ strainsVecLastConverged[npoint][i1] = strainsVec[npoint][i1]; }
			}

			foreach (var m in materialsAtGaussPoints) m.CreateState();
		}

		public IMatrix StiffnessMatrix() => DofEnumerator.GetTransformedMatrix(BuildStiffnessMatrix());

		public IMatrix PhysicsMatrix()
		{
			return StiffnessMatrix();
		}

		public (IReadOnlyList<double[]> strains, IReadOnlyList<double[]> stresses)
			UpdateStrainStressesAtGaussPoints(double[] localDisplacements)
		{
			int numberOfGPs = QuadratureForStiffness.IntegrationPoints.Count;
			var strains = new double[numberOfGPs][];
			var stresses = new double[numberOfGPs][];
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			var jacobianZero = new IsoparametricJacobian3D(Nodes,
							Interpolation.EvaluateNaturalGradientsAt(new NaturalPoint(xi: 0d, eta: 0d, zeta: 0d)));
			var detJ0 = jacobianZero.DirectDeterminant;
			var jacobianZeroInverse = jacobianZero.InverseMatrix;
			var transformationMat0 = CalculateTransformationMatrix(jacobianZeroInverse);
			var alphaVector = Vector.CreateFromArray(internalParamVectorEAS);
			for (int gp = 0; gp < numberOfGPs; gp++)
			{
				IMatrixView constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
				var localdN = CalculateShapeFunctionsLocalDerivatives(
					QuadratureForStiffness.IntegrationPoints[gp].Coordinates);
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				var jacobianMat = jacobian.DirectMatrix;
				var detJ = jacobian.DirectDeterminant;
				var jacobianInverse = jacobian.InverseMatrix;
				var transformationMatrix = CalculateTransformationMatrix(jacobianInverse);
				var deformation = CalculateDeformationMatrix(localdN, jacobianMat, QuadratureForStiffness.IntegrationPoints[gp].Coordinates,
					transformationMatrix);
				var Gamma = CalculateEnhancedStrainMatrixGamma(Matrix.CreateFromArray(transformationMat0), CreateEnhancedStrainsInterpolationMatrix(QuadratureForStiffness.IntegrationPoints[gp].Coordinates),
						detJ0, detJ);
				var strainsVecANS = deformation.Multiply(localDisplacements);
				var EnhStrainVector = Gamma * alphaVector;

				strains[gp] = (Vector.CreateFromArray(strainsVecANS) + EnhStrainVector).CopyToArray();
				stresses[gp] = constitutive.Multiply(strains[gp]);
			}

			return (strains, stresses);
		}

		private void UpdateElementEASParameters(double[] totalU)
		{
			Vector deltaU = Vector.CreateFromArray(totalU) - displacementVectorPreviousIncrement;
			internalParamVectorEAS = (Vector.CreateFromArray(internalParamVectorEAS) - internalEASDMatrix * (internalEASLMatrix * deltaU + Vector.CreateFromArray(internalForceVectorEAS))).CopyToArray();
			displacementVectorPreviousIncrement = Vector.CreateFromArray(totalU);
		}

		private double[] UpdateNodalCoordinates(double[] displacementVector)
		{
			double[] updatedCoor = new double[24];
			for (int i = 1; i <= 8; i++)
			{
				updatedCoor[3 * i - 3] = Nodes[i].X + displacementVector[3 * i - 3];
				updatedCoor[3 * i - 2] = Nodes[i].Y + displacementVector[3 * i - 2];
				updatedCoor[3 * i - 1] = Nodes[i].Z + displacementVector[3 * i - 1];
			}

			return updatedCoor;
		}

		private double[] InitialNodalCoordinates()
		{
			double[] initialCoordinates = new double[24];
			for (int i = 0; i <= 7; i++)
			{
				initialCoordinates[3 * i] = Nodes[i].X;
				initialCoordinates[3 * i + 1] = Nodes[i].Y;
				initialCoordinates[3 * i + 2] = Nodes[i].Z;
			}

			return initialCoordinates;
		}

		private double[,] CalculateTransformationMatrix(Matrix jacobianInverseMatrix)
		{
			double[,] T = new double[6, 6];
			T[0, 0] = jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[0, 0];
			T[0, 1] = jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[0, 1];
			T[0, 2] = jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[0, 2];
			T[0, 3] = jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[0, 1];
			T[0, 4] = jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[0, 2];
			T[0, 5] = jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[0, 0];

			T[1, 0] = jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[1, 0];
			T[1, 1] = jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[1, 1];
			T[1, 2] = jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[1, 2];
			T[1, 3] = jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[1, 1];
			T[1, 4] = jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[1, 2];
			T[1, 5] = jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[1, 0];

			T[2, 0] = jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[2, 0];
			T[2, 1] = jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[2, 1];
			T[2, 2] = jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[2, 2];
			T[2, 3] = jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[2, 1];
			T[2, 4] = jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[2, 2];
			T[2, 5] = jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[2, 0];

			T[3, 0] = 2 * jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[1, 0];
			T[3, 1] = 2 * jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[1, 1];
			T[3, 2] = 2 * jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[1, 2];

			T[3, 3] = jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[1, 1] +
								jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[1, 0];

			T[3, 4] = jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[1, 2] +
								jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[1, 1];

			T[3, 5] = jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[1, 0] +
								jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[2, 1];

			T[4, 0] = 2 * jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[2, 0];
			T[4, 1] = 2 * jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[2, 1];
			T[4, 2] = 2 * jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[2, 2];

			T[4, 3] = jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[2, 1] +
								jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[2, 0];

			T[4, 4] = jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[2, 2] +
								jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[2, 1];

			T[4, 5] = jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[2, 0] +
								jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[2, 2];

			T[5, 0] = 2 * jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[0, 0];
			T[5, 1] = 2 * jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[0, 1];
			T[5, 2] = 2 * jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[0, 2];

			T[5, 3] = jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[0, 1] +
								jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[0, 0];

			T[5, 4] = jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[0, 2] +
								jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[0, 1];

			T[5, 5] = jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[0, 0] +
								jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[0, 2];
			return T;
		}

		private Dictionary<string, double[]> CalculateShapeFunctionsLocalDerivatives(double[] naturalCoordinates)
		{
			double ksi = naturalCoordinates[0];
			double ihta = naturalCoordinates[1];
			double mhi = naturalCoordinates[2];

			double[] dN_ksi = new double[]
			{
				(1.0/8.0*(1+ihta)*(1+mhi)),
				(-1.0/8.0*(1+ihta)*(1+mhi)),
				(-1.0/8.0*(1-ihta)*(1+mhi)),
				(1.0/8.0*(1-ihta)*(1+mhi)),
				(1.0/8.0*(1+ihta)*(1-mhi)),
				(-1.0/8.0*(1+ihta)*(1-mhi)),
				(-1.0/8.0*(1-ihta)*(1-mhi)),
				(1.0/8.0*(1-ihta)*(1-mhi))
			};

			double[] dN_ihta = new double[]
			{
				(1.0/8.0*(1+ksi)*(1+mhi)),
				(1.0/8.0*(1-ksi)*(1+mhi)),
				(-1.0/8.0*(1-ksi)*(1+mhi)),
				(-1.0/8.0*(1+ksi)*(1+mhi)),
				(1.0/8.0*(1+ksi)*(1-mhi)),
				(1.0/8.0*(1-ksi)*(1-mhi)),
				(-1.0/8.0*(1-ksi)*(1-mhi)),
				(-1.0/8.0*(1+ksi)*(1-mhi))
			};

			double[] dN_mhi = new double[]
			{
				(1.0/8.0*(1+ksi)*(1+ihta)),
				(1.0/8.0*(1-ksi)*(1+ihta)),
				(1.0/8.0*(1-ksi)*(1-ihta)),
				(1.0/8.0*(1+ksi)*(1-ihta)),
				(-1.0/8.0*(1+ksi)*(1+ihta)),
				(-1.0/8.0*(1-ksi)*(1+ihta)),
				(-1.0/8.0*(1-ksi)*(1-ihta)),
				(-1.0/8.0*(1+ksi)*(1-ihta))
			};

			Dictionary<string, double[]> dN = new Dictionary<string, double[]>();
			dN.Add("ksi", dN_ksi);
			dN.Add("ihta", dN_ihta);
			dN.Add("mhi", dN_mhi);
			return dN;
		}

		private Matrix CreateEnhancedStrainsInterpolationMatrix(double[] ksi)
		{
			double[,] M = new double[6, 7];
			M[0, 0] = ksi[0];
			M[1, 1] = ksi[1];
			M[2, 2] = ksi[2];
			M[2, 5] = ksi[0] * ksi[2];
			M[2, 6] = ksi[1] * ksi[2];
			M[3, 3] = ksi[0];
			M[3, 4] = ksi[1];
			var result = Matrix.CreateFromArray(M);
			return result;
		}

		private Matrix CalculateEnhancedStrainMatrixGamma(Matrix transformationMat0, Matrix M,
		double detJ0, double detJ)
		{
			var scalar = detJ0 / detJ;
			var gamma = (transformationMat0 * M).Scale(scalar);
			return gamma;
		}

		private double[] CalculateJacobianPart1(double[] parametricCoordinates)
		{
			double[] jacobianPart1 = new double[3];
			double[] xNodalInitial = InitialNodalCoordinates();

			int k = 0;
			for (int i = 0; i < 8; i++)
			{
				jacobianPart1[0] = jacobianPart1[0] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 0, i);
				k = k + 3;
			}

			k = 1;
			for (int i = 0; i < 8; i++)
			{
				jacobianPart1[1] = jacobianPart1[1] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 0, i);
				k = k + 3;
			}

			k = 2;
			for (int i = 0; i < 8; i++)
			{
				jacobianPart1[2] = jacobianPart1[2] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 0, i);
				k = k + 3;
			}

			return jacobianPart1;
		}

		private double[] CalculateJacobianPart2(double[] parametricCoordinates)
		{
			double[] jacobianPart2 = new double[3];
			double[] xNodalInitial = InitialNodalCoordinates();

			int k = 0;
			for (int i = 0; i < 8; i++)
			{
				jacobianPart2[0] = jacobianPart2[0] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 1, i);
				k = k + 3;
			}

			k = 1;
			for (int i = 0; i < 8; i++)
			{
				jacobianPart2[1] = jacobianPart2[1] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 1, i);
				k = k + 3;
			}

			k = 2;
			for (int i = 0; i < 8; i++)
			{
				jacobianPart2[2] = jacobianPart2[2] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 1, i);
				k = k + 3;
			}

			return jacobianPart2;
		}

		private double[] CalculateJacobianPart3(double[] parametricCoordinates)
		{
			double[] jacobianPart3 = new double[3];
			double[] xNodalInitial = InitialNodalCoordinates();

			int k = 0;
			for (int i = 0; i < 8; i++)
			{
				jacobianPart3[0] = jacobianPart3[0] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 2, i);
				k = k + 3;
			}

			k = 1;
			for (int i = 0; i < 8; i++)
			{
				jacobianPart3[1] = jacobianPart3[1] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 2, i);
				k = k + 3;
			}

			k = 2;
			for (int i = 0; i < 8; i++)
			{
				jacobianPart3[2] = jacobianPart3[2] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 2, i);
				k = k + 3;
			}

			return jacobianPart3;
		}

		private Matrix CalculateDeformationMatrix(Dictionary<string, double[]> dNlocal, Matrix jacobianMatrix,
			double[] parametricCoordinates, double[,] transformationMatrix)
		{
			double[,] BmatrixLocal = new double[6, 24];
			//-----------------------------------------------------------------
			double[] a = new double[] { 1d, 1d, 0d };
			double[] b = new double[] { -1d, 1d, 0d };
			double[] c = new double[] { -1d, -1d, 0d };
			double[] d = new double[] { 1d, -1d, 0d };
			//-----------------------------------------------------------------
			double[] e = new double[] { 1d, 0d, 1d };
			double[] f = new double[] { -1d, 0d, 1d };
			double[] g = new double[] { -1d, 0d, -1d };
			double[] h = new double[] { 1d, 0d, -1d };
			//-----------------------------------------------------------------
			double[] j = new double[] { 0d, 1d, 1d };
			double[] k = new double[] { 0d, -1d, 1d };
			double[] l = new double[] { 0d, -1d, -1d };
			double[] m = new double[] { 0d, 1d, -1d };
			for (int i = 0; i < 8; i++)
			{
				BmatrixLocal[0, i * 3] = dNlocal["ksi"][i] * jacobianMatrix[0, 0];
				BmatrixLocal[0, i * 3 + 1] = dNlocal["ksi"][i] * jacobianMatrix[0, 1];
				BmatrixLocal[0, i * 3 + 2] = dNlocal["ksi"][i] * jacobianMatrix[0, 2];

				BmatrixLocal[1, i * 3] = dNlocal["ihta"][i] * jacobianMatrix[1, 0];
				BmatrixLocal[1, i * 3 + 1] = dNlocal["ihta"][i] * jacobianMatrix[1, 1];
				BmatrixLocal[1, i * 3 + 2] = dNlocal["ihta"][i] * jacobianMatrix[1, 2];

				//εζζ Assumed strain interpolation from points A, B, C, D
				BmatrixLocal[2, i * 3] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(a, 2, i) * CalculateJacobianPart3(a)[0] +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(b, 2, i) * CalculateJacobianPart3(b)[0] +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(c, 2, i) * CalculateJacobianPart3(c)[0] +
								(1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(d, 2, i) * CalculateJacobianPart3(d)[0];
				BmatrixLocal[2, i * 3 + 1] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(a, 2, i) * CalculateJacobianPart3(a)[1] +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(b, 2, i) * CalculateJacobianPart3(b)[1] +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(c, 2, i) * CalculateJacobianPart3(c)[1] +
								(1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(d, 2, i) * CalculateJacobianPart3(d)[1];
				BmatrixLocal[2, i * 3 + 2] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(a, 2, i) * CalculateJacobianPart3(a)[2] +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(b, 2, i) * CalculateJacobianPart3(b)[2] +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(c, 2, i) * CalculateJacobianPart3(c)[2] +
								(1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
								CalculateShapeFunctionsLocalDerivatives(d, 2, i) * CalculateJacobianPart3(d)[2];
				//-----------------------------------------------------------------
				BmatrixLocal[3, i * 3] = dNlocal["ihta"][i] * jacobianMatrix[0, 0] +
					dNlocal["ksi"][i] * jacobianMatrix[1, 0];
				BmatrixLocal[3, i * 3 + 1] = dNlocal["ihta"][i] * jacobianMatrix[0, 1] +
					dNlocal["ksi"][i] * jacobianMatrix[1, 1];
				BmatrixLocal[3, i * 3 + 2] = dNlocal["ihta"][i] * jacobianMatrix[0, 2] +
					dNlocal["ksi"][i] * jacobianMatrix[1, 2];
				//εηζ Assumed strain interpolation from points E,F,G,H
				BmatrixLocal[4, i * 3] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(e, 2, i) * CalculateJacobianPart2(e)[0] +
								CalculateShapeFunctionsLocalDerivatives(e, 1, i) * CalculateJacobianPart3(e)[0]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(f, 2, i) * CalculateJacobianPart2(f)[0] +
								CalculateShapeFunctionsLocalDerivatives(f, 1, i) * CalculateJacobianPart3(f)[0]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(g, 2, i) * CalculateJacobianPart2(g)[0] +
								CalculateShapeFunctionsLocalDerivatives(g, 1, i) * CalculateJacobianPart3(g)[0]) +
								(1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(h, 2, i) * CalculateJacobianPart2(h)[0] +
								CalculateShapeFunctionsLocalDerivatives(h, 1, i) * CalculateJacobianPart3(h)[0]);

				BmatrixLocal[4, i * 3 + 1] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(e, 2, i) * CalculateJacobianPart2(e)[1] +
								CalculateShapeFunctionsLocalDerivatives(e, 1, i) * CalculateJacobianPart3(e)[1]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(f, 2, i) * CalculateJacobianPart2(f)[1] +
								CalculateShapeFunctionsLocalDerivatives(f, 1, i) * CalculateJacobianPart3(f)[1]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(g, 2, i) * CalculateJacobianPart2(g)[1] +
								CalculateShapeFunctionsLocalDerivatives(g, 1, i) * CalculateJacobianPart3(g)[1]) +
								(1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(h, 2, i) * CalculateJacobianPart2(h)[1] +
								CalculateShapeFunctionsLocalDerivatives(h, 1, i) * CalculateJacobianPart3(h)[1]);

				BmatrixLocal[4, i * 3 + 2] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(e, 2, i) * CalculateJacobianPart2(e)[2] +
								CalculateShapeFunctionsLocalDerivatives(e, 1, i) * CalculateJacobianPart3(e)[2]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(f, 2, i) * CalculateJacobianPart2(f)[2] +
								CalculateShapeFunctionsLocalDerivatives(f, 1, i) * CalculateJacobianPart3(f)[2]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(g, 2, i) * CalculateJacobianPart2(g)[2] +
								CalculateShapeFunctionsLocalDerivatives(g, 1, i) * CalculateJacobianPart3(g)[2]) +
								(1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(h, 2, i) * CalculateJacobianPart2(h)[2] +
								CalculateShapeFunctionsLocalDerivatives(h, 1, i) * CalculateJacobianPart3(h)[2]);
				//-----------------------------------------------------------------
				//εζξ Assumed strain interpolation from points J, K, L, M
				BmatrixLocal[5, i * 3] = (1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(j, 0, i) * CalculateJacobianPart3(j)[0] +
								CalculateShapeFunctionsLocalDerivatives(j, 2, i) * CalculateJacobianPart1(j)[0]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(k, 0, i) * CalculateJacobianPart3(k)[0] +
								CalculateShapeFunctionsLocalDerivatives(k, 2, i) * CalculateJacobianPart1(k)[0]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(l, 0, i) * CalculateJacobianPart3(l)[0] +
								CalculateShapeFunctionsLocalDerivatives(l, 2, i) * CalculateJacobianPart1(l)[0]) +
								(1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(m, 0, i) * CalculateJacobianPart3(m)[0] +
								CalculateShapeFunctionsLocalDerivatives(m, 2, i) * CalculateJacobianPart1(m)[0]);

				BmatrixLocal[5, i * 3 + 1] = (1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(j, 0, i) * CalculateJacobianPart3(j)[1] +
								CalculateShapeFunctionsLocalDerivatives(j, 2, i) * CalculateJacobianPart1(j)[1]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(k, 0, i) * CalculateJacobianPart3(k)[1] +
								CalculateShapeFunctionsLocalDerivatives(k, 2, i) * CalculateJacobianPart1(k)[1]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(l, 0, i) * CalculateJacobianPart3(l)[1] +
								CalculateShapeFunctionsLocalDerivatives(l, 2, i) * CalculateJacobianPart1(l)[1]) +
								(1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(m, 0, i) * CalculateJacobianPart3(m)[1] +
								CalculateShapeFunctionsLocalDerivatives(m, 2, i) * CalculateJacobianPart1(m)[1]);

				BmatrixLocal[5, i * 3 + 2] = (1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(j, 0, i) * CalculateJacobianPart3(j)[2] +
								CalculateShapeFunctionsLocalDerivatives(j, 2, i) * CalculateJacobianPart1(j)[2]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(k, 0, i) * CalculateJacobianPart3(k)[2] +
								CalculateShapeFunctionsLocalDerivatives(k, 2, i) * CalculateJacobianPart1(k)[2]) +
								(1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(l, 0, i) * CalculateJacobianPart3(l)[2] +
								CalculateShapeFunctionsLocalDerivatives(l, 2, i) * CalculateJacobianPart1(l)[2]) +
								(1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
								(CalculateShapeFunctionsLocalDerivatives(m, 0, i) * CalculateJacobianPart3(m)[2] +
								CalculateShapeFunctionsLocalDerivatives(m, 2, i) * CalculateJacobianPart1(m)[2]);
				//-----------------------------------------------------------------
			}

			//var BmatrixGlobal = MatrixOperations.MatrixProduct(transformationMatrix, BmatrixLocal);
			var BmatrixGlobal = Matrix.CreateFromArray(transformationMatrix) * Matrix.CreateFromArray(BmatrixLocal);
			return BmatrixGlobal;
		}

		private double CalculateShapeFunctionsLocalDerivatives(double[] naturalCoordinates, int parametricAxis,
			int node)
		{
			double ksi = naturalCoordinates[0];
			double ihta = naturalCoordinates[1];
			double mhi = naturalCoordinates[2];
			if (parametricAxis == 0)
			{
				double dN_ksi = new double();
				switch (node)
				{
					case 0:
						dN_ksi = 1.0 / 8.0 * (1 + ihta) * (1 + mhi);
						break;
					case 1:
						dN_ksi = -1.0 / 8.0 * (1 + ihta) * (1 + mhi);
						break;
					case 2:
						dN_ksi = -1.0 / 8.0 * (1 - ihta) * (1 + mhi);
						break;
					case 3:
						dN_ksi = 1.0 / 8.0 * (1 - ihta) * (1 + mhi);
						break;
					case 4:
						dN_ksi = 1.0 / 8.0 * (1 + ihta) * (1 - mhi);
						break;
					case 5:
						dN_ksi = -1.0 / 8.0 * (1 + ihta) * (1 - mhi);
						break;
					case 6:
						dN_ksi = -1.0 / 8.0 * (1 - ihta) * (1 - mhi);
						break;
					case 7:
						dN_ksi = 1.0 / 8.0 * (1 - ihta) * (1 - mhi);
						break;
				}
				return dN_ksi;
			}
			else if (parametricAxis == 1)
			{
				double dN_ihta = new double();
				switch (node)
				{
					case 0:
						dN_ihta = 1.0 / 8.0 * (1 + ksi) * (1 + mhi);
						break;
					case 1:
						dN_ihta = 1.0 / 8.0 * (1 - ksi) * (1 + mhi);
						break;
					case 2:
						dN_ihta = -1.0 / 8.0 * (1 - ksi) * (1 + mhi);
						break;
					case 3:
						dN_ihta = -1.0 / 8.0 * (1 + ksi) * (1 + mhi);
						break;
					case 4:
						dN_ihta = 1.0 / 8.0 * (1 + ksi) * (1 - mhi);
						break;
					case 5:
						dN_ihta = 1.0 / 8.0 * (1 - ksi) * (1 - mhi);
						break;
					case 6:
						dN_ihta = -1.0 / 8.0 * (1 - ksi) * (1 - mhi);
						break;
					case 7:
						dN_ihta = -1.0 / 8.0 * (1 + ksi) * (1 - mhi);
						break;
				}
				return dN_ihta;
			}
			else
			{

				double dN_mhi = new double();
				switch (node)
				{
					case 0:
						dN_mhi = 1.0 / 8.0 * (1 + ksi) * (1 + ihta);
						break;
					case 1:
						dN_mhi = 1.0 / 8.0 * (1 - ksi) * (1 + ihta);
						break;
					case 2:
						dN_mhi = 1.0 / 8.0 * (1 - ksi) * (1 - ihta);
						break;
					case 3:
						dN_mhi = 1.0 / 8.0 * (1 + ksi) * (1 - ihta);
						break;
					case 4:
						dN_mhi = -1.0 / 8.0 * (1 + ksi) * (1 + ihta);
						break;
					case 5:
						dN_mhi = -1.0 / 8.0 * (1 - ksi) * (1 + ihta);
						break;
					case 6:
						dN_mhi = -1.0 / 8.0 * (1 - ksi) * (1 - ihta);
						break;
					case 7:
						dN_mhi = -1.0 / 8.0 * (1 + ksi) * (1 - ihta);
						break;
				}
				return dN_mhi;
			}
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
