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
using MGroup.MSolve.Discretization.Embedding;
using System.Linq;

namespace MGroup.FEM.Structural.Continuum
{
	/// <summary>
	/// Represents a continuum finite element for 3D problems. Specific elements (e.g. Hexa8, Hexa20, ...) can be created using
	/// the appropriate <see cref="IIsoparametricInterpolation3D_OLD"/>, <see cref="IQuadrature3D"/> etc. strategies. 
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class ContinuumElement3D : IStructuralElementType, ICell<INode>, IEmbeddedHostElement
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

		public ContinuumElement3D(IReadOnlyList<INode> nodes, IIsoparametricInterpolation3D interpolation,
			IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
			IGaussPointExtrapolation3D gaussPointExtrapolation,
			IReadOnlyList<IContinuumMaterial3D> materialsAtGaussPoints, ITransientAnalysisProperties dynamicProperties)
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
			lastStresses = new double[materialsAtGaussPoints.Count][];
			for (int gpoint = 0; gpoint < materialsAtGaussPoints.Count; gpoint++)
			{
				strainsVec[gpoint] = new double[6];
				strainsVecLastConverged[gpoint] = new double[6];
				lastStresses[gpoint] = new double[6];
			}
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

		//public bool ConstitutiveLawModified
		//{
		//	get
		//	{
		//		foreach (IContinuumMaterial3D material in materialsAtGaussPoints)
		//		{
		//			if (material.IsCurrentStateDifferent()) return true;
		//		}
		//		return false;
		//	}
		//}

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

		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads)
		//{
		//	int numberOfDofs = 3 * Nodes.Count;
		//	var accelerations = new double[numberOfDofs];
		//	IMatrix massMatrix = MassMatrix(element);

		//	foreach (var load in loads)
		//	{
		//		int index = 0;
		//		foreach (var nodalDOFTypes in dofTypes)
		//		{
		//			foreach (var dofType in nodalDOFTypes)
		//			{
		//				if (dofType == load.DOF) accelerations[index] += load.Amount;
		//				index++;
		//			}
		//		}
		//	}

		//	return massMatrix.Multiply(accelerations);
		//}

		public double[] CalculateResponseIntegral()
		{
			int numberOfDofs = 3 * Nodes.Count;
			var Forces = Vector.CreateZero(numberOfDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				Vector Stresses = Vector.CreateFromArray(lastStresses[gp]);
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

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
		{
			return CalculateResponseIntegral();
		}

		public Tuple<double[], double[]> CalculateResponse(double[] localDisplacements)
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
				lastStresses[gpo] = materialsAtGaussPoints[gpo].UpdateConstitutiveMatrixAndEvaluateResponse(strainsVecMinusLastConvergedValue);
				//To update with total strain simplY = materialsAtGaussPoints[npoint].UpdateMaterial(strainsVec[npoint]);
			}


			return new Tuple<double[], double[]>(strainsVec[materialsAtGaussPoints.Count-1], lastStresses[materialsAtGaussPoints.Count - 1]);
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

		public EmbeddedNode BuildHostElementEmbeddedNode(IElementType element, INode node,
			IEmbeddedDOFInHostTransformationVector transformation)
		{
			var points = GetNaturalCoordinates(element, (Node)node);
			if (points.Length == 0) return null;

			//element.EmbeddedNodes.Add(node);
			var embeddedNode = new EmbeddedNode(node, element, transformation.GetDependentDOFTypes);
			for (int i = 0; i < points.Length; i++) embeddedNode.Coordinates.Add(points[i]);
			return embeddedNode;
		}

		private double[] GetNaturalCoordinates(IElementType element, Node node)
		{
			double[] mins = new double[] { element.Nodes[0].X, element.Nodes[0].Y, element.Nodes[0].Z };
			double[] maxes = new double[] { element.Nodes[0].X, element.Nodes[0].Y, element.Nodes[0].Z };
			for (int i = 0; i < element.Nodes.Count; i++)
			{
				mins[0] = mins[0] > element.Nodes[i].X ? element.Nodes[i].X : mins[0];
				mins[1] = mins[1] > element.Nodes[i].Y ? element.Nodes[i].Y : mins[1];
				mins[2] = mins[2] > element.Nodes[i].Z ? element.Nodes[i].Z : mins[2];
				maxes[0] = maxes[0] < element.Nodes[i].X ? element.Nodes[i].X : maxes[0];
				maxes[1] = maxes[1] < element.Nodes[i].Y ? element.Nodes[i].Y : maxes[1];
				maxes[2] = maxes[2] < element.Nodes[i].Z ? element.Nodes[i].Z : maxes[2];
			}
			//return new double[] { (node.X - mins[0]) / ((maxes[0] - mins[0]) / 2) - 1,
			//    (node.Y - mins[1]) / ((maxes[1] - mins[1]) / 2) - 1,
			//    (node.Z - mins[2]) / ((maxes[2] - mins[2]) / 2) - 1 };

			bool maybeInsideElement = node.X <= maxes[0] && node.X >= mins[0] &&
				node.Y <= maxes[1] && node.Y >= mins[1] &&
				node.Z <= maxes[2] && node.Z >= mins[2];
			if (maybeInsideElement == false) return new double[0];

			const int jacobianSize = 3;
			const int maxIterations = 1000;
			const double tolerance = 1e-10;
			int iterations = 0;
			double deltaNaturalCoordinatesNormSquare = 100;
			double[] naturalCoordinates = new double[] { 0, 0, 0 };
			const double toleranceSquare = tolerance * tolerance;

			while (deltaNaturalCoordinatesNormSquare > toleranceSquare && iterations < maxIterations)
			{
				iterations++;
				var shapeFunctions = Interpolation.EvaluateFunctionsAt(new NaturalPoint(naturalCoordinates));
				double[] coordinateDifferences = new double[] { 0, 0, 0 };
				for (int i = 0; i < shapeFunctions.Length; i++)
				{
					coordinateDifferences[0] += shapeFunctions[i] * element.Nodes[i].X;
					coordinateDifferences[1] += shapeFunctions[i] * element.Nodes[i].Y;
					coordinateDifferences[2] += shapeFunctions[i] * element.Nodes[i].Z;
				}
				coordinateDifferences[0] = node.X - coordinateDifferences[0];
				coordinateDifferences[1] = node.Y - coordinateDifferences[1];
				coordinateDifferences[2] = node.Z - coordinateDifferences[2];

				double[,] faXYZ = GetCoordinatesTranspose(element);
				var nablaShapeFunctions = Interpolation.EvaluateNaturalGradientsAt(new NaturalPoint(naturalCoordinates));
				var inverseJacobian = new IsoparametricJacobian3D(element.Nodes.ToArray(), nablaShapeFunctions).InverseMatrix;

				double[] deltaNaturalCoordinates = new double[] { 0, 0, 0 };
				for (int i = 0; i < jacobianSize; i++)
					for (int j = 0; j < jacobianSize; j++)
						deltaNaturalCoordinates[i] += inverseJacobian[j, i] * coordinateDifferences[j];
				for (int i = 0; i < 3; i++)
					naturalCoordinates[i] += deltaNaturalCoordinates[i];

				deltaNaturalCoordinatesNormSquare = 0;
				for (int i = 0; i < 3; i++)
					deltaNaturalCoordinatesNormSquare += deltaNaturalCoordinates[i] * deltaNaturalCoordinates[i];
				//deltaNaturalCoordinatesNormSquare = Math.Sqrt(deltaNaturalCoordinatesNormSquare);
			}

			return naturalCoordinates.Count(x => Math.Abs(x) - 1.0 > tolerance) > 0 ? new double[0] : naturalCoordinates;
		}

		protected double[,] GetCoordinatesTranspose(IElementType element)
		{
			double[,] nodeCoordinatesXYZ = new double[3, dofTypes.Length];
			for (int i = 0; i < dofTypes.Length; i++)
			{
				nodeCoordinatesXYZ[0, i] = element.Nodes[i].X;
				nodeCoordinatesXYZ[1, i] = element.Nodes[i].Y;
				nodeCoordinatesXYZ[2, i] = element.Nodes[i].Z;
			}
			return nodeCoordinatesXYZ;
		}

		double[] IEmbeddedHostElement.GetShapeFunctionsForNode(IElementType element, EmbeddedNode node)
		{
			var naturalPoint = new NaturalPoint(node.Coordinates.ToArray());
			var shapeFunctions = Interpolation.EvaluateFunctionsAt(naturalPoint);
			var shapeFunctionGradients = Interpolation.EvaluateNaturalGradientsAt(naturalPoint);
			var jacobian = new IsoparametricJacobian3D(element.Nodes.ToArray(), shapeFunctionGradients);

			var shapeFunctionGradients2DArray = shapeFunctionGradients.CopyToArray2D();
			var shapeFunctionGradients1DArray = new double[shapeFunctionGradients2DArray.GetLength(0) * shapeFunctionGradients2DArray.GetLength(1)];
			for (int j = 0; j < shapeFunctionGradients2DArray.GetLength(1); j++)
			{
				for (int i = 0; i < shapeFunctionGradients2DArray.GetLength(0); i++)
				{
					shapeFunctionGradients1DArray[(shapeFunctionGradients2DArray.GetLength(0) * j) + i] = shapeFunctionGradients2DArray[i, j];
				}
			}

			var jacobianDirect2DArray = jacobian.DirectMatrix.CopyToArray2D();
			var jacobianInverse2DArray = jacobian.InverseMatrix.CopyToArray2D();
			var jacobianDirect1DArray = new double[jacobianDirect2DArray.GetLength(0) * jacobianDirect2DArray.GetLength(1)];
			var jacobianInverse1DArray = new double[jacobianInverse2DArray.GetLength(0) * jacobianInverse2DArray.GetLength(1)];
			for (int i = 0; i < jacobianDirect2DArray.GetLength(0); i++)
			{
				for (int j = 0; j < jacobianDirect2DArray.GetLength(1); j++)
				{
					jacobianDirect1DArray[(jacobianDirect2DArray.GetLength(1) * i) + j] = jacobianDirect2DArray[i, j];
					jacobianInverse1DArray[(jacobianInverse2DArray.GetLength(1) * i) + j] = jacobianInverse2DArray[i, j];
				}
			}
			var returnValueList = new List<double>();
			foreach (double shapeFunction in shapeFunctions)
			{
				returnValueList.Add(shapeFunction);
			}
			foreach (double shapeFunctionGradient in shapeFunctionGradients1DArray)
			{
				returnValueList.Add(shapeFunctionGradient);
			}
			foreach (double value in jacobianDirect1DArray)
			{
				returnValueList.Add(value);
			}
			foreach (double value in jacobianInverse1DArray)
			{
				returnValueList.Add(value);
			}

			return returnValueList.ToArray();
		}
	}
}
