using System;
using System.Collections.Generic;
using MGroup.FEM.Elements.SupportiveClasses;
using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.Geometry;

namespace MGroup.FEM.Structural.Elements
{
	public class Beam3DCorotationalQuaternion : Beam3DCorotationalAbstract
	{
		private readonly double[] lastDisplacements;
		private readonly double[] currentDisplacements;
		private readonly double[] displacementsOfCurrentIncrement;
		private readonly double[] initialAxisY;
		private readonly double[] initialAxisZ;
		private readonly double[] initialAxisX;
		private Quaternion quaternionCurrentNodeA;
		private Quaternion quaternionCurrentNodeB;
		private readonly double[] currentBeamAxis;
		private Quaternion quaternionLastNodeA;
		private Quaternion quaternionLastNodeB;

		/**
         * Creates a new instance of {@link Beam3DCorotationalIncremental} class.
         *
         * @param nodes
         *            The element nodes
         * @param material
         *            The element material
         * @param density
         *            The element density
         * @param beamSection
         *            The beam section.
         */
		public Beam3DCorotationalQuaternion(IList<Node> nodes, IIsotropicContinuumMaterial3D material, double density,
			BeamSection3D beamSection)
			: base(nodes, material, density, beamSection)
		{
			this.displacementsOfCurrentIncrement = new double[FREEDOM_DEGREE_COUNT];
			this.lastDisplacements = new double[FREEDOM_DEGREE_COUNT];
			this.currentDisplacements = new double[FREEDOM_DEGREE_COUNT];
			this.quaternionLastNodeA = Quaternion.OfZeroAngle();
			this.quaternionLastNodeB = Quaternion.OfZeroAngle();
			this.quaternionCurrentNodeA = Quaternion.OfZeroAngle();
			this.quaternionCurrentNodeB = Quaternion.OfZeroAngle();
			this.initialAxisX = new double[AXIS_COUNT];
			this.initialAxisY = new double[AXIS_COUNT];
			this.initialAxisZ = new double[AXIS_COUNT];
			this.currentBeamAxis = new double[AXIS_COUNT];
			this.InitializeElementAxes();
		}

		public override void SaveGeometryState()
		{
			displacementsOfCurrentIncrement.Clear();
			//displacementsOfCurrentIncrement.ScaleIntoThis(0d); 
			lastDisplacements.CopyFrom(currentDisplacements);

			double[] qA = quaternionCurrentNodeA.VectorPart.Copy();
			double[] qB = quaternionCurrentNodeB.VectorPart.Copy();
			this.quaternionLastNodeA = new Quaternion(quaternionCurrentNodeA.ScalarPart, qA);
			this.quaternionLastNodeB = new Quaternion(quaternionCurrentNodeB.ScalarPart, qB);
		}

		public override void UpdateState(double[] incrementalNodeDisplacements)
		{
			//displacementsOfCurrentIncrement.Add(new Vector(incrementalNodeDisplacements));
			//lastDisplacements.CopyTo(currentDisplacements.Data, 0);
			//currentDisplacements.Add(displacementsOfCurrentIncrement);

			//this.currentDisplacements.addIntoThis(this.lastDisplacements, incrementalNodeDisplacements);
			currentDisplacements.CopyFrom(lastDisplacements);
			currentDisplacements.AddIntoThis(incrementalNodeDisplacements);
			displacementsOfCurrentIncrement.CopyFrom(incrementalNodeDisplacements);

			var incrementalRotationsA = new double[AXIS_COUNT];
			var incrementalRotationsB = new double[AXIS_COUNT];

			incrementalRotationsA[0] = displacementsOfCurrentIncrement[3];
			incrementalRotationsA[1] = displacementsOfCurrentIncrement[4];
			incrementalRotationsA[2] = displacementsOfCurrentIncrement[5];
			incrementalRotationsB[0] = displacementsOfCurrentIncrement[9];
			incrementalRotationsB[1] = displacementsOfCurrentIncrement[10];
			incrementalRotationsB[2] = displacementsOfCurrentIncrement[11];

			double[] qA = quaternionLastNodeA.VectorPart.Copy();
			double[] qB = quaternionLastNodeB.VectorPart.Copy();
			this.quaternionCurrentNodeA = new Quaternion(quaternionLastNodeA.ScalarPart, qA);
			this.quaternionCurrentNodeB = new Quaternion(quaternionLastNodeB.ScalarPart, qB);
			this.quaternionCurrentNodeA.ApplyIncrementalRotation(incrementalRotationsA);
			this.quaternionCurrentNodeB.ApplyIncrementalRotation(incrementalRotationsB);

			double scalarA = this.quaternionCurrentNodeA.ScalarPart;
			double scalarB = this.quaternionCurrentNodeB.ScalarPart;
			var vectorPartA = this.quaternionCurrentNodeA.VectorPart;
			var vectorPartB = this.quaternionCurrentNodeB.VectorPart;
			double[] sumOfVectorParts = vectorPartA.Copy();
			sumOfVectorParts.AddIntoThis(vectorPartB);
			double sumOfVectorPartsNorm = sumOfVectorParts.Norm2();
			double scalarPartDifference = 0.5 * Math.Sqrt(((scalarA + scalarB) * (scalarA + scalarB)) + (sumOfVectorPartsNorm * sumOfVectorPartsNorm));
			double meanRotationScalarPart = (0.5 * (scalarA + scalarB)) / scalarPartDifference;
			double[] meanRotationVectorPart = vectorPartA.Copy();
			meanRotationVectorPart.AddIntoThis(vectorPartB);
			meanRotationVectorPart.ScaleIntoThis(1d / (2d * scalarPartDifference));
			double[] vectorPartDifference = vectorPartB.Copy();
			vectorPartDifference.ScaleIntoThis(scalarA);
			vectorPartDifference.AddIntoThis(vectorPartA.Scale(-scalarB));
			//vectorPartDifference.doPointwise(vectorPartA, DoubleBinaryOps.alphaPlusScaledBeta(-scalarB));
			vectorPartDifference.AddIntoThis(vectorPartA.CrossProduct(vectorPartB));
			vectorPartDifference.ScaleIntoThis(1d / (2d * scalarPartDifference));
			var meanRotationQuaternion = Quaternion.CreateFromIndividualParts(meanRotationScalarPart, meanRotationVectorPart);
			this.currentRotationMatrix = meanRotationQuaternion.GetRotationMatrix();
			this.CalculateUpdatedBeamAxis();
			this.UpdateRotationMatrix();
			this.UpdateNaturalDeformations(vectorPartDifference);

			//SaveGeometryState();
		}

		private void CalculateUpdatedBeamAxis()
		{
			currentRotationMatrix.MultiplyIntoResult(initialAxisX, beamAxisX);
			currentRotationMatrix.MultiplyIntoResult(initialAxisY, beamAxisY);
			currentRotationMatrix.MultiplyIntoResult(initialAxisZ, beamAxisZ);

			double dX = ((nodes[1].X - nodes[0].X) + this.currentDisplacements[6]) - this.currentDisplacements[0];
			double dY = ((nodes[1].Y - nodes[0].Y) + this.currentDisplacements[7]) - this.currentDisplacements[1];
			double dZ = ((nodes[1].Z - nodes[0].Z) + this.currentDisplacements[8]) - this.currentDisplacements[2];
			this.currentLength = Math.Sqrt((dX * dX) + (dY * dY) + (dZ * dZ));
			//var delta = new[] { dX, dY, dZ };
			this.currentBeamAxis[0] = dX / currentLength + beamAxisX[0];
			this.currentBeamAxis[1] = dY / currentLength + beamAxisX[1];
			this.currentBeamAxis[2] = dZ / currentLength + beamAxisX[2];
			this.currentBeamAxis.ScaleIntoThis(1d / this.currentBeamAxis.Norm2());

			//this.currentBeamAxis.divideAllIntoThis(delta, this.currentLength);
			//this.currentBeamAxis.add(this.beamAxisX);
			//this.currentBeamAxis.divideAll(this.currentBeamAxis.normEntrywise(EntrywiseNorms.L2));
		}

		private void InitializeElementAxes()
		{
			var globalVectorX = new double[AXIS_COUNT];
			var globalVectorY = new double[AXIS_COUNT];
			var globalVectorZ = new double[AXIS_COUNT];
			globalVectorX[0] = 1d;
			globalVectorY[1] = 1d;
			globalVectorZ[2] = 1d;
			double deltaX = nodes[1].X - nodes[0].X;
			double deltaY = nodes[1].Y - nodes[0].Y;
			double deltaZ = nodes[1].Z - nodes[0].Z;
			var elementAxis = new double[3];
			elementAxis[0] = deltaX;
			elementAxis[1] = deltaY;
			elementAxis[2] = deltaZ;
			elementAxis.ScaleIntoThis(1d / elementAxis.Norm2());

			currentRotationMatrix = RotationMatrix.CalculateRotationMatrix(globalVectorX, elementAxis);
			currentRotationMatrix.MultiplyIntoResult(globalVectorX, initialAxisX);
			currentRotationMatrix.MultiplyIntoResult(globalVectorY, initialAxisY);
			currentRotationMatrix.MultiplyIntoResult(globalVectorZ, initialAxisZ);
			beamAxisX.CopyFrom(initialAxisX);
			beamAxisY.CopyFrom(initialAxisY);
			beamAxisZ.CopyFrom(initialAxisZ);
		}

		private void UpdateNaturalDeformations(double[] vectorPartDifference)
		{
			double extension = this.currentLength - this.initialLength;
			double[] symmetricRotation = currentRotationMatrix.Multiply(vectorPartDifference, true);
			double[] antisymmetricRotation =
				currentRotationMatrix.Multiply(this.beamAxisX.CrossProduct(this.currentBeamAxis), true);

			//currentRotationMatrix.Transpose().Multiply(vectorPartDifference, symmetricRotation.Data);
			//currentRotationMatrix.Multiply(vectorPartDifference, symmetricRotation.Data);
			//currentRotationMatrix.Multiply(this.beamAxisX ^ this.currentBeamAxis, antisymmetricRotation.Data);
			//currentRotationMatrix.Transpose().Multiply(this.beamAxisX ^ this.currentBeamAxis, antisymmetricRotation.Data);
			//currentRotationMatrix.MultiplyTranspose2(this.beamAxisX ^ this.currentBeamAxis, antisymmetricRotation.Data);
			//final VectorView symmetricRotation = vectorPartDifference.leftMultiplyWithMatrix(this.currentRotationMatrix);
			//final VectorView antisymmetricRotation =
			//        GeometricUtils.crossProduct(this.beamAxisX, this.currentBeamAxis)
			//                        .leftMultiplyWithMatrix(this.currentRotationMatrix);

			this.naturalDeformations[0] = symmetricRotation[0] * 4.0;
			this.naturalDeformations[1] = symmetricRotation[1] * 4.0;
			this.naturalDeformations[2] = symmetricRotation[2] * 4.0;
			this.naturalDeformations[3] = extension;
			this.naturalDeformations[4] = antisymmetricRotation[1] * 4.0;
			this.naturalDeformations[5] = antisymmetricRotation[2] * 4.0;
		}

		private void UpdateRotationMatrix()
		{
			var rotationMatrixLocal = Matrix.CreateIdentity(AXIS_COUNT);
			rotationMatrixLocal.AxpyIntoThis(currentBeamAxis.TensorProduct(currentBeamAxis), -2d);
			//rotationMatrixLocal.LinearCombinationGOAT(new[] { -2d }, new[] { Matrix2D.FromVector(currentBeamAxis.Data) * Matrix2D.FromVectorTranspose(currentBeamAxis.Data) });

			double[] minusBeamAxisX = beamAxisX.Scale(-1.0);
			var tempBeamAxisY = beamAxisY.Copy();
			var tempBeamAxisZ = beamAxisZ.Copy();
			rotationMatrixLocal.MultiplyIntoResult(minusBeamAxisX, beamAxisX);
			rotationMatrixLocal.MultiplyIntoResult(tempBeamAxisY, beamAxisY);
			rotationMatrixLocal.MultiplyIntoResult(tempBeamAxisZ, beamAxisZ);
			this.currentRotationMatrix = RotationMatrix.CalculateFromOrthonormalVectors(this.beamAxisX, this.beamAxisY, this.beamAxisZ);
		}

	}
}
