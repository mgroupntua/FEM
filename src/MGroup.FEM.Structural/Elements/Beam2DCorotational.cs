using System;
using System.Collections.Generic;
using MGroup.FEM.Elements.SupportiveClasses;
using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.Geometry;

namespace MGroup.FEM.Structural.Elements
{
	public class Beam2DCorotational : Beam2DCorotationalAbstract
	{
		private readonly double[] lastDisplacements;
		private readonly double[] currentDisplacements;
		private readonly double[] displacementsOfCurrentIncrement;

		public Beam2DCorotational(IList<Node> nodes, IFiniteElementMaterial material, double density,
			BeamSection2D beamSection)
			: base(nodes, material, density, beamSection)
		{
			this.displacementsOfCurrentIncrement = new double[FREEDOM_DEGREE_COUNT];
			this.lastDisplacements = new double[FREEDOM_DEGREE_COUNT];
			this.currentDisplacements = new double[FREEDOM_DEGREE_COUNT];
			this.InitializeElementAxes();
		}

		public override void SaveGeometryState()
		{
			displacementsOfCurrentIncrement.Scale(0d);
			lastDisplacements.CopyFrom(currentDisplacements);
		}

		public override void UpdateState(double[] incrementalNodeDisplacements)
		{
			currentDisplacements.CopyFrom(lastDisplacements);
			currentDisplacements.AddIntoThis(incrementalNodeDisplacements);
			displacementsOfCurrentIncrement.CopyFrom(incrementalNodeDisplacements);

			double currentDisplacementX_A = currentDisplacements[0];
			double currentlDisplacementY_A = currentDisplacements[1];
			double currentRotationZ_A = currentDisplacements[2];
			double currentDisplacementX_B = currentDisplacements[3];
			double currentDisplacementY_B = currentDisplacements[4];
			double currentRotationZ_Β = currentDisplacements[5];

			double dX = ((nodes[1].X - nodes[0].X) + currentDisplacementX_B) - currentDisplacementX_A;
			double dY = ((nodes[1].Y - nodes[0].Y) + currentDisplacementY_B) - currentlDisplacementY_A;
			this.currentLength = Math.Sqrt((dX * dX) + (dY * dY));
			double axisAngle = 0d;
			if ((dY == 0) && (dX == currentLength)) { axisAngle = 0d; }
			else
				if ((dY == 0) && (dX == -currentLength)) { axisAngle = Math.PI; }
			else { axisAngle = 2.0 * Math.Atan((currentLength - dX) / dY); }

			double initialdX = nodes[1].X - nodes[0].X;
			double initialdY = nodes[1].Y - nodes[0].Y;
			double initialLength = Math.Sqrt((initialdX * initialdX) + (initialdY * initialdY));
			double axisAngleInitial = 0d;
			if ((initialdY == 0) && (initialdX == initialLength)) { axisAngleInitial = 0d; }
			else
				if ((initialdY == 0) && (initialdX == -initialLength)) { axisAngleInitial = Math.PI; }
			else { axisAngleInitial = 2.0 * Math.Atan((initialLength - dX) / dY); }

			double symmetricAngle = currentRotationZ_Β - currentRotationZ_A;
			double antiSymmetricAngle = currentRotationZ_Β + currentRotationZ_A - 2.0 * (axisAngle - axisAngleInitial);

			currentRotationMatrix = RotationMatrix.CalculateRotationMatrixBeam2D(axisAngle);
			double extension = this.currentLength - this.initialLength;
			this.naturalDeformations[0] = extension;
			this.naturalDeformations[1] = symmetricAngle;
			this.naturalDeformations[2] = antiSymmetricAngle;
		}

		private void InitializeElementAxes()
		{
			double deltaX = nodes[1].X - nodes[0].X;
			double deltaY = nodes[1].Y - nodes[0].Y;
			this.currentLength = Math.Sqrt((deltaX * deltaX) + (deltaY * deltaY));
			double axisAngle = 0d;
			if ((deltaY == 0) && (deltaX == currentLength)) { axisAngle = 0d; }
			else
				if ((deltaY == 0) && (deltaX == -currentLength)) { axisAngle = Math.PI; }
			else { axisAngle = 2.0 * Math.Atan((currentLength - deltaX) / deltaY); }

			currentRotationMatrix = RotationMatrix.CalculateRotationMatrixBeam2D(axisAngle);
		}
	}
}
