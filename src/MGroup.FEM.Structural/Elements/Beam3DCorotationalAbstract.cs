using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.FEM.Elements.SupportiveClasses;
using MGroup.FEM.Embedding;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Structural.Elements
{
	public abstract class Beam3DCorotationalAbstract : IFiniteElement, IEmbeddedElement
	{
		protected static readonly int NATURAL_DEFORMATION_COUNT = 6;
		protected static readonly int FREEDOM_DEGREE_COUNT = 12;
		protected static readonly int AXIS_COUNT = 3;
		protected static readonly int NODE_COUNT = 2;

		protected readonly static IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY, StructuralDof.RotationZ };
		protected readonly static IDofType[][] dofTypes = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		//protected static final List<Set<FreedomDegreeType>> FREEDOM_DEGREE_TYPES =
		//        Collections.nCopies(NODE_COUNT, FreedomDegreeTypeSets.X_Y_Z_ROTX_ROTY_ROTZ);
		protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		protected readonly IIsotropicContinuumMaterial3D material;
		protected readonly IList<Node> nodes;
		protected readonly double density;
		protected BeamSection3D beamSection;
		protected readonly double initialLength;
		protected double currentLength;
		protected Matrix currentRotationMatrix;
		protected double[] naturalDeformations;
		protected double[] beamAxisX;
		protected double[] beamAxisY;
		protected double[] beamAxisZ;
		private readonly List<EmbeddedNode> embeddedNodes = new List<EmbeddedNode>();


		public double RayleighAlpha { get; set; }
		public double RayleighBeta { get; set; }

		protected Beam3DCorotationalAbstract(IList<Node> nodes, IIsotropicContinuumMaterial3D material, double density,
			BeamSection3D beamSection)
		{
			this.nodes = nodes;
			this.material = material;
			this.density = density;
			this.beamSection = beamSection;
			this.initialLength = Math.Sqrt(Math.Pow(nodes[0].X - nodes[1].X, 2) + Math.Pow(nodes[0].Y - nodes[1].Y, 2)
				+ Math.Pow(nodes[0].Z - nodes[1].Z, 2));
			this.currentLength = this.initialLength;
			this.currentRotationMatrix = Matrix.CreateZero(AXIS_COUNT, AXIS_COUNT);
			this.naturalDeformations = new double[NATURAL_DEFORMATION_COUNT];
			this.beamAxisX = new double[AXIS_COUNT];
			this.beamAxisY = new double[AXIS_COUNT];
			this.beamAxisZ = new double[AXIS_COUNT];
		}

		public CellType CellType { get; } = CellType.Line;
		public int ID => 100;
		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;
		public IList<EmbeddedNode> EmbeddedNodes { get { return embeddedNodes; } }
		public bool MaterialModified => material.Modified;
		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		public abstract void SaveGeometryState();
		public abstract void UpdateState(double[] incrementalNodeDisplacements);

		private Matrix CalculateBlockRotationMatrix()
		{
			Matrix blockRotationMatrix = Matrix.CreateZero(FREEDOM_DEGREE_COUNT, FREEDOM_DEGREE_COUNT);
			int totalBlocks = 4;
			int blockSize = 3;
			double R11 = this.currentRotationMatrix[0, 0];
			double R12 = this.currentRotationMatrix[0, 1];
			double R13 = this.currentRotationMatrix[0, 2];
			double R21 = this.currentRotationMatrix[1, 0];
			double R22 = this.currentRotationMatrix[1, 1];
			double R23 = this.currentRotationMatrix[1, 2];
			double R31 = this.currentRotationMatrix[2, 0];
			double R32 = this.currentRotationMatrix[2, 1];
			double R33 = this.currentRotationMatrix[2, 2];

			for (int block = 0; block < totalBlocks; block++)
			{
				int cellDistanceCount = block * blockSize;

				blockRotationMatrix[cellDistanceCount, cellDistanceCount] = R11;
				blockRotationMatrix[cellDistanceCount, cellDistanceCount + 1] = R12;
				blockRotationMatrix[cellDistanceCount, cellDistanceCount + 2] = R13;

				blockRotationMatrix[cellDistanceCount + 1, cellDistanceCount] = R21;
				blockRotationMatrix[cellDistanceCount + 1, cellDistanceCount + 1] = R22;
				blockRotationMatrix[cellDistanceCount + 1, cellDistanceCount + 2] = R23;

				blockRotationMatrix[cellDistanceCount + 2, cellDistanceCount] = R31;
				blockRotationMatrix[cellDistanceCount + 2, cellDistanceCount + 1] = R32;
				blockRotationMatrix[cellDistanceCount + 2, cellDistanceCount + 2] = R33;
			}

			return blockRotationMatrix;
		}

		/**
	     * Calculates the constitutive stiffness of the element.
	     *
	     * @return The constitutive stiffness
	     */
		private Matrix CalculateConstitutiveStiffness()
		{
			var constitutiveStiffness = SymmetricMatrix.CreateZero(FREEDOM_DEGREE_COUNT);
			double E = this.material.YoungModulus;
			double G = E / (2d * (1d + this.material.PoissonRatio));
			double Iy = this.beamSection.InertiaY;
			double Iz = this.beamSection.InertiaZ;
			double J = this.beamSection.TorsionalInertia;
			double A = this.beamSection.Area;
			double Ay = this.beamSection.EffectiveAreaY;
			double Az = this.beamSection.EffectiveAreaZ;
			double l = this.currentLength;
			double lSqruared = l * l;
			double lCubed = l * l * l;
			double phiY = (12.0 * E * Iy) / (l * l * G * Az);
			double phiZ = (12.0 * E * Iz) / (l * l * G * Ay);
			double psiY = 1.0 / (1.0 + phiY);
			double psiZ = 1.0 / (1.0 + phiZ);
			double EAOverL = (E * A) / l;
			double GJOverL = (G * J) / l;
			double psiZ_E_Iz_12_OverlCubed = (12.0 * psiZ * E * Iz) / lCubed;
			double psiZ_E_Iz_6_OverlSquared = (6.0 * psiZ * E * Iz) / lSqruared;
			double psiY_E_12_Iy_OverlCubed = (12.0 * psiY * E * Iy) / lCubed;
			double psiY_E_6_Iy_OverlSquared = (6.0 * psiY * E * Iy) / lSqruared;
			double psiZ_3_Plus_1_E_Iz_Overl = (((3.0 * psiZ) + 1.0) * E * Iz) / l;
			double psiZ_3_Minus_1_E_Iz_Overl = (((3.0 * psiZ) - 1.0) * E * Iz) / l;
			double psiY_3_Plus_1_E_Iy_Overl = (((3.0 * psiY) + 1.0) * E * Iy) / l;
			double psiY_3_Minus_1_E_Iy_Overl = (((3.0 * psiY) - 1.0) * E * Iy) / l;

			constitutiveStiffness[0, 0] = EAOverL;
			constitutiveStiffness[0, 6] = -EAOverL;

			constitutiveStiffness[1, 1] = psiZ_E_Iz_12_OverlCubed;
			constitutiveStiffness[1, 5] = psiZ_E_Iz_6_OverlSquared;
			constitutiveStiffness[1, 7] = -psiZ_E_Iz_12_OverlCubed;
			constitutiveStiffness[1, 11] = psiZ_E_Iz_6_OverlSquared;

			constitutiveStiffness[2, 2] = psiY_E_12_Iy_OverlCubed;
			constitutiveStiffness[2, 4] = -psiY_E_6_Iy_OverlSquared;
			constitutiveStiffness[2, 8] = -psiY_E_12_Iy_OverlCubed;
			constitutiveStiffness[2, 10] = -psiY_E_6_Iy_OverlSquared;

			constitutiveStiffness[3, 3] = GJOverL;
			constitutiveStiffness[3, 9] = -GJOverL;

			constitutiveStiffness[4, 4] = psiY_3_Plus_1_E_Iy_Overl;
			constitutiveStiffness[4, 8] = psiZ_E_Iz_6_OverlSquared;
			constitutiveStiffness[4, 10] = psiY_3_Minus_1_E_Iy_Overl;

			constitutiveStiffness[5, 5] = psiZ_3_Plus_1_E_Iz_Overl;
			constitutiveStiffness[5, 7] = -psiY_E_6_Iy_OverlSquared;
			constitutiveStiffness[5, 11] = psiZ_3_Minus_1_E_Iz_Overl;

			constitutiveStiffness[6, 6] = EAOverL;

			constitutiveStiffness[7, 7] = psiZ_E_Iz_12_OverlCubed;
			constitutiveStiffness[7, 11] = -psiZ_E_Iz_6_OverlSquared;

			constitutiveStiffness[8, 8] = psiY_E_12_Iy_OverlCubed;
			constitutiveStiffness[8, 10] = psiY_E_6_Iy_OverlSquared;

			constitutiveStiffness[9, 9] = GJOverL;

			constitutiveStiffness[10, 10] = psiY_3_Plus_1_E_Iy_Overl;

			constitutiveStiffness[11, 11] = psiZ_3_Plus_1_E_Iz_Overl;

			return constitutiveStiffness.CopyToFullMatrix();
		}

		/**
	     * Calculates the forces in the global coordinate system.
	     *
	     * @return The forces in the global coordinate system
	     */
		private double[] CalculateForcesInGlobalSystem()
		{
			Matrix transformationMatrix = this.CalculateNaturalToGlobalTransormMatrix();
			double[] forcesNatural = this.CalculateForcesInNaturalSystem();
			double[] forcesGlobal = transformationMatrix.Multiply(forcesNatural);
			return forcesGlobal;
		}

		/**
	     * Calculates the forces in the local coordinate system.
	     *
	     * @return The forces in the local coordinate system
	     */
		private double[] CalculateForcesInLocalSystem()
		{
			Matrix naturalToLocal = this.CalculateNaturalToLocalTranformMatrix();
			double[] naturalForces = this.CalculateForcesInNaturalSystem();
			double[] forcesLocal = naturalToLocal.Multiply(naturalForces);
			return forcesLocal;
		}

		/**
	     * Calculates forces in the natural coordinate system.
	     *
	     * @return The forces in the natural coordinate system
	     */
		private double[] CalculateForcesInNaturalSystem()
		{
			var forcesNatural = new double[NATURAL_DEFORMATION_COUNT];
			double E = this.material.YoungModulus;
			double G = E / (2d * (1d + this.material.PoissonRatio));
			double Iy = this.beamSection.InertiaY;
			double Iz = this.beamSection.InertiaZ;
			double J = this.beamSection.TorsionalInertia;
			double A = this.beamSection.Area;
			double Ay = this.beamSection.EffectiveAreaY;
			double Az = this.beamSection.EffectiveAreaZ;
			double l = this.currentLength;
			double phiY = (12.0 * E * Iy) / (l * l * G * Az);
			double phiZ = (12.0 * E * Iz) / (l * l * G * Ay);
			double psiY = 1.0 / (1.0 + phiY);
			double psiZ = 1.0 / (1.0 + phiZ);
			forcesNatural[NaturalDeformationMode3D.TORSION] =
				(G * J * this.naturalDeformations[NaturalDeformationMode3D.TORSION]) / l;
			forcesNatural[NaturalDeformationMode3D.SYMMETRIC_BENDING_Y] =
				(E * Iy * this.naturalDeformations[NaturalDeformationMode3D.SYMMETRIC_BENDING_Y]) / l;
			forcesNatural[NaturalDeformationMode3D.SYMMETRIC_BENDING_Z] =
				(E * Iz * this.naturalDeformations[NaturalDeformationMode3D.SYMMETRIC_BENDING_Z]) / l;
			forcesNatural[NaturalDeformationMode3D.EXTENSION] =
				(E * A * this.naturalDeformations[NaturalDeformationMode3D.EXTENSION]) / l;
			forcesNatural[NaturalDeformationMode3D.ANTISYMMETRIC_BENDING_Y] =
				(3.0 * psiY * E * Iy * this.naturalDeformations[NaturalDeformationMode3D.ANTISYMMETRIC_BENDING_Y]) / l;
			forcesNatural[NaturalDeformationMode3D.ANTISYMMETRIC_BENDING_Z] =
				(3.0 * psiZ * E * Iz * this.naturalDeformations[NaturalDeformationMode3D.ANTISYMMETRIC_BENDING_Z]) / l;
			return forcesNatural;
		}

		/**
	     * Calculates the geometric stiffness of the element.
	     *
	     * @return The geometric stiffness
	     */
		private Matrix CalculateGeometricStiffness()
		{
			var geometricStiffness = SymmetricMatrix.CreateZero(FREEDOM_DEGREE_COUNT);
			var forcesInNaturalSystem = this.CalculateForcesInNaturalSystem();
			var forcesInLocalSystem = this.CalculateForcesInLocalSystem();
			double torsionalMoment = forcesInNaturalSystem[NaturalDeformationMode3D.TORSION];
			double axialForce = forcesInNaturalSystem[NaturalDeformationMode3D.EXTENSION];
			double momentY_A = forcesInLocalSystem[4];
			double momentZ_A = forcesInLocalSystem[5];
			double momentY_B = forcesInLocalSystem[10];
			double momentZ_B = forcesInLocalSystem[11];
			double length = this.currentLength;
			double Qy = -(momentZ_A + momentZ_B) / length;
			double Qz = (momentY_A + momentY_B) / length;
			double Qy_Over_L = Qy / length;
			double Qz_Over_L = Qz / length;
			double Qy_L_Over_6 = (Qy * length) / 6.0;
			double Qz_L_Over_6 = (Qz * length) / 6.0;
			double N_6_Over_5_L = (6.0 * axialForce) / (5.0 * length);
			double N_Over_10 = axialForce / 10.0;
			double N_L_4_Over_30 = (4.0 * axialForce * length) / 30.0;
			double N_L_Over_30 = (axialForce * length) / 30.0;
			double MyA_Over_L = momentY_A / length;
			double MyB_Over_L = momentY_B / length;
			double MzA_Over_L = momentZ_A / length;
			double MzB_Over_L = momentZ_B / length;
			double M_Over_L = torsionalMoment / length;
			double M_3_Over_6 = (3.0 * torsionalMoment) / 6.0;

			geometricStiffness[0, 1] = -Qy_Over_L;
			geometricStiffness[0, 2] = -Qz_Over_L;
			geometricStiffness[0, 7] = Qy_Over_L;
			geometricStiffness[0, 8] = Qz_Over_L;

			geometricStiffness[1, 1] = N_6_Over_5_L;
			geometricStiffness[1, 3] = MyA_Over_L;
			geometricStiffness[1, 4] = M_Over_L;
			geometricStiffness[1, 5] = N_Over_10;
			geometricStiffness[1, 6] = Qy_Over_L;
			geometricStiffness[1, 7] = -N_6_Over_5_L;
			geometricStiffness[1, 9] = MyB_Over_L;
			geometricStiffness[1, 10] = -M_Over_L;
			geometricStiffness[1, 11] = N_Over_10;

			geometricStiffness[2, 2] = N_6_Over_5_L;
			geometricStiffness[2, 3] = MzA_Over_L;
			geometricStiffness[2, 4] = -N_Over_10;
			geometricStiffness[2, 5] = M_Over_L;
			geometricStiffness[2, 6] = Qz_Over_L;
			geometricStiffness[2, 8] = -N_6_Over_5_L;
			geometricStiffness[2, 9] = MzB_Over_L;
			geometricStiffness[2, 10] = -N_Over_10;
			geometricStiffness[2, 11] = -M_Over_L;

			geometricStiffness[3, 4] = ((-2.0 * momentZ_A) + momentZ_B) / 6.0;
			geometricStiffness[3, 5] = ((2.0 * momentY_A) - momentY_B) / 6.0;
			geometricStiffness[3, 7] = -MyA_Over_L;
			geometricStiffness[3, 8] = -MzA_Over_L;
			geometricStiffness[3, 10] = Qy_L_Over_6;
			geometricStiffness[3, 11] = Qz_L_Over_6;

			geometricStiffness[4, 4] = N_L_4_Over_30;
			geometricStiffness[4, 7] = -M_Over_L;
			geometricStiffness[4, 8] = +N_Over_10;
			geometricStiffness[4, 9] = Qy_L_Over_6;
			geometricStiffness[4, 10] = -N_L_Over_30;
			geometricStiffness[4, 11] = M_3_Over_6;

			geometricStiffness[5, 5] = N_L_4_Over_30;
			geometricStiffness[5, 7] = -N_Over_10;
			geometricStiffness[5, 8] = -M_Over_L;
			geometricStiffness[5, 9] = Qz_L_Over_6;
			geometricStiffness[5, 10] = -M_3_Over_6;
			geometricStiffness[5, 11] = -N_L_Over_30;

			geometricStiffness[6, 7] = -Qy_Over_L;
			geometricStiffness[6, 8] = -Qz_Over_L;

			geometricStiffness[7, 7] = N_6_Over_5_L;
			geometricStiffness[7, 9] = -MyB_Over_L;
			geometricStiffness[7, 10] = M_Over_L;
			geometricStiffness[7, 11] = -N_Over_10;

			geometricStiffness[8, 8] = N_6_Over_5_L;
			geometricStiffness[8, 9] = -MzB_Over_L;
			geometricStiffness[8, 10] = N_Over_10;
			geometricStiffness[8, 11] = M_Over_L;

			geometricStiffness[9, 10] = ((-2.0 * momentZ_B) + momentZ_A) / 6.0;
			geometricStiffness[9, 11] = ((+2.0 * momentY_B) - momentY_A) / 6.0;

			geometricStiffness[10, 10] = N_L_4_Over_30;

			geometricStiffness[11, 11] = N_L_4_Over_30;

			return geometricStiffness.CopyToFullMatrix();
		}

		/**
	     * Calculates the stiffness matrix in the local coordinate system.
	     *
	     * @return The stiffness matrix in the local coordinate system.
	     */
		private Matrix CalculateLocalStiffnessMatrix()
		{
			Matrix constitutivePart = this.CalculateConstitutiveStiffness();
			Matrix geometricPart = this.CalculateGeometricStiffness();
			constitutivePart.AddIntoThis(geometricPart);
			return constitutivePart;
		}

		/**
	     * Calculates the transformation matrix from natural to local coordinate system.
	     *
	     * @return The natural to local transformation matrix
	     */
		private Matrix CalculateNaturalToGlobalTransormMatrix()
		{
			var transformMatrix = Matrix.CreateZero(FREEDOM_DEGREE_COUNT, NATURAL_DEFORMATION_COUNT);
			double L = this.currentLength;
			double nx1 = this.beamAxisX[0];
			double nx2 = this.beamAxisX[1];
			double nx3 = this.beamAxisX[2];
			double ny1 = this.beamAxisY[0];
			double ny2 = this.beamAxisY[1];
			double ny3 = this.beamAxisY[2];
			double nz1 = this.beamAxisZ[0];
			double nz2 = this.beamAxisZ[1];
			double nz3 = this.beamAxisZ[2];

			transformMatrix[0, 3] = -nx1;
			transformMatrix[0, 4] = (-2.0 * nz1) / L;
			transformMatrix[0, 5] = (2.0 * ny1) / L;

			transformMatrix[1, 3] = -nx2;
			transformMatrix[1, 4] = (-2.0 * nz2) / L;
			transformMatrix[1, 5] = (2.0 * ny2) / L;

			transformMatrix[2, 3] = -nx3;
			transformMatrix[2, 4] = (-2.0 * nz3) / L;
			transformMatrix[2, 5] = (2.0 * ny3) / L;

			transformMatrix[3, 0] = -nx1;
			transformMatrix[3, 1] = -ny1;
			transformMatrix[3, 2] = -nz1;
			transformMatrix[3, 4] = ny1;
			transformMatrix[3, 5] = nz1;

			transformMatrix[4, 0] = -nx2;
			transformMatrix[4, 1] = -ny2;
			transformMatrix[4, 2] = -nz2;
			transformMatrix[4, 4] = ny2;
			transformMatrix[4, 5] = nz2;

			transformMatrix[5, 0] = -nx3;
			transformMatrix[5, 1] = -ny3;
			transformMatrix[5, 2] = -nz3;
			transformMatrix[5, 4] = ny3;
			transformMatrix[5, 5] = nz3;

			transformMatrix[6, 3] = nx1;
			transformMatrix[6, 4] = (2.0 * nz1) / L;
			transformMatrix[6, 5] = (-2.0 * ny1) / L;

			transformMatrix[7, 3] = nx2;
			transformMatrix[7, 4] = (2.0 * nz2) / L;
			transformMatrix[7, 5] = (-2.0 * ny2) / L;

			transformMatrix[8, 3] = nx3;
			transformMatrix[8, 4] = (2.0 * nz3) / L;
			transformMatrix[8, 5] = (-2.0 * ny3) / L;

			transformMatrix[9, 0] = nx1;
			transformMatrix[9, 1] = ny1;
			transformMatrix[9, 2] = nz1;
			transformMatrix[9, 4] = ny1;
			transformMatrix[9, 5] = nz1;

			transformMatrix[10, 0] = nx2;
			transformMatrix[10, 1] = ny2;
			transformMatrix[10, 2] = nz2;
			transformMatrix[10, 4] = ny2;
			transformMatrix[10, 5] = nz2;

			transformMatrix[11, 0] = nx3;
			transformMatrix[11, 1] = ny3;
			transformMatrix[11, 2] = nz3;
			transformMatrix[11, 4] = ny3;
			transformMatrix[11, 5] = nz3;

			return transformMatrix;
		}

		/**
	     * Calculates the transformation matrix from natural to local coordinate system.
	     *
	     * @return The natural to local transformation matrix
	     */
		private Matrix CalculateNaturalToLocalTranformMatrix()
		{
			var transformMatrix = Matrix.CreateZero(FREEDOM_DEGREE_COUNT, NATURAL_DEFORMATION_COUNT);
			double L = this.currentLength;
			double nx1 = 1.0;
			double ny2 = 1.0;
			double nz3 = 1.0;

			transformMatrix[0, 3] = -nx1;

			transformMatrix[1, 5] = (2.0 * ny2) / L;

			transformMatrix[2, 4] = (-2.0 * nz3) / L;

			transformMatrix[3, 0] = -nx1;

			transformMatrix[4, 1] = -ny2;
			transformMatrix[4, 4] = ny2;

			transformMatrix[5, 2] = -nz3;
			transformMatrix[5, 5] = nz3;

			transformMatrix[6, 3] = nx1;

			transformMatrix[7, 5] = (-2.0 * ny2) / L;

			transformMatrix[8, 4] = (2.0 * nz3) / L;

			transformMatrix[9, 0] = nx1;

			transformMatrix[10, 1] = ny2;
			transformMatrix[10, 4] = ny2;

			transformMatrix[11, 2] = nz3;
			transformMatrix[11, 5] = nz3;

			return transformMatrix;
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;

		public IMatrix StiffnessMatrix(IElement element)
		{
			Matrix rotationMatrixBlock = this.CalculateBlockRotationMatrix();
			Matrix localStiffnessMatrix = this.CalculateLocalStiffnessMatrix();
			Matrix s = rotationMatrixBlock.MultiplyRight(localStiffnessMatrix).MultiplyRight(rotationMatrixBlock, false, true);
			return dofEnumerator.GetTransformedMatrix(s);
		}

		public IMatrix MassMatrix(IElement element)
		{
			//throw new NotImplementedException();
			double area = beamSection.Area;
			double inertiaY = beamSection.InertiaY;
			double inertiaZ = beamSection.InertiaZ;
			double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
			double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
			double z2 = Math.Pow(element.Nodes[1].Z - element.Nodes[0].Z, 2);
			double L = Math.Sqrt(x2 + y2 + z2);
			double fullMass = density * area * L;

			var massMatrix = Matrix.CreateZero(FREEDOM_DEGREE_COUNT, FREEDOM_DEGREE_COUNT);
			massMatrix[0, 0] = (1.0 / 3.0) * fullMass;
			massMatrix[0, 6] = (1.0 / 6.0) * fullMass;

			massMatrix[1, 1] = (13.0 / 35.0) * fullMass;
			massMatrix[1, 5] = (11.0 * L / 210.0) * fullMass;
			massMatrix[1, 7] = (9.0 / 70.0) * fullMass;
			massMatrix[1, 11] = -(13.0 * L / 420.0) * fullMass;

			massMatrix[2, 2] = (13.0 / 35.0) * fullMass;
			massMatrix[2, 4] = -(11.0 * L / 210.0) * fullMass;
			massMatrix[2, 8] = (9.0 / 70.0) * fullMass;
			massMatrix[2, 10] = -(13.0 * L / 420.0) * fullMass;

			massMatrix[3, 3] = ((inertiaY + inertiaZ) / (3.0 * area)) * fullMass;
			massMatrix[3, 9] = ((inertiaY + inertiaZ) / (6.0 * area)) * fullMass;

			massMatrix[4, 4] = ((L * L) / 105.0) * fullMass;
			massMatrix[4, 8] = -(13.0 * L / 420.0) * fullMass;
			massMatrix[4, 10] = -((L * L) / 105.0) * fullMass;

			massMatrix[5, 5] = ((L * L) / 105.0) * fullMass;
			massMatrix[5, 7] = (13.0 * L / 420.0) * fullMass;
			massMatrix[5, 11] = -((L * L) / 105.0) * fullMass;

			massMatrix[6, 6] = (1.0 / 3.0) * fullMass;

			massMatrix[7, 7] = (13.0 / 35.0) * fullMass;
			massMatrix[7, 11] = -(11.0 * L / 210.0) * fullMass;

			massMatrix[8, 8] = (13.0 / 35.0) * fullMass;
			massMatrix[8, 10] = (11.0 * L / 210.0) * fullMass;

			massMatrix[9, 9] = ((inertiaY + inertiaZ) / (3.0 * area)) * fullMass;

			massMatrix[10, 10] = ((L * L) / 105.0) * fullMass;

			massMatrix[11, 11] = ((L * L) / 105.0) * fullMass;

			massMatrix[6, 0] = (1.0 / 6.0) * fullMass;
			massMatrix[5, 1] = (11.0 * L / 210.0) * fullMass;
			massMatrix[7, 1] = (9.0 / 70.0) * fullMass;
			massMatrix[11, 1] = -(13.0 * L / 420.0) * fullMass;
			massMatrix[4, 2] = -(11.0 * L / 210.0) * fullMass;
			massMatrix[8, 2] = (9.0 / 70.0) * fullMass;
			massMatrix[10, 2] = -(13.0 * L / 420.0) * fullMass;
			massMatrix[9, 3] = ((inertiaY + inertiaZ) / (6.0 * area)) * fullMass;
			massMatrix[8, 4] = -(13.0 * L / 420.0) * fullMass;
			massMatrix[10, 4] = -((L * L) / 105.0) * fullMass;
			massMatrix[7, 5] = (13.0 * L / 420.0) * fullMass;
			massMatrix[11, 5] = -((L * L) / 105.0) * fullMass;
			massMatrix[11, 7] = -(11.0 * L / 210.0) * fullMass;
			massMatrix[10, 8] = (11.0 * L / 210.0) * fullMass;

			return dofEnumerator.GetTransformedMatrix(massMatrix);
		}

		public IMatrix DampingMatrix(IElement element)
		{
			IMatrix k = StiffnessMatrix(element);
			IMatrix m = MassMatrix(element);
			k.LinearCombinationIntoThis(RayleighBeta, m, RayleighAlpha);
			return dofEnumerator.GetTransformedMatrix(k);
		}

		public void ResetMaterialModified() => material.ResetModified();

		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			UpdateState(dofEnumerator.GetTransformedDisplacementsVector(localdDisplacements));
			//TODO: Should calculate strains and update material as well
			//material.UpdateMaterial(strains);
			//TODO: Should calculate stresses as well
			return new Tuple<double[], double[]>(new double[FREEDOM_DEGREE_COUNT], new double[FREEDOM_DEGREE_COUNT]);
		}

		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			double[] internalForces = this.CalculateForcesInGlobalSystem();
			return dofEnumerator.GetTransformedForcesVector(internalForces);
		}

		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			var accelerations = new double[6];
			IMatrix massMatrix = MassMatrix(element);

			int index = 0;
			foreach (MassAccelerationLoad load in loads)
				foreach (IDofType[] nodalDOFTypes in dofs)
					foreach (IDofType dofType in nodalDOFTypes)
					{
						if (dofType == load.DOF) accelerations[index] += load.Amount;
						index++;
					}

			return massMatrix.Multiply(accelerations);
		}

		public void SaveMaterialState()
		{
			SaveGeometryState();
			material.SaveState();
		}

		public void ClearMaterialState() => material.ClearState();

		public void ClearMaterialStresses() => material.ClearStresses();

		public Dictionary<IDofType, int> GetInternalNodalDOFs(Element element, Node node)
		{
			throw new NotImplementedException();
		}

		public double[] GetLocalDOFValues(Element hostElement, double[] hostDOFValues)
		{
			throw new NotImplementedException();
		}
	}
}
