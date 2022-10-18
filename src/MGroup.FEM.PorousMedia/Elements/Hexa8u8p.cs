using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using MGroup.Constitutive.PorousMedia;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.DataStructures;

namespace MGroup.FEM.PorousMedia.Elements
{
	public class Hexa8u8p : IPorousElementType
	{
		protected static int iInt = 2;
		protected static int iInt2 = iInt * iInt;
		protected static int iInt3 = iInt * iInt * iInt;
		private readonly static IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, PorousMediaDof.Pressure };
		private readonly static IDofType[][] dofTypes = new IDofType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
			nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
		protected IIsotropicContinuumMaterial3D[] materialsAtGaussPoints;
		protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private double[] fluidDisplacements;
		private double[][] lastStresses;
		private double youngModulus, poissonRatio;

		#region Fortran imports
		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH8GAUSSMATRICES",
		//	CallingConvention = CallingConvention.Cdecl)]
		//protected static extern void CalcH8GaussMatrices(ref int iInt, [MarshalAs(UnmanagedType.LPArray)] double[,] faXYZ,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faWeight, [MarshalAs(UnmanagedType.LPArray)] double[,] faS,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,] faDS, [MarshalAs(UnmanagedType.LPArray)] double[,,] faJ,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faDetJ, [MarshalAs(UnmanagedType.LPArray)] double[,,] faB);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH8STRAINS",
		//	CallingConvention = CallingConvention.Cdecl)]
		//private static extern void CalcH8Strains(ref int iInt,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,,] faB, [MarshalAs(UnmanagedType.LPArray)] double[] fau,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,] faStrains);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH8FORCES",
		//	CallingConvention = CallingConvention.Cdecl)]
		//private static extern void CalcH8Forces(ref int iInt,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,,] faB, [MarshalAs(UnmanagedType.LPArray)] double[] faWeight,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,] faStresses,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faForces);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH20U8PFORCESWATERACC",
		//	CallingConvention = CallingConvention.Cdecl)]
		//private static extern void CalcH20u8pForcesWaterAcc(ref int iInt,
		//	[MarshalAs(UnmanagedType.LPArray)] bool[] alImpermeable, ref double ffDensity,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] afPermeability,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] afXw,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] afSw,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] afAcc,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,] afS,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,,] afB, [MarshalAs(UnmanagedType.LPArray)] double[] afWeights,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] afLocalForces);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH8K",
		//	CallingConvention = CallingConvention.Cdecl)]
		//protected static extern void CalcH8K(ref int iInt, [MarshalAs(UnmanagedType.LPArray)] double[,,] faE,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,,] faB, [MarshalAs(UnmanagedType.LPArray)] double[] faWeight,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faK);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH8MLUMPED",
		//	CallingConvention = CallingConvention.Cdecl)]
		//private static extern void CalcH8MLumped(ref int iInt, ref double fDensity,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faWeight, [MarshalAs(UnmanagedType.LPArray)] double[] faM);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH20U8PH",
		//	CallingConvention = CallingConvention.Cdecl)]
		//private static extern void CalcH20u8pH(ref int iInt, [MarshalAs(UnmanagedType.LPArray)] double[] faPermeability,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,] faS, [MarshalAs(UnmanagedType.LPArray)] double[,,] faB,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faWeight, [MarshalAs(UnmanagedType.LPArray)] double[] faH);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH20U8PS",
		//	CallingConvention = CallingConvention.Cdecl)]
		//private static extern void CalcH20u8pS(ref int iInt, [MarshalAs(UnmanagedType.LPArray)] double[] faXwDivQ,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,] faS,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faWeight, [MarshalAs(UnmanagedType.LPArray)] double[] faM);

		//[DllImport("femelements.dll",
		//	EntryPoint = "CALCH8U8PQMINUS",
		//	CallingConvention = CallingConvention.Cdecl)]
		//private static extern void CalcH8u8pQMinus(ref int iInt, ref double fPoreA,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faXw,
		//	[MarshalAs(UnmanagedType.LPArray)] double[,,] faB, [MarshalAs(UnmanagedType.LPArray)] double[,] faS,
		//	[MarshalAs(UnmanagedType.LPArray)] double[] faWeight, [MarshalAs(UnmanagedType.LPArray)] double[,] faQ);

		private void CalcH8K(ref int iInt, double[,,] afE, double[,,] faB, double[] faWeight, double[] faK)
		{
			throw new NotImplementedException();
		}

		private void CalcH8GaussMatrices(ref int iInt, double[,] faXYZ, double[] faWeight, double[,] faS, double[,] faDS, double[,,] faJ, double[] faDetJ, double[,,] faB)
		{
			throw new NotImplementedException();
		}

		private void CalcH8MLumped(ref int iInt, ref double fDensity, double[] faWeight, double[] faM)
		{
			throw new NotImplementedException();
		}
		private void CalcH8Strains(ref int iInt, double[,,] faB, double[] solidDisplacements, double[,] faStrains)
		{
			throw new NotImplementedException();
		}

		private void CalcH20u8pH(ref int iInt, double[] faPermeability, double[,] faS, double[,,] faB, double[] faWeight, double[] faH)
		{
			throw new NotImplementedException();
		}
		private void CalcH20u8pForcesWaterAcc(ref int iInt, bool[] impermeableDOFs, ref double fD, double[] faPermeability, double[] faXw, double[] faSw, double[] waterAcc, double[,] faS, double[,,] faB, double[] faWeight, double[] fluidDrags)
		{
			throw new NotImplementedException();
		}
		private void CalcH8Forces(ref int iInt, double[,,] faB, double[] faWeight, double[,] faStresses, double[] solidForces)
		{
			throw new NotImplementedException();
		}
		private void CalcH8u8pQMinus(ref int iInt, ref double fPoreA, double[] faXw, double[,,] faB, double[,] faS, double[] faWeight, double[,] faQ)
		{
			throw new NotImplementedException();
		}
		private void CalcH20u8pS(ref int iInt, double[] faXwDivQ, double[,] faS, double[] faWeight, double[] faSaturation)
		{
			throw new NotImplementedException();
		}
		#endregion

		protected Hexa8u8p()
		{
		}

		public Hexa8u8p(IIsotropicContinuumMaterial3D material, double youngModulus, double poissonRatio)
		{
			throw new NotImplementedException();

			this.youngModulus = youngModulus;
			this.poissonRatio = poissonRatio;
			//materialsAtGaussPoints = new IIsotropicContinuumMaterial3D[Hexa8u8p.iInt3];
			//lastStresses = new double[Hexa8u8p.iInt3][];
			//for (int i = 0; i < Hexa8u8p.iInt3; i++)
			//{
			//    materialsAtGaussPoints[i] = (IIsotropicContinuumMaterial3D)material.Clone();
			//    lastStresses[i] = new double[6];
			//}
		}

		public Hexa8u8p(double youngModulus, double poissonRatio, IIsotropicContinuumMaterial3D material, 
			IElementDofEnumerator dofEnumerator) : this(material, youngModulus, poissonRatio)
		{
			this.dofEnumerator = dofEnumerator;
		}

		public CellType CellType { get; } = CellType.Hexa8;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		public double FluidBulkModulus { get; set; }
		public double SolidDensity { get; set; }
		public double FluidDensity { get; set; }
		public double Permeability { get; set; }
		public double Porosity { get; set; }
		public double Saturation { get; set; }
		public double PoreA { get; set; }
		public double Xw { get; set; }
		public double Cs { get; set; }
		public double RayleighAlpha { get; set; }
		public double RayleighBeta { get; set; }

		public double SolidBulkModulus => youngModulus / (3 - 6 * poissonRatio);

		public double Density => Porosity * Saturation * FluidDensity + (1 - Porosity) * SolidDensity;

		public double QInv => Cs + Porosity * Saturation / FluidBulkModulus + (PoreA - Porosity) * Saturation / SolidBulkModulus;

		protected double[,] GetCoordinates()
		{
			double[,] faXYZ = new double[dofTypes.Length, 3];
			for (int i = 0; i < dofTypes.Length; i++)
			{
				faXYZ[i, 0] = Nodes[i].X;
				faXYZ[i, 1] = Nodes[i].Y;
				faXYZ[i, 2] = Nodes[i].Z;
			}
			return faXYZ;
		}

		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofTypes;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		public virtual IMatrix StiffnessMatrix()
		{
			double[,,] afE = new double[iInt3, 6, 6];

			for (int i = 0; i < iInt3; i++)
			{
				IMatrixView constitutive = materialsAtGaussPoints[i].ConstitutiveMatrix;
				for (int j = 0; j < 6; j++)
				{
					for (int k = 0; k < 6; k++) afE[i, j, k] = constitutive[j, k];
				}
			}

			double[,] faXYZ = GetCoordinates();
			double[,] faDS = new double[iInt3, 24];
			double[,] faS = new double[iInt3, 8];
			double[,,] faB = new double[iInt3, 24, 6];
			double[] faDetJ = new double[iInt3];
			double[,,] faJ = new double[iInt3, 3, 3];
			double[] faWeight = new double[iInt3];
			double[] faK = new double[300];
			CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
			CalcH8K(ref iInt, afE, faB, faWeight, faK);
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromPackedRowMajorArray(faK));
		}
		public IMatrix PhysicsMatrix()
		{
			return StiffnessMatrix();
		}
		public IMatrix MassMatrix()
		{
			double[,] faXYZ = GetCoordinates();
			double[,] faDS = new double[iInt3, 24];
			double[,] faS = new double[iInt3, 8];
			double[,,] faB = new double[iInt3, 24, 6];
			double[] faDetJ = new double[iInt3];
			double[,,] faJ = new double[iInt3, 3, 3];
			double[] faWeight = new double[iInt3];
			double[] faM = new double[300];
			double fDensity = Density;
			CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
			CalcH8MLumped(ref iInt, ref fDensity, faWeight, faM);
			return SymmetricMatrix.CreateFromPackedRowMajorArray(faM);
		}


		public IMatrix DampingMatrix()
		{
			IMatrix k = StiffnessMatrix();
			IMatrix m = MassMatrix();
			k.LinearCombinationIntoThis(RayleighBeta, m, RayleighAlpha);
			return k;
		}

		private static int GetSolidDOFs()
		{
			int totalDisplacementDOFs = 0;
			foreach (IDofType[] nodalDOFs in dofTypes)
				foreach (IDofType dofType in nodalDOFs)
					if (dofType != PorousMediaDof.Pressure) totalDisplacementDOFs++;
			return totalDisplacementDOFs;
		}

		private static int GetAllDOFs()
		{
			int totalDisplacementDOFs = 0;
			foreach (IDofType[] nodalDOFs in dofTypes)
				foreach (IDofType dofType in nodalDOFs)
					totalDisplacementDOFs++;
			return totalDisplacementDOFs;
		}

		private static double[] ExtractSolidVector(double[] localVector)
		{
			int localPos = 0;
			int solidPos = 0;
			double[] solidVector = new double[GetSolidDOFs()];
			foreach (IDofType[] nodalDOFs in dofTypes)
				foreach (IDofType dofType in nodalDOFs)
				{
					if (dofType != PorousMediaDof.Pressure)
					{
						solidVector[solidPos] = localVector[localPos];
						solidPos++;
					}
					localPos++;
				}

			return solidVector;
		}

		private static double[] ExtractFluidVector(double[] localVector)
		{
			int localPos = 0;
			int fluidPos = 0;
			double[] fluidVector = new double[GetAllDOFs() - GetSolidDOFs()];
			foreach (IDofType[] nodalDOFs in dofTypes)
				foreach (IDofType dofType in nodalDOFs)
				{
					if (dofType == PorousMediaDof.Pressure)
					{
						fluidVector[fluidPos] = localVector[localPos];
						fluidPos++;
					}
					localPos++;
				}

			return fluidVector;
		}

		private static void ScatterFromSolidVector(double[] solidVector, double[] totalVector)
		{
			int localPos = 0;
			int solidPos = 0;
			foreach (IDofType[] nodalDOFs in dofTypes)
				foreach (IDofType dofType in nodalDOFs)
				{
					if (dofType != PorousMediaDof.Pressure)
					{
						totalVector[localPos] = solidVector[solidPos];
						solidPos++;
					}
					localPos++;
				}
		}

		private static void ScatterFromFluidVector(double[] fluidVector, double[] totalVector)
		{
			int localPos = 0;
			int fluidPos = 0;
			foreach (IDofType[] nodalDOFs in dofTypes)
				foreach (IDofType dofType in nodalDOFs)
				{
					if (dofType == PorousMediaDof.Pressure)
					{
						totalVector[localPos] = fluidVector[fluidPos];
						fluidPos++;
					}
					localPos++;
				}
		}

		// TODO: Replaced dstrains with strains (goat: this might go wrong)
		public Tuple<double[], double[]> CalculateResponse(double[] localDisplacements)
		{
			fluidDisplacements = ExtractFluidVector(localDisplacements);
			// H correction
			fluidDisplacements.ScaleIntoThis(-1.0);

			double[,] faXYZ = GetCoordinates();
			double[,] faDS = new double[iInt3, 24];
			double[,] faS = new double[iInt3, 8];
			double[,,] faB = new double[iInt3, 24, 6];
			double[] faDetJ = new double[iInt3];
			double[,,] faJ = new double[iInt3, 3, 3];
			double[] faWeight = new double[iInt3];
			double[,] faStrains = new double[iInt3, 6];
			double[,] fadStrains = new double[iInt3, 6];
			CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);

			//double[] soliddDisplacements = ExtractSolidVector(localdDisplacements);
			double[] solidDisplacements = ExtractSolidVector(localDisplacements);
			//CalcH8Strains(ref iInt, faB, soliddDisplacements, fadStrains);
			CalcH8Strains(ref iInt, faB, solidDisplacements, faStrains);

			double[] dStrains = new double[6];
			double[] strains = new double[6];
			for (int i = 0; i < materialsAtGaussPoints.Length; i++)
			{
				//for (int j = 0; j < 6; j++) dStrains[j] = fadStrains[i, j];
				for (int j = 0; j < 6; j++) strains[j] = faStrains[i, j];
				//materialsAtGaussPoints[i].UpdateConstitutiveMatrixAndEvaluateResponse(dStrains);
				lastStresses[i] = materialsAtGaussPoints[i].UpdateConstitutiveMatrixAndEvaluateResponse(strains);
			}

			return new Tuple<double[], double[]>(strains, lastStresses[materialsAtGaussPoints.Length - 1]);
		}



		public double[] CalculateResponseIntegralForLogging( double[] localDisplacements)
		{
			return CalculateResponseIntegral();
		}

		public double[] CalculateResponseIntegral()
		{
			double[,] faStresses = new double[iInt3, 6];
			for (int i = 0; i < materialsAtGaussPoints.Length; i++)
				for (int j = 0; j < 6; j++) faStresses[i, j] = lastStresses[i][j];

			double[,] faXYZ = GetCoordinates();
			double[,] faDS = new double[iInt3, 24];
			double[,] faS = new double[iInt3, 8];
			double[,,] faB = new double[iInt3, 24, 6];
			double[] faDetJ = new double[iInt3];
			double[,,] faJ = new double[iInt3, 3, 3];
			double[] faWeight = new double[iInt3];
			double[] solidForces = new double[24];
			CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
			CalcH8Forces(ref iInt, faB, faWeight, faStresses, solidForces);

			double[] faH = new double[36];
			double[] faPermeability = new double[] { Permeability, Permeability, Permeability, Permeability,
				Permeability, Permeability, Permeability, Permeability };
			CalcH20u8pH(ref iInt, faPermeability, faS, faB, faWeight, faH);

			//double[,] faQ = new double[8, 24];
			//double fPoreA = PoreA;
			//double[] faXw = new double[] { Xw, Xw, Xw, Xw, Xw, Xw, Xw, Xw };
			//CalcH8u8pQMinus(ref iInt, ref fPoreA, faXw, faB, faS, faWeight, faQ);
			//Matrix<double> q = (new Matrix<double>(faQ)).Transpose();

			//// Changed! Check for errors...
			////Vector<double> fluidDisplacements = new Vector<double>(ExtractFluidVector(localDisplacements));
			//var fluidDisplacements = ExtractFluidVector(localTotalDisplacements);

			//double[] solidAndFluidForces = solidForces;
			////double[] solidAndFluidForces = q * fluidDisplacements + (new Vector<double>(solidForces));
			//// H correction
			//fluidDisplacements.ScaleIntoThis(-1.0);
			double[] fluidDrags = SymmetricMatrix.CreateFromPackedRowMajorArray(faH).Multiply(fluidDisplacements);

			double[] totalForces = new double[GetAllDOFs()];
			ScatterFromFluidVector(fluidDrags, totalForces);
			ScatterFromSolidVector(solidForces, totalForces);
			return totalForces;
		}


		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads)
		//{
		//	var accelerations = new double[24];
		//	int index = 0;
		//	foreach (MassAccelerationLoad load in loads)
		//	{
		//		if (load.DOF == StructuralDof.TranslationX)
		//			for (int i = 0; i < 8; i++) accelerations[i * 3] += load.Amount;
		//		else if (load.DOF == StructuralDof.TranslationY)
		//			for (int i = 0; i < 8; i++) accelerations[i * 3 + 1] += load.Amount;
		//		else if (load.DOF == StructuralDof.TranslationZ)
		//			for (int i = 0; i < 8; i++) accelerations[i * 3 + 2] += load.Amount;
		//		else
		//			throw new InvalidOperationException("Cannot handle global acceleration for water pore when NOT translational.");
		//	}
		//	double[] solidForces = MassMatrix(element).Multiply(accelerations);

		//	double[,] faXYZ = GetCoordinates(element);
		//	double[,] faDS = new double[iInt3, 24];
		//	double[,] faS = new double[iInt3, 8];
		//	double[,,] faB = new double[iInt3, 24, 6];
		//	double[] faDetJ = new double[iInt3];
		//	double[,,] faJ = new double[iInt3, 3, 3];
		//	double[] faWeight = new double[iInt3];
		//	CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);

		//	double[] waterAcc = new double[3];
		//	foreach (MassAccelerationLoad load in loads)
		//	{
		//		if (load.DOF == StructuralDof.TranslationX) waterAcc[0] = load.Amount;
		//		else if (load.DOF == StructuralDof.TranslationY) waterAcc[1] = load.Amount;
		//		else if (load.DOF == StructuralDof.TranslationZ) waterAcc[2] = load.Amount;
		//		else throw new InvalidOperationException(
		//			"Cannot handle global acceleration for water pore when NOT translational.");
		//	}

		//	bool[] impermeableDOFs = new bool[8];
		//	for (int n = 0; n < element.Nodes.Count; ++n)
		//	{
		//		foreach (Constraint constraint in element.Nodes[n].Constraints)
		//		{
		//			if (constraint.DOF == PorousMediaDof.Pressure)
		//			{
		//				impermeableDOFs[n] = true;
		//			}
		//		}
		//	}

		//	double fD = FluidDensity;
		//	double[] faPermeability = new double[] { Permeability, Permeability, Permeability, Permeability,
		//		Permeability, Permeability, Permeability, Permeability };
		//	double[] faXw = new double[] { Xw, Xw, Xw, Xw, Xw, Xw, Xw, Xw };
		//	double[] faSw = new double[] { Saturation, Saturation, Saturation, Saturation, Saturation,
		//		Saturation, Saturation, Saturation };
		//	double[] fluidDrags = new double[8];
		//	CalcH20u8pForcesWaterAcc(ref iInt, impermeableDOFs, ref fD, faPermeability,
		//		faXw, faSw, waterAcc, faS, faB, faWeight, fluidDrags);

		//	// H correction
		//	for (int i = 0; i < 8; i++)
		//		fluidDrags[i] = -fluidDrags[i];
		//	double[] totalForces = new double[GetAllDOFs()];
		//	ScatterFromFluidVector(fluidDrags, totalForces);
		//	ScatterFromSolidVector(solidForces, totalForces);
		//	return totalForces;
		//}

		public double[] CalculateSolidForcesFromPorePressures(IElementType element, double[] porePressures)
		{
			IMatrix Q = CouplingMatrix();
			return Q.Multiply(porePressures, true);
		}

		//public bool ConstitutiveLawModified
		//{
		//	get
		//	{
		//		foreach (IContinuumMaterial3D material in materialsAtGaussPoints)
		//			if (material.IsCurrentStateDifferent()) return true;
		//		return false;
		//	}
		//}

		//public void ResetConstitutiveLawModified()
		//{
		//	foreach (IContinuumMaterial3D material in materialsAtGaussPoints) material.ResetModified();
		//}

		//public void ClearMaterialState()
		//{
		//	foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearState();
		//}

		public void SaveConstitutiveLawState(IHaveState externalState)
		{
			foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.CreateState();
		}

		//public void ClearMaterialStresses()
		//{
		//	foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
		//}

		#endregion

		#region IPorousFiniteElement Members

		public IMatrix PermeabilityMatrix()
		{
			double[,] faXYZ = GetCoordinates();
			double[,] faDS = new double[iInt3, 24];
			double[,] faS = new double[iInt3, 8];
			double[,,] faB = new double[iInt3, 24, 6];
			double[] faDetJ = new double[iInt3];
			double[,,] faJ = new double[iInt3, 3, 3];
			double[] faWeight = new double[iInt3];
			double[] faPermeability = new double[] { Permeability, Permeability, Permeability, Permeability,
				Permeability, Permeability, Permeability, Permeability };
			double[] faH = new double[36];
			CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
			CalcH20u8pH(ref iInt, faPermeability, faS, faB, faWeight, faH);
			return SymmetricMatrix.CreateFromPackedRowMajorArray(faH);
		}

		// Rows are fluid DOFs and columns are solid DOFs
		public IMatrix CouplingMatrix()
		{
			double[,] faXYZ = GetCoordinates();
			double[,] faDS = new double[iInt3, 24];
			double[,] faS = new double[iInt3, 8];
			double[,,] faB = new double[iInt3, 24, 6];
			double[] faDetJ = new double[iInt3];
			double[,,] faJ = new double[iInt3, 3, 3];
			double[] faWeight = new double[iInt3];
			//double fDensity = Density;
			double fPoreA = PoreA;
			double[] faXw = new double[] { Xw, Xw, Xw, Xw, Xw, Xw, Xw, Xw };
			double[,] faQ = new double[8, 24];
			CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
			CalcH8u8pQMinus(ref iInt, ref fPoreA, faXw, faB, faS, faWeight, faQ);
			return Matrix.CreateFromArray(faQ);
		}


		public IMatrix SaturationMatrix()
		{
			double[,] faXYZ = GetCoordinates();
			double[,] faDS = new double[iInt3, 24];
			double[,] faS = new double[iInt3, 8];
			double[,,] faB = new double[iInt3, 24, 6];
			double[] faDetJ = new double[iInt3];
			double[,,] faJ = new double[iInt3, 3, 3];
			double[] faWeight = new double[iInt3];
			double[] faXwDivQ = new double[] { Xw * QInv, Xw * QInv, Xw * QInv, Xw * QInv, Xw * QInv,
				Xw * QInv, Xw * QInv, Xw * QInv};
			double[] faSaturation = new double[36];
			CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
			CalcH20u8pS(ref iInt, faXwDivQ, faS, faWeight, faSaturation);
			return SymmetricMatrix.CreateFromPackedRowMajorArray(faSaturation);
		}


		#endregion
	}
}
