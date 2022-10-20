// --------------------------------------------------------------------------------------------------------------------
// <copyright file="VonMisesMaterial3DVisco.cs" company="National Technical University of Athens">
//   To be decided
// </copyright>
// <summary>
//   Class for 3D Von Mises Viscoplastic materials.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

using System;

using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;

namespace MGroup.Constitutive.Structural.Continuum
{
	/// <summary>
	///   Class for 3D Elastic Visco materials.
	/// </summary>
	/// <a href = "http://en.wikipedia.org/wiki/Von_Mises_yield_criterion">Wikipedia -Von Mises yield criterion</a>
	public class ElasticMaterial3DVisco : IIsotropicContinuumMaterial3D, ITransientConstitutiveLaw
	{
		private const string STRESS_X = "Stress X";
		private const string STRESS_Y = "Stress Y";
		private const string STRESS_Z = "Stress Z";
		private const string STRESS_XY = "Stress XY";
		private const string STRESS_XZ = "Stress XZ";
		private const string STRESS_YZ = "Stress YZ";

		private const string STRAIN_X = "Strain X";
		private const string STRAIN_Y = "Strain Y";
		private const string STRAIN_Z = "Strain Z";
		private const string STRAIN_XY = "Strain XY";
		private const string STRAIN_XZ = "Strain XZ";
		private const string STRAIN_YZ = "Strain YZ";

		private const string VISCOUS_STRAIN_X_0 = "Viscous Strain X of maxwell element 0";
		private const string VISCOUS_STRAIN_Y_0 = "Viscous Strain Y of maxwell element 0";
		private const string VISCOUS_STRAIN_Z_0 = "Viscous Strain Z of maxwell element 0";
		private const string VISCOUS_STRAIN_XY_0 = "Viscous Strain XY of maxwell element 0";
		private const string VISCOUS_STRAIN_XZ_0 = "Viscous Strain XZ of maxwell element 0";
		private const string VISCOUS_STRAIN_YZ_0 = "Viscous Strain YZ of maxwell element 0";
		private const string VOLUMETRIC_VISCOUS_STRAIN_0 = "Volumetric Viscous Strain of maxwell element 0";

		private const string VISCOUS_STRAIN_X_1 = "Viscous Strain X of maxwell element 1";
		private const string VISCOUS_STRAIN_Y_1 = "Viscous Strain Y of maxwell element 1";
		private const string VISCOUS_STRAIN_Z_1 = "Viscous Strain Z of maxwell element 1";
		private const string VISCOUS_STRAIN_XY_1 = "Viscous Strain XY of maxwell element 1";
		private const string VISCOUS_STRAIN_XZ_1 = "Viscous Strain XZ of maxwell element 1";
		private const string VISCOUS_STRAIN_YZ_1 = "Viscous Strain YZ of maxwell element 1";
		private const string VOLUMETRIC_VISCOUS_STRAIN_1 = "Volumetric Viscous Strain of maxwell element 1";

		private const string VISCOUS_STRAIN_X_2 = "Viscous Strain X of maxwell element 2";
		private const string VISCOUS_STRAIN_Y_2 = "Viscous Strain Y of maxwell element 2";
		private const string VISCOUS_STRAIN_Z_2 = "Viscous Strain Z of maxwell element 2";
		private const string VISCOUS_STRAIN_XY_2 = "Viscous Strain XY of maxwell element 2";
		private const string VISCOUS_STRAIN_XZ_2 = "Viscous Strain XZ of maxwell element 2";
		private const string VISCOUS_STRAIN_YZ_2 = "Viscous Strain YZ of maxwell element 2";
		private const string VOLUMETRIC_VISCOUS_STRAIN_2 = "Volumetric Viscous Strain of maxwell element 2";

		private const string VISCOUS_STRAIN_X_3 = "Viscous Strain X of maxwell element 3";
		private const string VISCOUS_STRAIN_Y_3 = "Viscous Strain Y of maxwell element 3";
		private const string VISCOUS_STRAIN_Z_3 = "Viscous Strain Z of maxwell element 3";
		private const string VISCOUS_STRAIN_XY_3 = "Viscous Strain XY of maxwell element 3";
		private const string VISCOUS_STRAIN_XZ_3 = "Viscous Strain XZ of maxwell element 3";
		private const string VISCOUS_STRAIN_YZ_3 = "Viscous Strain YZ of maxwell element 3";
		private const string VOLUMETRIC_VISCOUS_STRAIN_3 = "Volumetric Viscous Strain of maxwell element 3";

		private const string VISCOUS_STRAIN_X_4 = "Viscous Strain X of maxwell element 4";
		private const string VISCOUS_STRAIN_Y_4 = "Viscous Strain Y of maxwell element 4";
		private const string VISCOUS_STRAIN_Z_4 = "Viscous Strain Z of maxwell element 4";
		private const string VISCOUS_STRAIN_XY_4 = "Viscous Strain XY of maxwell element 4";
		private const string VISCOUS_STRAIN_XZ_4 = "Viscous Strain XZ of maxwell element 4";
		private const string VISCOUS_STRAIN_YZ_4 = "Viscous Strain YZ of maxwell element 4";
		private const string VOLUMETRIC_VISCOUS_STRAIN_4 = "Volumetric Viscous Strain of maxwell element 4";

		private const string VISCOUS_STRAIN_X_5 = "Viscous Strain X of maxwell element 5";
		private const string VISCOUS_STRAIN_Y_5 = "Viscous Strain Y of maxwell element 5";
		private const string VISCOUS_STRAIN_Z_5 = "Viscous Strain Z of maxwell element 5";
		private const string VISCOUS_STRAIN_XY_5 = "Viscous Strain XY of maxwell element 5";
		private const string VISCOUS_STRAIN_XZ_5 = "Viscous Strain XZ of maxwell element 5";
		private const string VISCOUS_STRAIN_YZ_5 = "Viscous Strain YZ of maxwell element 5";
		private const string VOLUMETRIC_VISCOUS_STRAIN_5 = "Volumetric Viscous Strain of maxwell element 5";

		private const string VISCOUS_STRAIN_X_6 = "Viscous Strain X of maxwell element 6";
		private const string VISCOUS_STRAIN_Y_6 = "Viscous Strain Y of maxwell element 6";
		private const string VISCOUS_STRAIN_Z_6 = "Viscous Strain Z of maxwell element 6";
		private const string VISCOUS_STRAIN_XY_6 = "Viscous Strain XY of maxwell element 6";
		private const string VISCOUS_STRAIN_XZ_6 = "Viscous Strain XZ of maxwell element 6";
		private const string VISCOUS_STRAIN_YZ_6 = "Viscous Strain YZ of maxwell element 6";
		private const string VOLUMETRIC_VISCOUS_STRAIN_6 = "Volumetric Viscous Strain of maxwell element 6";

		private const string VISCOUS_STRAIN_X_7 = "Viscous Strain X of maxwell element 7";
		private const string VISCOUS_STRAIN_Y_7 = "Viscous Strain Y of maxwell element 7";
		private const string VISCOUS_STRAIN_Z_7 = "Viscous Strain Z of maxwell element 7";
		private const string VISCOUS_STRAIN_XY_7 = "Viscous Strain XY of maxwell element 7";
		private const string VISCOUS_STRAIN_XZ_7 = "Viscous Strain XZ of maxwell element 7";
		private const string VISCOUS_STRAIN_YZ_7 = "Viscous Strain YZ of maxwell element 7";
		private const string VOLUMETRIC_VISCOUS_STRAIN_7 = "Volumetric Viscous Strain of maxwell element 7";

		private const string VISCOUS_STRAIN_X_8 = "Viscous Strain X of maxwell element 8";
		private const string VISCOUS_STRAIN_Y_8 = "Viscous Strain Y of maxwell element 8";
		private const string VISCOUS_STRAIN_Z_8 = "Viscous Strain Z of maxwell element 8";
		private const string VISCOUS_STRAIN_XY_8 = "Viscous Strain XY of maxwell element 8";
		private const string VISCOUS_STRAIN_XZ_8 = "Viscous Strain XZ of maxwell element 8";
		private const string VISCOUS_STRAIN_YZ_8 = "Viscous Strain YZ of maxwell element 8";
		private const string VOLUMETRIC_VISCOUS_STRAIN_8 = "Volumetric Viscous Strain of maxwell element 8";

		private const string VISCOUS_STRAIN_X_9 = "Viscous Strain X of maxwell element 9";
		private const string VISCOUS_STRAIN_Y_9 = "Viscous Strain Y of maxwell element 9";
		private const string VISCOUS_STRAIN_Z_9 = "Viscous Strain Z of maxwell element 9";
		private const string VISCOUS_STRAIN_XY_9 = "Viscous Strain XY of maxwell element 9";
		private const string VISCOUS_STRAIN_XZ_9 = "Viscous Strain XZ of maxwell element 9";
		private const string VISCOUS_STRAIN_YZ_9 = "Viscous Strain YZ of maxwell element 9";
		private const string VOLUMETRIC_VISCOUS_STRAIN_9 = "Volumetric Viscous Strain of maxwell element 9";

		private const string VISCOUS_STRAIN_X_10 = "Viscous Strain X of maxwell element 10";
		private const string VISCOUS_STRAIN_Y_10 = "Viscous Strain Y of maxwell element 10";
		private const string VISCOUS_STRAIN_Z_10 = "Viscous Strain Z of maxwell element 10";
		private const string VISCOUS_STRAIN_XY_10 = "Viscous Strain XY of maxwell element 10";
		private const string VISCOUS_STRAIN_XZ_10 = "Viscous Strain XZ of maxwell element 10";
		private const string VISCOUS_STRAIN_YZ_10 = "Viscous Strain YZ of maxwell element 10";
		private const string VOLUMETRIC_VISCOUS_STRAIN_10 = "Volumetric Viscous Strain of maxwell element 10";

		private const string VISCOUS_STRAIN_X_11 = "Viscous Strain X of maxwell element 11";
		private const string VISCOUS_STRAIN_Y_11 = "Viscous Strain Y of maxwell element 11";
		private const string VISCOUS_STRAIN_Z_11 = "Viscous Strain Z of maxwell element 11";
		private const string VISCOUS_STRAIN_XY_11 = "Viscous Strain XY of maxwell element 11";
		private const string VISCOUS_STRAIN_XZ_11 = "Viscous Strain XZ of maxwell element 11";
		private const string VISCOUS_STRAIN_YZ_11 = "Viscous Strain YZ of maxwell element 11";
		private const string VOLUMETRIC_VISCOUS_STRAIN_11 = "Volumetric Viscous Strain of maxwell element 11";

		private const string VISCOUS_STRAIN_X_12 = "Viscous Strain X of maxwell element 12";
		private const string VISCOUS_STRAIN_Y_12 = "Viscous Strain Y of maxwell element 12";
		private const string VISCOUS_STRAIN_Z_12 = "Viscous Strain Z of maxwell element 12";
		private const string VISCOUS_STRAIN_XY_12 = "Viscous Strain XY of maxwell element 12";
		private const string VISCOUS_STRAIN_XZ_12 = "Viscous Strain XZ of maxwell element 12";
		private const string VISCOUS_STRAIN_YZ_12 = "Viscous Strain YZ of maxwell element 12";
		private const string VOLUMETRIC_VISCOUS_STRAIN_12 = "Volumetric Viscous Strain of maxwell element 12";

		private GenericConstitutiveLawState currentState;

		private double timePrevious = double.NaN, timeCurrent;

		/// <summary>
		///   The Poisson ratio value of an incompressible solid.
		/// </summary>
		private const double PoissonRatioForIncompressibleSolid = 0.5;

		/// <summary>
		///   The total number of strains.
		/// </summary>
		private const int TotalStrains = 6;

		/// <summary>
		///   The total number of stresses.
		/// </summary>
		private const int TotalStresses = TotalStrains;

		/// <summary>
		///   An array needed for the formulation of the consistent constitutive matrix.
		/// </summary>
		private static readonly double[,] SupportiveMatrixForConsistentConstitutiveMatrix = new[,]
			{
				{  0, -1, -1, 0,   0,   0   },
				{ -1,  0, -1, 0,   0,   0   },
				{ -1, -1,  0, 0,   0,   0   },
				{  0,  0,  0, 0.5, 0,   0   },
				{  0,  0,  0, 0,   0.5, 0   },
				{  0,  0,  0, 0,   0,   0.5 }
			};

		/// <summary>
		///   An array needed for the formulation of the consistent constitutive matrix.
		/// </summary>
		private static readonly double[,] SupportiveMatrixForConsistentConstitutiveMatrixAIdM = new[,]
			{
			{ 2.0 / 3.0, - 1.0 / 3.0, - 1.0 / 3.0, 0.0, 0.0, 0.0 },
			{-1.0/3.0, 2.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0 },
			{-1.0/3.0, -1.0/3.0, 2.0/3.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.5, 0.0, 0.0},
			{ 0.0, 0.0, 0.0, 0.0, 0.5, 0.0},
			{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.5}
			};

		/// <summary>
		///   An array needed for the formulation of the consistent constitutive matrix.
		/// </summary>
		private static readonly double[,] SupportiveMatrixForConsistentConstitutiveMatrixIIM = new[,]
			{
				{1.0, 1.0, 1.0, 0, 0, 0},
				{1.0, 1.0, 1.0, 0, 0, 0},
				{1.0, 1.0, 1.0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0}
			};
		/// <summary>
		///   An array needed for the formulation of the consistent constitutive matrix.
		/// </summary>
		private static readonly double[,] SupportiveMatrixForConsistentConstitutiveMatrixI4M = new[,]
			{
				{1.0, 0.0, 0.0, 0, 0, 0},
				{0.0, 1.0, 0.0, 0, 0, 0},
				{0.0, 0.0, 1.0, 0, 0, 0},
				{0, 0, 0, 0.5, 0, 0},
				{0, 0, 0, 0, 0.5, 0},
				{0, 0, 0, 0, 0, 0.5}
			};

		/// <summary>
		///   An array needed for the formulation of deviatoric stress or strain arrays.
		/// </summary>
		private static readonly double[] SupportiveArrayForDeviatoricFormulation = { 1d, 1d, 1d, 0d, 0d, 0d };

		/// <summary>
		///   An array needed for the formulation of engineering strains array.
		/// </summary>
		private static readonly double[] SupportiveArrayForEngineeringStrains = { 1d, 1d, 1d, 0.5, 0.5, 0.5 };


		/// <summary>
		///   The constitutive matrix of the material while still in the elastic region.
		/// </summary>
		private readonly Matrix elasticConstitutiveMatrix;

		/// <summary>
		///   The Poisson ratio.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Poisson%27s_ratio">Wikipedia - Poisson's Ratio</a>
		/// </remarks>
		private double poissonRatio;

		/// <summary>
		///   The instantaneous shear modulus.
		/// </summary>
		private readonly double instantaneousShearModulus;

		/// <summary>
		///   The instantaneous bulk modulus.
		/// </summary>
		private readonly double instantaneousBulkModulus;

		/// <summary>
		///   The young modulus.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Young%27s_modulus">Wikipedia - Young's Modulus</a>
		/// </remarks>
		private double youngModulus;

		/// <summary>
		///   Moduli time scale parameter:
		///   "LONG_TERM" or "INSTANTANEOUS".
		///   Default is "LONG_TERM"
		/// </summary>
		private string moduliTimeScale;

		/// <summary>
		///   Viscocity data:
		///   Prony series parameters
		///   [gnormi, knormi,taui].
		/// </summary>
		private double[,] viscoData;

		/// <summary>
		///   Viscocity data:
		///   Prony series parameters
		///   numberOfMaxwelElements.
		///   Number of Maxwell elements
		/// </summary>
		private int numberOfMaxwelElements;

		/// <summary>
		///   Viscocity data:
		///   Prony series parameters of Maxwell elements
		///   gnormi : normalized shear modulus 
		///   knormi : normalized bulk modulus 
		///   taui   : normalized relaxation time
		/// </summary>
		private double[] gnormi, knormi, taui;

		/// <summary>
		///   The increment of time.
		///   This is related to the history of incemental load. It is set 1.0 here for quick testing
		/// </summary>
		private double dtime { get => double.IsNaN(timePrevious)  ? 1d : timeCurrent - timePrevious; }

		/// <summary>
		///   The constitutive matrix of the material.
		/// </summary>
		private Matrix constitutiveMatrix;

		/// <summary>
		///   The array of incremental strains.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Deformation_%28engineering%29">Deformation (engineering)</a>
		/// </remarks>
		private double[] incrementalStrains = new double[6];

		/// <summary>
		///   Indicates whether this <see cref = "IStructuralMaterial" /> is modified.
		/// </summary>
		private bool modified;

		/// <summary>
		///   The array of stresses.
		/// </summary>
		private double[] stresses = new double[6];

		/// <summary>
		///   The array of new stresses.
		/// </summary>
		private double[] stressesNew = new double[6];

		/// <summary>
		///   The array of strains.
		/// </summary>
		private double[] strains = new double[6];

		/// <summary>
		///   The array of new strains.
		/// </summary>
		private double[] strainsNew = new double[6];

		/// <summary>
		///   The array of state variables
		///   maximum number Of Maxwel Elements = 13
		/// </summary>
		private double[] stateVariables = new double[13 * 7];

		/// <summary>
		///   The array of new state variables
		///   maximum number Of Maxwel Elements = 13
		/// </summary>
		private double[] stateVariablesNew = new double[13 * 7];

		/// <summary>
		///   Initializes a new instance of the <see cref = "VonMisesMaterial3DVisco" /> class.
		/// </summary>
		/// <param name = "youngModulus">
		///   The young modulus.
		/// </param>
		/// <param name = "poissonRatio">
		///   The Poisson ratio.
		/// </param>
		/// <exception cref = "ArgumentException"> When Poisson ratio is equal to 0.5.</exception>
		public ElasticMaterial3DVisco(double youngModulus, double poissonRatio, string moduliTimeScale, double[,] viscoData)
		{
			this.youngModulus = youngModulus;

			if (poissonRatio == PoissonRatioForIncompressibleSolid)
			{
				throw new ArgumentException(
					"Poisson ratio cannot be" + PoissonRatioForIncompressibleSolid + "(incompressible solid)");
			}

			this.poissonRatio = poissonRatio;

			var shearModulusTemp = this.youngModulus / (2.0 * (1.0 + this.poissonRatio));
			var bulkModulusTemp = this.youngModulus / (3.0 * (1.0 - (2.0 * this.poissonRatio)));

			this.moduliTimeScale = moduliTimeScale;
			this.viscoData = viscoData;
			this.numberOfMaxwelElements = this.viscoData.GetLength(0);

			this.gnormi = new double[numberOfMaxwelElements];
			this.knormi = new double[numberOfMaxwelElements];
			this.taui = new double[numberOfMaxwelElements];

			double gnormiSum = 0.0;
			double knormiSum = 0.0;

			for (int iNME = 0; iNME < numberOfMaxwelElements; iNME++)
			{
				this.gnormi[iNME] = viscoData[iNME, 0];
				this.knormi[iNME] = viscoData[iNME, 1];
				this.taui[iNME] = viscoData[iNME, 2];

				gnormiSum += this.gnormi[iNME];
				knormiSum += this.knormi[iNME];
			}

			//Static Analysis based in instantaneous material response so:
			if (moduliTimeScale.Equals("LONG_TERM"))
			{
				shearModulusTemp = shearModulusTemp / (1.0 - gnormiSum);        //instant Shear modulus
				bulkModulusTemp = bulkModulusTemp / (1.0 - knormiSum); //instant Bulk modulus
			}
			instantaneousShearModulus = shearModulusTemp;
			instantaneousBulkModulus = bulkModulusTemp;

			//TODO: Initialization is done prior to setting current time
			var shearModulus = 0.0;
			var bulkModulus = 0.0;
			for (int iNME = 0; iNME < numberOfMaxwelElements; iNME++)
			{
				var dT = dtime / taui[iNME];
				if (dT > 1E-7)
				{
					shearModulus += gnormi[iNME] * (1.0 / dT) * (dT + Math.Exp(-dT) - 1.0);
					bulkModulus += knormi[iNME] * (1.0 / dT) * (dT + Math.Exp(-dT) - 1.0);
				}
				else
				{
					shearModulus += (1.0 / 2.0) * gnormi[iNME] * dT;
					bulkModulus += (1.0 / 2.0) * knormi[iNME] * dT;
				}

			}

			shearModulus = instantaneousShearModulus * (1.0 - shearModulus);
			bulkModulus = instantaneousBulkModulus * (1.0 - bulkModulus);

			this.elasticConstitutiveMatrix = Matrix.CreateZero(TotalStresses, TotalStrains);
			for (int k1 = 0; k1 < TotalStresses; k1++)
			{
				for (int k2 = 0; k2 < TotalStresses; k2++)
				{
					this.elasticConstitutiveMatrix[k2, k1] = (2d * shearModulus * SupportiveMatrixForConsistentConstitutiveMatrixI4M[k2, k1]) + (bulkModulus - ((2d / 3d) * shearModulus)) * SupportiveMatrixForConsistentConstitutiveMatrixIIM[k2, k1];
				}
			}

			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(STRESS_X, 0d),
				(STRESS_Y, 0d),
				(STRESS_Z, 0d),
				(STRESS_XY, 0d),
				(STRESS_XZ, 0d),
				(STRESS_YZ, 0d),
				(STRAIN_X, 0d),
				(STRAIN_Y, 0d),
				(STRAIN_Z, 0d),
				(STRAIN_XY, 0d),
				(STRAIN_XZ, 0d),
				(STRAIN_YZ, 0d),
				(VISCOUS_STRAIN_X_0, 0d),
				(VISCOUS_STRAIN_Y_0, 0d),
				(VISCOUS_STRAIN_Z_0, 0d),
				(VISCOUS_STRAIN_XY_0, 0d),
				(VISCOUS_STRAIN_XZ_0, 0d),
				(VISCOUS_STRAIN_YZ_0, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_0, 0d),
				(VISCOUS_STRAIN_X_1, 0d),
				(VISCOUS_STRAIN_Y_1, 0d),
				(VISCOUS_STRAIN_Z_1, 0d),
				(VISCOUS_STRAIN_XY_1, 0d),
				(VISCOUS_STRAIN_XZ_1, 0d),
				(VISCOUS_STRAIN_YZ_1, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_1, 0d),
				(VISCOUS_STRAIN_X_2, 0d),
				(VISCOUS_STRAIN_Y_2, 0d),
				(VISCOUS_STRAIN_Z_2, 0d),
				(VISCOUS_STRAIN_XY_2, 0d),
				(VISCOUS_STRAIN_XZ_2, 0d),
				(VISCOUS_STRAIN_YZ_2, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_2, 0d),
				(VISCOUS_STRAIN_X_3, 0d),
				(VISCOUS_STRAIN_Y_3, 0d),
				(VISCOUS_STRAIN_Z_3, 0d),
				(VISCOUS_STRAIN_XY_3, 0d),
				(VISCOUS_STRAIN_XZ_3, 0d),
				(VISCOUS_STRAIN_YZ_3, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_3, 0d),
				(VISCOUS_STRAIN_X_4, 0d),
				(VISCOUS_STRAIN_Y_4, 0d),
				(VISCOUS_STRAIN_Z_4, 0d),
				(VISCOUS_STRAIN_XY_4, 0d),
				(VISCOUS_STRAIN_XZ_4, 0d),
				(VISCOUS_STRAIN_YZ_4, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_4, 0d),
				(VISCOUS_STRAIN_X_5, 0d),
				(VISCOUS_STRAIN_Y_5, 0d),
				(VISCOUS_STRAIN_Z_5, 0d),
				(VISCOUS_STRAIN_XY_5, 0d),
				(VISCOUS_STRAIN_XZ_5, 0d),
				(VISCOUS_STRAIN_YZ_5, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_5, 0d),
				(VISCOUS_STRAIN_X_6, 0d),
				(VISCOUS_STRAIN_Y_6, 0d),
				(VISCOUS_STRAIN_Z_6, 0d),
				(VISCOUS_STRAIN_XY_6, 0d),
				(VISCOUS_STRAIN_XZ_6, 0d),
				(VISCOUS_STRAIN_YZ_6, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_6, 0d),
				(VISCOUS_STRAIN_X_7, 0d),
				(VISCOUS_STRAIN_Y_7, 0d),
				(VISCOUS_STRAIN_Z_7, 0d),
				(VISCOUS_STRAIN_XY_7, 0d),
				(VISCOUS_STRAIN_XZ_7, 0d),
				(VISCOUS_STRAIN_YZ_7, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_7, 0d),
				(VISCOUS_STRAIN_X_8, 0d),
				(VISCOUS_STRAIN_Y_8, 0d),
				(VISCOUS_STRAIN_Z_8, 0d),
				(VISCOUS_STRAIN_XY_8, 0d),
				(VISCOUS_STRAIN_XZ_8, 0d),
				(VISCOUS_STRAIN_YZ_8, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_8, 0d),
				(VISCOUS_STRAIN_X_9, 0d),
				(VISCOUS_STRAIN_Y_9, 0d),
				(VISCOUS_STRAIN_Z_9, 0d),
				(VISCOUS_STRAIN_XY_9, 0d),
				(VISCOUS_STRAIN_XZ_9, 0d),
				(VISCOUS_STRAIN_YZ_9, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_9, 0d),
				(VISCOUS_STRAIN_X_10, 0d),
				(VISCOUS_STRAIN_Y_10, 0d),
				(VISCOUS_STRAIN_Z_10, 0d),
				(VISCOUS_STRAIN_XY_10, 0d),
				(VISCOUS_STRAIN_XZ_10, 0d),
				(VISCOUS_STRAIN_YZ_10, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_10, 0d),
				(VISCOUS_STRAIN_X_11, 0d),
				(VISCOUS_STRAIN_Y_11, 0d),
				(VISCOUS_STRAIN_Z_11, 0d),
				(VISCOUS_STRAIN_XY_11, 0d),
				(VISCOUS_STRAIN_XZ_11, 0d),
				(VISCOUS_STRAIN_YZ_11, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_11, 0d),
				(VISCOUS_STRAIN_X_12, 0d),
				(VISCOUS_STRAIN_Y_12, 0d),
				(VISCOUS_STRAIN_Z_12, 0d),
				(VISCOUS_STRAIN_XY_12, 0d),
				(VISCOUS_STRAIN_XZ_12, 0d),
				(VISCOUS_STRAIN_YZ_12, 0d),
				(VOLUMETRIC_VISCOUS_STRAIN_12, 0d),
			});
		}

		public double[] Coordinates { get; set; }

		/// <summary>
		///   Gets the constitutive matrix.
		/// </summary>
		/// <value>
		///   The constitutive matrix.
		/// </value>
		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (this.constitutiveMatrix == null) UpdateConstitutiveMatrixAndEvaluateResponse(new double[6]);
				return constitutiveMatrix;
			}
		}

		/// <summary>
		///   Gets the ID of the material.
		/// </summary>
		/// <value>
		///   The id.
		/// </value>
		public int ID => 1;

		/// <summary>
		///   Gets the incremental strains of the finite element's material.
		/// </summary>
		/// <value>
		///   The incremental strains.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Deformation_%28engineering%29">Deformation (engineering)</a>
		/// </remarks>
		public double[] IncrementalStrains => this.incrementalStrains;

		/// <summary>
		///   Gets a value indicating whether this <see cref = "IStructuralMaterial" /> is modified.
		/// </summary>
		/// <value>
		///   <c>true</c> if modified; otherwise, <c>false</c>.
		/// </value>
		public bool IsCurrentStateDifferent() => modified;

		/// <summary>
		///   Gets the Poisson ratio.
		/// </summary>
		/// <value>
		///   The Poisson ratio.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Poisson%27s_ratio">Wikipedia - Poisson's Ratio</a>
		/// </remarks>
		public double PoissonRatio
		{
			get
			{
				return this.poissonRatio;
			}
			set
			{
				this.poissonRatio = value;
			}
		}

		/// <summary>
		///   Gets the stresses of the finite element's material.
		/// </summary>
		/// <value>
		///   The stresses.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Stress_%28mechanics%29">Stress (mechanics)</a>
		/// </remarks>
		public double[] Stresses => this.stressesNew;

		/// <summary>
		///   Gets the strains of the finite element's material.
		/// </summary>
		/// <value>
		///   The stains.
		/// </value>
		public double[] Strains => this.strainsNew;

		/// <summary>
		///   Gets the Young's Modulus.
		/// </summary>
		/// <value>
		///   The young modulus.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Young%27s_modulus">Wikipedia - Young's Modulus</a>
		/// </remarks>
		public double YoungModulus
		{
			get => this.youngModulus;
			set => this.youngModulus = value;
		}

		/// <summary>
		///   Creates a new object that is a copy of the current instance.
		/// </summary>
		/// <returns>
		///   A new object that is a copy of this instance.
		/// </returns>
		public object Clone()
		{
			return new ElasticMaterial3DVisco(this.youngModulus, this.poissonRatio, this.moduliTimeScale, this.viscoData)
			{
				modified = this.IsCurrentStateDifferent(),
				incrementalStrains = incrementalStrains.Copy(),
				stresses = stresses.Copy(),
				strains = strains.Copy(),
				stateVariables = stateVariables.Copy(),
			};
		}

		/// <summary>
		///   Resets the indicator of whether the material is modified.
		/// </summary>
		public void ResetModified() => this.modified = false;

		/// <summary>
		///   Clears the stresses of the element's material.
		/// </summary>
		public void ClearStresses()
		{
			stresses.Clear();
			stressesNew.Clear();
		}

		public void ClearState()
		{
			modified = false;
			constitutiveMatrix.Clear();
			incrementalStrains.Clear();
			//strains.Clear();
			//strainsNew.Clear();
			stresses.Clear();
			stressesNew.Clear();
			stateVariables.Clear();
			stateVariablesNew.Clear();
		}

		/// <summary>
		///   Saves the state of the element's material.
		/// </summary>
		public GenericConstitutiveLawState CreateState()
		{
			stresses.CopyFrom(stressesNew);
			strains.CopyFrom(strainsNew);
			stateVariables.CopyFrom(stateVariablesNew);

			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(STRESS_X, stresses[0]),
				(STRESS_Y, stresses[1]),
				(STRESS_Z, stresses[2]),
				(STRESS_XY, stresses[3]),
				(STRESS_XZ, stresses[4]),
				(STRESS_YZ, stresses[5]),
				(STRAIN_X, strains[0]),
				(STRAIN_Y, strains[1]),
				(STRAIN_Z, strains[2]),
				(STRAIN_XY, strains[3]),
				(STRAIN_XZ, strains[4]),
				(STRAIN_YZ, strains[5]),
				(VISCOUS_STRAIN_X_0, stateVariables[0 + 7*0]),
				(VISCOUS_STRAIN_Y_0, stateVariables[1 + 7*0]),
				(VISCOUS_STRAIN_Z_0, stateVariables[2 + 7*0]),
				(VISCOUS_STRAIN_XY_0, stateVariables[3 + 7*0]),
				(VISCOUS_STRAIN_XZ_0, stateVariables[4 + 7*0]),
				(VISCOUS_STRAIN_YZ_0, stateVariables[5 + 7*0]),
				(VOLUMETRIC_VISCOUS_STRAIN_0, stateVariables[6 + 7*0]),
				(VISCOUS_STRAIN_X_1, stateVariables[0 + 7*1]),
				(VISCOUS_STRAIN_Y_1, stateVariables[1 + 7*1]),
				(VISCOUS_STRAIN_Z_1, stateVariables[2 + 7*1]),
				(VISCOUS_STRAIN_XY_1, stateVariables[3 + 7*1]),
				(VISCOUS_STRAIN_XZ_1, stateVariables[4 + 7*1]),
				(VISCOUS_STRAIN_YZ_1, stateVariables[5 + 7*1]),
				(VOLUMETRIC_VISCOUS_STRAIN_1, stateVariables[6 + 7*1]),
				(VISCOUS_STRAIN_X_2, stateVariables[0 + 7*2]),
				(VISCOUS_STRAIN_Y_2, stateVariables[1 + 7*2]),
				(VISCOUS_STRAIN_Z_2, stateVariables[2 + 7*2]),
				(VISCOUS_STRAIN_XY_2, stateVariables[3 + 7*2]),
				(VISCOUS_STRAIN_XZ_2, stateVariables[4 + 7*2]),
				(VISCOUS_STRAIN_YZ_2, stateVariables[5 + 7*2]),
				(VOLUMETRIC_VISCOUS_STRAIN_2, stateVariables[6 + 7*2]),
				(VISCOUS_STRAIN_X_3, stateVariables[0 + 7*3]),
				(VISCOUS_STRAIN_Y_3, stateVariables[1 + 7*3]),
				(VISCOUS_STRAIN_Z_3, stateVariables[2 + 7*3]),
				(VISCOUS_STRAIN_XY_3, stateVariables[3 + 7*3]),
				(VISCOUS_STRAIN_XZ_3, stateVariables[4 + 7*3]),
				(VISCOUS_STRAIN_YZ_3, stateVariables[5 + 7*3]),
				(VOLUMETRIC_VISCOUS_STRAIN_3, stateVariables[6 + 7*3]),
				(VISCOUS_STRAIN_X_4, stateVariables[0 + 7*4]),
				(VISCOUS_STRAIN_Y_4, stateVariables[1 + 7*4]),
				(VISCOUS_STRAIN_Z_4, stateVariables[2 + 7*4]),
				(VISCOUS_STRAIN_XY_4, stateVariables[3 + 7*4]),
				(VISCOUS_STRAIN_XZ_4, stateVariables[4 + 7*4]),
				(VISCOUS_STRAIN_YZ_4, stateVariables[5 + 7*4]),
				(VOLUMETRIC_VISCOUS_STRAIN_4, stateVariables[6 + 7*4]),
				(VISCOUS_STRAIN_X_5, stateVariables[0 + 7*5]),
				(VISCOUS_STRAIN_Y_5, stateVariables[1 + 7*5]),
				(VISCOUS_STRAIN_Z_5, stateVariables[2 + 7*5]),
				(VISCOUS_STRAIN_XY_5, stateVariables[3 + 7*5]),
				(VISCOUS_STRAIN_XZ_5, stateVariables[4 + 7*5]),
				(VISCOUS_STRAIN_YZ_5, stateVariables[5 + 7*5]),
				(VOLUMETRIC_VISCOUS_STRAIN_5, stateVariables[6 + 7*5]),
				(VISCOUS_STRAIN_X_6, stateVariables[0 + 7*6]),
				(VISCOUS_STRAIN_Y_6, stateVariables[1 + 7*6]),
				(VISCOUS_STRAIN_Z_6, stateVariables[2 + 7*6]),
				(VISCOUS_STRAIN_XY_6, stateVariables[3 + 7*6]),
				(VISCOUS_STRAIN_XZ_6, stateVariables[4 + 7*6]),
				(VISCOUS_STRAIN_YZ_6, stateVariables[5 + 7*6]),
				(VOLUMETRIC_VISCOUS_STRAIN_6, stateVariables[6 + 7*6]),
				(VISCOUS_STRAIN_X_7, stateVariables[0 + 7*7]),
				(VISCOUS_STRAIN_Y_7, stateVariables[1 + 7*7]),
				(VISCOUS_STRAIN_Z_7, stateVariables[2 + 7*7]),
				(VISCOUS_STRAIN_XY_7, stateVariables[3 + 7*7]),
				(VISCOUS_STRAIN_XZ_7, stateVariables[4 + 7*7]),
				(VISCOUS_STRAIN_YZ_7, stateVariables[5 + 7*7]),
				(VOLUMETRIC_VISCOUS_STRAIN_7, stateVariables[6 + 7*7]),
				(VISCOUS_STRAIN_X_8, stateVariables[0 + 7*8]),
				(VISCOUS_STRAIN_Y_8, stateVariables[1 + 7*8]),
				(VISCOUS_STRAIN_Z_8, stateVariables[2 + 7*8]),
				(VISCOUS_STRAIN_XY_8, stateVariables[3 + 7*8]),
				(VISCOUS_STRAIN_XZ_8, stateVariables[4 + 7*8]),
				(VISCOUS_STRAIN_YZ_8, stateVariables[5 + 7*8]),
				(VOLUMETRIC_VISCOUS_STRAIN_8, stateVariables[6 + 7*8]),
				(VISCOUS_STRAIN_X_9, stateVariables[0 + 7*9]),
				(VISCOUS_STRAIN_Y_9, stateVariables[1 + 7*9]),
				(VISCOUS_STRAIN_Z_9, stateVariables[2 + 7*9]),
				(VISCOUS_STRAIN_XY_9, stateVariables[3 + 7*9]),
				(VISCOUS_STRAIN_XZ_9, stateVariables[4 + 7*9]),
				(VISCOUS_STRAIN_YZ_9, stateVariables[5 + 7*9]),
				(VOLUMETRIC_VISCOUS_STRAIN_9, stateVariables[6 + 7*9]),
				(VISCOUS_STRAIN_X_10, stateVariables[0 + 7*10]),
				(VISCOUS_STRAIN_Y_10, stateVariables[1 + 7*10]),
				(VISCOUS_STRAIN_Z_10, stateVariables[2 + 7*10]),
				(VISCOUS_STRAIN_XY_10, stateVariables[3 + 7*10]),
				(VISCOUS_STRAIN_XZ_10, stateVariables[4 + 7*10]),
				(VISCOUS_STRAIN_YZ_10, stateVariables[5 + 7*10]),
				(VOLUMETRIC_VISCOUS_STRAIN_10, stateVariables[6 + 7*10]),
				(VISCOUS_STRAIN_X_11, stateVariables[0 + 7*11]),
				(VISCOUS_STRAIN_Y_11, stateVariables[1 + 7*11]),
				(VISCOUS_STRAIN_Z_11, stateVariables[2 + 7*11]),
				(VISCOUS_STRAIN_XY_11, stateVariables[3 + 7*11]),
				(VISCOUS_STRAIN_XZ_11, stateVariables[4 + 7*11]),
				(VISCOUS_STRAIN_YZ_11, stateVariables[5 + 7*11]),
				(VOLUMETRIC_VISCOUS_STRAIN_11, stateVariables[6 + 7*11]),
				(VISCOUS_STRAIN_X_12, stateVariables[0 + 7*12]),
				(VISCOUS_STRAIN_Y_12, stateVariables[1 + 7*12]),
				(VISCOUS_STRAIN_Z_12, stateVariables[2 + 7*12]),
				(VISCOUS_STRAIN_XY_12, stateVariables[3 + 7*12]),
				(VISCOUS_STRAIN_XZ_12, stateVariables[4 + 7*12]),
				(VISCOUS_STRAIN_YZ_12, stateVariables[5 + 7*12]),
				(VOLUMETRIC_VISCOUS_STRAIN_12, stateVariables[6 + 7*12]),
			});
			return currentState;
		}

		IHaveState ICreateState.CreateState() => CreateState();
		public GenericConstitutiveLawState CurrentState
		{
			get => currentState;
			set
			{
				currentState = value;
				stresses[0] = currentState.StateValues[STRESS_X];
				stresses[1] = currentState.StateValues[STRESS_Y];
				stresses[2] = currentState.StateValues[STRESS_Z];
				stresses[3] = currentState.StateValues[STRESS_XY];
				stresses[4] = currentState.StateValues[STRESS_XZ];
				stresses[5] = currentState.StateValues[STRESS_YZ];
				strains[0] = currentState.StateValues[STRAIN_X];
				strains[1] = currentState.StateValues[STRAIN_Y];
				strains[2] = currentState.StateValues[STRAIN_Z];
				strains[3] = currentState.StateValues[STRAIN_XY];
				strains[4] = currentState.StateValues[STRAIN_XZ];
				strains[5] = currentState.StateValues[STRAIN_YZ];
				stateVariables[0 + 7 * 0] = currentState.StateValues[VISCOUS_STRAIN_X_0];
				stateVariables[1 + 7 * 0] = currentState.StateValues[VISCOUS_STRAIN_Y_0];
				stateVariables[2 + 7 * 0] = currentState.StateValues[VISCOUS_STRAIN_Z_0];
				stateVariables[3 + 7 * 0] = currentState.StateValues[VISCOUS_STRAIN_XY_0];
				stateVariables[4 + 7 * 0] = currentState.StateValues[VISCOUS_STRAIN_XZ_0];
				stateVariables[5 + 7 * 0] = currentState.StateValues[VISCOUS_STRAIN_YZ_0];
				stateVariables[6 + 7 * 0] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_0];
				stateVariables[0 + 7 * 1] = currentState.StateValues[VISCOUS_STRAIN_X_1];
				stateVariables[1 + 7 * 1] = currentState.StateValues[VISCOUS_STRAIN_Y_1];
				stateVariables[2 + 7 * 1] = currentState.StateValues[VISCOUS_STRAIN_Z_1];
				stateVariables[3 + 7 * 1] = currentState.StateValues[VISCOUS_STRAIN_XY_1];
				stateVariables[4 + 7 * 1] = currentState.StateValues[VISCOUS_STRAIN_XZ_1];
				stateVariables[5 + 7 * 1] = currentState.StateValues[VISCOUS_STRAIN_YZ_1];
				stateVariables[6 + 7 * 1] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_1];
				stateVariables[0 + 7 * 2] = currentState.StateValues[VISCOUS_STRAIN_X_2];
				stateVariables[1 + 7 * 2] = currentState.StateValues[VISCOUS_STRAIN_Y_2];
				stateVariables[2 + 7 * 2] = currentState.StateValues[VISCOUS_STRAIN_Z_2];
				stateVariables[3 + 7 * 2] = currentState.StateValues[VISCOUS_STRAIN_XY_2];
				stateVariables[4 + 7 * 2] = currentState.StateValues[VISCOUS_STRAIN_XZ_2];
				stateVariables[5 + 7 * 2] = currentState.StateValues[VISCOUS_STRAIN_YZ_2];
				stateVariables[6 + 7 * 2] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_2];
				stateVariables[0 + 7 * 3] = currentState.StateValues[VISCOUS_STRAIN_X_3];
				stateVariables[1 + 7 * 3] = currentState.StateValues[VISCOUS_STRAIN_Y_3];
				stateVariables[2 + 7 * 3] = currentState.StateValues[VISCOUS_STRAIN_Z_3];
				stateVariables[3 + 7 * 3] = currentState.StateValues[VISCOUS_STRAIN_XY_3];
				stateVariables[4 + 7 * 3] = currentState.StateValues[VISCOUS_STRAIN_XZ_3];
				stateVariables[5 + 7 * 3] = currentState.StateValues[VISCOUS_STRAIN_YZ_3];
				stateVariables[6 + 7 * 3] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_3];
				stateVariables[0 + 7 * 4] = currentState.StateValues[VISCOUS_STRAIN_X_4];
				stateVariables[1 + 7 * 4] = currentState.StateValues[VISCOUS_STRAIN_Y_4];
				stateVariables[2 + 7 * 4] = currentState.StateValues[VISCOUS_STRAIN_Z_4];
				stateVariables[3 + 7 * 4] = currentState.StateValues[VISCOUS_STRAIN_XY_4];
				stateVariables[4 + 7 * 4] = currentState.StateValues[VISCOUS_STRAIN_XZ_4];
				stateVariables[5 + 7 * 4] = currentState.StateValues[VISCOUS_STRAIN_YZ_4];
				stateVariables[6 + 7 * 4] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_4];
				stateVariables[0 + 7 * 5] = currentState.StateValues[VISCOUS_STRAIN_X_5];
				stateVariables[1 + 7 * 5] = currentState.StateValues[VISCOUS_STRAIN_Y_5];
				stateVariables[2 + 7 * 5] = currentState.StateValues[VISCOUS_STRAIN_Z_5];
				stateVariables[3 + 7 * 5] = currentState.StateValues[VISCOUS_STRAIN_XY_5];
				stateVariables[4 + 7 * 5] = currentState.StateValues[VISCOUS_STRAIN_XZ_5];
				stateVariables[5 + 7 * 5] = currentState.StateValues[VISCOUS_STRAIN_YZ_5];
				stateVariables[6 + 7 * 5] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_5];
				stateVariables[0 + 7 * 6] = currentState.StateValues[VISCOUS_STRAIN_X_6];
				stateVariables[1 + 7 * 6] = currentState.StateValues[VISCOUS_STRAIN_Y_6];
				stateVariables[2 + 7 * 6] = currentState.StateValues[VISCOUS_STRAIN_Z_6];
				stateVariables[3 + 7 * 6] = currentState.StateValues[VISCOUS_STRAIN_XY_6];
				stateVariables[4 + 7 * 6] = currentState.StateValues[VISCOUS_STRAIN_XZ_6];
				stateVariables[5 + 7 * 6] = currentState.StateValues[VISCOUS_STRAIN_YZ_6];
				stateVariables[6 + 7 * 6] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_6];
				stateVariables[0 + 7 * 7] = currentState.StateValues[VISCOUS_STRAIN_X_7];
				stateVariables[1 + 7 * 7] = currentState.StateValues[VISCOUS_STRAIN_Y_7];
				stateVariables[2 + 7 * 7] = currentState.StateValues[VISCOUS_STRAIN_Z_7];
				stateVariables[3 + 7 * 7] = currentState.StateValues[VISCOUS_STRAIN_XY_7];
				stateVariables[4 + 7 * 7] = currentState.StateValues[VISCOUS_STRAIN_XZ_7];
				stateVariables[5 + 7 * 7] = currentState.StateValues[VISCOUS_STRAIN_YZ_7];
				stateVariables[6 + 7 * 7] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_7];
				stateVariables[0 + 7 * 8] = currentState.StateValues[VISCOUS_STRAIN_X_8];
				stateVariables[1 + 7 * 8] = currentState.StateValues[VISCOUS_STRAIN_Y_8];
				stateVariables[2 + 7 * 8] = currentState.StateValues[VISCOUS_STRAIN_Z_8];
				stateVariables[3 + 7 * 8] = currentState.StateValues[VISCOUS_STRAIN_XY_8];
				stateVariables[4 + 7 * 8] = currentState.StateValues[VISCOUS_STRAIN_XZ_8];
				stateVariables[5 + 7 * 8] = currentState.StateValues[VISCOUS_STRAIN_YZ_8];
				stateVariables[6 + 7 * 8] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_8];
				stateVariables[0 + 7 * 9] = currentState.StateValues[VISCOUS_STRAIN_X_9];
				stateVariables[1 + 7 * 9] = currentState.StateValues[VISCOUS_STRAIN_Y_9];
				stateVariables[2 + 7 * 9] = currentState.StateValues[VISCOUS_STRAIN_Z_9];
				stateVariables[3 + 7 * 9] = currentState.StateValues[VISCOUS_STRAIN_XY_9];
				stateVariables[4 + 7 * 9] = currentState.StateValues[VISCOUS_STRAIN_XZ_9];
				stateVariables[5 + 7 * 9] = currentState.StateValues[VISCOUS_STRAIN_YZ_9];
				stateVariables[6 + 7 * 9] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_9];
				stateVariables[0 + 7 * 10] = currentState.StateValues[VISCOUS_STRAIN_X_10];
				stateVariables[1 + 7 * 10] = currentState.StateValues[VISCOUS_STRAIN_Y_10];
				stateVariables[2 + 7 * 10] = currentState.StateValues[VISCOUS_STRAIN_Z_10];
				stateVariables[3 + 7 * 10] = currentState.StateValues[VISCOUS_STRAIN_XY_10];
				stateVariables[4 + 7 * 10] = currentState.StateValues[VISCOUS_STRAIN_XZ_10];
				stateVariables[5 + 7 * 10] = currentState.StateValues[VISCOUS_STRAIN_YZ_10];
				stateVariables[6 + 7 * 10] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_10];
				stateVariables[0 + 7 * 11] = currentState.StateValues[VISCOUS_STRAIN_X_11];
				stateVariables[1 + 7 * 11] = currentState.StateValues[VISCOUS_STRAIN_Y_11];
				stateVariables[2 + 7 * 11] = currentState.StateValues[VISCOUS_STRAIN_Z_11];
				stateVariables[3 + 7 * 11] = currentState.StateValues[VISCOUS_STRAIN_XY_11];
				stateVariables[4 + 7 * 11] = currentState.StateValues[VISCOUS_STRAIN_XZ_11];
				stateVariables[5 + 7 * 11] = currentState.StateValues[VISCOUS_STRAIN_YZ_11];
				stateVariables[6 + 7 * 11] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_11];
				stateVariables[0 + 7 * 12] = currentState.StateValues[VISCOUS_STRAIN_X_12];
				stateVariables[1 + 7 * 12] = currentState.StateValues[VISCOUS_STRAIN_Y_12];
				stateVariables[2 + 7 * 12] = currentState.StateValues[VISCOUS_STRAIN_Z_12];
				stateVariables[3 + 7 * 12] = currentState.StateValues[VISCOUS_STRAIN_XY_12];
				stateVariables[4 + 7 * 12] = currentState.StateValues[VISCOUS_STRAIN_XZ_12];
				stateVariables[5 + 7 * 12] = currentState.StateValues[VISCOUS_STRAIN_YZ_12];
				stateVariables[6 + 7 * 12] = currentState.StateValues[VOLUMETRIC_VISCOUS_STRAIN_12];
			}
		}

		/// <summary>
		///   Updates the element's material with the provided incremental strains.
		/// </summary>
		/// <param name = "strainsIncrement">The incremental strains to use for the next step.</param>
		public double[] UpdateConstitutiveMatrixAndEvaluateResponse(double[] strainsIncrement)
		{
			incrementalStrains.CopyFrom(strainsIncrement);
			this.CalculateNextStressStrainPoint();

			return stressesNew;
		}

		/// <summary>
		///   Builds the consistent tangential constitutive matrix.
		/// </summary>
		/// <param name = "value1"> This is a constant already calculated in the calling method. </param>
		/// <remarks>
		///   We need an additional constant here equal to: 2*G*G*deltaPlasticStrain*sqrt(3/J2elastic)
		///   Since value1 is already calculated, the additional constant can be calculated through it:
		///   value3 = 2 * G * value1;
		/// </remarks>
		private void BuildConsistentTangentialConstitutiveMatrix()
		{
			this.constitutiveMatrix = Matrix.CreateZero(TotalStresses, TotalStrains);

			var shearModulus = 0.0;
			var bulkModulus = 0.0;
			for (int iNME = 0; iNME < numberOfMaxwelElements; iNME++)
			{
				var dT = dtime / taui[iNME];
				if (dT > 1E-7)
				{
					shearModulus += gnormi[iNME] * (1.0 / dT) * (dT + Math.Exp(-dT) - 1.0);
					bulkModulus += knormi[iNME] * (1.0 / dT) * (dT + Math.Exp(-dT) - 1.0);
				}
				else
				{
					shearModulus += (1.0 / 2.0) * gnormi[iNME] * dT;
					bulkModulus += (1.0 / 2.0) * knormi[iNME] * dT;
				}

			}

			shearModulus = instantaneousShearModulus * (1.0 - shearModulus);
			bulkModulus = instantaneousBulkModulus * (1.0 - bulkModulus);

			for (int k1 = 0; k1 < TotalStresses; k1++)
			{
				for (int k2 = 0; k2 < TotalStresses; k2++)
				{
					this.constitutiveMatrix[k2, k1] = (2d * shearModulus * SupportiveMatrixForConsistentConstitutiveMatrixI4M[k2, k1]) + ((bulkModulus - ((2d / 3d) * shearModulus)) * SupportiveMatrixForConsistentConstitutiveMatrixIIM[k2, k1]);
				}
			}
		}

		/// <summary>
		///   Calculates the next stress-strain point.
		/// </summary>
		/// <exception cref = "InvalidOperationException"> When the new plastic strain is less than the previous one.</exception>
		private void CalculateNextStressStrainPoint()
		{

			// EVALUATE DEVIATORIC AND VOLUMETRIC STRAIN AT START OF INCREMENT t_n
			// and ASSIGN STORED DATA FROM PREVIOUS TIME
			var f_n = 0d;
			var Df = 0d;
			for (int i = 0; i < 3; i++)
			{
				f_n += strains[i];
				Df += incrementalStrains[i];
			}

			var e_n = new double[6];
			var De = new double[6];
			for (int i = 0; i < 6; i++)
			{
				e_n[i] = strains[i] - (1d / 3d) * f_n * SupportiveArrayForDeviatoricFormulation[i];
				De[i] = incrementalStrains[i] - (1d / 3d) * Df * SupportiveArrayForDeviatoricFormulation[i];
			}

			var ei_n = new double[6, numberOfMaxwelElements];
			var fi_n = new double[numberOfMaxwelElements];
			for (int iNME = 0; iNME < numberOfMaxwelElements; iNME++)
			{
				for (int i = 0; i < 6; i++)
				{
					ei_n[i, iNME] = stateVariables[i + (7 * iNME)];
				}
				fi_n[iNME] = stateVariables[6 + (7 * iNME)];
			}


			// EVALUATE VISCOUS (CREEP) STRAIN IN EACH TERM OF THE SERIES
			var Dei = new double[6, numberOfMaxwelElements];
			var Dfi = new double[numberOfMaxwelElements];
			for (int iNME = 0; iNME < numberOfMaxwelElements; iNME++)
			{
				var dT = dtime / taui[iNME];
				for (int i = 0; i < 6; i++)
				{
					if (dT > 1E-7)
					{
						Dei[i, iNME] = ((1d / dT) * (dT + Math.Exp(-dT) - 1d) * De[i]) + ((1d - Math.Exp(-dT)) * (e_n[i] - ei_n[i, iNME]));
					}
					else
					{
						Dei[i, iNME] = dT * ((0.5 * De[i]) + e_n[i] - ei_n[i, iNME]);
					}
				}
				if (dT > 1E-7)
				{
					Dfi[iNME] = ((1.0 / dT) * (dT + Math.Exp(-dT) - 1.0) * Df) + ((1.0 - Math.Exp(-dT)) * (f_n - fi_n[iNME]));
				}
				else
				{
					Dfi[iNME] = dT * ((0.5 * Df) + f_n - fi_n[iNME]);
				}
			}

			// Calculate viscous strains
			var De_viscous = new double[6];
			var Df_viscous = 0d;
			var ei_n1 = new double[6, numberOfMaxwelElements];
			var fi_n1 = new double[numberOfMaxwelElements];
			for (int iNME = 0; iNME < numberOfMaxwelElements; iNME++)
			{
				for (int i = 0; i < 6; i++)
				{
					De_viscous[i] += gnormi[iNME] * Dei[i, iNME];
					ei_n1[i, iNME] = ei_n[i, iNME] + Dei[i, iNME];
				}
				Df_viscous += knormi[iNME] * Dfi[iNME];
				fi_n1[iNME] = fi_n[iNME] + Dfi[iNME];
			}

			for (int i = 0; i < 6; i++)
			{
				this.stressesNew[i] = this.stresses[i] + (2.0 * instantaneousShearModulus * (De[i] - De_viscous[i]) * SupportiveArrayForEngineeringStrains[i]) + (instantaneousBulkModulus * Df * SupportiveArrayForDeviatoricFormulation[i]);   // PROSOXI(+BULK0 * Df)
			}

			this.BuildConsistentTangentialConstitutiveMatrix();

			// Update strains
			for (int i = 0; i < 6; i++)
			{
				strainsNew[i] = strains[i] + incrementalStrains[i];
			}
			// Update state variables
			for (int iNME = 0; iNME < numberOfMaxwelElements; iNME++)
			{
				for (int i = 0; i < 6; i++)
				{
					stateVariablesNew[(iNME * 7) + i] = ei_n1[i, iNME];
				}
				stateVariablesNew[(iNME * 7) + 6] = fi_n1[iNME];
			}

			this.modified = true;

		}

		/// <summary>
		///   Calculates and returns the first stress invariant (I1).
		/// </summary>
		/// <returns> The first stress invariant (I1).</returns>
		public double GetFirstStressInvariant(double[] stresses) => stresses[0] + stresses[1] + stresses[2];

		/// <summary>
		///   Calculates and returns the mean hydrostatic stress.
		/// </summary>
		/// <returns> The mean hydrostatic stress.</returns>
		public double GetMeanStress(double[] stresses) => GetFirstStressInvariant(stresses) / 3.0;

		/// <summary>
		///   Calculates and returns the second stress invariant (I2).
		/// </summary>
		/// <returns> The second stress invariant (I2).</returns>
		public double GetSecondStressInvariant(double[] stresses)
			=> (stresses[0] * stresses[1]) + (stresses[1] * stresses[2]) + (stresses[0] * stresses[2])
			- Math.Pow(stresses[5], 2) - Math.Pow(stresses[3], 2) - Math.Pow(stresses[4], 2);

		/// <summary>
		///   Calculates and returns the stress deviator tensor in vector form.
		/// </summary>
		/// <returns> The stress deviator tensor in vector form.</returns>
		public double[] GetStressDeviator(double[] stresses)
		{
			var hydrostaticStress = this.GetMeanStress(stresses);
			var stressDeviator = new double[]
			{
				stresses[0] - hydrostaticStress,
				stresses[1] - hydrostaticStress,
				stresses[2] - hydrostaticStress,
				stresses[3],
				stresses[4],
				stresses[5]
			};

			return stressDeviator;
		}

		/// <summary>
		///   Calculates and returns the third stress invariant (I3).
		/// </summary>
		/// <returns> The third stress invariant (I3). </returns>
		public double GetThirdStressInvariant(double[] stresses)
			=> (stresses[0] * stresses[1] * stresses[2]) + (2 * stresses[5] * stresses[3] * stresses[4])
			- (Math.Pow(stresses[5], 2) * stresses[2]) - (Math.Pow(stresses[3], 2) * stresses[0])
			- (Math.Pow(stresses[4], 2) * stresses[1]);

		/// <summary>
		///   Returns the first stress invariant of the stress deviator tensor (J1), which is zero.
		/// </summary>
		/// <returns> The first stress invariant of the stress deviator tensor (J1). </returns>
		public double GetDeviatorFirstStressInvariant(double[] stresses) => 0;

		/// <summary>
		///   Calculates and returns the second stress invariant of the stress deviator tensor (J2).
		/// </summary>
		/// <returns> The second stress invariant of the stress deviator tensor (J2). </returns>
		public double GetDeviatorSecondStressInvariant(double[] stresses)
		{
			double i1 = this.GetFirstStressInvariant(stresses);
			double i2 = this.GetSecondStressInvariant(stresses);

			double j2 = (1 / 3d * Math.Pow(i1, 2)) - i2;
			return j2;
		}

		/// <summary>
		///   Calculates and returns the third stress invariant of the stress deviator tensor (J3).
		/// </summary>
		/// <returns> The third deviator stress invariant (J3). </returns>
		public double GetDeviatorThirdStressInvariant(double[] stresses)
		{
			double i1 = this.GetFirstStressInvariant(stresses);
			double i2 = this.GetSecondStressInvariant(stresses);
			double i3 = this.GetThirdStressInvariant(stresses);

			double j3 = (2 / 27 * Math.Pow(i1, 3)) - (1 / 3 * i1 * i2) + i3;
			return j3;
		}

		public void SetCurrentTime(double time)
		{
			timePrevious = timeCurrent;
			timeCurrent = time;
		}
	}
}
