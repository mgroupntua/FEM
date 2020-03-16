using System;
using MGroup.MSolve.Discretization.Commons;
using Xunit;

namespace MGroup.FEM.Tests.IntegrationTests
{
	using Constitutive.Structural;
	using MSolve.Constitutive;

	public static class BondSlipTest
	{
		[Fact]
		public static void CheckBondSlipMaterial()
		{
			(double[][] stressHistory, double[][,] constitutiveMatrixHistory) = CheckStressStrainBonSlipMaterial();
			int[] positions = new int[8] { 25, 50, 75, 100, 125, 150, 175, 200 };
			var calculated_stresses = new double[positions.GetLength(0), 2];
			var calculated_Const = new double[2 * positions.GetLength(0), 2];

			for (int i1 = 0; i1 < positions.GetLength(0); i1++)
			{
				calculated_stresses[i1, 0] = stressHistory[positions[i1] - 1][0];
				calculated_stresses[i1, 1] = stressHistory[positions[i1] - 1][1];

				calculated_Const[2 * (i1) + 0, 0] = constitutiveMatrixHistory[positions[i1] - 1][0, 0];
				calculated_Const[2 * (i1) + 0, 1] = constitutiveMatrixHistory[positions[i1] - 1][0, 1];
				calculated_Const[2 * (i1) + 1, 0] = constitutiveMatrixHistory[positions[i1] - 1][1, 0];
				calculated_Const[2 * (i1) + 1, 1] = constitutiveMatrixHistory[positions[i1] - 1][1, 1];
			}




			double[,] stress_data = new double[8, 2] {{50.518148554092285,29.166666666666675},
													  {88.045916051417962,50.833333333333265},
													  {93.097730906827223,53.749999999999886},
													  {98.149545762236471,56.666666666666501},
													  {87.035553080336157,50.249999999999829},
													  {11.258330249197762,6.499999999999830},
													  {-64.518892581940605,-37.250000000000163},
													  {-77.942286340599466,-45.000000000000057},};

			double[,] const_data = new double[16, 2] {{100.000000000000000,0.000000000000000},
													  {0.000000000000000,100.000000000000000},
													  {32.499999999999943,-38.971143170299698},
													  {-38.971143170299705,77.500000000000057},
													  {32.499999999999908,-38.971143170299676},
													  {-38.971143170299669,77.500000000000099},
													  {32.499999999999851,-38.971143170299662},
													  {-38.971143170299662,77.500000000000142},
													  {100.000000000000000,0.000000000000000},
													  {0.000000000000000,100.000000000000000},
													  {100.000000000000000,0.000000000000000},
													  {0.000000000000000,100.000000000000000},
													  {100.000000000000000,0.000000000000000},
													  {0.000000000000000,100.000000000000000},
													  {32.500000000000028,-38.971143170299747},
													  {-38.971143170299754,77.499999999999972},};

			Assert.True(AreDisplacementsSame(calculated_stresses, stress_data));
			Assert.True(AreDisplacementsSame(calculated_Const, const_data));
		}

		public static (double[][], double[][,]) CheckStressStrainBonSlipMaterial()
		{
			//VectorExtensions.AssignTotalAffinityCount();
			BondSlipCohMat material1 = new BondSlipCohMat(100, 10, 100, 100, new double[2], new double[2], 1e-10);
			int loadsteps_2 = 120;
			double[][] DeltaEhist = new double[2 * loadsteps_2][];
			double phi_metakinhshs = ((double)30 / (double)360) * 2 * Math.PI;
			double[] epsilon_max = new double[3] { 2.8 * Math.Cos(phi_metakinhshs), 2.8 * Math.Sin(phi_metakinhshs), 2.8 };
			for (int i1 = 0; i1 < loadsteps_2; i1++)
			{ DeltaEhist[i1] = new double[3] { (1 / (double)loadsteps_2) * epsilon_max[0], (1 / (double)loadsteps_2) * epsilon_max[1], (1 / (double)loadsteps_2) * epsilon_max[2] }; }
			for (int i1 = loadsteps_2; i1 < 2 * loadsteps_2; i1++)
			{ DeltaEhist[i1] = new double[3] { -1.5 * (1 / (double)loadsteps_2) * epsilon_max[0], -1.5 * (1 / (double)loadsteps_2) * epsilon_max[1], -1.5 * (1 / (double)loadsteps_2) * epsilon_max[2] }; }
			double[][] Ehist = new double[2 * loadsteps_2][];
			Ehist[0] = new double[3] { DeltaEhist[0][0], DeltaEhist[0][1], DeltaEhist[0][2] };
			for (int i1 = 1; i1 < 2 * loadsteps_2; i1++)
			{ Ehist[i1] = new double[3] { Ehist[i1 - 1][0] + DeltaEhist[i1][0], Ehist[i1 - 1][1] + DeltaEhist[i1][1], Ehist[i1 - 1][2] + DeltaEhist[i1][2] }; }


			(double[][] stressHistory, double[][,] constitutiveMatrixHistory) = StressStrainHistory(Ehist, material1);

			return (stressHistory, constitutiveMatrixHistory);

		}

		public static bool AreDisplacementsSame(double[,] expectedValues,
			double[,] computedValues)
		{
			var comparer = new ValueComparer(1E-8);
			for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
			{
				for (int i2 = 0; i2 < expectedValues.GetLength(1); i2++)
				{
					if (!comparer.AreEqual(expectedValues[i1, i2], computedValues[i1, i2]))
					{
						return false;
					}
				}
			}
			return true;
		}

		private static (double[][] stressHistory, double[][,] constitutiveMatrixHistory) StressStrainHistory(double[][] strainHistory, ICohesiveZoneMaterial3D testedMaterial)
		{
			double[][] stressHistory = new double[strainHistory.GetLength(0)][];
			double[][,] constitutiveMatrixHistory = new double[strainHistory.GetLength(0)][,];

			for (int l = 0; l < strainHistory.GetLength(0); l++)
			{
				if (l == 42)
				{
					Console.Write("breakPointIsHere");
				}
				testedMaterial.UpdateMaterial(strainHistory[l]);
				testedMaterial.SaveState();
				stressHistory[l] = new double[testedMaterial.Tractions.Length];
				testedMaterial.Tractions.CopyTo(stressHistory[l], 0);
				constitutiveMatrixHistory[l] = new double[testedMaterial.ConstitutiveMatrix.NumColumns, testedMaterial.ConstitutiveMatrix.NumRows];

				for (int m = 0; m < testedMaterial.ConstitutiveMatrix.NumColumns; m++)
				{
					for (int n = 0; n < testedMaterial.ConstitutiveMatrix.NumRows; n++)
					{ constitutiveMatrixHistory[l][m, n] = testedMaterial.ConstitutiveMatrix[m, n]; }
				}
			}

			return (stressHistory, constitutiveMatrixHistory);
		}
	}
}
