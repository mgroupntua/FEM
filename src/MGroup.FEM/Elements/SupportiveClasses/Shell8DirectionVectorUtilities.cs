using System;

namespace MGroup.FEM.Elements.SupportiveClasses
{
	public class Shell8DirectionVectorUtilities
	{
		public static (double[][] tU, double[][] tUvec) GetInitialDirectionVectorValues(double[][] oVn_i)
		{
			double[][] tU = new double[8][];
			double[][] tUvec = new double[8][];
			for (int j = 0; j < 8; j++)
			{
				tU[j] = new double[6];
				tUvec[j] = new double[6];
				for (int k = 0; k < 3; k++) { tU[j][3 + k] = oVn_i[j][k]; }


				tUvec[j][0] = tU[j][5];
				tUvec[j][1] = 0;
				tUvec[j][2] = -tU[j][3];

				double tV1norm = Math.Sqrt(tUvec[j][0] * tUvec[j][0] + tUvec[j][1] * tUvec[j][1] + tUvec[j][2] * tUvec[j][2]);

				tUvec[j][0] = tUvec[j][0] / tV1norm;
				tUvec[j][1] = tUvec[j][1] / tV1norm;
				tUvec[j][2] = tUvec[j][2] / tV1norm;

				tUvec[j][3] = tU[j][3 + 1] * tUvec[j][2] - tU[j][3 + 2] * tUvec[j][1];
				tUvec[j][4] = tU[j][3 + 2] * tUvec[j][0] - tU[j][3 + 0] * tUvec[j][2];
				tUvec[j][5] = tU[j][3 + 0] * tUvec[j][1] - tU[j][3 + 1] * tUvec[j][0];
			}

			return (tU, tUvec);
		}


		/// <summary>
		/// Used by cohesive shell elements for now
		/// </summary>
		/// <param name="ak"></param>
		/// <param name="bk"></param>
		/// <param name="n_vector"></param>
		/// <param name="tU">this will be updated</param>
		/// <param name="tUvec">this will be updated</param>
		/// 
		public static void RotateNodalDirectionVectors(double ak, double bk, int n_vector, double[][] tU, double[][] tUvec)
		{
			double[,] Q = new double[3, 3];
			double[,] Q2 = new double[3, 3];

			double[] tdtVn = new double[3];
			double[] tdtV1 = new double[3];
			double[] tdtV2 = new double[3];
			double[] theta_vec = new double[3];
			double[,] s_k = new double[3, 3];

			for (int j = 0; j < 3; j++)
			{
				theta_vec[j] = ak * tUvec[n_vector][j] + bk * tUvec[n_vector][3 + j];
			}
			double theta = Math.Sqrt((theta_vec[0] * theta_vec[0]) + (theta_vec[1] * theta_vec[1]) + (theta_vec[2] * theta_vec[2]));
			if (theta > 0)
			{
				s_k[0, 1] = -theta_vec[2];
				s_k[0, 2] = theta_vec[1];
				s_k[1, 0] = theta_vec[2];
				s_k[1, 2] = -theta_vec[0];
				s_k[2, 0] = -theta_vec[1];
				s_k[2, 1] = theta_vec[0];

				for (int j = 0; j < 3; j++)
				{
					for (int m = 0; m < 3; m++)
					{
						Q[j, m] = (Math.Sin(theta) / theta) * s_k[j, m];
					}
				}

				for (int m = 0; m < 3; m++)
				{
					Q[m, m] += 1;
				}
				double gk1 = 0.5 * ((Math.Sin(0.5 * theta) / (0.5 * theta)) * (Math.Sin(0.5 * theta) / (0.5 * theta)));
				for (int j = 0; j < 3; j++)
				{
					for (int m = 0; m < 3; m++)
					{
						for (int n = 0; n < 3; n++)
						{ Q2[j, m] += gk1 * s_k[j, n] * s_k[n, m]; }
					}
				}
				for (int j = 0; j < 3; j++)
				{
					for (int m = 0; m < 3; m++)
					{
						Q[j, m] += Q2[j, m];
					}
				}
				//
				for (int j = 0; j < 3; j++)
				{
					tdtVn[j] = 0;
					for (int m = 0; m < 3; m++)
					{
						tdtVn[j] += Q[j, m] * tU[n_vector][3 + m];
					}
				}

				for (int j = 0; j < 3; j++)
				{
					tU[n_vector][3 + j] = tdtVn[j];
				}
				//
				for (int j = 0; j < 3; j++)
				{
					tdtV1[j] = 0;
					for (int m = 0; m < 3; m++)
					{
						tdtV1[j] += Q[j, m] * tUvec[n_vector][m];
					}
				}

				for (int j = 0; j < 3; j++)
				{
					tUvec[n_vector][j] = tdtV1[j];
				}
				//
				for (int j = 0; j < 3; j++)
				{
					tdtV2[j] = 0;
					for (int m = 0; m < 3; m++)
					{
						tdtV2[j] += Q[j, m] * tUvec[n_vector][3 + m];
					}
				}

				for (int j = 0; j < 3; j++)
				{
					tUvec[n_vector][3 + j] = tdtV2[j];
				}
			}
		}
	}
}
