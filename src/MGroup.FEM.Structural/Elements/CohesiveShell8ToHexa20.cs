using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.FEM.Elements.SupportiveClasses;
using MGroup.FEM.Embedding;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;
using MGroup.FEM.Interpolation;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;

//TODO: current nodal coordinates should be managed by the analyzer, instead of each element calculting and storing them independently. The same applies for direction vectors of shells. 
//TODO: direction vectors creation and update could be handled by a dedicated class that will be composed into this element. Which element would update them then?
//TODO: perhaps separate the cases fot when the shell is over or under the cohesive element ibto 2 subclasses
namespace MGroup.FEM.Structural.Elements
{
	/// <summary>
	/// Cohesive element for modeling of delamination and debonding effects between parts modeled with shell and hexa20 elements.
	/// Authors: Gerasimos Sotiropoulos
	/// </summary>
	public class CohesiveShell8ToHexa20 : IStructuralFiniteElement, IEmbeddedElement
	{
		protected readonly static IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		protected readonly static IDofType[] nodalDOFTypes2 = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };
		protected readonly static IDofType[][] dofTypes = new IDofType[][] { nodalDOFTypes2, nodalDOFTypes2, nodalDOFTypes2,
			nodalDOFTypes2, nodalDOFTypes2, nodalDOFTypes2, nodalDOFTypes2, nodalDOFTypes2,nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
			nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
		protected readonly ICohesiveZoneMaterial3D[] materialsAtGaussPoints;
		private readonly InterpolationShell8Cohesive interpolation = InterpolationShell8Cohesive.UniqueInstance;
		private readonly List<EmbeddedNode> embeddedNodes = new List<EmbeddedNode>();
		protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		private int nGaussPoints;
		public double[][] oVn_i { get; set; }
		public double[] tk { get; set; }

		/// <summary>
		/// Zero denotes bottom surface of the shell element
		/// </summary>
		public int ShellElementSide { get; set; }

		private double[][] tU;   // dimensions: 8x6
		private double[][] tUvec;// dimensions: 8x6

		/// <summary>
		/// Initial coordinates of shell midsurface nodes
		/// </summary>
		private double[][] ox_i_shell_midsurface;

		/// <summary>
		/// Initial nodel coordinates of 16 node inner cohesive element
		/// </summary>
		private double[][] ox_i;

		/// <summary>
		/// Unrolled current nodal coordinates of 16 node inner cohesive element
		/// </summary>
		private double[] x_local;

		/// <summary>
		/// Auxiliary variables for calculating rotations and updating the element's direction vectors.
		/// </summary>      
		private double[] ak_total = new double[8];
		private double[] bk_total = new double[8];
		public bool MatrixIsNotInitialized = true;

		public CohesiveShell8ToHexa20(ICohesiveZoneMaterial3D material, IQuadrature2D quadratureForStiffness)
		{
			this.QuadratureForStiffness = quadratureForStiffness;
			this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;
			materialsAtGaussPoints = new ICohesiveZoneMaterial3D[nGaussPoints];
			for (int i = 0; i < nGaussPoints; i++) materialsAtGaussPoints[i] = material.Clone();
		}

		public int ID => 13;
		public CellType CellType { get; } = CellType.Unknown;

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		public IList<EmbeddedNode> EmbeddedNodes => embeddedNodes;

		public IQuadrature2D QuadratureForStiffness { get; }

		private void GetInitialGeometricDataForMidsurface(IElement element)
		{
			ox_i_shell_midsurface = new double[8][];
			for (int j = 0; j < 8; j++)
			{
				ox_i_shell_midsurface[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
			}
		}

		private void UpdateCoordinateDataForDirectVectorsAndMidsurface(double[] localdisplacements, out double[][] tx_i_shell_midsurface)
		{
			double ak;
			double bk;
			tx_i_shell_midsurface = new double[8][];
			for (int j = 0; j < 8; j++)
			{
				tx_i_shell_midsurface[j] = new double[3];
			}

			for (int k = 0; k < 8; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					tx_i_shell_midsurface[k][l] = ox_i_shell_midsurface[k][l] + localdisplacements[5 * k + l];
				}

				//update tU & tUvec 
				tU[k][0] = localdisplacements[5 * k + 0];
				tU[k][1] = localdisplacements[5 * k + 1];
				tU[k][2] = localdisplacements[5 * k + 2];
				ak = localdisplacements[5 * k + 3] - ak_total[k];
				ak_total[k] = localdisplacements[5 * k + 3];
				bk = localdisplacements[5 * k + 4] - bk_total[k];
				bk_total[k] = localdisplacements[5 * k + 4];
				Shell8DirectionVectorUtilities.RotateNodalDirectionVectors(ak, bk, k, tU, tUvec);
			}
		}

		private void GetInitialGeometricDataAndInitializeMatrices(IElement element)
		{
			(tU, tUvec) = Shell8DirectionVectorUtilities.GetInitialDirectionVectorValues(oVn_i);
			this.GetInitialGeometricDataForMidsurface(element);
			ox_i = new double[16][];
			if (ShellElementSide == 0)
			{
				for (int j = 0; j < 8; j++)
				{
					ox_i[j] = new double[] { ox_i_shell_midsurface[j][0]-0.5* tk[j] * tU[j][3], ox_i_shell_midsurface[j][1] - 0.5 * tk[j] * tU[j][4],
											 ox_i_shell_midsurface[j][2]-0.5* tk[j] * tU[j][5], };
				}
				for (int j = 8; j < 16; j++)
				{ ox_i[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, }; }
			}
			else
			{
				for (int j = 0; j < 8; j++)
				{ ox_i[j] = new double[] { element.Nodes[j + 8].X, element.Nodes[j + 8].Y, element.Nodes[j + 8].Z, }; }
				for (int j = 8; j < 16; j++)
				{
					ox_i[j] = new double[] { ox_i_shell_midsurface[j-8][0]+0.5* tk[j-8] * tU[j-8][3], ox_i_shell_midsurface[j-8][1] + 0.5 * tk[j-8] * tU[j-8][4],
											 ox_i_shell_midsurface[j-8][2]+0.5* tk[j-8] * tU[j-8][5], };
				}
			}

			x_local = new double[48];
		}

		private double[][] UpdateCoordinateDataAndCalculateDisplacementVector(double[] localdisplacements)
		{
			IReadOnlyList<Matrix> shapeFunctionDerivatives = interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			IReadOnlyList<double[]> N1 = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
			IReadOnlyList<Matrix> N3 = interpolation.EvaluateN3ShapeFunctionsReorganized(QuadratureForStiffness);

			double[,] u_prok = new double[3, 8];
			double[,] x_bar = new double[3, 8];

			double[] e_ksi = new double[3];
			double e_ksi_norm;
			double[] e_heta = new double[3];
			double[] e_1 = new double[3];
			double[] e_2 = new double[3];
			double[] e_3 = new double[3];
			double e_3_norm;
			double[] u = new double[3];

			double[] coh_det_J_t = new double[nGaussPoints];

			double[][] Delta = new double[nGaussPoints][];
			double[][] c_1 = new double[nGaussPoints][];
			for (int j = 0; j < nGaussPoints; j++)
			{
				Delta[j] = new double[3];
				c_1[j] = new double[3];
			}

			double[][,] R = new double[nGaussPoints][,]; //TODO: maybe cache R
			for (int j = 0; j < nGaussPoints; j++)
			{
				R[j] = new double[3, 3];
			}

			this.UpdateCoordinateDataForDirectVectorsAndMidsurface(localdisplacements, out double[][] tx_i_shell_midsurface);

			// Update x_local
			if (ShellElementSide == 0)
			{
				for (int j = 0; j < 8; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						x_local[3 * j + k] = tx_i_shell_midsurface[j][k] - 0.5 * tk[j] * tU[j][3 + k];
					}
				}
				for (int j = 8; j < 16; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						x_local[3 * j + k] = ox_i[j][k] + localdisplacements[40 + 3 * (j - 8) + k];
					}
				}
			}
			else
			{
				for (int j = 0; j < 8; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						x_local[3 * j + k] = ox_i[j][k] + localdisplacements[40 + 3 * j + k];
					}
				}
				for (int j = 8; j < 16; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						x_local[3 * j + k] = tx_i_shell_midsurface[j - 8][k] + 0.5 * tk[j - 8] * tU[j - 8][3 + k];
					}
				}
			}

			for (int j = 0; j < 8; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					u_prok[k, j] = x_local[k + 3 * j] - x_local[24 + k + 3 * j];
					x_bar[k, j] = x_local[k + 3 * j] + x_local[24 + k + 3 * j];
				}
			}

			//Calculate Delta for all GPs
			for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
			{
				for (int l = 0; l < 3; l++)
				{
					e_ksi[l] = 0;
					e_heta[l] = 0;
					for (int m = 0; m < 8; m++) // must be 4 in cohesive 8-node
					{
						e_ksi[l] += shapeFunctionDerivatives[npoint1][0, m] * x_bar[l, m];
						e_heta[l] += shapeFunctionDerivatives[npoint1][1, m] * x_bar[l, m];
					}
					e_ksi[l] = 0.5 * e_ksi[l];
					e_heta[l] = 0.5 * e_heta[l];
				}
				this.Cross(e_ksi, e_heta, e_3);
				e_3_norm = Math.Sqrt(e_3[0] * e_3[0] + e_3[1] * e_3[1] + e_3[2] * e_3[2]);
				e_ksi_norm = Math.Sqrt(e_ksi[0] * e_ksi[0] + e_ksi[1] * e_ksi[1] + e_ksi[2] * e_ksi[2]);
				for (int l = 0; l < 3; l++)
				{
					e_3[l] = e_3[l] / e_3_norm;
					e_1[l] = e_ksi[l] / e_ksi_norm;
				}
				this.Cross(e_1, e_3, e_2);
				for (int l = 0; l < 3; l++)
				{
					R[npoint1][l, 0] = e_1[l];
					R[npoint1][l, 1] = e_2[l];
					R[npoint1][l, 2] = e_3[l];
				}
				for (int l = 0; l < 3; l++)
				{ u[l] = 0; }
				for (int l = 0; l < 3; l++)
				{
					for (int m = 0; m < 8; m++)  // must be changed for cohesive 8-nodes
					{
						u[l] += u_prok[l, m] * N1[npoint1][m];
					}
				}
				for (int l = 0; l < 3; l++)

				{
					for (int m = 0; m < 3; m++)
					{
						Delta[npoint1][l] += R[npoint1][m, l] * u[m];
					}
				}
			}
			return Delta;
		}

		private Tuple<Matrix[], double[]> CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations()
		{
			IReadOnlyList<double[]> N1 = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
			IReadOnlyList<Matrix> N3 = interpolation.EvaluateN3ShapeFunctionsReorganized(QuadratureForStiffness);
			IReadOnlyList<Matrix> shapeFunctionDerivatives = interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			double[] integrationsCoeffs = new double[nGaussPoints];
			Matrix[] RtN3 = new Matrix[nGaussPoints];
			double[,] x_bar = new double[3, 8];

			double[] e_1 = new double[3];
			double[] e_2 = new double[3];
			double[] e_3 = new double[3];
			double e_3_norm;

			double[] coh_det_J_t = new double[nGaussPoints];

			double[][] c_1 = new double[nGaussPoints][];
			for (int j = 0; j < nGaussPoints; j++)
			{
				c_1[j] = new double[3];
			}

			Matrix[] R = new Matrix[nGaussPoints]; //TODO: perhaps cache matrices in InitializeMatrices() where RtN3 is calculated
			for (int j = 0; j < nGaussPoints; j++)
			{ R[j] = Matrix.CreateZero(3, 3); }

			for (int j = 0; j < 8; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					x_bar[k, j] = x_local[k + 3 * j] + x_local[24 + k + 3 * j];
				}
			}

			// Calculate Delta for all GPs
			for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
			{
				double[] e_ksi = new double[3];
				double[] e_heta = new double[3];
				for (int l = 0; l < 3; l++)
				{
					for (int m = 0; m < 8; m++) // must be 4 in cohesive 8-nodes
					{
						e_ksi[l] += shapeFunctionDerivatives[npoint1][0, m] * x_bar[l, m];
						e_heta[l] += shapeFunctionDerivatives[npoint1][1, m] * x_bar[l, m];
					}
					e_ksi[l] = 0.5 * e_ksi[l];
					e_heta[l] = 0.5 * e_heta[l];
				}
				this.Cross(e_ksi, e_heta, e_3);
				e_3_norm = Math.Sqrt(e_3[0] * e_3[0] + e_3[1] * e_3[1] + e_3[2] * e_3[2]);
				double e_ksi_norm = Math.Sqrt(e_ksi[0] * e_ksi[0] + e_ksi[1] * e_ksi[1] + e_ksi[2] * e_ksi[2]);
				for (int l = 0; l < 3; l++)
				{
					e_3[l] = e_3[l] / e_3_norm;
					e_1[l] = e_ksi[l] / e_ksi_norm;
				}
				this.Cross(e_1, e_3, e_2);
				for (int l = 0; l < 3; l++)
				{
					R[npoint1][l, 0] = e_1[l];
					R[npoint1][l, 1] = e_2[l];
					R[npoint1][l, 2] = e_3[l];

				}

				this.Cross(e_ksi, e_heta, c_1[npoint1]);
				coh_det_J_t[npoint1] = Math.Sqrt(c_1[npoint1][0] * c_1[npoint1][0] + c_1[npoint1][1] * c_1[npoint1][1] + c_1[npoint1][2] * c_1[npoint1][2]);
				integrationsCoeffs[npoint1] = coh_det_J_t[npoint1] * QuadratureForStiffness.IntegrationPoints[npoint1].Weight;

				// Calculate RtN3 here instead of in InitializeRN3() and then in UpdateForces()
				RtN3[npoint1] = R[npoint1].Transpose() * N3[npoint1];
			}
			return new Tuple<Matrix[], double[]>(RtN3, integrationsCoeffs);
		}

		private void Cross(double[] A, double[] B, double[] C) //TODO: replace with linear algebra method
		{
			C[0] = A[1] * B[2] - A[2] * B[1];
			C[1] = A[2] * B[0] - A[0] * B[2];
			C[2] = A[0] * B[1] - A[1] * B[0];
		}

		private double[,] CalculateTMatrix(IElement element)
		{
			double[,] T;
			T = new double[24, 40];
			double[,] eye3 = new double[3, 3];
			for (int m = 0; m < 3; m++)
			{
				eye3[m, m] = 1;
			}

			for (int m = 0; m < 8; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					for (int l = 0; l < 3; l++)
					{
						T[3 * m + n, 5 * m + l] = eye3[n, l];
					}
				}
			}

			if (ShellElementSide == 0)
			{
				for (int m = 0; m < 8; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						T[3 * m + n, 5 * m + 3] = 0.5 * tk[m] * tUvec[m][3 + n];
						T[3 * m + n, 5 * m + 4] = -0.5 * tk[m] * tUvec[m][n];
					}
				}
			}
			else
			{
				for (int m = 0; m < 8; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						T[3 * m + n, 5 * m + 3] = -0.5 * tk[m] * tUvec[m][3 + n];
						T[3 * m + n, 5 * m + 4] = +0.5 * tk[m] * tUvec[m][n];
					}
				}
			}

			return T;
		}

		private double[] MultiplyForcesForEmbedding(double[] fxk1_coh, IElement element)
		{
			double[,] T = CalculateTMatrix(element);
			double[] fxk2_coh = new double[64];
			if (ShellElementSide == 0)
			{
				for (int n = 0; n < 40; n++)
				{
					for (int p = 0; p < 24; p++)
					{
						fxk2_coh[n] += T[p, n] * fxk1_coh[p];
					}
				}
				for (int n = 40; n < 64; n++) //24 dof values will be updated starting from index n-16 = 40-16
				{
					fxk2_coh[n] = fxk1_coh[n - 16];
				}

			}
			else
			{
				for (int n = 0; n < 24; n++)
				{
					fxk2_coh[40 + n] = fxk1_coh[n];
				}
				for (int n = 0; n < 40; n++)
				{
					for (int p = 0; p < 24; p++)
					{
						fxk2_coh[n] += T[p, n] * fxk1_coh[p + 24];
					}
				}
			}

			return fxk2_coh;
		}

		private double[,] MultiplyStifnessMatrixForEmbedding(double[,] k_element_coh, IElement element)
		{
			double[,] T = CalculateTMatrix(element);
			double[,] k_element_coh2 = new double[64, 64];
			double[,] Kii_A = new double[24, 40];
			if (ShellElementSide == 0)
			{
				for (int n = 0; n < 24; n++)
				{
					for (int p = 0; p < 40; p++)
					{
						for (int k = 0; k < 24; k++)
						{
							Kii_A[n, p] += k_element_coh[n, k] * T[k, p];
						}
					}
				}

				for (int n = 0; n < 40; n++)
				{
					for (int p = 0; p < 40; p++)
					{
						for (int k = 0; k < 24; k++)
						{
							k_element_coh2[n, p] += T[k, n] * Kii_A[k, p];
						}
					}
				}

				for (int n = 0; n < 40; n++)
				{
					for (int p = 0; p < 24; p++)
					{
						for (int k = 0; k < 24; k++)
						{
							k_element_coh2[n, 40 + p] += T[k, n] * k_element_coh[k, 24 + p];
						}
					}
				}

				for (int n = 0; n < 24; n++)
				{
					for (int p = 0; p < 40; p++)
					{
						for (int k = 0; k < 24; k++)
						{
							k_element_coh2[40 + n, p] += k_element_coh[24 + n, k] * T[k, p];
						}
					}
				}

				// Bottom right submatrix of Tt * K * T
				for (int n = 0; n < 24; n++)
				{
					for (int p = 0; p < 24; p++)
					{
						k_element_coh2[40 + n, 40 + p] = k_element_coh[24 + n, 24 + p];
					}
				}
			}
			else
			{
				for (int n = 0; n < 24; n++)
				{
					for (int p = 0; p < 40; p++)
					{
						for (int k = 0; k < 24; k++)
						{
							Kii_A[n, p] += k_element_coh[24 + n, 24 + k] * T[k, p];
						}
					}
				}

				// Copy upper left submatrix of Tt*K*T to bottom right
				for (int n = 0; n < 24; n++)
				{
					for (int p = 0; p < 24; p++)
					{
						k_element_coh2[40 + n, 40 + p] = k_element_coh[n, p];
					}
				}

				// Copy upper right submatrix of Tt*K*T to bottom left
				for (int n = 0; n < 24; n++)
				{
					for (int p = 0; p < 40; p++)
					{
						for (int k = 0; k < 24; k++)
						{
							k_element_coh2[40 + n, p] += k_element_coh[n, 24 + k] * T[k, p];
						}
					}
				}

				// submatrix 21 Tt_K_T -->12
				for (int n = 0; n < 40; n++)
				{
					for (int p = 0; p < 24; p++)
					{
						for (int k = 0; k < 24; k++)
						{
							k_element_coh2[n, 40 + p] += T[k, n] * k_element_coh[24 + k, p];
						}
					}
				}

				// submatrix 22 of Tt_K_T  -->11
				for (int n = 0; n < 40; n++)
				{
					for (int p = 0; p < 40; p++)
					{
						for (int k = 0; k < 24; k++)
						{
							k_element_coh2[n, p] += T[k, n] * Kii_A[k, p];
						}
					}
				}

			}

			return k_element_coh2;
		}

		private double[] UpdateForces(IElement element, Matrix[] RtN3, double[] integrationCoeffs)
		{
			double[] fxk2_coh = new double[64];
			double[] fxk1_coh = new double[48]; // TODO: must be 24 in cohesive 8 node

			for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
			{
				double[] T_int_integration_coeffs = new double[3];
				for (int l = 0; l < 3; l++)
				{
					T_int_integration_coeffs[l] = materialsAtGaussPoints[npoint1].Tractions[l] * integrationCoeffs[npoint1];
				}

				double[] r_int_1 = new double[24];
				for (int l = 0; l < 24; l++)
				{
					for (int m = 0; m < 3; m++)
					{ r_int_1[l] += RtN3[npoint1][m, l] * T_int_integration_coeffs[m]; }
				}
				for (int l = 0; l < 24; l++)
				{
					fxk1_coh[l] += r_int_1[l];
					fxk1_coh[24 + l] += (-r_int_1[l]);
				}
			}

			fxk2_coh = this.MultiplyForcesForEmbedding(fxk1_coh, element);
			return fxk2_coh;
		}

		private double[,] UpdateKmatrices(IElement element, Matrix[] RtN3, double[] integrationCoeffs)
		{
			double[,] k_element_coh2 = new double[64, 64];
			double[,] k_stoixeiou_coh = new double[48, 48];


			for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
			{
				Matrix D_tan_sunt_ol = Matrix.CreateZero(3, 3);
				for (int l = 0; l < 3; l++)
				{
					for (int m = 0; m < 3; m++)
					{
						D_tan_sunt_ol[l, m] = materialsAtGaussPoints[npoint1].ConstitutiveMatrix[l, m] * integrationCoeffs[npoint1];// D_tan[npoint1][l, m] * integrationCoeffs[npoint1];
					}
				}

				Matrix D_RtN3_sunt_ol = D_tan_sunt_ol * RtN3[npoint1];
				Matrix M = RtN3[npoint1].Transpose() * D_RtN3_sunt_ol;

				for (int l = 0; l < 24; l++)
				{
					for (int m = 0; m < 24; m++)
					{
						k_stoixeiou_coh[l, m] += M[l, m];
						k_stoixeiou_coh[l, 24 + m] += -M[l, m];
						k_stoixeiou_coh[24 + l, m] += -M[l, m];
						k_stoixeiou_coh[24 + l, 24 + m] += M[l, m];
					}
				}
			}

			k_element_coh2 = this.MultiplyStifnessMatrixForEmbedding(k_stoixeiou_coh, element);
			return k_element_coh2;
		}

		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localTotalDisplacementsSuperElement, double[] localdDisplacementsSuperElement)
		{
			double[][] Delta = new double[nGaussPoints][];
			double[] localTotalDisplacements = dofEnumerator.GetTransformedDisplacementsVector(localTotalDisplacementsSuperElement); // embedding
																																	 //double[] localDisplacements = dofEnumerator.GetTransformedVector(localdDisplacementsSuperElement); // embedding

			Delta = this.UpdateCoordinateDataAndCalculateDisplacementVector(localTotalDisplacements);
			for (int i = 0; i < materialsAtGaussPoints.Length; i++)
			{
				materialsAtGaussPoints[i].UpdateMaterial(Delta[i]);
			}
			return new Tuple<double[], double[]>(Delta[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Tractions);
		}

		public double[] CalculateForces(IElement element, double[] localTotalDisplacementsSuperElement, double[] localdDisplacementsSuperelement)
		{
			double[] fxk2_coh = new double[64];
			Tuple<Matrix[], double[]> RtN3AndIntegrationCoeffs;
			RtN3AndIntegrationCoeffs = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations();
			Matrix[] RtN3;
			RtN3 = RtN3AndIntegrationCoeffs.Item1;
			double[] integrationCoeffs;
			integrationCoeffs = RtN3AndIntegrationCoeffs.Item2;

			fxk2_coh = this.UpdateForces(element, RtN3, integrationCoeffs);
			return dofEnumerator.GetTransformedForcesVector(fxk2_coh);// embedding
		}

		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
			=> CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

		public virtual IMatrix StiffnessMatrix(IElement element)
		{
			double[,] k_stoixeiou_coh2 = new double[64, 64];
			if (MatrixIsNotInitialized)
			{
				this.GetInitialGeometricDataAndInitializeMatrices(element);
				this.UpdateCoordinateDataAndCalculateDisplacementVector(new double[64]); //returns Delta that can't be used for the initial material state
				MatrixIsNotInitialized = false;
			}

			Tuple<Matrix[], double[]> RtN3AndIntegrationCoeffs;
			RtN3AndIntegrationCoeffs = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations();
			Matrix[] RtN3;
			RtN3 = RtN3AndIntegrationCoeffs.Item1;
			double[] integrationCoeffs;
			integrationCoeffs = RtN3AndIntegrationCoeffs.Item2;

			k_stoixeiou_coh2 = this.UpdateKmatrices(element, RtN3, integrationCoeffs);
			IMatrix element_stiffnessMatrix = Matrix.CreateFromArray(k_stoixeiou_coh2);
			return dofEnumerator.GetTransformedMatrix(element_stiffnessMatrix); // embedding
		}

		public void ResetMaterialModified()
		{
			foreach (ICohesiveZoneMaterial3D material in materialsAtGaussPoints) material.ResetModified();
		}

		public void ClearMaterialState()
		{
			foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.ClearState();
		}

		public void SaveMaterialState()
		{
			foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.SaveState();
		}

		public void ClearMaterialStresses()
		{
			foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.ClearTractions();
		}

		public bool MaterialModified
		{
			get
			{
				foreach (ICohesiveZoneMaterial3D material in materialsAtGaussPoints)
					if (material.Modified) return true;
				return false;
			}
		}

		public virtual IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;

		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			return new double[64];
		}

		public virtual IMatrix MassMatrix(IElement element)
		{
			return Matrix.CreateZero(64, 64);
		}

		public virtual IMatrix DampingMatrix(IElement element)
		{
			return Matrix.CreateZero(64, 64);
		}

		#region EMBEDDED

		public Dictionary<IDofType, int> GetInternalNodalDOFs(Element element, Node node)//
		{
			int index = 0;
			foreach (var elementNode in element.Nodes)
			{
				if (node.ID == elementNode.ID)
					break;
				index++;
			}
			if (index >= 16)
				throw new ArgumentException(String.Format("GetInternalNodalDOFs: Node {0} not found in element {1}.", node.ID, element.ID));

			if (index >= 8)
			{
				int index2 = index - 8;
				return new Dictionary<IDofType, int>() { { StructuralDof.TranslationX, 39 + 3 * index2 + 1 }, { StructuralDof.TranslationY, 39 + 3 * index2 + 2 }, { StructuralDof.TranslationZ, 39 + 3 * index2 + 3 } };
			}
			else
			{
				return new Dictionary<IDofType, int>() { { StructuralDof.TranslationX, + 5 * index + 0 }, { StructuralDof.TranslationY, + 5 * index + 1 }, { StructuralDof.TranslationZ, + 5 * index + 2 },
														{ StructuralDof.RotationX, + 5 * index + 3 }, { StructuralDof.RotationY, + 5 * index + 4 }};
			}
		}

		public double[] GetLocalDOFValues(Element hostElement, double[] hostDOFValues) // omoiws Beam3D
		{
			//if (transformation == null)
			//    throw new InvalidOperationException("Requested embedded node values for element that has no embedded nodes.");
			//if (hostElementList == null)
			//    throw new InvalidOperationException("Requested host element list for element that has no embedded nodes.");
			//int index = hostElementList.IndexOf(hostElement);
			//if (index < 0)
			//    throw new ArgumentException("Requested host element is not inside host element list.");

			//double[] values = new double[transformation.Columns];
			//int multiplier = hostElement.ElementType.DofEnumerator.GetDOFTypes(hostElement).SelectMany(d => d).Count();
			//int vectorIndex = 0;
			//for (int i = 0; i < index; i++)
			//    vectorIndex += isNodeEmbedded[i] ? 3 : multiplier;
			//Array.Copy(hostDOFValues, 0, values, vectorIndex, multiplier);

			//return (transformation * new Vector<double>(values)).Data;

			return dofEnumerator.GetTransformedDisplacementsVector(hostDOFValues);
		}
		#endregion
	}
}
