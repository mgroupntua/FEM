using System;
using System.Collections.Generic;
using System.Linq;
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
using MGroup.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.FEM.Elements
{
	public class cohesiveElement : IStructuralFiniteElement
	{
		protected readonly static IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		protected readonly IDofType[][] dofTypes;
		private readonly IIsoparametricInterpolation2D interpolation;
		protected readonly ICohesiveZoneMaterial3D[] materialsAtGaussPoints;
		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();



		private int nGaussPoints;

		public bool MatrixIsNotInitialized = true;

		protected cohesiveElement()
		{
		}

		public cohesiveElement(ICohesiveZoneMaterial3D material, IQuadrature2D quadratureForStiffness, IIsoparametricInterpolation2D interpolation)
		{
			this.QuadratureForStiffness = quadratureForStiffness;
			this.interpolation = interpolation;
			this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;

			materialsAtGaussPoints = new ICohesiveZoneMaterial3D[nGaussPoints];
			for (int i = 0; i < nGaussPoints; i++)
				materialsAtGaussPoints[i] = material.Clone();

			dofTypes = new IDofType[interpolation.NumFunctions * 2][];
			for (int i = 0; i < interpolation.NumFunctions * 2; i++)
			{
				dofTypes[i] = new IDofType[]
				{
					StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
				};
			}
		}

		//public IReadOnlyList<IFiniteElementMaterial> Materials => materialsAtGaussPoints;

		public CellType CellType { get; } = CellType.Unknown;

		public IQuadrature2D QuadratureForStiffness { get; }


		public int endeixiShapeFunctionAndGaussPointData = 1;




		private double[][] intialCoordinates; //den einai apo afta pou orizei o xrhsths
		private double[] deformedCoordinates; // to dianusma x ths matlab sunarthshs pou einai apo t_x_global_pr
		private void GetInitialGeometricDataAndInitializeMatrices(IElement element)
		{
			intialCoordinates = new double[2 * interpolation.NumFunctions][];
			for (int j = 0; j < 2 * interpolation.NumFunctions; j++)
			{
				intialCoordinates[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
			}

			deformedCoordinates = new double[2 * 3 * interpolation.NumFunctions];


		}

		private double[][] UpdateCoordinateData(double[] localdisplacements) // sto shell8disp sto calculate forces kaleitai me this.UpdateCoordinateData(localTotalDisplacements);
		{
			IReadOnlyList<Matrix> shapeFunctionDerivatives = interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			IReadOnlyList<double[]> shapeFunctionsAtGpoints = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);

			double[,] deltaUCartesian = new double[3, 2 * interpolation.NumFunctions];
			double[,] xMidsurface = new double[3, 2 * interpolation.NumFunctions];

			double[] tangentVectorKsi = new double[3];
			double tangentVectorKsiNorm;
			double[] tangentVectorHeta = new double[3];
			double[] tangentVectorNormalised1 = new double[3];
			double[] tangentVectorNormalised2 = new double[3];
			double[] normalVectorE3 = new double[3];
			double normalVectorNOrm;
			double[] deltaUGpoint = new double[3];


			double[][] deltaUlocal = new double[nGaussPoints][];
			double[][] auxC1 = new double[nGaussPoints][];
			for (int j = 0; j < nGaussPoints; j++)
			{
				deltaUlocal[j] = new double[3];
				auxC1[j] = new double[3];
			}
			double[][,] R = new double[nGaussPoints][,]; //TODO: maybe cache R
			for (int j = 0; j < nGaussPoints; j++)
			{
				R[j] = new double[3, 3];
			}

			for (int j = 0; j < 2 * interpolation.NumFunctions; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					deformedCoordinates[3 * j + k] = intialCoordinates[j][k] + localdisplacements[3 * j + k];
				}
			}

			for (int j = 0; j < interpolation.NumFunctions; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					deltaUCartesian[k, j] = deformedCoordinates[k + 3 * j] - deformedCoordinates[3 * interpolation.NumFunctions + k + 3 * j];
					xMidsurface[k, j] = deformedCoordinates[k + 3 * j] + deformedCoordinates[3 * interpolation.NumFunctions + k + 3 * j];
				}
			}
			// sunexeia ews upologismou tou Delta gia ola ta gp

			for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
			{
				for (int l = 0; l < 3; l++)
				{
					tangentVectorKsi[l] = 0;
					tangentVectorHeta[l] = 0;
					for (int m = 0; m < interpolation.NumFunctions; m++) // tha ginei 4 sto cohesive 8 node
					{
						tangentVectorKsi[l] += shapeFunctionDerivatives[npoint1][m, 0] * xMidsurface[l, m];
						tangentVectorHeta[l] += shapeFunctionDerivatives[npoint1][m, 1] * xMidsurface[l, m];
					}
					tangentVectorKsi[l] = 0.5 * tangentVectorKsi[l];
					tangentVectorHeta[l] = 0.5 * tangentVectorHeta[l];
				}
				normalVectorE3 = Vector.CreateFromArray(tangentVectorKsi).CrossProduct(Vector.CreateFromArray(tangentVectorHeta)).CopyToArray();
				//this.Cross(tangentVectorKsi, tangentVectorHeta, normalVectorE3);
				normalVectorNOrm = Math.Sqrt(normalVectorE3[0] * normalVectorE3[0] + normalVectorE3[1] * normalVectorE3[1] + normalVectorE3[2] * normalVectorE3[2]);
				tangentVectorKsiNorm = Math.Sqrt(tangentVectorKsi[0] * tangentVectorKsi[0] + tangentVectorKsi[1] * tangentVectorKsi[1] + tangentVectorKsi[2] * tangentVectorKsi[2]);
				for (int l = 0; l < 3; l++)
				{
					normalVectorE3[l] = normalVectorE3[l] / normalVectorNOrm;
					tangentVectorNormalised1[l] = tangentVectorKsi[l] / tangentVectorKsiNorm;
				}
				tangentVectorNormalised2 = Vector.CreateFromArray(tangentVectorNormalised1).CrossProduct(Vector.CreateFromArray(normalVectorE3)).CopyToArray();
				//this.Cross(tangentVectorNormalised1, normalVectorE3, tangentVectorNormalised2);
				for (int l = 0; l < 3; l++)
				{
					R[npoint1][l, 0] = tangentVectorNormalised1[l];
					R[npoint1][l, 1] = tangentVectorNormalised2[l];
					R[npoint1][l, 2] = normalVectorE3[l];

				}
				for (int l = 0; l < 3; l++)
				{ deltaUGpoint[l] = 0; }
				for (int l = 0; l < 3; l++)
				{
					for (int m = 0; m < interpolation.NumFunctions; m++)  // pithanws gia to cohesive 8 node na gineiews 4 to m
					{
						deltaUGpoint[l] += deltaUCartesian[l, m] * shapeFunctionsAtGpoints[npoint1][m];
					}
				}
				for (int l = 0; l < 3; l++)
				{ deltaUlocal[npoint1][l] = 0; }
				for (int l = 0; l < 3; l++)
				{
					for (int m = 0; m < 3; m++)
					{
						deltaUlocal[npoint1][l] += R[npoint1][m, l] * deltaUGpoint[m];
					}
				}

				auxC1[npoint1] = Vector.CreateFromArray(tangentVectorKsi).CrossProduct(Vector.CreateFromArray(tangentVectorHeta)).CopyToArray();
				//this.Cross(tangentVectorKsi, tangentVectorHeta, auxC1[npoint1]);

			}

			return deltaUlocal;

		}

		private double[,] ReShapeShapeFunctionDataMatrix(double[] N1)
		{
			var auxillaryN3 = new double[3, 3 * interpolation.NumFunctions];
			for (int l = 0; l < 3; l++)
			{
				for (int m = 0; m < interpolation.NumFunctions; m++)
				{ auxillaryN3[l, l + 3 * m] = N1[m]; }
			}
			return auxillaryN3;
		}
		private Tuple<Matrix[], double[]> CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations()
		{
			IReadOnlyList<double[]> N1 = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
			IReadOnlyList<Matrix> N3 = N1.Select(x => Matrix.CreateFromArray(ReShapeShapeFunctionDataMatrix(x))).ToList();
			IReadOnlyList<Matrix> shapeFunctionDerivatives = interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			double[] integrationsCoeffs = new double[nGaussPoints];
			Matrix[] auxRtTimesN3 = new Matrix[nGaussPoints];
			double[,] auxMidsurfaceDeformedSum = new double[3, 8];

			double[] tangentVector1normalised = new double[3];
			double[] tangentVector2normalised = new double[3];
			double[] normalVector3normalised = new double[3];
			double normalVectorNorm;

			double[] auxDeformedDeterminant = new double[nGaussPoints];

			double[][] auxC1 = new double[nGaussPoints][];
			for (int j = 0; j < nGaussPoints; j++)
			{
				auxC1[j] = new double[3];
			}

			Matrix[] R = new Matrix[nGaussPoints]; //TODO: perhaps cache matrices in InitializeMatrices() where RtN3 is calculated
			for (int j = 0; j < nGaussPoints; j++)
			{ R[j] = Matrix.CreateZero(3, 3); }

			for (int j = 0; j < interpolation.NumFunctions; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					auxMidsurfaceDeformedSum[k, j] = deformedCoordinates[k + 3 * j] + deformedCoordinates[3 * interpolation.NumFunctions + k + 3 * j];
				}
			}

			// Calculate Delta for all GPs
			for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
			{
				double[] tangentVectorKsi = new double[3];
				double[] tangentVectorHEta = new double[3];
				for (int l = 0; l < 3; l++)
				{
					for (int m = 0; m < interpolation.NumFunctions; m++) // must be 4 in cohesive 8-nodes
					{
						tangentVectorKsi[l] += shapeFunctionDerivatives[npoint1][m, 0] * auxMidsurfaceDeformedSum[l, m];
						tangentVectorHEta[l] += shapeFunctionDerivatives[npoint1][m, 1] * auxMidsurfaceDeformedSum[l, m];
					}
					tangentVectorKsi[l] = 0.5 * tangentVectorKsi[l];
					tangentVectorHEta[l] = 0.5 * tangentVectorHEta[l];
				}
				normalVector3normalised = Vector.CreateFromArray(tangentVectorKsi).CrossProduct(Vector.CreateFromArray(tangentVectorHEta)).CopyToArray();
				//this.Cross(tangentVectorKsi, tangentVectorHEta, normalVector3normalised);
				normalVectorNorm = Math.Sqrt(normalVector3normalised[0] * normalVector3normalised[0] + normalVector3normalised[1] * normalVector3normalised[1] + normalVector3normalised[2] * normalVector3normalised[2]);
				double tangentVectorKsiNorm = Math.Sqrt(tangentVectorKsi[0] * tangentVectorKsi[0] + tangentVectorKsi[1] * tangentVectorKsi[1] + tangentVectorKsi[2] * tangentVectorKsi[2]);
				for (int l = 0; l < 3; l++)
				{
					normalVector3normalised[l] = normalVector3normalised[l] / normalVectorNorm;
					tangentVector1normalised[l] = tangentVectorKsi[l] / tangentVectorKsiNorm;
				}
				tangentVector2normalised = Vector.CreateFromArray(tangentVector1normalised).CrossProduct(Vector.CreateFromArray(normalVector3normalised)).CopyToArray();
				//this.Cross(tangentVector1normalised, normalVector3normalised, tangentVector2normalised);
				for (int l = 0; l < 3; l++)
				{
					R[npoint1][l, 0] = tangentVector1normalised[l];
					R[npoint1][l, 1] = tangentVector2normalised[l];
					R[npoint1][l, 2] = normalVector3normalised[l];

				}

				auxC1[npoint1] = Vector.CreateFromArray(tangentVectorKsi).CrossProduct(Vector.CreateFromArray(tangentVectorHEta)).CopyToArray();
				//.this.Cross(tangentVectorKsi, tangentVectorHEta, auxC1[npoint1]);
				auxDeformedDeterminant[npoint1] = Math.Sqrt(auxC1[npoint1][0] * auxC1[npoint1][0] + auxC1[npoint1][1] * auxC1[npoint1][1] + auxC1[npoint1][2] * auxC1[npoint1][2]);
				integrationsCoeffs[npoint1] = auxDeformedDeterminant[npoint1] * QuadratureForStiffness.IntegrationPoints[npoint1].Weight;

				// Calculate RtN3 here instead of in InitializeRN3() and then in UpdateForces()
				auxRtTimesN3[npoint1] = R[npoint1].Transpose() * N3[npoint1];
			}
			return new Tuple<Matrix[], double[]>(auxRtTimesN3, integrationsCoeffs);
		}




		private double[] UpdateForces(IElement element, Matrix[] RtN3, double[] integrationCoeffs)
		{
			double[] elementForces = new double[3 * 2 * interpolation.NumFunctions]; // TODO: must be 24 in cohesive 8 node

			for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
			{
				double[] tractionsTimeIntegrationCoeefs = new double[3];
				for (int l = 0; l < 3; l++)
				{
					tractionsTimeIntegrationCoeefs[l] = materialsAtGaussPoints[npoint1].Tractions[l] * integrationCoeffs[npoint1];
				}

				double[] rVector = new double[3 * interpolation.NumFunctions];
				for (int l = 0; l < 3 * interpolation.NumFunctions; l++)
				{
					for (int m = 0; m < 3; m++)
					{ rVector[l] += RtN3[npoint1][m, l] * tractionsTimeIntegrationCoeefs[m]; }
				}
				for (int l = 0; l < 3 * interpolation.NumFunctions; l++)
				{
					elementForces[l] += rVector[l];
					elementForces[3 * interpolation.NumFunctions + l] += (-rVector[l]);
				}
			}

			return elementForces;
		}

		private double[,] UpdateKmatrices(IElement element, Matrix[] auxRtTimesN3, double[] integrationCoeffs)
		{
			double[,] elementStiffness = new double[3 * 2 * interpolation.NumFunctions, 3 * 2 * interpolation.NumFunctions];


			for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
			{
				Matrix DtangentTimesIntegrationCoeffs = Matrix.CreateZero(3, 3);
				for (int l = 0; l < 3; l++)
				{
					for (int m = 0; m < 3; m++)
					{
						DtangentTimesIntegrationCoeffs[l, m] = materialsAtGaussPoints[npoint1].ConstitutiveMatrix[l, m] * integrationCoeffs[npoint1];// D_tan[npoint1][l, m] * integrationCoeffs[npoint1];
					}
				}

				Matrix auxDtanTimesRtTimesN3TimesIntegrcoef = DtangentTimesIntegrationCoeffs * auxRtTimesN3[npoint1];
				Matrix M = auxRtTimesN3[npoint1].Transpose() * auxDtanTimesRtTimesN3TimesIntegrcoef;

				for (int l = 0; l < 3 * interpolation.NumFunctions; l++)
				{
					for (int m = 0; m < 3 * interpolation.NumFunctions; m++)
					{
						elementStiffness[l, m] += M[l, m];
						elementStiffness[l, 3 * interpolation.NumFunctions + m] += -M[l, m];
						elementStiffness[3 * interpolation.NumFunctions + l, m] += -M[l, m];
						elementStiffness[3 * interpolation.NumFunctions + l, 3 * interpolation.NumFunctions + m] += M[l, m];
					}
				}
			}

			return elementStiffness;
		}

		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localTotalDisplacements, double[] localdDisplacements)
		{
			double[][] Delta = UpdateCoordinateData(localTotalDisplacements);
			for (int i = 0; i < materialsAtGaussPoints.Length; i++)
			{
				materialsAtGaussPoints[i].UpdateMaterial(Delta[i]);
			}
			return new Tuple<double[], double[]>(Delta[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Tractions);
			//TODO mono to teleftaio dianusma tha epistrefei?
		}

		public double[] CalculateForces(IElement element, double[] localTotalDisplacements, double[] localdDisplacements)
		{
			Tuple<Matrix[], double[]> auxRtTimesN3AndIntegrationCoeffs;
			auxRtTimesN3AndIntegrationCoeffs = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations();
			Matrix[] auxRtTimesN3;
			auxRtTimesN3 = auxRtTimesN3AndIntegrationCoeffs.Item1;
			double[] integrationCoeffs;
			integrationCoeffs = auxRtTimesN3AndIntegrationCoeffs.Item2;

			double[] elementForces = UpdateForces(element, auxRtTimesN3, integrationCoeffs);
			return elementForces;
		}

		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
		{
			return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
		}

		public virtual IMatrix StiffnessMatrix(IElement element)
		{
			double[,] elementStiffness;
			if (MatrixIsNotInitialized)
			{
				this.GetInitialGeometricDataAndInitializeMatrices(element);
				this.UpdateCoordinateData(new double[2 * 3 * interpolation.NumFunctions]); //returns Delta that can't be used for the initial material state
				MatrixIsNotInitialized = false;
			}

			Tuple<Matrix[], double[]> RtN3AndIntegrationCoeffs;
			RtN3AndIntegrationCoeffs = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations();
			Matrix[] auxRtTimesN3;
			auxRtTimesN3 = RtN3AndIntegrationCoeffs.Item1;
			double[] integrationCoeffs;
			integrationCoeffs = RtN3AndIntegrationCoeffs.Item2;

			elementStiffness = this.UpdateKmatrices(element, auxRtTimesN3, integrationCoeffs);
			IMatrix elementStiffnessMatrix = Matrix.CreateFromArray(elementStiffness);
			return elementStiffnessMatrix; // embedding
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
			////temporary
			//n_incr += 1;
			//if (n_incr == 17)
			//{ n_incr += 0; }
		}

		public void ClearMaterialStresses()
		{
			foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.ClearTractions();
		}

		// omoiws me hexa 8 shell8disp implemented
		public int ID
		{
			get { return 14; }
		}
		public ElementDimensions ElementDimensions
		{
			get { return ElementDimensions.ThreeD; }
		}

		public IElementDofEnumerator DOFEnumerator
		{
			get { return DofEnumerator; }
			set { DofEnumerator = value; }
		}

		public virtual IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;

		// aplopoihtika implemented mhdenikes masses gia cohesive - not implemented
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			return new double[2 * 3 * interpolation.NumFunctions];
		}

		public virtual IMatrix MassMatrix(IElement element)
		{
			return Matrix.CreateZero(3 * 2 * interpolation.NumFunctions, 3 * 2 * interpolation.NumFunctions);
		}

		public virtual IMatrix DampingMatrix(IElement element)
		{

			return Matrix.CreateZero(3 * 2 * interpolation.NumFunctions, 3 * 2 * interpolation.NumFunctions);
		}

	}

}
