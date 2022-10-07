using System;
using System.Collections.Generic;
using System.Xml.Linq;

using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.Cohesive;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Constitutive.Structural.Line;
using MGroup.FEM.Structural.Line;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Embedding;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Numerics.Interpolation;

namespace MGroup.FEM.Structural.Embedding
{
    public class CohesiveBeam3DToBeam3D : IStructuralElementType, IEmbeddedElement
	{
		private List<INode> nodes = new List<INode>();
        protected readonly static IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY, StructuralDof.RotationZ };
        protected readonly static IDofType[][] dofTypes = new IDofType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
            nodalDOFTypes };
        protected readonly ICohesiveZoneMaterial[] materialsAtGaussPoints;
        private readonly InterpolationTruss1D interpolation = InterpolationTruss1D.UniqueInstance;
        private int nGaussPoints;
        protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
        private SupportiveCohesiveBeam3DCorotationalQuaternion supportive_beam;
        private SupportiveCohesiveBeam3DCorotationalQuaternion supportive_clone;
        private readonly double perimeter;

        /// <summary>
        /// Initial nodel coordinates of 4 node inner cohesive element
        /// </summary>
        private double[][] initialNodalCoordinates; // ox_i

        /// <summary>
        /// Unrolled current nodal coordinates of 4 node inner cohesive element
        /// </summary>
        private double[] currentNodalCoordinates; // x_local

        public bool MatrixIsNotInitialized = true;
        public double[,] k_cohesive_element_total {get; private set;}
		private double[][] tractions;
		private double[] localDisplacements;
		private double[] localDisplacementsLastConverged;
		//protected CohesiveBeam3DToBeam3D()
		//{
		//}

		public CohesiveBeam3DToBeam3D(List<INode> nodes, ICohesiveZoneMaterial material, IQuadrature1D quadratureForStiffness, IList<INode> nodes_beam,
			IList<INode> nodes_clone, IIsotropicContinuumMaterial3D material_beam, double density, BeamSection3D beamSection, double perimeter)
        {
            this.QuadratureForStiffness = quadratureForStiffness;
            this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;
            materialsAtGaussPoints = new ICohesiveZoneMaterial[nGaussPoints];
			for (int i = 0; i < nGaussPoints; i++)
			{
				var materialCopy = (ICloneable)material;
				materialsAtGaussPoints[i] = (ICohesiveZoneMaterial)materialCopy.Clone();
			}
            this.supportive_beam = new SupportiveCohesiveBeam3DCorotationalQuaternion(nodes_beam, material_beam, density, beamSection);
            this.supportive_clone = new SupportiveCohesiveBeam3DCorotationalQuaternion(nodes_clone, material_beam, density, beamSection);
            this.perimeter = perimeter;
			this.nodes = nodes;
        }

        public IQuadrature1D QuadratureForStiffness { get; }

        private void GetInitialGeometricDataAndInitializeMatrices(IElementType element)
        {
            initialNodalCoordinates = new double[4][];

            //TODOgsoim: why 0 to 4??
            for (int j = 0; j < 4; j++)
            { initialNodalCoordinates[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, }; }

            currentNodalCoordinates = new double[12]; //12:without rotations, 24:with rotations
        }

        private double[][] UpdateCoordinateDataAndCalculateDisplacementVector(double[] localdisplacements)
        {
            //Matrix shapeFunctionDerivatives = interpolation.EvaluateGradientsAt();
            IReadOnlyList<Vector> N1 = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            //IReadOnlyList<Matrix> N3 = interpolation.EvaluateN3ShapeFunctionsReorganized(QuadratureForStiffness); //Shape functions matrix [N_beam]

            double[,] u_prok = new double[3, 2];// [3d-axes, nodes] - of the mid-surface(#i_m, #j_m) 
            double[,] x_bar = new double[3, 2];
            double[] u = new double[3];

            double[][] Delta = new double[nGaussPoints][];
            for (int j = 0; j < nGaussPoints; j++)
            {
                Delta[j] = new double[3];
            }

            //double[][,] R = new double[nGaussPoints][,]; //TODO: maybe cache R
            //for (int j = 0; j < nGaussPoints; j++)
            //{
            //    R[j] = new double[3, 3];
            //    R[j][0, 0] = 0.5 * (supportive_beam.currentRotationMatrix[0, 0] + supportive_clone.currentRotationMatrix[0, 0]);
            //    R[j][0, 1] = 0.5 * (supportive_beam.currentRotationMatrix[0, 1] + supportive_clone.currentRotationMatrix[0, 1]);
            //    R[j][0, 2] = 0.5 * (supportive_beam.currentRotationMatrix[0, 2] + supportive_clone.currentRotationMatrix[0, 2]);
            //    R[j][1, 0] = 0.5 * (supportive_beam.currentRotationMatrix[1, 0] + supportive_clone.currentRotationMatrix[1, 0]);
            //    R[j][1, 1] = 0.5 * (supportive_beam.currentRotationMatrix[1, 1] + supportive_clone.currentRotationMatrix[1, 1]);
            //    R[j][1, 2] = 0.5 * (supportive_beam.currentRotationMatrix[1, 2] + supportive_clone.currentRotationMatrix[1, 2]);
            //    R[j][2, 0] = 0.5 * (supportive_beam.currentRotationMatrix[2, 0] + supportive_clone.currentRotationMatrix[2, 0]);
            //    R[j][2, 1] = 0.5 * (supportive_beam.currentRotationMatrix[2, 1] + supportive_clone.currentRotationMatrix[2, 1]);
            //    R[j][2, 2] = 0.5 * (supportive_beam.currentRotationMatrix[2, 2] + supportive_clone.currentRotationMatrix[2, 2]);
            //}

            Matrix R = CalculateRotationMatrix();

            // Update x_local
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    currentNodalCoordinates[3 * j + k] = initialNodalCoordinates[j][k] + localdisplacements[6 * j + k];
                }
            }
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    u_prok[k, j] = currentNodalCoordinates[k + 3 * j] - currentNodalCoordinates[6 + k + 3 * j];
                    x_bar[k, j] = currentNodalCoordinates[k + 3 * j] + currentNodalCoordinates[6 + k + 3 * j];
                }
            }

            //Calculate Delta for all GPs
            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                for (int l = 0; l < 3; l++)
                { u[l] = 0; }
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 2; m++)  // must be 2 for beam cohesive 8-nodes
                    {
                        u[l] += u_prok[l, m] * N1[npoint1][m];
                    }
                }
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        Delta[npoint1][l] += R[m, l] * u[m];
                    }
                }
            }
            return Delta;
        }

        private Tuple<Matrix[], double[]> CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations()
        {
            //IReadOnlyList<double[]> N1 = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Matrix> N3 = EvaluateN3ShapeFunctionsReorganized(QuadratureForStiffness); //Shape functions matrix [N_beam]
                                                                                                                  //Matrix shapeFunctionDerivatives = interpolation.EvaluateGradientsAt();

            double[] integrationsCoeffs = new double[nGaussPoints];
            Matrix[] RtN3 = new Matrix[nGaussPoints];
            //Matrix[] R = new Matrix[nGaussPoints]; //TODO: perhaps cache matrices in InitializeMatrices() where RtN3 is calculated
            //for (int j = 0; j < nGaussPoints; j++)
            //{ R[j] = Matrix.CreateZero(3, 3); }

            // Calculate Delta for all GPs
            Matrix R = CalculateRotationMatrix();
            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                //double[,] u_prok = new double[3, 2];// [3d-axes, nodes] - of the mid-surface(#i_m, #j_m) 
                //double[] u = new double[3];

                double[][] Delta = new double[nGaussPoints][];
                for (int j = 0; j < nGaussPoints; j++)
                {
                    Delta[j] = new double[3];
                }

                //R[npoint1][0, 0] = 0.5 * (supportive_beam.currentRotationMatrix[0, 0] + supportive_clone.currentRotationMatrix[0, 0]);
                //R[npoint1][0, 1] = 0.5 * (supportive_beam.currentRotationMatrix[0, 1] + supportive_clone.currentRotationMatrix[0, 1]);
                //R[npoint1][0, 2] = 0.5 * (supportive_beam.currentRotationMatrix[0, 2] + supportive_clone.currentRotationMatrix[0, 2]);
                //R[npoint1][1, 0] = 0.5 * (supportive_beam.currentRotationMatrix[1, 0] + supportive_clone.currentRotationMatrix[1, 0]);
                //R[npoint1][1, 1] = 0.5 * (supportive_beam.currentRotationMatrix[1, 1] + supportive_clone.currentRotationMatrix[1, 1]);
                //R[npoint1][1, 2] = 0.5 * (supportive_beam.currentRotationMatrix[1, 2] + supportive_clone.currentRotationMatrix[1, 2]);
                //R[npoint1][2, 0] = 0.5 * (supportive_beam.currentRotationMatrix[2, 0] + supportive_clone.currentRotationMatrix[2, 0]);
                //R[npoint1][2, 1] = 0.5 * (supportive_beam.currentRotationMatrix[2, 1] + supportive_clone.currentRotationMatrix[2, 1]);
                //R[npoint1][2, 2] = 0.5 * (supportive_beam.currentRotationMatrix[2, 2] + supportive_clone.currentRotationMatrix[2, 2]);

                integrationsCoeffs[npoint1] = ((supportive_beam.currentLength + supportive_clone.currentLength) / 4.0) * QuadratureForStiffness.IntegrationPoints[npoint1].Weight;


                // Calculate RtN3 here instead of in InitializeRN3() and then in UpdateForces()
                RtN3[npoint1] = R.Transpose() * N3[npoint1];
            }
            return new Tuple<Matrix[], double[]>(RtN3, integrationsCoeffs);
        }

        private double[] UpdateForces(/*IElementType element, */Matrix[] RtN3, double[] integrationCoeffs, double[] localTotalDisplacements)
        {
            double[] fxk1_coh = new double[24]; // Beam3D: 4 nodes, 6 dofs/node (translations + rotations)

            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                double[] T_int_integration_coeffs = new double[3];
                for (int l = 0; l < 3; l++)
                {
                    T_int_integration_coeffs[l] = /*materialsAtGaussPoints[npoint1].Tractions[l]*/Tractions[npoint1][l] * integrationCoeffs[npoint1] * perimeter;
                }
                double[] r_int_1 = new double[6];
                for (int l = 0; l < 6; l++)
                {
                    for (int m = 0; m < 3; m++)
                    { r_int_1[l] += RtN3[npoint1][m, l] * T_int_integration_coeffs[m]; }
                }
                for (int l = 0; l < 3; l++)
                {
                    fxk1_coh[l] += r_int_1[l];
                    fxk1_coh[12 + l] += (-r_int_1[l]);
                }
                for (int l = 6; l < 9; l++)
                {
                    fxk1_coh[l] += r_int_1[l - 3];
                    fxk1_coh[12 + l] += (-r_int_1[l - 3]);
                }

                Matrix R = CalculateRotationMatrix();
                Matrix ConstaintsRotations = Matrix.CreateZero(3, 3);
                ConstaintsRotations[0, 0] =  1000.0;
                ConstaintsRotations[1, 1] =  1000.0;
                ConstaintsRotations[2, 2] =  1000.0;
                Matrix Constr_R = ConstaintsRotations * R;
                Matrix M2 = R.Transpose() * Constr_R;
                var r_int_2a = M2 * Vector.CreateFromArray(new double[3] {localTotalDisplacements[3]-localTotalDisplacements[15],
                localTotalDisplacements[4] - localTotalDisplacements[16], localTotalDisplacements[5]-localTotalDisplacements[17] }) * integrationCoeffs[npoint1];
                var r_int_2b = M2 * Vector.CreateFromArray(new double[3] {localTotalDisplacements[9]-localTotalDisplacements[21],
                localTotalDisplacements[10] - localTotalDisplacements[22], localTotalDisplacements[11]-localTotalDisplacements[23] }) * integrationCoeffs[npoint1];

                for (int i = 0; i < 3; i++)
                {
                    fxk1_coh[i + 3] += r_int_2a[i];
                    fxk1_coh[i + 9] += r_int_2b[i];

                    fxk1_coh[i + 15] += -r_int_2a[i];
                    fxk1_coh[i + 21] += -r_int_2b[i];
                }
            }

            //Matrix R = CalculateRotationMatrix();
            //Matrix ConstaintsRotations = Matrix.CreateZero(3, 3);
            //ConstaintsRotations[0, 0] = 1.0; // 10.0; //
            //Matrix Constr_R = ConstaintsRotations * R;
            //Matrix M2 = R.Transpose() * Constr_R;
            //var r_int_2a = M2 * Vector.CreateFromArray(new double[3] {localTotalDisplacements[3]-localTotalDisplacements[15],
            //    localTotalDisplacements[4] - localTotalDisplacements[16], localTotalDisplacements[5]-localTotalDisplacements[17] });
            //var r_int_2b = M2 * Vector.CreateFromArray(new double[3] {localTotalDisplacements[9]-localTotalDisplacements[21],
            //    localTotalDisplacements[10] - localTotalDisplacements[22], localTotalDisplacements[11]-localTotalDisplacements[23] });

            //for (int i = 0; i < 3; i++)
            //{
            //    fxk1_coh[i + 3] = r_int_2a[i];
            //    fxk1_coh[i + 9] = r_int_2b[i];

            //    fxk1_coh[i + 15] = -r_int_2a[i];
            //    fxk1_coh[i + 21] = -r_int_2b[i];
            //}


            return fxk1_coh;
        }

        private double[,] UpdateKmatrices(IElementType element, Matrix[] RtN3, double[] integrationCoeffs)
        {
            k_cohesive_element_total = new double[24, 24];
            //double[,] k_cohesive_element = new double[12, 12];

            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                Matrix D_tan_sunt_ol = Matrix.CreateZero(3, 3);
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        D_tan_sunt_ol[l, m] = materialsAtGaussPoints[npoint1].ConstitutiveMatrix[l, m] * integrationCoeffs[npoint1];
                    }
                }
                Matrix D_RtN3_sunt_ol = D_tan_sunt_ol * RtN3[npoint1];
                Matrix M = RtN3[npoint1].Transpose() * D_RtN3_sunt_ol * perimeter; //TODOgsoim: is perimeter wright here??

                //k_cohesive_element_total - Contrains rotations, that create zero-elements in the diagonal of the total stiffness matrix.
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        k_cohesive_element_total[l, m] += M[l, m];
                        k_cohesive_element_total[l, 12 + m] += -M[l, m];
                        k_cohesive_element_total[12 + l, m] += -M[l, m];
                        k_cohesive_element_total[12 + l, 12 + m] += M[l, m];
                    }
                }
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        k_cohesive_element_total[l + 6, m] += M[l + 3, m];
                        k_cohesive_element_total[l + 6, 12 + m] += -M[l + 3, m];
                        k_cohesive_element_total[l + 18, m] += -M[l + 3, m];
                        k_cohesive_element_total[l + 18, 12 + m] += M[l + 3, m];
                    }
                }
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        k_cohesive_element_total[l, m + 6] += M[l, m + 3];
                        k_cohesive_element_total[l, m + 18] += -M[l, m + 3];
                        k_cohesive_element_total[l + 12, m + 6] += -M[l, m + 3];
                        k_cohesive_element_total[l + 12, m + 18] += M[l, m + 3];
                    }
                }
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        k_cohesive_element_total[l + 6, m + 6] += M[l + 3, m + 3];
                        k_cohesive_element_total[l + 6, m + 18] += -M[l + 3, m + 3];
                        k_cohesive_element_total[l + 18, m + 6] += -M[l + 3, m + 3];
                        k_cohesive_element_total[l + 18, m + 18] += M[l + 3, m + 3];
                    }
                }

                // ***NEW***
                Matrix R = CalculateRotationMatrix();
                Matrix ConstaintsRotations = Matrix.CreateZero(3, 3);
                ConstaintsRotations[0, 0] =  1000.0;
                ConstaintsRotations[1, 1] =  1000.0;
                ConstaintsRotations[2, 2] =  1000.0;
                Matrix Constr_R = ConstaintsRotations * R;
                Matrix M2 = R.Transpose() * Constr_R * integrationCoeffs[npoint1] * perimeter; //TODOgsoim: should perimeter be here??

                for (int ii = 0; ii < 3; ii++)
                {
                    for (int jj = 0; jj < 3; jj++)
                    {
                        k_cohesive_element_total[ii + 3, jj + 3] += M2[ii, jj]; //
                        k_cohesive_element_total[ii + 3, jj + 15] += -M2[ii, jj];
                        k_cohesive_element_total[ii + 9, jj + 9] += M2[ii, jj]; //
                        k_cohesive_element_total[ii + 9, jj + 21] += -M2[ii, jj];
                        k_cohesive_element_total[ii + 15, jj + 3] += -M2[ii, jj];
                        k_cohesive_element_total[ii + 15, jj + 15] += M2[ii, jj]; //
                        k_cohesive_element_total[ii + 21, jj + 9] += -M2[ii, jj];
                        k_cohesive_element_total[ii + 21, jj + 21] += M2[ii, jj]; //
                    }
                }
            }
            return k_cohesive_element_total; //k_cohesive_element; //   
        }

        public Tuple<double[], double[]> CalculateResponse(double[] localTotalDisplacementsSuperElement)
        {
			if (this.localDisplacements == null)
			{
				this.localDisplacements = new double[localTotalDisplacementsSuperElement.Length];
			}
			if (localDisplacementsLastConverged == null)
			{
				localDisplacementsLastConverged = new double[localTotalDisplacementsSuperElement.Length];
			}
			double[] localdDisplacementsSuperElement = new double[this.localDisplacements.Length];
			
			for (int i = 0; i < this.localDisplacements.Length; i++)
			{
				localdDisplacementsSuperElement[i] = localTotalDisplacementsSuperElement[i] - localDisplacementsLastConverged[i];
			}
			this.localDisplacements = localTotalDisplacementsSuperElement;
            double[][] Delta = new double[nGaussPoints][];
            double[] localTotalDisplacements = dofEnumerator.GetTransformedDisplacementsVector(localTotalDisplacementsSuperElement);
			double[] localTotaldDisplacements = dofEnumerator.GetTransformedDisplacementsVector(localdDisplacementsSuperElement);
			double[] localTotaldDisplacements_beam = new double[12];
			double[] localTotaldDisplacements_clone = new double[12];

			for (int i1 = 0; i1 < 12; i1++)
			{
				localTotaldDisplacements_beam[i1] = localTotaldDisplacements[12 + i1];
				localTotaldDisplacements_clone[i1] = localTotaldDisplacements[i1];
			}

			supportive_beam.CalculateStresses(localTotaldDisplacements_beam);
			supportive_clone.CalculateStresses(localTotaldDisplacements_clone);
			Delta = this.UpdateCoordinateDataAndCalculateDisplacementVector(localTotalDisplacements);
			
			double[][] tractions = new double[materialsAtGaussPoints.Length][];
            for (int i = 0; i < materialsAtGaussPoints.Length; i++)
            {
				tractions[i] = materialsAtGaussPoints[i].UpdateConstitutiveMatrixAndEvaluateResponse(Delta[i]);
            }
			this.tractions = tractions;
            return new Tuple<double[], double[]>(Delta[materialsAtGaussPoints.Length - 1], tractions[materialsAtGaussPoints.Length - 1]);
        }

        public double[] CalculateResponseIntegral()
        {
			double[] localTotalDisplacementsSuperElement = this.localDisplacements;
            double[] localTotalDisplacements = dofEnumerator.GetTransformedDisplacementsVector(localTotalDisplacementsSuperElement);
            double[] fxk2_coh;
            Tuple<Matrix[], double[]> RtN3AndIntegrationCoeffs;
            RtN3AndIntegrationCoeffs = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations(); // Rt * Nbeam
            Matrix[] RtN3;
            RtN3 = RtN3AndIntegrationCoeffs.Item1;
            double[] integrationCoeffs;
            integrationCoeffs = RtN3AndIntegrationCoeffs.Item2;
            fxk2_coh = this.UpdateForces(/*element, */RtN3, integrationCoeffs, localTotalDisplacements); // sxesh 18 k 19
            return dofEnumerator.GetTransformedForcesVector(fxk2_coh);// embedding
        }

        public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
        {
            return CalculateResponseIntegral();
        }

        public IMatrix StiffnessMatrix(IElementType element)
        {
            double[,] k_stoixeiou_coh2;
            if (MatrixIsNotInitialized)
            {
                this.GetInitialGeometricDataAndInitializeMatrices(element);
                this.UpdateCoordinateDataAndCalculateDisplacementVector(new double[24]); //returns Delta that can't be used for the initial material state
                MatrixIsNotInitialized = false;
            }

            Tuple<Matrix[], double[]> RtN3AndIntegrationCoeffs;
            RtN3AndIntegrationCoeffs = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations(); //Rt*Nbeam
            Matrix[] RtN3;
            RtN3 = RtN3AndIntegrationCoeffs.Item1;
            double[] integrationCoeffs;
            integrationCoeffs = RtN3AndIntegrationCoeffs.Item2;
            k_stoixeiou_coh2 = this.UpdateKmatrices(element, RtN3, integrationCoeffs);
            IMatrix element_stiffnessMatrix = Matrix.CreateFromArray(k_stoixeiou_coh2);
            return dofEnumerator.GetTransformedMatrix(element_stiffnessMatrix);
        }

        //public bool MaterialModified
        //{
        //    get
        //    {
        //        foreach (ICohesiveZoneMaterial3D material in materialsAtGaussPoints)
        //            if (material.Modified) return true;
        //        return false;
        //    }
        //}

        //public void ResetMaterialModified()
        //{
        //    foreach (ICohesiveZoneMaterial3D material in materialsAtGaussPoints) material.ResetModified();
        //}

        //public void ClearMaterialState()
        //{
        //    foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.ClearState();
        //}

        //public void SaveMaterialState()
        //{
        //    foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.SaveState();
        //    supportive_beam.SaveMaterialState();
        //    supportive_clone.SaveMaterialState();
        //}

        //public void ClearMaterialStresses()
        //{
        //    foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.ClearTractions();
        //}

        //public int ID
        //{
        //    get { return 13; }
        //}
        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public IElementDofEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

		//public virtual IList<IList<IDofType>> GetElementDOFTypes(IElementType element)
		//{
		//	return dofTypes;
		//}

		//public double[] CalculateAccelerationForces(IElementType element, IList<MassAccelerationLoad> loads)
		//{
		//    return new double[64];
		//}

		public virtual IMatrix MassMatrix(IElementType element)
        {
            return Matrix.CreateZero(64, 64);
        }

        public virtual IMatrix DampingMatrix(IElementType element)
        {

            return Matrix.CreateZero(64, 64);
        }

        #region EMBEDDED
        private readonly IReadOnlyList<EmbeddedNode> embeddedNodes = new List<EmbeddedNode>();
        public IList<EmbeddedNode> EmbeddedNodes { get { return (IList<EmbeddedNode>)embeddedNodes; } }

        public IElementDofEnumerator DofEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public Dictionary<IDofType, int> GetInternalNodalDOFs(IElementType element, INode node)//
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

        public double[] GetLocalDOFValues(IElementType hostElement, double[] hostDOFValues) // omoiws Beam3D
        {
            //if (transformation == null)
            //    throw new InvalidOperationException("Requested embedded node values for element that has no embedded nodes.");
            //if (hostElementList == null)
            //    throw new InvalidOperationException("Requested host element list for element that has no embedded nodes.");
            //int index = hostElementList.IndexOf(hostElement);
            //if (index < 0)
            //    throw new ArgumentException("Requested host element is not inside host element list.");

            //double[] values = new double[transformation.Columns];
            //int multiplier = hostElement.ElementType.DOFEnumerator.GetDOFTypes(hostElement).SelectMany(d => d).Count();
            //int vectorIndex = 0;
            //for (int i = 0; i < index; i++)
            //    vectorIndex += isNodeEmbedded[i] ? 3 : multiplier;
            //Array.Copy(hostDOFValues, 0, values, vectorIndex, multiplier);

            //return (transformation * new Vector<double>(values)).Data;
            return dofEnumerator.GetTransformedDisplacementsVector(hostDOFValues);
        }

        //public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        //{
        //    throw new NotImplementedException();
        //}

        //IMatrix IElementType.StiffnessMatrix(IElement element)
        //{
        //    throw new NotImplementedException();
        //}

        //public IMatrix MassMatrix(IElement element)
        //{
        //    throw new NotImplementedException();
        //}

        //public IMatrix DampingMatrix(IElement element)
        //{
        //    throw new NotImplementedException();
        //}

        //public virtual IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElementType element) => dofTypes;

        //Dictionary<IDofType, int> IEmbeddedElement.GetInternalNodalDOFs(Element element, Node node)
        //{
        //    throw new NotImplementedException();
        //}

        //double[] IEmbeddedElement.GetLocalDOFValues(Element hostElement, double[] hostDOFValues)
        //{
        //    throw new NotImplementedException();
        //}
        #endregion

        public Matrix CalculateRotationMatrix()
        {
            Matrix R = Matrix.CreateZero(3, 3);
            // R = 0.5*(supportive_beam.currentRotationMatrix + supportive_clone.currentRotationMatrix)
            R[0, 1] = 0.5 * (supportive_beam.currentRotationMatrix[0, 1] + supportive_clone.currentRotationMatrix[0, 1]);
            R[0, 2] = 0.5 * (supportive_beam.currentRotationMatrix[0, 2] + supportive_clone.currentRotationMatrix[0, 2]);
            R[1, 0] = 0.5 * (supportive_beam.currentRotationMatrix[1, 0] + supportive_clone.currentRotationMatrix[1, 0]);
            R[0, 0] = 0.5 * (supportive_beam.currentRotationMatrix[0, 0] + supportive_clone.currentRotationMatrix[0, 0]);
            R[1, 1] = 0.5 * (supportive_beam.currentRotationMatrix[1, 1] + supportive_clone.currentRotationMatrix[1, 1]);
            R[1, 2] = 0.5 * (supportive_beam.currentRotationMatrix[1, 2] + supportive_clone.currentRotationMatrix[1, 2]);
            R[2, 0] = 0.5 * (supportive_beam.currentRotationMatrix[2, 0] + supportive_clone.currentRotationMatrix[2, 0]);
            R[2, 1] = 0.5 * (supportive_beam.currentRotationMatrix[2, 1] + supportive_clone.currentRotationMatrix[2, 1]);
            R[2, 2] = 0.5 * (supportive_beam.currentRotationMatrix[2, 2] + supportive_clone.currentRotationMatrix[2, 2]);
            return R;
		}

		public double[][] Tractions => tractions;
		public IMatrix StiffnessMatrix() => StiffnessMatrix(this);
		public IMatrix MassMatrix() => throw new NotImplementedException();
		public IMatrix DampingMatrix() => throw new NotImplementedException();
		public void SaveConstitutiveLawState()
		{
			localDisplacementsLastConverged = new double[localDisplacements.Length];
			for (int i = 0; i < localDisplacements.Length; i++)
			{
				localDisplacementsLastConverged[i] = localDisplacements[i];
			}	
			foreach (ICohesiveZoneMaterial m in materialsAtGaussPoints) m.CreateState();
			supportive_beam.SaveMaterialState();
			supportive_clone.SaveMaterialState();
		}
		//public Tuple<double[], double[]> CalculateResponse(double[] localDisplacements) => throw new NotImplementedException();
		//public double[] CalculateResponseIntegral() => throw new NotImplementedException();
		//public double[] CalculateResponseIntegralForLogging(double[] localDisplacements) => throw new NotImplementedException();
		public IMatrix PhysicsMatrix()
		{
			return StiffnessMatrix();
		}
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofTypes;

		//public IReadOnlyList<IFiniteElementMaterial> Materials { get; }

		public CellType CellType { get; }
		int IElementType.ID { get; set; }

		public IReadOnlyList<INode> Nodes{ get => nodes; }

		public int SubdomainID { get; set; }

		public IReadOnlyList<Matrix> EvaluateN3ShapeFunctionsReorganized(IQuadrature1D quadrature)
		{
			var cachedN3AtGPs = new Dictionary<IQuadrature1D, IReadOnlyList<Matrix>>();
			bool isCached = cachedN3AtGPs.TryGetValue(quadrature,
				out IReadOnlyList<Matrix> N3AtGPs);
			if (isCached)
			{
				return N3AtGPs;
			}
			else
			{
				IReadOnlyList<Vector> N1 = this.interpolation.EvaluateFunctionsAtGaussPoints(quadrature);
				N3AtGPs = GetN3ShapeFunctionsReorganized(quadrature, N1);
				cachedN3AtGPs.Add(quadrature, N3AtGPs);
				return N3AtGPs;
			}
		}

		private IReadOnlyList<Matrix> GetN3ShapeFunctionsReorganized(IQuadrature1D quadrature, IReadOnlyList<Vector> N1)
		{
			//TODO reorganize cohesive shell  to use only N1 (not reorganised)

			int nGaussPoints = quadrature.IntegrationPoints.Count;
			var N3 = new Matrix[nGaussPoints]; // shapeFunctionsgpData
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				double ksi = quadrature.IntegrationPoints[npoint].Xi;
				var N3gp = Matrix.CreateZero(3, 6); //8=nShapeFunctions;
				for (int l = 0; l < 3; l++)
				{
					for (int m = 0; m < 2; m++) N3gp[l, l + 3 * m] = N1[npoint][m];
				}
				N3[npoint] = N3gp;
			}
			return N3;
		}
	}
}
