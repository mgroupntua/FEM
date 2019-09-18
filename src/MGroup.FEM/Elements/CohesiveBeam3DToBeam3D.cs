//using ISAAR.MSolve.Discretization;
//using ISAAR.MSolve.Discretization.Integration.Quadratures;
//using ISAAR.MSolve.Discretization.Interfaces;
//using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
//using ISAAR.MSolve.FEM.Embedding;
//using ISAAR.MSolve.FEM.Entities;
//using ISAAR.MSolve.FEM.Interfaces;
//using ISAAR.MSolve.FEM.Interpolation;
//using ISAAR.MSolve.Materials.Interfaces;
//using ISAAR.MSolve.Numerical.LinearAlgebra;
//using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
//using System;
//using System.Collections.Generic;
//using IEmbeddedElement = ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;

//namespace ISAAR.MSolve.FEM.Elements
//{
//    class CohesiveBeam3DToBeam3D
//    {
//        protected readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
//        protected readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
//            nodalDOFTypes };
//        protected readonly ICohesiveZoneMaterial3D[] materialsAtGaussPoints;
//        private readonly InterpolationTruss1D interpolation = InterpolationTruss1D.UniqueInstance;
//        protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

//        private int nGaussPoints;

//        /// <summary>
//        /// Initial nodel coordinates of 4 node inner cohesive element
//        /// </summary>
//        private double[][] ox_i; // initialNodalCoordinates

//        /// <summary>
//        /// Unrolled current nodal coordinates of 4 node inner cohesive element
//        /// </summary>
//        private double[] x_local; // currentNodalCoordinates

//        public bool MatrixIsNotInitialized = true;

//        protected CohesiveBeam3DToBeam3D()
//        {
//        }

//        public CohesiveBeam3DToBeam3D(ICohesiveZoneMaterial3D material, IQuadrature1D quadratureForStiffness)
//        {
//            this.QuadratureForStiffness = quadratureForStiffness;
//            this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;
//            materialsAtGaussPoints = new ICohesiveZoneMaterial3D[nGaussPoints];
//            for (int i = 0; i < nGaussPoints; i++) materialsAtGaussPoints[i] = material.Clone();
//        }

//        public IQuadrature1D QuadratureForStiffness { get; }

//        private void GetInitialGeometricDataAndInitializeMatrices(IElement element)
//        {
//            ox_i = new double[4][];

//            for (int j = 0; j < 4; j++)
//            { ox_i[j] = new double[] { element.INodes[j].X, element.INodes[j].Y, element.INodes[j].Z, }; }

//            x_local = new double[12];
//        }

//        private double[][] UpdateCoordinateDataAndCalculateDisplacementVector(double[] localdisplacements)
//        {
//            //IReadOnlyList<Matrix2D> shapeFunctionDerivatives = interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
//            double[,] shapeFunctionDerivatives = interpolation.EvaluateGradientsAt();
            
//            IReadOnlyList<Vector> N1 = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            
//            IReadOnlyList<Matrix2D> N3 = interpolation.EvaluateN3ShapeFunctionsReorganized(QuadratureForStiffness); //Shape functions matrix [N_beam]

//            double[,] u_prok = new double[3, 2];
//            double[,] x_bar = new double[3, 2];

//            double[] e_ksi = new double[3];
//            double e_ksi_norm;
//            double[] e_heta = new double[3];
//            double[] e_1 = new double[3];
//            double[] e_2 = new double[3];
//            double[] e_3 = new double[3];
//            double e_3_norm;
//            double[] u = new double[3];

//            double[] coh_det_J_t = new double[nGaussPoints];

//            double[][] Delta = new double[nGaussPoints][];
//            double[][] c_1 = new double[nGaussPoints][];
//            for (int j = 0; j < nGaussPoints; j++)
//            {
//                Delta[j] = new double[3];
//                c_1[j] = new double[3];
//            }

//            double[][,] R = new double[nGaussPoints][,]; //TODO: maybe cache R
//            for (int j = 0; j < nGaussPoints; j++)
//            {
//                R[j] = new double[3, 3];
//            }


//            // Update x_local
//            for (int j = 0; j < 4; j++)
//            {
//                for (int k = 0; k < 3; k++)
//                {
//                    x_local[3 * j + k] = ox_i[j][k] + localdisplacements[3 * j + k];
//                }
//            }


//            for (int j = 0; j < 2; j++)
//            {
//                for (int k = 0; k < 3; k++)
//                {
//                    u_prok[k, j] = x_local[k + 3 * j] - x_local[6 + k + 3 * j];
//                    x_bar[k, j] = x_local[k + 3 * j] + x_local[6 + k + 3 * j];
//                }
//            }

//            //Calculate Delta for all GPs
//            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
//            {
//                for (int l = 0; l < 3; l++)
//                {
//                    e_ksi[l] = 0;
//                    e_heta[l] = 0;
//                    for (int m = 0; m < 2; m++) // must be 4 in cohesive 8-node
//                    {
//                        e_ksi[l] += shapeFunctionDerivatives[npoint1][0, m] * x_bar[l, m];
//                        e_heta[l] += shapeFunctionDerivatives[npoint1][0, m] * x_bar[l, m];
//                    }
//                    e_ksi[l] = 0.5 * e_ksi[l];
//                    e_heta[l] = 0.5 * e_heta[l];
//                }
//                this.Cross(e_ksi, e_heta, e_3);
//                e_3_norm = Math.Sqrt(e_3[0] * e_3[0] + e_3[1] * e_3[1] + e_3[2] * e_3[2]);
//                e_ksi_norm = Math.Sqrt(e_ksi[0] * e_ksi[0] + e_ksi[1] * e_ksi[1] + e_ksi[2] * e_ksi[2]);
//                for (int l = 0; l < 3; l++)
//                {
//                    e_3[l] = e_3[l] / e_3_norm;
//                    e_1[l] = e_ksi[l] / e_ksi_norm;
//                }
//                this.Cross(e_1, e_3, e_2);
//                for (int l = 0; l < 3; l++)
//                {
//                    R[npoint1][l, 0] = e_1[l];
//                    R[npoint1][l, 1] = e_2[l];
//                    R[npoint1][l, 2] = e_3[l];
//                }
//                for (int l = 0; l < 3; l++)
//                { u[l] = 0; }
//                for (int l = 0; l < 3; l++)
//                {
//                    for (int m = 0; m < 8; m++)  // must be changed for cohesive 8-nodes
//                    {
//                        u[l] += u_prok[l, m] * N1[npoint1][m];
//                    }
//                }
//                for (int l = 0; l < 3; l++)

//                {
//                    for (int m = 0; m < 3; m++)
//                    {
//                        Delta[npoint1][l] += R[npoint1][m, l] * u[m];
//                    }
//                }
//            }
//            return Delta;
//        }

//        private Tuple<Matrix2D[], double[]> CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations()
//        {
//            IReadOnlyList<Vector> N1 = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
//            IReadOnlyList<Matrix2D> N3 = interpolation.EvaluateN3ShapeFunctionsReorganized(QuadratureForStiffness);
//            IReadOnlyList<Matrix2D> shapeFunctionDerivatives = interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

//            double[] integrationsCoeffs = new double[nGaussPoints];
//            Matrix2D[] RtN3 = new Matrix2D[nGaussPoints];
//            double[,] x_bar = new double[3, 8];

//            double[] e_1 = new double[3];
//            double[] e_2 = new double[3];
//            double[] e_3 = new double[3];
//            double e_3_norm;

//            double[] coh_det_J_t = new double[nGaussPoints];

//            double[][] c_1 = new double[nGaussPoints][];
//            for (int j = 0; j < nGaussPoints; j++)
//            {
//                c_1[j] = new double[3];
//            }

//            Matrix2D[] R = new Matrix2D[nGaussPoints]; //TODO: perhaps cache matrices in InitializeMatrices() where RtN3 is calculated
//            for (int j = 0; j < nGaussPoints; j++)
//            { R[j] = new Matrix2D(3, 3); }

//            for (int j = 0; j < 8; j++)
//            {
//                for (int k = 0; k < 3; k++)
//                {
//                    x_bar[k, j] = x_local[k + 3 * j] + x_local[24 + k + 3 * j];
//                }
//            }

//            // Calculate Delta for all GPs
//            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
//            {
//                double[] e_ksi = new double[3];
//                double[] e_heta = new double[3];
//                for (int l = 0; l < 3; l++)
//                {
//                    for (int m = 0; m < 8; m++) // must be 4 in cohesive 8-nodes
//                    {
//                        e_ksi[l] += shapeFunctionDerivatives[npoint1][0, m] * x_bar[l, m];
//                        e_heta[l] += shapeFunctionDerivatives[npoint1][1, m] * x_bar[l, m];
//                    }
//                    e_ksi[l] = 0.5 * e_ksi[l];
//                    e_heta[l] = 0.5 * e_heta[l];
//                }
//                this.Cross(e_ksi, e_heta, e_3);
//                e_3_norm = Math.Sqrt(e_3[0] * e_3[0] + e_3[1] * e_3[1] + e_3[2] * e_3[2]);
//                double e_ksi_norm = Math.Sqrt(e_ksi[0] * e_ksi[0] + e_ksi[1] * e_ksi[1] + e_ksi[2] * e_ksi[2]);
//                for (int l = 0; l < 3; l++)
//                {
//                    e_3[l] = e_3[l] / e_3_norm;
//                    e_1[l] = e_ksi[l] / e_ksi_norm;
//                }
//                this.Cross(e_1, e_3, e_2);
//                for (int l = 0; l < 3; l++)
//                {
//                    R[npoint1][l, 0] = e_1[l];
//                    R[npoint1][l, 1] = e_2[l];
//                    R[npoint1][l, 2] = e_3[l];

//                }

//                this.Cross(e_ksi, e_heta, c_1[npoint1]);
//                coh_det_J_t[npoint1] = Math.Sqrt(c_1[npoint1][0] * c_1[npoint1][0] + c_1[npoint1][1] * c_1[npoint1][1] + c_1[npoint1][2] * c_1[npoint1][2]);
//                integrationsCoeffs[npoint1] = coh_det_J_t[npoint1] * QuadratureForStiffness.IntegrationPoints[npoint1].Weight;

//                // Calculate RtN3 here instead of in InitializeRN3() and then in UpdateForces()
//                RtN3[npoint1] = R[npoint1].Transpose() * N3[npoint1];
//            }
//            return new Tuple<Matrix2D[], double[]>(RtN3, integrationsCoeffs);
//        }

//        private double[] UpdateForces(Element element, Matrix2D[] RtN3, double[] integrationCoeffs)
//        {
//            double[] fxk2_coh = new double[64];
//            double[] fxk1_coh = new double[48]; // TODO: must be 24 in cohesive 8 node

//            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
//            {
//                double[] T_int_integration_coeffs = new double[3];
//                for (int l = 0; l < 3; l++)
//                {
//                    T_int_integration_coeffs[l] = materialsAtGaussPoints[npoint1].Tractions[l] * integrationCoeffs[npoint1];
//                }

//                double[] r_int_1 = new double[24];
//                for (int l = 0; l < 24; l++)
//                {
//                    for (int m = 0; m < 3; m++)
//                    { r_int_1[l] += RtN3[npoint1][m, l] * T_int_integration_coeffs[m]; }
//                }
//                for (int l = 0; l < 24; l++)
//                {
//                    fxk1_coh[l] += r_int_1[l];
//                    fxk1_coh[24 + l] += (-r_int_1[l]);
//                }
//            }

//            fxk2_coh = this.MultiplyForcesForEmbedding(fxk1_coh, element);
//            return fxk2_coh;
//        }

//        private double[,] UpdateKmatrices(IElement element, Matrix2D[] RtN3, double[] integrationCoeffs)
//        {
//            double[,] k_element_coh2 = new double[64, 64];
//            double[,] k_stoixeiou_coh = new double[48, 48];


//            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
//            {
//                Matrix2D D_tan_sunt_ol = new Matrix2D(3, 3);
//                for (int l = 0; l < 3; l++)
//                {
//                    for (int m = 0; m < 3; m++)
//                    {
//                        D_tan_sunt_ol[l, m] = materialsAtGaussPoints[npoint1].ConstitutiveMatrix[l, m] * integrationCoeffs[npoint1];// D_tan[npoint1][l, m] * integrationCoeffs[npoint1];
//                    }
//                }

//                Matrix2D D_RtN3_sunt_ol = D_tan_sunt_ol * RtN3[npoint1];
//                Matrix2D M = RtN3[npoint1].Transpose() * D_RtN3_sunt_ol;

//                for (int l = 0; l < 24; l++)
//                {
//                    for (int m = 0; m < 24; m++)
//                    {
//                        k_stoixeiou_coh[l, m] += M[l, m];
//                        k_stoixeiou_coh[l, 24 + m] += -M[l, m];
//                        k_stoixeiou_coh[24 + l, m] += -M[l, m];
//                        k_stoixeiou_coh[24 + l, 24 + m] += M[l, m];
//                    }
//                }
//            }

//            k_element_coh2 = this.MultiplyStifnessMatrixForEmbedding(k_stoixeiou_coh, element);
//            return k_element_coh2;
//        }

//        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localTotalDisplacementsSuperElement, double[] localdDisplacementsSuperElement)
//        {
//            double[][] Delta = new double[nGaussPoints][];
//            double[] localTotalDisplacements = dofEnumerator.GetTransformedDisplacementsVector(localTotalDisplacementsSuperElement); // embedding
//                                                                                                                                     //double[] localDisplacements = dofEnumerator.GetTransformedVector(localdDisplacementsSuperElement); // embedding

//            Delta = this.UpdateCoordinateDataAndCalculateDisplacementVector(localTotalDisplacements);
//            for (int i = 0; i < materialsAtGaussPoints.Length; i++)
//            {
//                materialsAtGaussPoints[i].UpdateMaterial(Delta[i]);
//            }
//            return new Tuple<double[], double[]>(Delta[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Tractions);
//        }

//        public double[] CalculateForces(Element element, double[] localTotalDisplacementsSuperElement, double[] localdDisplacementsSuperelement)
//        {
//            double[] fxk2_coh = new double[64];
//            Tuple<Matrix2D[], double[]> RtN3AndIntegrationCoeffs;
//            RtN3AndIntegrationCoeffs = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations(); // Rt * Nbeam
//            Matrix2D[] RtN3;
//            RtN3 = RtN3AndIntegrationCoeffs.Item1;
//            double[] integrationCoeffs;
//            integrationCoeffs = RtN3AndIntegrationCoeffs.Item2;

//            fxk2_coh = this.UpdateForces(element, RtN3, integrationCoeffs); // sxesh 18 k 19
//            return dofEnumerator.GetTransformedForcesVector(fxk2_coh);// embedding
//        }

//        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
//        {
//            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
//        }

//        public virtual IMatrix2D StiffnessMatrix(IElement element)
//        {
//            double[,] k_stoixeiou_coh2 = new double[64, 64];
//            if (MatrixIsNotInitialized)
//            {
//                this.GetInitialGeometricDataAndInitializeMatrices(element);
//                this.UpdateCoordinateDataAndCalculateDisplacementVector(new double[64]); //returns Delta that can't be used for the initial material state
//                MatrixIsNotInitialized = false;
//            }

//            Tuple<Matrix2D[], double[]> RtN3AndIntegrationCoeffs;
//            RtN3AndIntegrationCoeffs = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations();//Rt *Nbeam
//            Matrix2D[] RtN3;
//            RtN3 = RtN3AndIntegrationCoeffs.Item1;
//            double[] integrationCoeffs;
//            integrationCoeffs = RtN3AndIntegrationCoeffs.Item2;

//            k_stoixeiou_coh2 = this.UpdateKmatrices(element, RtN3, integrationCoeffs); //sxesh 8
//            IMatrix2D element_stiffnessMatrix = new Matrix2D(k_stoixeiou_coh2);
//            return dofEnumerator.GetTransformedMatrix(element_stiffnessMatrix); // embedding
//        }


//        public bool MaterialModified
//        {
//            get
//            {
//                foreach (ICohesiveZoneMaterial3D material in materialsAtGaussPoints)
//                    if (material.Modified) return true;
//                return false;
//            }
//        }

//        public void ResetMaterialModified()
//        {
//            foreach (ICohesiveZoneMaterial3D material in materialsAtGaussPoints) material.ResetModified();
//        }

//        public void ClearMaterialState()
//        {
//            foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.ClearState();
//        }

//        public void SaveMaterialState()
//        {
//            foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.SaveState();
//        }

//        public void ClearMaterialStresses()
//        {
//            foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.ClearTractions();
//        }

//        public int ID
//        {
//            get { return 13; }
//        }
//        public ElementDimensions ElementDimensions
//        {
//            get { return ElementDimensions.ThreeD; }
//        }

//        public IElementDOFEnumerator DOFEnumerator
//        {
//            get { return dofEnumerator; }
//            set { dofEnumerator = value; }
//        }

//        public virtual IList<IList<DOFType>> GetElementDOFTypes(IElement element)
//        {
//            return dofTypes;
//        }

//        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
//        {
//            return new double[64];
//        }

//        public virtual IMatrix2D MassMatrix(IElement element)
//        {
//            return new Matrix2D(64, 64);
//        }

//        public virtual IMatrix2D DampingMatrix(IElement element)
//        {

//            return new Matrix2D(64, 64);
//        }

//        #region EMBEDDED
//        private readonly List<EmbeddedNode> embeddedNodes = new List<EmbeddedNode>();
//        public IList<EmbeddedNode> EmbeddedNodes { get { return embeddedNodes; } }

//        public Dictionary<DOFType, int> GetInternalNodalDOFs(Element element, Node node)//
//        {
//            int index = 0;
//            foreach (var elementNode in element.Nodes)
//            {
//                if (node.ID == elementNode.ID)
//                    break;
//                index++;
//            }
//            if (index >= 16)
//                throw new ArgumentException(String.Format("GetInternalNodalDOFs: Node {0} not found in element {1}.", node.ID, element.ID));

//            if (index >= 8)
//            {
//                int index2 = index - 8;
//                return new Dictionary<DOFType, int>() { { DOFType.X, 39 + 3 * index2 + 1 }, { DOFType.Y, 39 + 3 * index2 + 2 }, { DOFType.Z, 39 + 3 * index2 + 3 } };
//            }
//            else
//            {
//                return new Dictionary<DOFType, int>() { { DOFType.X, + 5 * index + 0 }, { DOFType.Y, + 5 * index + 1 }, { DOFType.Z, + 5 * index + 2 },
//                                                        { DOFType.RotX, + 5 * index + 3 }, { DOFType.RotY, + 5 * index + 4 }};
//            }
//        }

//        public double[] GetLocalDOFValues(Element hostElement, double[] hostDOFValues) // omoiws Beam3D
//        {
//            //if (transformation == null)
//            //    throw new InvalidOperationException("Requested embedded node values for element that has no embedded nodes.");
//            //if (hostElementList == null)
//            //    throw new InvalidOperationException("Requested host element list for element that has no embedded nodes.");
//            //int index = hostElementList.IndexOf(hostElement);
//            //if (index < 0)
//            //    throw new ArgumentException("Requested host element is not inside host element list.");

//            //double[] values = new double[transformation.Columns];
//            //int multiplier = hostElement.ElementType.DOFEnumerator.GetDOFTypes(hostElement).SelectMany(d => d).Count();
//            //int vectorIndex = 0;
//            //for (int i = 0; i < index; i++)
//            //    vectorIndex += isNodeEmbedded[i] ? 3 : multiplier;
//            //Array.Copy(hostDOFValues, 0, values, vectorIndex, multiplier);

//            //return (transformation * new Vector<double>(values)).Data;

//            return dofEnumerator.GetTransformedDisplacementsVector(hostDOFValues);
//        }
//        #endregion

//    }
//}
