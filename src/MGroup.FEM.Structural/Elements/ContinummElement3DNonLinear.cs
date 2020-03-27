using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;
using MGroup.FEM.Interpolation;
using MGroup.FEM.Interpolation.GaussPointExtrapolation;
using MGroup.FEM.Interpolation.Jacobians;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;
using MGroup.MSolve.Geometry.Coordinates;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.FEM.Embedding;
using System.Linq;


namespace ISAAR.MSolve.FEM.Elements
{
    /// <summary>
    /// Continuum finite Element for 3d problems with material and geometric nonlinearities
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class ContinummElement3DNonLinear : IStructuralFiniteElement, IEmbeddedHostElement
    {
        protected readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
        protected readonly IDofType[][] dofTypes;
        protected readonly IContinuumMaterial3D[] materialsAtGaussPoints;
        protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
        
        private readonly int nGaussPoints;
        private bool isInitialized = false;

        private double[][] ox_i; //not defined by user. 8 arrays of 3 elements
        private double[][] tu_i;
        private double[] integrationCoeffs;

        private double[][] GLvec;
        private double[][] GLvec_last_converged;

        protected ContinummElement3DNonLinear()
        {
        }

        public ContinummElement3DNonLinear(IReadOnlyList<Node> nodes,IContinuumMaterial3D material, IQuadrature3D quadratureForStiffness,
             IIsoparametricInterpolation3D interpolation)
        {
            this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;
            this.QuadratureForStiffness = quadratureForStiffness;
            this.Interpolation = interpolation;

            
            materialsAtGaussPoints = new IContinuumMaterial3D[nGaussPoints];
            for (int i = 0; i < nGaussPoints; i++)
                materialsAtGaussPoints[i] = (IContinuumMaterial3D)material.Clone();

            dofTypes = new IDofType[nodes.Count][];
            for (int i = 0; i < nodes.Count; i++)
            {
                dofTypes[i] = new IDofType[]
                {
                    StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
                };
            }
        }

        public IIsoparametricInterpolation3D Interpolation { get; }
        public IQuadrature3D QuadratureForStiffness { get; }

        public int ID => 13;
        public CellType CellType => Interpolation.CellType;

        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

        public IElementDofEnumerator DofEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public IReadOnlyList<IFiniteElementMaterial> Materials => materialsAtGaussPoints;

        public bool MaterialModified
        {
            get
            {
                foreach (IContinuumMaterial3D material in materialsAtGaussPoints)
                    if (material.Modified) return true;
                return false;
            }
        }

        private Matrix[] GetBL13Hexa(IReadOnlyList<Matrix> ll1_hexa)
        {
            Matrix[] BL13_hexa;
            BL13_hexa = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                BL13_hexa[npoint] = Matrix.CreateZero(9, 3*ll1_hexa[npoint].NumRows);
                for (int m = 0; m < ll1_hexa[npoint].NumRows; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL13_hexa[npoint][n, 3 * m + 0] = ll1_hexa[npoint][ m,n];
                        BL13_hexa[npoint][n + 3, 3 * m + 1] = ll1_hexa[npoint][m,n];
                        BL13_hexa[npoint][n + 6, 3 * m + 2] = ll1_hexa[npoint][m, n];
                    }
                }
            }
            return BL13_hexa;
        }

        private Matrix[] GetBL11a_hexa(Matrix[] J_0inv_hexa)
        {
            Matrix[] BL11a_hexa = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                BL11a_hexa[gpoint] = Matrix.CreateZero(6, 9);
                for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL11a_hexa[gpoint][m, 3 * m + n] = J_0inv_hexa[gpoint][m, n];
                    }
                }
                for (int n = 0; n < 3; n++)
                {
                    BL11a_hexa[gpoint][3, n] = J_0inv_hexa[gpoint][1, n]; // calculate 4th data line
                    BL11a_hexa[gpoint][3, 3 + n] = J_0inv_hexa[gpoint][0, n];
                    BL11a_hexa[gpoint][4, 3 + n] = J_0inv_hexa[gpoint][2, n]; // calculate 5th data line
                    BL11a_hexa[gpoint][4, 6 + n] = J_0inv_hexa[gpoint][1, n];
                    BL11a_hexa[gpoint][5, 0 + n] = J_0inv_hexa[gpoint][2, n]; // calculate 6th data line
                    BL11a_hexa[gpoint][5, 6 + n] = J_0inv_hexa[gpoint][0, n];
                }
            }

            return BL11a_hexa;
        }

        private Matrix[] GetBL12_hexa(Matrix[] J_0inv_hexa)
        {
            Matrix[] BL12_hexa = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                BL12_hexa[gpoint] = Matrix.CreateZero(9, 9);
                for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL12_hexa[gpoint][m, 3 * m + n] = J_0inv_hexa[gpoint][0, n];
                    }
                }
                for (int m = 0; m < 3; m++) // calculate  data lines 4:6
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL12_hexa[gpoint][3 + m, 3 * m + n] = J_0inv_hexa[gpoint][1, n];
                    }
                }
                for (int m = 0; m < 3; m++) // calculate  data lines 7:8
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL12_hexa[gpoint][6 + m, 3 * m + n] = J_0inv_hexa[gpoint][2, n];
                    }
                }

            }

            return BL12_hexa;
        }

        private Matrix[] GetBL01_hexa(Matrix[] J_0inv_hexa)
        {
            Matrix[] BL01_hexa = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                BL01_hexa[gpoint] = Matrix.CreateZero(6, 9);
                for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL01_hexa[gpoint][m, 3 * m + n] = J_0inv_hexa[gpoint][m, n];
                    }
                }
                for (int n = 0; n < 3; n++)
                {
                    BL01_hexa[gpoint][3, n] = J_0inv_hexa[gpoint][1, n]; // calculate 4th data line
                    BL01_hexa[gpoint][3, 3 + n] = J_0inv_hexa[gpoint][0, n];
                    BL01_hexa[gpoint][4, 3 + n] = J_0inv_hexa[gpoint][2, n]; // calculate 5th data line
                    BL01_hexa[gpoint][4, 6 + n] = J_0inv_hexa[gpoint][1, n];
                    BL01_hexa[gpoint][5, 0 + n] = J_0inv_hexa[gpoint][2, n]; // calculate 6th data line
                    BL01_hexa[gpoint][5, 6 + n] = J_0inv_hexa[gpoint][0, n];
                }
            }
            return BL01_hexa;
        }

        private Matrix[] GetBNL1_hexa(Matrix[] J_0inv_hexa)
        {
            Matrix[] BNL1_hexa = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                BNL1_hexa[gpoint] = Matrix.CreateZero(9, 9);
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            BNL1_hexa[gpoint][3 * m + n, 3 * m + p] = J_0inv_hexa[gpoint][n, p];
                        }
                    }
                }
            }
            return BNL1_hexa;
        }

        private void CalculateInitialConfigurationData(IElement element)
        {
            int numNodes = element.Nodes.Count;
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            Matrix[] BL13_hexa;
            BL13_hexa = GetBL13Hexa(shapeFunctionNaturalDerivatives);
            

            Matrix[] BNL1_hexa;

            ox_i = new double[numNodes][];
            tu_i = new double[numNodes][];

            var jacobians = shapeFunctionNaturalDerivatives.Select(x=> new IsoparametricJacobian3D(element.Nodes,x));

            Matrix[] J_0inv_hexa = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] detJ_0 = jacobians.Select(x => x.DirectDeterminant).ToArray();

            



            integrationCoeffs = new double[nGaussPoints];

            BNL1_hexa = GetBNL1_hexa(J_0inv_hexa);

            for (int j = 0; j < numNodes; j++)
            {
                ox_i[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
                }

            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                integrationCoeffs[gpoint] = detJ_0[gpoint] * QuadratureForStiffness.IntegrationPoints[gpoint].Weight;

            }

            tu_i = new double[numNodes][];
            GLvec = new double[nGaussPoints][];
            GLvec_last_converged = new double[nGaussPoints][];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                GLvec[gpoint] = new double[6];
                GLvec_last_converged[gpoint] = new double[6];
            }
            for (int k = 0; k < numNodes; k++)
            {
                tu_i[k] = new double[3];
            }
            isInitialized = true;

        }

        private void UpdateCoordinateData(double[] localdisplacements, out double[][] tx_i)
        {
            int numNodes = localdisplacements.Length / 3;
            tx_i = new double[numNodes][];
            for (int j = 0; j < numNodes; j++)
            {
                tx_i[j] = new double[3];
                for (int k = 0; k < 3; k++)
                {
                    tu_i[j][k] = localdisplacements[3 * j + k];
                    tx_i[j][k] = ox_i[j][k] + tu_i[j][k];
                }
            }
        }

        private void CalculateStrains(double[] localdisplacements, IElement element, double[][] tx_i) 
        {
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));
            Matrix[] J_0inv_hexa = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] detJ_0 = jacobians.Select(x => x.DirectDeterminant).ToArray();
            //TODO: possibility of caching ll1_hexa or J_0inv

            Matrix[] DGtr = new Matrix[nGaussPoints];
            Matrix[] GL = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                DGtr[npoint] = Matrix.CreateZero(3, 3);
                GL[npoint] = Matrix.CreateZero(3, 3);
            }

            var jacobiansDeformed = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(tx_i, x,false)).ToArray();
            Matrix[] J_1 = jacobiansDeformed.Select(x => x.DirectMatrix).ToArray();


            //
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                //
                DGtr[npoint] = J_0inv_hexa[npoint] * J_1[npoint];

                //
                GL[npoint] = DGtr[npoint] * DGtr[npoint].Transpose();
                for (int m = 0; m < 3; m++)
                {
                    GL[npoint][m, m] += -1;
                }
                GL[npoint].ScaleIntoThis(0.5);
                
                //
                for (int m = 0; m < 3; m++)
                {
                    GLvec[npoint][m] = GL[npoint][m, m];
                }
                GLvec[npoint][3] = 2 * GL[npoint][0, 1];
                GLvec[npoint][4] = 2 * GL[npoint][1, 2];
                GLvec[npoint][5] = 2 * GL[npoint][2, 0];
            }

        }

        private double[] UpdateForces(IElement element)
        {
            //TODO: the gauss point loop should be the outer one


            // Matrices that are not currently cached are calculated here.
            int numNodes = element.Nodes.Count();
            Matrix ll2 = Matrix.CreateZero(numNodes, 3);
            for (int m = 0; m < numNodes; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    ll2[m, n] = tu_i[m][n];
                }
            }
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));
            Matrix[] J_0inv_hexa = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] detJ_0 = jacobians.Select(x => x.DirectDeterminant).ToArray();

            Matrix[] BL13_hexa;
            BL13_hexa = GetBL13Hexa(shapeFunctionNaturalDerivatives);
            Matrix[] BL11a_hexa; // dimension number of gpoints
            Matrix[] BL12_hexa;
            Matrix[] BL01_hexa;
            BL11a_hexa = GetBL11a_hexa(J_0inv_hexa);
            BL12_hexa = GetBL12_hexa(J_0inv_hexa);
            BL01_hexa = GetBL01_hexa(J_0inv_hexa);

            //INITIALIZATION of MAtrixes that are currently not cached
            double[][] integrCoeff_Spkvec = new double[nGaussPoints][];
            Matrix[] BL = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                integrCoeff_Spkvec[gpoint] = new double[6];
                BL[gpoint] = Matrix.CreateZero(6, 3*numNodes);
            }

            double[][] fxk1 = new double[nGaussPoints + 1][];
            for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
            {
                fxk1[npoint] = new double[3*numNodes];
            }

            Matrix[] BL11 = new Matrix[nGaussPoints];
            Matrix[] BL1112sun01_hexa = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                BL11[npoint] = Matrix.CreateZero(6, 9);
                BL1112sun01_hexa[npoint] = Matrix.CreateZero(6, 9);
            }

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {

                integrCoeff_Spkvec[npoint] = materialsAtGaussPoints[npoint].Stresses.Scale(integrationCoeffs[npoint]);

                //
                Matrix l_cyrcumflex;//= Matrix.CreateZero(3, 3);
                l_cyrcumflex = shapeFunctionNaturalDerivatives[npoint].Transpose() * ll2;

                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            BL11[npoint][m, n] += BL11a_hexa[npoint][m, p] * l_cyrcumflex[p, n];
                            BL11[npoint][m, 3 + n] += BL11a_hexa[npoint][m, 3 + p] * l_cyrcumflex[p, n];
                            BL11[npoint][m, 6 + n] += BL11a_hexa[npoint][m, 6 + p] * l_cyrcumflex[p, n];
                        }
                    }
                }

                //
                BL1112sun01_hexa[npoint] = BL11[npoint] * BL12_hexa[npoint];
                BL1112sun01_hexa[npoint].AddIntoThis(BL01_hexa[npoint]);
                
                // 
                BL[npoint] = BL1112sun01_hexa[npoint] * BL13_hexa[npoint];
                
                //              
                fxk1[npoint] = BL[npoint].Multiply(integrCoeff_Spkvec[npoint], true);
            }

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                fxk1[nGaussPoints].AddIntoThis(fxk1[npoint]);                
            }

            return fxk1[nGaussPoints];
        }

        private Matrix UpdateKmatrices(IElement element)
        {
            int numNodes = element.Nodes.Count();
            Matrix k_element = Matrix.CreateZero(3*numNodes, 3*numNodes);


            // initialization of matrices that are not cached currently
            double[][] integrCoeff_Spkvec = new double[nGaussPoints][];
            Matrix[] BL = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                integrCoeff_Spkvec[gpoint] = new double[6];
                BL[gpoint] = Matrix.CreateZero(6, 3*numNodes);

            }
            Matrix ll2 = Matrix.CreateZero(numNodes, 3);
            for (int m = 0; m < numNodes; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    ll2[m, n] = tu_i[m][n];
                }
            }
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));
            Matrix[] J_0inv_hexa = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] detJ_0 = jacobians.Select(x => x.DirectDeterminant).ToArray();

            Matrix[] BL13_hexa;
            BL13_hexa = GetBL13Hexa(shapeFunctionNaturalDerivatives);
            Matrix[] BL11a_hexa; // dimension: gpoints
            Matrix[] BL12_hexa;
            Matrix[] BL01_hexa;
            BL11a_hexa = GetBL11a_hexa(J_0inv_hexa);
            BL12_hexa = GetBL12_hexa(J_0inv_hexa);
            BL01_hexa = GetBL01_hexa(J_0inv_hexa);

            Matrix[] BL11 = new Matrix[nGaussPoints];
            Matrix[] BL1112sun01_hexa = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                BL11[npoint] = Matrix.CreateZero(6, 9);
                BL1112sun01_hexa[npoint] = Matrix.CreateZero(6, 9); //TODO this may be unnescessary
            }



            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {

                // 
                integrCoeff_Spkvec[npoint] = materialsAtGaussPoints[npoint].Stresses.Scale(integrationCoeffs[npoint]);
                
                //
                Matrix l_cyrcumflex = Matrix.CreateZero(3, 3);
                l_cyrcumflex = shapeFunctionNaturalDerivatives[npoint].Transpose() * ll2;
                 
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            BL11[npoint][m, n] += BL11a_hexa[npoint][m, p] * l_cyrcumflex[p, n];
                            BL11[npoint][m, 3 + n] += BL11a_hexa[npoint][m, 3 + p] * l_cyrcumflex[p, n];
                            BL11[npoint][m, 6 + n] += BL11a_hexa[npoint][m, 6 + p] * l_cyrcumflex[p, n];
                        }
                    }
                }

                // 
                BL1112sun01_hexa[npoint] = BL11[npoint] * BL12_hexa[npoint];
                BL1112sun01_hexa[npoint].AddIntoThis(BL01_hexa[npoint]);
                
                //
                BL[npoint] = BL1112sun01_hexa[npoint] * BL13_hexa[npoint];
                
            }
            // TODO: BL and above calculations can cached from calculate forces method

            Matrix[] BNL1_hexa;
            Matrix[] BNL_hexa;
            BNL1_hexa = GetBNL1_hexa(J_0inv_hexa);
            BNL_hexa = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
               //BNL_hexa[gpoint] = Matrix.CreateZero(9, 3*numNodes); //todo this may be unnescessary

                BNL_hexa[gpoint] = BNL1_hexa[gpoint] * BL13_hexa[gpoint];
                
            }


            Matrix[] integrCoeff_Spk = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                integrCoeff_Spk[npoint] = Matrix.CreateZero(3, 3);
            }

            Matrix[] kl_ = new Matrix[nGaussPoints + 1];
            Matrix[] knl_ = new Matrix[nGaussPoints + 1];
            for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
            {
                kl_[npoint] = Matrix.CreateZero(3*numNodes, 3*numNodes);
                knl_[npoint] = Matrix.CreateZero(3 * numNodes, 3 * numNodes);
            }



            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                Matrix integrCoeff_SPK_epi_BNL_hexa = Matrix.CreateZero(9, 3*numNodes); //TODO
                Matrix integrCoeff_cons_disp = Matrix.CreateZero(6, 6); //TODO
                Matrix integrCoeff_cons_disp_epi_BL = Matrix.CreateZero(6, 3*numNodes);//TODO

                //
                integrCoeff_Spk[npoint][0, 0] = integrCoeff_Spkvec[npoint][0];
                integrCoeff_Spk[npoint][0, 1] = integrCoeff_Spkvec[npoint][3];
                integrCoeff_Spk[npoint][0, 2] = integrCoeff_Spkvec[npoint][5];
                integrCoeff_Spk[npoint][1, 0] = integrCoeff_Spkvec[npoint][3];
                integrCoeff_Spk[npoint][1, 1] = integrCoeff_Spkvec[npoint][1];
                integrCoeff_Spk[npoint][1, 2] = integrCoeff_Spkvec[npoint][4];
                integrCoeff_Spk[npoint][2, 0] = integrCoeff_Spkvec[npoint][5];
                integrCoeff_Spk[npoint][2, 1] = integrCoeff_Spkvec[npoint][4];
                integrCoeff_Spk[npoint][2, 2] = integrCoeff_Spkvec[npoint][2];

                //
                IMatrixView consDisp = materialsAtGaussPoints[npoint].ConstitutiveMatrix;

                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 6; n++)
                    {
                        integrCoeff_cons_disp[m, n] = integrationCoeffs[npoint] * consDisp[m, n];
                    }
                }

                //
                integrCoeff_cons_disp_epi_BL = integrCoeff_cons_disp * BL[npoint];
                
                //
                kl_[npoint] = BL[npoint].Transpose() * integrCoeff_cons_disp_epi_BL;                

                //
                for (int m = 0; m < 3; m++) // 3x24 dimensions
                {
                    for (int n = 0; n < 3*numNodes; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            integrCoeff_SPK_epi_BNL_hexa[m, n] += integrCoeff_Spk[npoint][m, p] * BNL_hexa[npoint][p, n];
                            integrCoeff_SPK_epi_BNL_hexa[3 + m, n] += integrCoeff_Spk[npoint][m, p] * BNL_hexa[npoint][3 + p, n];
                            integrCoeff_SPK_epi_BNL_hexa[6 + m, n] += integrCoeff_Spk[npoint][m, p] * BNL_hexa[npoint][6 + p, n];
                        }
                    }
                }

                //
                knl_[npoint] = BNL_hexa[npoint].Transpose() * integrCoeff_SPK_epi_BNL_hexa;                
            }

            // Add contributions of each gp on the total element stiffness matrix k_element            
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                for (int m = 0; m < 3*numNodes; m++)
                {
                    for (int n = 0; n < 3*numNodes; n++)
                    {
                        kl_[nGaussPoints][m, n] += kl_[npoint][m, n];
                        knl_[nGaussPoints][m, n] += knl_[npoint][m, n];
                    }
                }
            }
            for (int m = 0; m < 3 * numNodes; m++)
            {
                for (int n = 0; n < 3 * numNodes; n++)
                {
                    k_element[m, n] = kl_[nGaussPoints][m, n] + knl_[nGaussPoints][m, n];
                }
            }

            return k_element;
        }
        
        public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            this.UpdateCoordinateData(localTotalDisplacements, out double[][] tx_i);
            this.CalculateStrains(localTotalDisplacements, element, tx_i);
            double[] GLvec_strain_minus_last_converged_value = new double[6];
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                GLvec_strain_minus_last_converged_value = new double[6] 
                {
                    GLvec[npoint][0]- GLvec_last_converged[npoint][0],
                    GLvec[npoint][1] - GLvec_last_converged[npoint][1],
                    GLvec[npoint][2] - GLvec_last_converged[npoint][2],
                    GLvec[npoint][3]- GLvec_last_converged[npoint][3],
                    GLvec[npoint][4]- GLvec_last_converged[npoint][4],
                    GLvec[npoint][5]- GLvec_last_converged[npoint][5]
                };
                materialsAtGaussPoints[npoint].UpdateMaterial(GLvec_strain_minus_last_converged_value); 
                //To update with total strain simplY = materialsAtGaussPoints[npoint].UpdateMaterial(GLvec[npoint]);
            }
            return new Tuple<double[], double[]>(GLvec_strain_minus_last_converged_value, 
                materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);
            //TODO return data with total strains data would be:
            //return new Tuple<double[], double[]>(GLvec[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);
            //TODO: why return only the strain- stress of the gausspoint that is last on the array, Where is it needed?
        }

        public double[] CalculateForces(IElement element, double[] localTotalDisplacements, double[] localdDisplacements) 
            => this.UpdateForces(element);

        public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements) 
            => CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

        public virtual IMatrix StiffnessMatrix(IElement element)
        {
            if (!isInitialized)
            {
                int numNodes = element.Nodes.Count();
                this.CalculateInitialConfigurationData(element);
                var localTotalDisplacements = new double[3*numNodes];
                this.UpdateCoordinateData(localTotalDisplacements, out double[][] tx_i);
                this.CalculateStrains(localTotalDisplacements, element, tx_i);
            }
            Matrix elementStiffness = this.UpdateKmatrices(element);
            //It doesn't implement Iembedded to return dof.Enumerator.GetTransformedMatrix
            return elementStiffness;
        }

        public void ResetMaterialModified()
        {
            foreach (IContinuumMaterial3D material in materialsAtGaussPoints) material.ResetModified();
        }

        public void ClearMaterialState()
        {
            //TODO: the next throws an exception. Investigate. Possible changes in Analyzers may be the cause.
            //foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearState();
        }

        public void SaveMaterialState()
        {
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                for (int i1 = 0; i1 < 6; i1++)
                { GLvec_last_converged[npoint][i1] = GLvec[npoint][i1]; }
            }

            foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.SaveState();
        }

        public void ClearMaterialStresses()
        {
            foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
        }

        public virtual IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;

        #region not implemented
        public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
        {
            throw new NotImplementedException();
        }

        public virtual IMatrix MassMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public virtual IMatrix DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }
        #endregion


        #region IEmbeddedHostElement
        protected double[,] GetCoordinatesTranspose(IElement element)
        {
            double[,] faXYZ = new double[3, dofTypes.Length];
            for (int i = 0; i < dofTypes.Length; i++)
            {
                faXYZ[0, i] = element.Nodes[i].X;
                faXYZ[1, i] = element.Nodes[i].Y;
                faXYZ[2, i] = element.Nodes[i].Z;
            }
            return faXYZ;
        }

        //TODO: This should be handled by InterpolationHexa8Reverse
        private double[] CalcH8Shape(double fXi, double fEta, double fZeta)
        {
            const double fSqC125 = 0.5;
            double fXiP = (1.0 + fXi) * fSqC125;
            double fEtaP = (1.0 + fEta) * fSqC125;
            double fZetaP = (1.0 + fZeta) * fSqC125;
            double fXiM = (1.0 - fXi) * fSqC125;
            double fEtaM = (1.0 - fEta) * fSqC125;
            double fZetaM = (1.0 - fZeta) * fSqC125;

            double[] H8Shape = new double[8]; // Warning: shape function data not in hexa8fixed order.

            H8Shape[6] = fXiM * fEtaM * fZetaM;
            H8Shape[7] = fXiP * fEtaM * fZetaM;
            H8Shape[4] = fXiP * fEtaP * fZetaM;
            H8Shape[5] = fXiM * fEtaP * fZetaM;
            H8Shape[2] = fXiM * fEtaM * fZetaP;
            H8Shape[3] = fXiP * fEtaM * fZetaP;
            H8Shape[0] = fXiP * fEtaP * fZetaP;
            H8Shape[1] = fXiM * fEtaP * fZetaP;

            return H8Shape;
        }

        //TODO: This should be handled by InterpolationHexa8Reverse
        private double[] CalcH8NablaShape(double fXi, double fEta, double fZeta)
        {
            const double fSq125 = 0.35355339059327376220042218105242;
            double fXiP = (1.0 + fXi) * fSq125;
            double fEtaP = (1.0 + fEta) * fSq125;
            double fZetaP = (1.0 + fZeta) * fSq125;
            double fXiM = (1.0 - fXi) * fSq125;
            double fEtaM = (1.0 - fEta) * fSq125;
            double fZetaM = (1.0 - fZeta) * fSq125;

            double[] faDS = new double[24];
            faDS[6] = -fEtaM * fZetaM;
            faDS[4] = fEtaP * fZetaM;
            faDS[2] = -fEtaM * fZetaP;
            faDS[0] = fEtaP * fZetaP;
            faDS[7] = -faDS[6];
            faDS[5] = -faDS[4];
            faDS[3] = -faDS[2];
            faDS[1] = -faDS[0];


            faDS[14] = -fXiM * fZetaM;
            faDS[15] = -fXiP * fZetaM;
            faDS[10] = -fXiM * fZetaP;
            faDS[11] = -fXiP * fZetaP;
            faDS[12] = -faDS[15];
            faDS[13] = -faDS[14];
            faDS[8] = -faDS[11];
            faDS[9] = -faDS[10];


            faDS[22] = -fXiM * fEtaM;
            faDS[23] = -fXiP * fEtaM;
            faDS[20] = -fXiP * fEtaP;
            faDS[21] = -fXiM * fEtaP;
            faDS[18] = -faDS[22];
            faDS[19] = -faDS[23];
            faDS[16] = -faDS[20];
            faDS[17] = -faDS[21];

            return faDS;
        }

        protected static double determinantTolerance = 0.00000001;
        //TODO: This should be handled by JacobianHexa8Reverse
        private Tuple<double[,], double[,], double> CalcH8JDetJ(double[,] faXYZ, double[] faDS)
        {
            double[,] faJ = new double[3, 3];
            faJ[0, 0] = faDS[0] * faXYZ[0, 0] + faDS[1] * faXYZ[0, 1] + faDS[2] * faXYZ[0, 2] + faDS[3] * faXYZ[0, 3] + faDS[4] * faXYZ[0, 4] + faDS[5] * faXYZ[0, 5] + faDS[6] * faXYZ[0, 6] + faDS[7] * faXYZ[0, 7];
            faJ[0, 1] = faDS[0] * faXYZ[1, 0] + faDS[1] * faXYZ[1, 1] + faDS[2] * faXYZ[1, 2] + faDS[3] * faXYZ[1, 3] + faDS[4] * faXYZ[1, 4] + faDS[5] * faXYZ[1, 5] + faDS[6] * faXYZ[1, 6] + faDS[7] * faXYZ[1, 7];
            faJ[0, 2] = faDS[0] * faXYZ[2, 0] + faDS[1] * faXYZ[2, 1] + faDS[2] * faXYZ[2, 2] + faDS[3] * faXYZ[2, 3] + faDS[4] * faXYZ[2, 4] + faDS[5] * faXYZ[2, 5] + faDS[6] * faXYZ[2, 6] + faDS[7] * faXYZ[2, 7];
            faJ[1, 0] = faDS[8] * faXYZ[0, 0] + faDS[9] * faXYZ[0, 1] + faDS[10] * faXYZ[0, 2] + faDS[11] * faXYZ[0, 3] + faDS[12] * faXYZ[0, 4] + faDS[13] * faXYZ[0, 5] + faDS[14] * faXYZ[0, 6] + faDS[15] * faXYZ[0, 7];
            faJ[1, 1] = faDS[8] * faXYZ[1, 0] + faDS[9] * faXYZ[1, 1] + faDS[10] * faXYZ[1, 2] + faDS[11] * faXYZ[1, 3] + faDS[12] * faXYZ[1, 4] + faDS[13] * faXYZ[1, 5] + faDS[14] * faXYZ[1, 6] + faDS[15] * faXYZ[1, 7];
            faJ[1, 2] = faDS[8] * faXYZ[2, 0] + faDS[9] * faXYZ[2, 1] + faDS[10] * faXYZ[2, 2] + faDS[11] * faXYZ[2, 3] + faDS[12] * faXYZ[2, 4] + faDS[13] * faXYZ[2, 5] + faDS[14] * faXYZ[2, 6] + faDS[15] * faXYZ[2, 7];
            faJ[2, 0] = faDS[16] * faXYZ[0, 0] + faDS[17] * faXYZ[0, 1] + faDS[18] * faXYZ[0, 2] + faDS[19] * faXYZ[0, 3] + faDS[20] * faXYZ[0, 4] + faDS[21] * faXYZ[0, 5] + faDS[22] * faXYZ[0, 6] + faDS[23] * faXYZ[0, 7];
            faJ[2, 1] = faDS[16] * faXYZ[1, 0] + faDS[17] * faXYZ[1, 1] + faDS[18] * faXYZ[1, 2] + faDS[19] * faXYZ[1, 3] + faDS[20] * faXYZ[1, 4] + faDS[21] * faXYZ[1, 5] + faDS[22] * faXYZ[1, 6] + faDS[23] * faXYZ[1, 7];
            faJ[2, 2] = faDS[16] * faXYZ[2, 0] + faDS[17] * faXYZ[2, 1] + faDS[18] * faXYZ[2, 2] + faDS[19] * faXYZ[2, 3] + faDS[20] * faXYZ[2, 4] + faDS[21] * faXYZ[2, 5] + faDS[22] * faXYZ[2, 6] + faDS[23] * faXYZ[2, 7];

            double fDet1 = faJ[0, 0] * (faJ[1, 1] * faJ[2, 2] - faJ[2, 1] * faJ[1, 2]);
            double fDet2 = -faJ[0, 1] * (faJ[1, 0] * faJ[2, 2] - faJ[2, 0] * faJ[1, 2]);
            double fDet3 = faJ[0, 2] * (faJ[1, 0] * faJ[2, 1] - faJ[2, 0] * faJ[1, 1]);
            double fDetJ = fDet1 + fDet2 + fDet3;
            if (fDetJ < determinantTolerance)
            {
                throw new ArgumentException(
                    $"Jacobian determinant is negative or under tolerance ({fDetJ} < {determinantTolerance})." 
                     + " Check the order of nodes or the element geometry.");
            }

            double fDetInv = 1.0 / fDetJ;
            double[,] faJInv = new double[3, 3];
            faJInv[0, 0] = (faJ[1, 1] * faJ[2, 2] - faJ[2, 1] * faJ[1, 2]) * fDetInv;
            faJInv[1, 0] = (faJ[2, 0] * faJ[1, 2] - faJ[1, 0] * faJ[2, 2]) * fDetInv;
            faJInv[2, 0] = (faJ[1, 0] * faJ[2, 1] - faJ[2, 0] * faJ[1, 1]) * fDetInv;
            faJInv[0, 1] = (faJ[2, 1] * faJ[0, 2] - faJ[0, 1] * faJ[2, 2]) * fDetInv;
            faJInv[1, 1] = (faJ[0, 0] * faJ[2, 2] - faJ[2, 0] * faJ[0, 2]) * fDetInv;
            faJInv[2, 1] = (faJ[2, 0] * faJ[0, 1] - faJ[2, 1] * faJ[0, 0]) * fDetInv;
            faJInv[0, 2] = (faJ[0, 1] * faJ[1, 2] - faJ[1, 1] * faJ[0, 2]) * fDetInv;
            faJInv[1, 2] = (faJ[1, 0] * faJ[0, 2] - faJ[0, 0] * faJ[1, 2]) * fDetInv;
            faJInv[2, 2] = (faJ[0, 0] * faJ[1, 1] - faJ[1, 0] * faJ[0, 1]) * fDetInv;

            return new Tuple<double[,], double[,], double>(faJ, faJInv, fDetJ);
        }

        public EmbeddedNode BuildHostElementEmbeddedNode(Element element, Node node,
            IEmbeddedDOFInHostTransformationVector transformationVector)
        {
            if (!(Interpolation.CellType == CellType.Hexa8))
            {
                throw new NotImplementedException();
            }

            var points = GetNaturalCoordinates(element, node);
            if (points.Length == 0) return null;

            element.EmbeddedNodes.Add(node);
            var embeddedNode = new EmbeddedNode(node, element, transformationVector.GetDependentDOFTypes);
            for (int i = 0; i < points.Length; i++) embeddedNode.Coordinates.Add(points[i]);
            return embeddedNode;
        }

        private double[] GetNaturalCoordinates(IElement element, Node node)
        {
            if (!(Interpolation.CellType == CellType.Hexa8))
            {
                throw new NotImplementedException();
            }

            double[] mins = new double[] { element.Nodes[0].X, element.Nodes[0].Y, element.Nodes[0].Z };
            double[] maxes = new double[] { element.Nodes[0].X, element.Nodes[0].Y, element.Nodes[0].Z };
            for (int i = 0; i < element.Nodes.Count; i++)
            {
                mins[0] = mins[0] > element.Nodes[i].X ? element.Nodes[i].X : mins[0];
                mins[1] = mins[1] > element.Nodes[i].Y ? element.Nodes[i].Y : mins[1];
                mins[2] = mins[2] > element.Nodes[i].Z ? element.Nodes[i].Z : mins[2];
                maxes[0] = maxes[0] < element.Nodes[i].X ? element.Nodes[i].X : maxes[0];
                maxes[1] = maxes[1] < element.Nodes[i].Y ? element.Nodes[i].Y : maxes[1];
                maxes[2] = maxes[2] < element.Nodes[i].Z ? element.Nodes[i].Z : maxes[2];
            }
            //return new double[] { (node.X - mins[0]) / ((maxes[0] - mins[0]) / 2) - 1,
            //    (node.Y - mins[1]) / ((maxes[1] - mins[1]) / 2) - 1,
            //    (node.Z - mins[2]) / ((maxes[2] - mins[2]) / 2) - 1 };

            bool maybeInsideElement = node.X <= maxes[0] && node.X >= mins[0] &&
                node.Y <= maxes[1] && node.Y >= mins[1] &&
                node.Z <= maxes[2] && node.Z >= mins[2];
            if (maybeInsideElement == false) return new double[0];

            const int jacobianSize = 3;
            const int maxIterations = 1000;
            const double tolerance = 1e-10;
            int iterations = 0;
            double deltaNaturalCoordinatesNormSquare = 100;
            double[] naturalCoordinates = new double[] { 0, 0, 0 };
            const double toleranceSquare = tolerance * tolerance;

            while (deltaNaturalCoordinatesNormSquare > toleranceSquare && iterations < maxIterations)
            {
                iterations++;
                var shapeFunctions = CalcH8Shape(naturalCoordinates[0], naturalCoordinates[1], naturalCoordinates[2]);
                double[] coordinateDifferences = new double[] { 0, 0, 0 };
                for (int i = 0; i < shapeFunctions.Length; i++)
                {
                    coordinateDifferences[0] += shapeFunctions[i] * element.Nodes[i].X;
                    coordinateDifferences[1] += shapeFunctions[i] * element.Nodes[i].Y;
                    coordinateDifferences[2] += shapeFunctions[i] * element.Nodes[i].Z;
                }
                coordinateDifferences[0] = node.X - coordinateDifferences[0];
                coordinateDifferences[1] = node.Y - coordinateDifferences[1];
                coordinateDifferences[2] = node.Z - coordinateDifferences[2];

                double[,] faXYZ = GetCoordinatesTranspose(element);
                double[] nablaShapeFunctions = CalcH8NablaShape(naturalCoordinates[0], naturalCoordinates[1], naturalCoordinates[2]);
                var inverseJacobian = CalcH8JDetJ(faXYZ, nablaShapeFunctions).Item2;

                double[] deltaNaturalCoordinates = new double[] { 0, 0, 0 };
                for (int i = 0; i < jacobianSize; i++)
                    for (int j = 0; j < jacobianSize; j++)
                        deltaNaturalCoordinates[i] += inverseJacobian[j, i] * coordinateDifferences[j];
                for (int i = 0; i < 3; i++)
                    naturalCoordinates[i] += deltaNaturalCoordinates[i];

                deltaNaturalCoordinatesNormSquare = 0;
                for (int i = 0; i < 3; i++)
                    deltaNaturalCoordinatesNormSquare += deltaNaturalCoordinates[i] * deltaNaturalCoordinates[i];
                //deltaNaturalCoordinatesNormSquare = Math.Sqrt(deltaNaturalCoordinatesNormSquare);
            }

            return naturalCoordinates.Count(x => Math.Abs(x) - 1.0 > tolerance) > 0 ? new double[0] : naturalCoordinates;
        }

        public double[] GetShapeFunctionsForNode(Element element, EmbeddedNode node)
        {
            double[,] elementCoordinates = GetCoordinatesTranspose(element);
            var shapeFunctions = CalcH8Shape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
            var nablaShapeFunctions = CalcH8NablaShape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
            var jacobian = CalcH8JDetJ(elementCoordinates, nablaShapeFunctions);

            return new double[]
            {
                shapeFunctions[0], shapeFunctions[1], shapeFunctions[2], shapeFunctions[3], shapeFunctions[4], shapeFunctions[5], shapeFunctions[6], shapeFunctions[7],
                nablaShapeFunctions[0], nablaShapeFunctions[1], nablaShapeFunctions[2], nablaShapeFunctions[3], nablaShapeFunctions[4], nablaShapeFunctions[5], nablaShapeFunctions[6], nablaShapeFunctions[7],
                nablaShapeFunctions[8], nablaShapeFunctions[9], nablaShapeFunctions[10], nablaShapeFunctions[11], nablaShapeFunctions[12], nablaShapeFunctions[13], nablaShapeFunctions[14], nablaShapeFunctions[15],
                nablaShapeFunctions[16], nablaShapeFunctions[17], nablaShapeFunctions[18], nablaShapeFunctions[19], nablaShapeFunctions[20], nablaShapeFunctions[21], nablaShapeFunctions[22], nablaShapeFunctions[23],
                jacobian.Item1[0, 0], jacobian.Item1[0, 1], jacobian.Item1[0, 2], jacobian.Item1[1, 0], jacobian.Item1[1, 1], jacobian.Item1[1, 2], jacobian.Item1[2, 0], jacobian.Item1[2, 1], jacobian.Item1[2, 2],
                jacobian.Item2[0, 0], jacobian.Item2[0, 1], jacobian.Item2[0, 2], jacobian.Item2[1, 0], jacobian.Item2[1, 1], jacobian.Item2[1, 2], jacobian.Item2[2, 0], jacobian.Item2[2, 1], jacobian.Item2[2, 2]
            };
        }

        #endregion

    }


}
