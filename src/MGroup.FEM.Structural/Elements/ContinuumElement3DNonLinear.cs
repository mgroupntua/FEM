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
    public class ContinuumElement3DNonLinear : IStructuralFiniteElement//, IEmbeddedHostElement
    {
        protected readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
        protected readonly IDofType[][] dofTypes;
        protected readonly IContinuumMaterial3D[] materialsAtGaussPoints;
        protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
        
        private readonly int nGaussPoints;
        private bool isInitialized = false;

        private double[][] intialCoordinates; //not defined by user. 8 arrays of 3 elements
        private double[][] totalDisplacements;
        private double[] integrationCoeffs;

        private double[][] strainsVec;
        private double[][] strainsVecLastConverged;

        protected ContinuumElement3DNonLinear()
        {
        }

        public ContinuumElement3DNonLinear(IReadOnlyList<Node> nodes,IContinuumMaterial3D material, IQuadrature3D quadratureForStiffness,
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

        private Matrix[] Getbl13DeformationMatrices(IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives)
        {
            Matrix[] bl13Matrices;
            bl13Matrices = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                bl13Matrices[npoint] = Matrix.CreateZero(9, 3*shapeFunctionNaturalDerivatives[npoint].NumRows);
                for (int m = 0; m < shapeFunctionNaturalDerivatives[npoint].NumRows; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl13Matrices[npoint][n, 3 * m + 0] = shapeFunctionNaturalDerivatives[npoint][ m,n];
                        bl13Matrices[npoint][n + 3, 3 * m + 1] = shapeFunctionNaturalDerivatives[npoint][m,n];
                        bl13Matrices[npoint][n + 6, 3 * m + 2] = shapeFunctionNaturalDerivatives[npoint][m, n];
                    }
                }
            }
            return bl13Matrices;
        }

        private Matrix[] Getbl11aDeformationMatrices(Matrix[] jacobianInverse)
        {
            Matrix[] bl11aMatrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                bl11aMatrices[gpoint] = Matrix.CreateZero(6, 9);
                for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl11aMatrices[gpoint][m, 3 * m + n] = jacobianInverse[gpoint][m, n];
                    }
                }
                for (int n = 0; n < 3; n++)
                {
                    bl11aMatrices[gpoint][3, n] = jacobianInverse[gpoint][1, n]; // calculate 4th data line
                    bl11aMatrices[gpoint][3, 3 + n] = jacobianInverse[gpoint][0, n];
                    bl11aMatrices[gpoint][4, 3 + n] = jacobianInverse[gpoint][2, n]; // calculate 5th data line
                    bl11aMatrices[gpoint][4, 6 + n] = jacobianInverse[gpoint][1, n];
                    bl11aMatrices[gpoint][5, 0 + n] = jacobianInverse[gpoint][2, n]; // calculate 6th data line
                    bl11aMatrices[gpoint][5, 6 + n] = jacobianInverse[gpoint][0, n];
                }
            }

            return bl11aMatrices;
        }

        private Matrix[] GetBL12DeformationMatrices(Matrix[] jacobianInverse)
        {
            Matrix[] bl12Marices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                bl12Marices[gpoint] = Matrix.CreateZero(9, 9);
                for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl12Marices[gpoint][m, 3 * m + n] = jacobianInverse[gpoint][0, n];
                    }
                }
                for (int m = 0; m < 3; m++) // calculate  data lines 4:6
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl12Marices[gpoint][3 + m, 3 * m + n] = jacobianInverse[gpoint][1, n];
                    }
                }
                for (int m = 0; m < 3; m++) // calculate  data lines 7:8
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl12Marices[gpoint][6 + m, 3 * m + n] = jacobianInverse[gpoint][2, n];
                    }
                }

            }

            return bl12Marices;
        }

        private Matrix[] Getbl01MDeformationMatrices(Matrix[] jacobianInverse)
        {
            Matrix[] bl01Matrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                bl01Matrices[gpoint] = Matrix.CreateZero(6, 9);
                for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl01Matrices[gpoint][m, 3 * m + n] = jacobianInverse[gpoint][m, n];
                    }
                }
                for (int n = 0; n < 3; n++)
                {
                    bl01Matrices[gpoint][3, n] = jacobianInverse[gpoint][1, n]; // calculate 4th data line
                    bl01Matrices[gpoint][3, 3 + n] = jacobianInverse[gpoint][0, n];
                    bl01Matrices[gpoint][4, 3 + n] = jacobianInverse[gpoint][2, n]; // calculate 5th data line
                    bl01Matrices[gpoint][4, 6 + n] = jacobianInverse[gpoint][1, n];
                    bl01Matrices[gpoint][5, 0 + n] = jacobianInverse[gpoint][2, n]; // calculate 6th data line
                    bl01Matrices[gpoint][5, 6 + n] = jacobianInverse[gpoint][0, n];
                }
            }
            return bl01Matrices;
        }

        private Matrix[] GetAuxilliaryDeformationbnl1Matrices(Matrix[] jacobianInverse)
        {
            Matrix[] bnl1Matrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                bnl1Matrices[gpoint] = Matrix.CreateZero(9, 9);
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            bnl1Matrices[gpoint][3 * m + n, 3 * m + p] = jacobianInverse[gpoint][n, p];
                        }
                    }
                }
            }
            return bnl1Matrices;
        }

        private void CalculateInitialConfigurationData(IElement element)
        {
            int numNodes = element.Nodes.Count;
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            Matrix[] bl13Matrices;
            bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);
            

            Matrix[] bnl1Matrices;

            intialCoordinates = new double[numNodes][];
            totalDisplacements = new double[numNodes][];

            var jacobians = shapeFunctionNaturalDerivatives.Select(x=> new IsoparametricJacobian3D(element.Nodes,x));

            Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();

            



            integrationCoeffs = new double[nGaussPoints];

            bnl1Matrices = GetAuxilliaryDeformationbnl1Matrices(jacobianInverse);

            for (int j = 0; j < numNodes; j++)
            {
                intialCoordinates[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
                }

            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                integrationCoeffs[gpoint] = jacobianDeterminants[gpoint] * QuadratureForStiffness.IntegrationPoints[gpoint].Weight;

            }

            totalDisplacements = new double[numNodes][];
            strainsVec = new double[nGaussPoints][];
            strainsVecLastConverged = new double[nGaussPoints][];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                strainsVec[gpoint] = new double[6];
                strainsVecLastConverged[gpoint] = new double[6];
            }
            for (int k = 0; k < numNodes; k++)
            {
                totalDisplacements[k] = new double[3];
            }
            isInitialized = true;

        }

        private void UpdateCoordinateData(double[] localdisplacements, out double[][] deformedCoordinates)
        {
            int numNodes = localdisplacements.Length / 3;
            deformedCoordinates = new double[numNodes][];
            for (int j = 0; j < numNodes; j++)
            {
                deformedCoordinates[j] = new double[3];
                for (int k = 0; k < 3; k++)
                {
                    totalDisplacements[j][k] = localdisplacements[3 * j + k];
                    deformedCoordinates[j][k] = intialCoordinates[j][k] + totalDisplacements[j][k];
                }
            }
        }

        private void CalculateStrains(double[] localdisplacements, IElement element, double[][] deformedCoordinates) 
        {
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));
            Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();
            //TODO: possibility of caching shapeFunctionNaturalDerivatives or J_0inv

            Matrix[] deformationGradients = new Matrix[nGaussPoints];
            Matrix[] strains = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                deformationGradients[npoint] = Matrix.CreateZero(3, 3);
                strains[npoint] = Matrix.CreateZero(3, 3);
            }

            var jacobiansDeformed = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(deformedCoordinates, x,false)).ToArray();
            Matrix[] jacobiansDeformedMatrices = jacobiansDeformed.Select(x => x.DirectMatrix).ToArray();


            //
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                //
                deformationGradients[npoint] = jacobianInverse[npoint] * jacobiansDeformedMatrices[npoint];

                //
                strains[npoint] = deformationGradients[npoint] * deformationGradients[npoint].Transpose();
                for (int m = 0; m < 3; m++)
                {
                    strains[npoint][m, m] += -1;
                }
                strains[npoint].ScaleIntoThis(0.5);
                
                //
                for (int m = 0; m < 3; m++)
                {
                    strainsVec[npoint][m] = strains[npoint][m, m];
                }
                strainsVec[npoint][3] = 2 * strains[npoint][0, 1];
                strainsVec[npoint][4] = 2 * strains[npoint][1, 2];
                strainsVec[npoint][5] = 2 * strains[npoint][2, 0];
            }

        }

        private double[] UpdateForces(IElement element)
        {
            //TODO: the gauss point loop should be the outer one


            // Matrices that are not currently cached are calculated here.
            int numNodes = element.Nodes.Count();
            Matrix totalDisplacementsMatrixReordered = Matrix.CreateZero(numNodes, 3);
            for (int m = 0; m < numNodes; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    totalDisplacementsMatrixReordered[m, n] = totalDisplacements[m][n];
                }
            }
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));
            Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();

            Matrix[] bl13Matrices;
            bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);
            Matrix[] bl11aMatrices; // dimension number of gpoints
            Matrix[] BL12Matrices;
            Matrix[] bl01Matrices;
            bl11aMatrices = Getbl11aDeformationMatrices(jacobianInverse);
            BL12Matrices = GetBL12DeformationMatrices(jacobianInverse);
            bl01Matrices = Getbl01MDeformationMatrices(jacobianInverse);

            //INITIALIZATION of MAtrixes that are currently not cached
            double[][] integrCoeffsTimesStresses = new double[nGaussPoints][];
            Matrix[] blMatrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                integrCoeffsTimesStresses[gpoint] = new double[6];
                blMatrices[gpoint] = Matrix.CreateZero(6, 3*numNodes);
            }

            double[][] forces = new double[nGaussPoints + 1][];
            for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
            {
                forces[npoint] = new double[3*numNodes];
            }

            Matrix[] bl11Matrices = new Matrix[nGaussPoints];
            Matrix[] bL1112Plus01Matrices = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                bl11Matrices[npoint] = Matrix.CreateZero(6, 9);
                bL1112Plus01Matrices[npoint] = Matrix.CreateZero(6, 9);
            }

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {

                integrCoeffsTimesStresses[npoint] = materialsAtGaussPoints[npoint].Stresses.Scale(integrationCoeffs[npoint]);

                //
                Matrix lcyrcumflex;//= Matrix.CreateZero(3, 3);
                lcyrcumflex = shapeFunctionNaturalDerivatives[npoint].Transpose() * totalDisplacementsMatrixReordered;

                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            bl11Matrices[npoint][m, n] += bl11aMatrices[npoint][m, p] * lcyrcumflex[p, n];
                            bl11Matrices[npoint][m, 3 + n] += bl11aMatrices[npoint][m, 3 + p] * lcyrcumflex[p, n];
                            bl11Matrices[npoint][m, 6 + n] += bl11aMatrices[npoint][m, 6 + p] * lcyrcumflex[p, n];
                        }
                    }
                }

                //
                bL1112Plus01Matrices[npoint] = bl11Matrices[npoint] * BL12Matrices[npoint];
                bL1112Plus01Matrices[npoint].AddIntoThis(bl01Matrices[npoint]);
                
                // 
                blMatrices[npoint] = bL1112Plus01Matrices[npoint] * bl13Matrices[npoint];
                
                //              
                forces[npoint] = blMatrices[npoint].Multiply(integrCoeffsTimesStresses[npoint], true);
            }

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                forces[nGaussPoints].AddIntoThis(forces[npoint]);                
            }

            return forces[nGaussPoints];
        }

        private Matrix UpdateKmatrices(IElement element)
        {
            int numNodes = element.Nodes.Count();
            Matrix elementStiffnessMatrix = Matrix.CreateZero(3*numNodes, 3*numNodes);


            // initialization of matrices that are not cached currently
            double[][] integrCoeffsTimesSpkvec = new double[nGaussPoints][];
            Matrix[] BLMatrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                integrCoeffsTimesSpkvec[gpoint] = new double[6];
                BLMatrices[gpoint] = Matrix.CreateZero(6, 3*numNodes);

            }
            Matrix totalDisplacementsMatrixReordered = Matrix.CreateZero(numNodes, 3);
            for (int m = 0; m < numNodes; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    totalDisplacementsMatrixReordered[m, n] = totalDisplacements[m][n];
                }
            }
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));
            Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();

            Matrix[] bl13Matrices;
            bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);
            Matrix[] bl11aMatrices; // dimension: gpoints
            Matrix[] BL12Matrices;
            Matrix[] bl01Matrices;
            bl11aMatrices = Getbl11aDeformationMatrices(jacobianInverse);
            BL12Matrices = GetBL12DeformationMatrices(jacobianInverse);
            bl01Matrices = Getbl01MDeformationMatrices(jacobianInverse);

            Matrix[] BL11Matrices = new Matrix[nGaussPoints];
            Matrix[] bl1112Plus01Mtrices = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                BL11Matrices[npoint] = Matrix.CreateZero(6, 9);
                bl1112Plus01Mtrices[npoint] = Matrix.CreateZero(6, 9); //TODO this may be unnescessary
            }



            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {

                // 
                integrCoeffsTimesSpkvec[npoint] = materialsAtGaussPoints[npoint].Stresses.Scale(integrationCoeffs[npoint]);
                
                //
                Matrix lcyrcumflex = Matrix.CreateZero(3, 3);
                lcyrcumflex = shapeFunctionNaturalDerivatives[npoint].Transpose() * totalDisplacementsMatrixReordered;
                 
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            BL11Matrices[npoint][m, n] += bl11aMatrices[npoint][m, p] * lcyrcumflex[p, n];
                            BL11Matrices[npoint][m, 3 + n] += bl11aMatrices[npoint][m, 3 + p] * lcyrcumflex[p, n];
                            BL11Matrices[npoint][m, 6 + n] += bl11aMatrices[npoint][m, 6 + p] * lcyrcumflex[p, n];
                        }
                    }
                }

                // 
                bl1112Plus01Mtrices[npoint] = BL11Matrices[npoint] * BL12Matrices[npoint];
                bl1112Plus01Mtrices[npoint].AddIntoThis(bl01Matrices[npoint]);
                
                //
                BLMatrices[npoint] = bl1112Plus01Mtrices[npoint] * bl13Matrices[npoint];
                
            }
            // TODO: BL and above calculations can cached from calculate forces method

            Matrix[] bnl1Matrices;
            Matrix[] bnlMatrices;
            bnl1Matrices = GetAuxilliaryDeformationbnl1Matrices(jacobianInverse);
            bnlMatrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
               //bnlMatrices[gpoint] = Matrix.CreateZero(9, 3*numNodes); //todo this may be unnescessary

                bnlMatrices[gpoint] = bnl1Matrices[gpoint] * bl13Matrices[gpoint];
                
            }


            Matrix[] integrCoeffsTimesStresses = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                integrCoeffsTimesStresses[npoint] = Matrix.CreateZero(3, 3);
            }

            Matrix[] klStiffnessMatrixContributions = new Matrix[nGaussPoints + 1];
            Matrix[] knlStiffnessMatrixContributions = new Matrix[nGaussPoints + 1];
            for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
            {
                klStiffnessMatrixContributions[npoint] = Matrix.CreateZero(3*numNodes, 3*numNodes);
                knlStiffnessMatrixContributions[npoint] = Matrix.CreateZero(3 * numNodes, 3 * numNodes);
            }



            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                Matrix integrCoeffsTimesStressesTimesbnlMatrices = Matrix.CreateZero(9, 3*numNodes); //TODO
                Matrix integrCoeffsTimesConsMatrix = Matrix.CreateZero(6, 6); //TODO
                Matrix integrCoeffTimesConsMatrixTimesBLMatrices = Matrix.CreateZero(6, 3*numNodes);//TODO

                //
                integrCoeffsTimesStresses[npoint][0, 0] = integrCoeffsTimesSpkvec[npoint][0];
                integrCoeffsTimesStresses[npoint][0, 1] = integrCoeffsTimesSpkvec[npoint][3];
                integrCoeffsTimesStresses[npoint][0, 2] = integrCoeffsTimesSpkvec[npoint][5];
                integrCoeffsTimesStresses[npoint][1, 0] = integrCoeffsTimesSpkvec[npoint][3];
                integrCoeffsTimesStresses[npoint][1, 1] = integrCoeffsTimesSpkvec[npoint][1];
                integrCoeffsTimesStresses[npoint][1, 2] = integrCoeffsTimesSpkvec[npoint][4];
                integrCoeffsTimesStresses[npoint][2, 0] = integrCoeffsTimesSpkvec[npoint][5];
                integrCoeffsTimesStresses[npoint][2, 1] = integrCoeffsTimesSpkvec[npoint][4];
                integrCoeffsTimesStresses[npoint][2, 2] = integrCoeffsTimesSpkvec[npoint][2];

                //
                IMatrixView consDisp = materialsAtGaussPoints[npoint].ConstitutiveMatrix;

                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 6; n++)
                    {
                        integrCoeffsTimesConsMatrix[m, n] = integrationCoeffs[npoint] * consDisp[m, n];
                    }
                }

                //
                integrCoeffTimesConsMatrixTimesBLMatrices = integrCoeffsTimesConsMatrix * BLMatrices[npoint];
                
                //
                klStiffnessMatrixContributions[npoint] = BLMatrices[npoint].Transpose() * integrCoeffTimesConsMatrixTimesBLMatrices;                

                //
                for (int m = 0; m < 3; m++) // 3x24 dimensions
                {
                    for (int n = 0; n < 3*numNodes; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            integrCoeffsTimesStressesTimesbnlMatrices[m, n] += integrCoeffsTimesStresses[npoint][m, p] * bnlMatrices[npoint][p, n];
                            integrCoeffsTimesStressesTimesbnlMatrices[3 + m, n] += integrCoeffsTimesStresses[npoint][m, p] * bnlMatrices[npoint][3 + p, n];
                            integrCoeffsTimesStressesTimesbnlMatrices[6 + m, n] += integrCoeffsTimesStresses[npoint][m, p] * bnlMatrices[npoint][6 + p, n];
                        }
                    }
                }

                //
                knlStiffnessMatrixContributions[npoint] = bnlMatrices[npoint].Transpose() * integrCoeffsTimesStressesTimesbnlMatrices;                
            }

            // Add contributions of each gp on the total element stiffness matrix elementStiffnessMatrix            
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                for (int m = 0; m < 3*numNodes; m++)
                {
                    for (int n = 0; n < 3*numNodes; n++)
                    {
                        klStiffnessMatrixContributions[nGaussPoints][m, n] += klStiffnessMatrixContributions[npoint][m, n];
                        knlStiffnessMatrixContributions[nGaussPoints][m, n] += knlStiffnessMatrixContributions[npoint][m, n];
                    }
                }
            }
            for (int m = 0; m < 3 * numNodes; m++)
            {
                for (int n = 0; n < 3 * numNodes; n++)
                {
                    elementStiffnessMatrix[m, n] = klStiffnessMatrixContributions[nGaussPoints][m, n] + knlStiffnessMatrixContributions[nGaussPoints][m, n];
                }
            }

            return elementStiffnessMatrix;
        }
        
        public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            this.UpdateCoordinateData(localTotalDisplacements, out double[][] deformedCoordinates);
            this.CalculateStrains(localTotalDisplacements, element, deformedCoordinates);
            double[] strainsVecMinusLastConvergedValues = new double[6];
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                strainsVecMinusLastConvergedValues = new double[6] 
                {
                    strainsVec[npoint][0]- strainsVecLastConverged[npoint][0],
                    strainsVec[npoint][1] - strainsVecLastConverged[npoint][1],
                    strainsVec[npoint][2] - strainsVecLastConverged[npoint][2],
                    strainsVec[npoint][3]- strainsVecLastConverged[npoint][3],
                    strainsVec[npoint][4]- strainsVecLastConverged[npoint][4],
                    strainsVec[npoint][5]- strainsVecLastConverged[npoint][5]
                };
                materialsAtGaussPoints[npoint].UpdateMaterial(strainsVecMinusLastConvergedValues); 
                //To update with total strain simplY = materialsAtGaussPoints[npoint].UpdateMaterial(strainsvec[npoint]);
            }
            return new Tuple<double[], double[]>(strainsVecMinusLastConvergedValues, 
                materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);
            //TODO return data with total strains data would be:
            //return new Tuple<double[], double[]>(strainsvec[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);
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
                this.UpdateCoordinateData(localTotalDisplacements, out double[][] deformedCoordinates);
                this.CalculateStrains(localTotalDisplacements, element, deformedCoordinates);
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
                { strainsVecLastConverged[npoint][i1] = strainsVec[npoint][i1]; }
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

		
    }


}
