using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.FEM.Structural.Helpers
{
    public static class TransformationMethods
    {
        public static double[,] Calculate_dSdE_from_dSdEe(double[,] Cijrs_el, double[,] F,
            double lamda_g, double[,] SpkMat, double[,] F_el, double[,] SPK_el, double[,] Fg_inv, double detFg)
        {
            double[,] Cinpk_el = ExpandCijrsMatrixToTensor(Cijrs_el);
            double[,] Aijkl_el = Transform3DCinpkToAijkl(Cinpk_el, F_el, SPK_el);

            double[,] dFel_dF = new double[9, 9];

            dFel_dF[0, 0] = ((double)1) / lamda_g;
            dFel_dF[1, 1] = dFel_dF[0, 0];
            dFel_dF[2, 2] = dFel_dF[0, 0];

            //1
            //dFel_dF[3, 7] = dFel_dF[0, 0];
            //dFel_dF[4, 8] = dFel_dF[0, 0];
            //dFel_dF[5, 6] = dFel_dF[0, 0];

            //dFel_dF[6, 5] = dFel_dF[0, 0];
            //dFel_dF[7, 3] = dFel_dF[0, 0];
            //dFel_dF[8, 4] = dFel_dF[0, 0];

            //2
            for (int i1 = 3; i1 < 9; i1++)
            {
                dFel_dF[i1, i1] = dFel_dF[0, 0];
            }

            var dFel_dFmat = Matrix.CreateFromArray(dFel_dF);
            var dWel_dFelmat = Matrix.CreateFromArray(Aijkl_el);

            var dPel_dF = dWel_dFelmat * dFel_dFmat;

            //paper
            double[,] fg_inv = new double[9, 9];

            fg_inv[0, 0] = Fg_inv[0, 0];
            fg_inv[0, 5] = Fg_inv[2, 0];
            fg_inv[0, 7] = Fg_inv[1, 0];

            fg_inv[1, 1] = Fg_inv[1, 1];
            fg_inv[1, 3] = Fg_inv[0, 1];
            fg_inv[1, 8] = Fg_inv[2, 1];

            fg_inv[2, 2] = Fg_inv[2, 2];
            fg_inv[2, 4] = Fg_inv[1, 2];
            fg_inv[2, 6] = Fg_inv[0, 2];

            fg_inv[3, 1] = Fg_inv[1, 0];
            fg_inv[3, 3] = Fg_inv[0, 0];
            fg_inv[3, 8] = Fg_inv[2, 0];

            fg_inv[4, 2] = Fg_inv[2, 1];
            fg_inv[4, 4] = Fg_inv[1, 1];
            fg_inv[4, 6] = Fg_inv[0, 1];

            fg_inv[5, 0] = Fg_inv[0, 2];
            fg_inv[5, 5] = Fg_inv[2, 2];
            fg_inv[5, 7] = Fg_inv[1, 2];

            fg_inv[6, 2] = Fg_inv[2, 0];
            fg_inv[6, 4] = Fg_inv[1, 0];
            fg_inv[6, 6] = Fg_inv[0, 0];

            fg_inv[7, 0] = Fg_inv[0, 1];
            fg_inv[7, 5] = Fg_inv[2, 1];
            fg_inv[7, 7] = Fg_inv[1, 1];

            fg_inv[8, 1] = Fg_inv[1, 2];
            fg_inv[8, 3] = Fg_inv[0, 2];
            fg_inv[8, 8] = Fg_inv[2, 2];

            var Fg_inv_mat = Matrix.CreateFromArray(fg_inv);

            var dP_dF = Fg_inv_mat * dPel_dF * detFg;

            var Cinpk = Transform_d2WdFtrdFtr_to_Cijrs(dP_dF.CopyToArray2D(), SpkMat, F);

            double[,] Cijrs = CombineCinpkTensorTermsIntoMatrix(Cinpk);


            #region elegxoi variation
            //if (ElementStiffnessesReused.performCalculations)
            //{
            //    if ((ElementStiffnessesReused.gpNumber == ElementStiffnessesReused.gpNumberToCheck) && ElementStiffnessesReused.saveSt-iffnessMatrixState)
            //    {
            //        ElementStiffnessesReused.saveOriginalState = true;

            //        ElementStiffnessesReused.ProccessVariable(25, da3_dksidrds[0, 0], true, 3 * k + 0);
            //        ElementStiffnessesReused.ProccessVariable(25, da3_dksidrds[0, 1], true, 3 * k + 1);
            //        ElementStiffnessesReused.ProccessVariable(25, da3_dksidrds[0, 2], true, 3 * k + 2);
            //        ElementStiffnessesReused.ProccessVariable(25, da3_dksidr[0], false);
            //              ElementStiffnessesReused.saveOriginalState = false;

            //    }


            //}
            #endregion

            return Cijrs;
        }

        private static double[,] ExpandCijrsMatrixToTensor(double[,] Cijrs)
        {
            double[,] Cijrs_columns = new double[6, 9];
            for (int i1 = 0; i1 < 6; i1++)
            {
                Cijrs_columns[i1, 0] = Cijrs[i1, 0];
                Cijrs_columns[i1, 1] = Cijrs[i1, 1];
                Cijrs_columns[i1, 2] = Cijrs[i1, 2];
                Cijrs_columns[i1, 3] = Cijrs[i1, 3];
                Cijrs_columns[i1, 4] = Cijrs[i1, 4];
                Cijrs_columns[i1, 5] = Cijrs[i1, 5];
                Cijrs_columns[i1, 7] = Cijrs[i1, 3];
                Cijrs_columns[i1, 8] = Cijrs[i1, 4];
                Cijrs_columns[i1, 6] = Cijrs[i1, 5];
            }

            double[,] Cinpk = new double[9, 9];
            for (int i1 = 0; i1 < 9; i1++)
            {
                Cinpk[0, i1] = Cijrs_columns[0, i1];
                Cinpk[1, i1] = Cijrs_columns[1, i1];
                Cinpk[2, i1] = Cijrs_columns[2, i1];
                Cinpk[3, i1] = Cijrs_columns[3, i1];
                Cinpk[4, i1] = Cijrs_columns[4, i1];
                Cinpk[5, i1] = Cijrs_columns[5, i1];
                Cinpk[7, i1] = Cijrs_columns[3, i1];
                Cinpk[8, i1] = Cijrs_columns[4, i1];
                Cinpk[6, i1] = Cijrs_columns[5, i1];
            }


            return Cinpk;
        }

        private static double[,] Transform3DCinpkToAijkl(double[,] Cinpk, double[,] F_rve, double[,] SPKMat)
        {
            int[,] thesi_ijkl = new int[3, 3] { { 0, 3, 6 }, { 7, 1, 4 }, { 5, 8, 2 } }; // 11 22 33 12 23 31 13 21 32
            int row_ij;//anaferetai sto Aijkl_3D mhtrwo
            int col_kl;// anaferetai sto Aijkl_3D mhtrwo

            double[,] Aijkl = new double[9, 9];
            int[,] dlj = new int[3, 3] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

            // exwterika loop einai oi theseis pou gemizoun sto Aijkl
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    row_ij = thesi_ijkl[i, j];
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            col_kl = thesi_ijkl[k, l];

                            Aijkl[row_ij, col_kl] += SPKMat[i, k] * dlj[l, j];

                            //eswterika loop einai to athroisma logw n kai p
                            for (int n = 0; n < 3; n++)
                            {
                                int rowC = thesi_ijkl[i, n];
                                for (int p = 0; p < 3; p++)
                                {
                                    int colC = thesi_ijkl[p, k];

                                    Aijkl[row_ij, col_kl] += Cinpk[rowC, colC] * F_rve[j, n] * F_rve[l, p];
                                    //Aijkl[row, col] += Cinpk[rowC, colC] * F_rve[j, n] * F_rve[l, p] + SPKMat[i, k] * dlj[l, j];
                                }
                            }

                        }
                    }
                }
            }

            return Aijkl;
        }

        private static double[,] Transform_d2WdFtrdFtr_to_Cijrs(double[,] Aijkl, double[,] SPK, double[,] F)
        {
            int[,] i_seira = { { 1, 2, 3 }, { 3, 1, 2 }, { 2, 3, 1 } };
            int[,] k_seira = { { 1, 2, 3 }, { 3, 1, 2 }, { 2, 3, 1 } };

            double[,] Cinpk = new double[9, 9];

            double[,] F__F__ = new double[9, 9];


            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            F__F__[3 * i1 + k, 3 * j1 + l] = F[k, j1] * F[i1, l];
                        }
                    }

                }
            }

            Matrix F__F__Mat = Matrix.CreateFromArray(F__F__);

            double[,] multipleRHSs = new double[9, 9];
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int k1 = 0; k1 < 3; k1++)
                {
                    double[] A_j_l = new double[9] { Aijkl[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[k1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[k1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[k1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[k1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[k1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[k1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[k1, 2] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[k1, 2] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[k1, 2] - 1) + k1]};

                    double[] sec_term = new double[9] { -SPK[i1, k1], 0, 0, 0, -SPK[i1, k1], 0, 0, 0, -SPK[i1, k1] };

                    int RHScolumn = i1 * 3 + k1;
                    for (int a1 = 0; a1 < 9; a1++)
                    {
                        multipleRHSs[a1, RHScolumn] = A_j_l[a1] + sec_term[a1];
                    }

                }
            }

            //TODO use solution multiple RHSs when pavailable :Matrix2D MultipleSolutions = F__F__Mat.SolveLU(new Matrix2D(multipleRHSs), true);
            Matrix inverse = F__F__Mat.Invert();
            Matrix MultipleSolutions = Matrix.CreateZero(9, 9);
            for (int i1 = 0; i1 < 9; i1++)
            {
                double[] RHS = new double[9];
                for (int i2 = 0; i2 < 9; i2++)
                {
                    RHS[i2] = multipleRHSs[i2, i1];
                }
                Vector solution = inverse * Vector.CreateFromArray(RHS);
                for (int i2 = 0; i2 < 9; i2++)
                {
                    MultipleSolutions[i2, i1] = solution[i2];
                }

            }


            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int k1 = 0; k1 < 3; k1++)
                {
                    int RHScolumn = i1 * 3 + k1;

                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 0] = MultipleSolutions[0, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 1] = MultipleSolutions[1, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 2] = MultipleSolutions[2, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 0] = MultipleSolutions[3, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 1] = MultipleSolutions[4, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 2] = MultipleSolutions[5, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 0] = MultipleSolutions[6, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 1] = MultipleSolutions[7, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 2] = MultipleSolutions[8, RHScolumn];

                }
            }


            return Cinpk;
        }

        private static double[,] CombineCinpkTensorTermsIntoMatrix(double[,] Cinpk)
        {
            // transformation se 6x6 se 2 vhmata

            double[,] Cijrs_columns = new double[9, 6];
            for (int i1 = 0; i1 < 9; i1++)
            {
                Cijrs_columns[i1, 0] = Cinpk[i1, 0];
                Cijrs_columns[i1, 1] = Cinpk[i1, 1];
                Cijrs_columns[i1, 2] = Cinpk[i1, 2];
                Cijrs_columns[i1, 3] = 0.5 * (Cinpk[i1, 3] + Cinpk[i1, 7]);
                Cijrs_columns[i1, 4] = 0.5 * (Cinpk[i1, 4] + Cinpk[i1, 8]);
                Cijrs_columns[i1, 5] = 0.5 * (Cinpk[i1, 5] + Cinpk[i1, 6]);
            }

            double[,] Cijrs = new double[6, 6];

            for (int j1 = 0; j1 < 6; j1++)
            {
                Cijrs[0, j1] = Cijrs_columns[0, j1];
                Cijrs[1, j1] = Cijrs_columns[1, j1];
                Cijrs[2, j1] = Cijrs_columns[2, j1];
                Cijrs[3, j1] = 0.5 * (Cijrs_columns[3, j1] + Cijrs_columns[7, j1]);
                Cijrs[4, j1] = 0.5 * (Cijrs_columns[4, j1] + Cijrs_columns[8, j1]);
                Cijrs[5, j1] = 0.5 * (Cijrs_columns[5, j1] + Cijrs_columns[6, j1]);
            }

            return Cijrs;
        }



    }
}