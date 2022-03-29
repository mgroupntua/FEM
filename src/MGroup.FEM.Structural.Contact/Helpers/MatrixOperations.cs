using System;
using System.IO;
using System.Threading.Tasks;

namespace MSolve.FEM.Structural.Contact.Helpers
{
    static class MatrixOperations
    {
        public static double[,] TempVariable;
        public static double[,] TempVariable2;
        public static double[,] TempVariable3;
        public static double[,] TempVariable4;
        public static bool ParallelCalculations { get; set; } = false;
        public static void PrintMatrix(double[,] matrix)
        {
            int matrixRows = matrix.GetLength(0);
            int matrixCols = matrix.GetLength(1);
            for (int row = 0; row < matrixRows; row++)
            {
                for (int col = 0; col < matrixCols; col++)
                {
                    Console.Write(String.Format("{0}\t", matrix[row, col]));
                }

                Console.WriteLine();
            }

        }

        public static double[,] CreateRandomMatrix(int rows, int columns)
        {
            double[,] randomMatrix = new double[rows, columns];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    randomMatrix[i, j] = new Random().NextDouble();
                }
            }
            return randomMatrix;
        }

        public static double[,] CreateDiagonalMatrix(int dimension, double diagonalNumber)
        {
            double[,] diagMatrix = new double[dimension, dimension];
            for (int i = 0; i < dimension; i++)
            {
                diagMatrix[i, i] = diagonalNumber;
            }
            return diagMatrix;
        }

        public static double[,] Transpose(double[,] matrix)
        {
            int matrixRows = matrix.GetLength(0);
            int matrixCols = matrix.GetLength(1);
            double[,] m = new double[matrixCols, matrixRows];
            for (int row = 0; row < matrixRows; row++)
            {
                for (int col = 0; col < matrixCols; col++)
                {
                    m[col, row] = matrix[row, col];
                }
            }
            return m;
        }

        public static double[,] MatrixProduct(double[,] matrix1, double[,] matrix2)
        {

            if (ParallelCalculations == true)
            {
                TempVariable2 = matrix1;
                TempVariable3 = matrix2;
                TempVariable4 = new double[matrix1.GetLength(0), matrix2.GetLength(1)];
                int rowsForEachThread = 400;
                Task first = Task.Run(() => ParallelMatrixProductCalculations(0, rowsForEachThread));
                Task second = Task.Run(() => ParallelMatrixProductCalculations(rowsForEachThread, rowsForEachThread * 2));
                Task third = Task.Run(() => ParallelMatrixProductCalculations(rowsForEachThread * 2, rowsForEachThread * 3));
                Task fourth = Task.Run(() => ParallelMatrixProductCalculations(rowsForEachThread * 3, rowsForEachThread * 4));
                Task fifth = Task.Run(() => ParallelMatrixProductCalculations(rowsForEachThread * 4, rowsForEachThread * 5));
                Task.WaitAll(first, second, third, fourth, fifth);
                //Task.WaitAll(first);
                return TempVariable4;
            }
            else
            {
                int matrix1rows = matrix1.GetLength(0);
                int matrix2cols = matrix2.GetLength(1);
                double[,] productMatrix = new double[matrix1rows, matrix2cols];
                for (int i = 0; i < matrix1rows; i++)
                {
                    for (int j = 0; j < matrix2cols; j++)
                    {
                        double sum = 0;
                        for (int k = 0; k < matrix2.GetLength(0); k++)
                        {
                            sum = sum + matrix1[i, k] * matrix2[k, j];
                        }
                        productMatrix[i, j] = sum;
                    }
                }
                return productMatrix;
            }
        }

        private static void ParallelMatrixProductCalculations(int min, int max)
        {
            for (int i = min; i < max; i++)
            {
                for (int j = 0; j < TempVariable3.GetLength(1); j++)
                {
                    double sum = 0;
                    for (int k = 0; k < TempVariable3.GetLength(0); k++)
                    {
                        sum = sum + TempVariable2[i, k] * TempVariable3[k, j];
                    }
                    TempVariable4[i, j] = sum;
                }
            }
        }

        public static double[,] MatrixAddition(double[,] matrix1, double[,] matrix2)
        {
            int threads, threadsAsync;
            threads = 4;
            //ThreadPool.GetAvailableThreads(out threads, out threadsAsync);
            if (matrix1.GetLength(0) != matrix2.GetLength(0) || matrix1.GetLength(1) != matrix2.GetLength(1))
            {
                throw new Exception("Not equally sized matrices");
            }
            if (ParallelCalculations == true)
            {
                TempVariable = (double[,])matrix1.Clone();

                int totalRows = matrix1.GetLength(0);
                int rowsForEachThread = 400;//totalRows / threads;
                Task[] tasks = new Task[threads];
                int k = 0;
                //for (int i = 0; i < threads; i++)
                //{
                //    tasks[i] = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, i * rowsForEachThread, i * rowsForEachThread + rowsForEachThread));
                //    //k = k + rowsForEachThread;
                //}

                Task first = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 0, rowsForEachThread));
                Task second = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, rowsForEachThread, rowsForEachThread * 2));
                Task third = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, rowsForEachThread * 2, rowsForEachThread * 3));
                Task fourth = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, rowsForEachThread * 3, rowsForEachThread * 4));
                Task fifth = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, rowsForEachThread * 4, rowsForEachThread * 5));
                Task.WaitAll(first, second, third, fourth, fifth);
                //Task.WaitAll(tasks);
                return TempVariable;
            }
            else
            {
                int matrixrows = matrix1.GetLength(0);
                int matrixcols = matrix1.GetLength(1);
                double[,] sumMatrix = new double[matrixrows, matrixcols];

                for (int row = 0; row < matrixrows; row++)
                {
                    for (int col = 0; col < matrixcols; col++)
                    {
                        sumMatrix[row, col] = matrix1[row, col] + matrix2[row, col];
                    }
                }
                return sumMatrix;
            }

        }

        private static void DoAdditionCalculations(int row, double[,] matrix1, double[,] matrix2)
        {
            if (matrix1.GetLength(0) != matrix2.GetLength(0) || matrix1.GetLength(1) != matrix2.GetLength(1))
            {
                throw new Exception("Not equally sized matrices");
            }
            int matrixrows = matrix1.GetLength(0);
            int matrixcols = matrix1.GetLength(1);
            double[,] sumMatrix = new double[matrixrows, matrixcols];

            for (int col = 0; col < matrixcols; col++)
            {
                sumMatrix[row, col] = matrix1[row, col] + matrix2[row, col];
            }

            TempVariable = sumMatrix;
        }

        public static void MatrixAdditionParallel(double[,] matrix1, double[,] matrix2)
        {
            Parallel.For(0, 2000, row => DoAdditionCalculations(row, matrix1, matrix2));
        }

        public static double[,] MatrixSubtraction(double[,] matrix1, double[,] matrix2)
        {
            int matrixrows = matrix1.GetLength(0);
            int matrixcols = matrix1.GetLength(1);
            double[,] sumMatrix = new double[matrixrows, matrixcols];

            for (int row = 0; row < matrixrows; row++)
            {
                for (int col = 0; col < matrixcols; col++)
                {
                    sumMatrix[row, col] = matrix1[row, col] - matrix2[row, col];
                }
            }
            return sumMatrix;
        }

        public static double[,] MatrixAdditionParallel2(double[,] matrix1, double[,] matrix2)
        {
            int totalRows = 2000;
            TempVariable = matrix1;
            Task first = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 0, 500));
            Task second = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 500, 1000));
            Task third = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 1000, 1500));
            Task fourth = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 1500, 2000));
            Task.WaitAll(first, second, third, fourth);
            return TempVariable;
        }
        private static void MatrixAdditionParallel2Calculations(double[,] matrix2, int min, int max)
        {
            //if (matrix1.GetLength(0) != matrix2.GetLength(0) || matrix1.GetLength(1) != matrix2.GetLength(1))
            //{
            //    throw new Exception("Not equally sized matrices");
            //}
            //int matrixrows = matrix1.GetLength(0);
            int matrixcols = matrix2.GetLength(1);
            //double[,] sumMatrix = new double[matrixrows, matrixcols];

            for (int row = min; row < max; row++)
            {
                for (int col = 0; col < matrixcols; col++)
                {
                    TempVariable[row, col] = TempVariable[row, col] + matrix2[row, col];
                }
            }

        }

        public static double[,] ScalarMatrixProduct(double scalar, double[,] matrix)
        {
            int matrixrows = matrix.GetLength(0);
            int matrixcols = matrix.GetLength(1);
            for (int row = 0; row < matrixrows; row++)
            {
                for (int col = 0; col < matrixcols; col++)
                {
                    matrix[row, col] = scalar * matrix[row, col];
                }
            }
            return matrix;
        }

        public static double[,] ScalarMatrixProductNew(double scalar, double[,] matrix)
        {
            int matrixrows = matrix.GetLength(0);
            int matrixcols = matrix.GetLength(1);
            double[,] resultMatrix = new double[matrixrows, matrixcols];
            for (int row = 0; row < matrixrows; row++)
            {
                for (int col = 0; col < matrixcols; col++)
                {
                    resultMatrix[row, col] = scalar * matrix[row, col];
                }
            }
            return resultMatrix;
        }

        public static double[] ScalarByVectorProduct(double scalarFactor, double[] initialVector)
        {
            int dimension = initialVector.GetLength(0);
            double[] finalVector = new double[dimension];

            for (int row = 0; row < dimension; row++)
            {
                finalVector[row] = scalarFactor * initialVector[row];
            }
            return finalVector;
        }

        public static double[,] PutZerosInDiag(double[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                matrix[i, i] = 0;
            }
            return matrix;
        }

        public static double[,] GetDiagMatrix(double[,] matrix)
        {
            double[,] diagMatrix = new double[matrix.GetLength(0), matrix.GetLength(1)];
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                diagMatrix[i, i] = matrix[i, i];
            }
            return diagMatrix;
        }

        public static double[,] InvertDiagMatrix(double[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                matrix[i, i] = 1 / matrix[i, i];
            }
            return matrix;
        }

        public static void PrintMatrixToFile(double[,] matrix, string path)
        {
            int rows = matrix.GetLength(0);
            int columns = matrix.GetLength(1);
            string[] dataToPrint = new string[rows];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    dataToPrint[i] = dataToPrint[i] + "\t" + matrix[i, j];
                }
            }
            File.WriteAllLines(path, dataToPrint);
        }

    }
}
