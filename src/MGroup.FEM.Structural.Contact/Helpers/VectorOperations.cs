using System;
using System.IO;

namespace MSolve.FEM.Structural.Contact.Helpers
{
    static class VectorOperations
    {
        public static void PrintVector(double[] vector)
        {
            int vectorRows = vector.GetLength(0);
            for (int row = 0; row < vectorRows; row++)
            {
                Console.WriteLine(String.Format("{0}", vector[row]));
            }
        }

        public static void PrintVectorToFile(double[] vector, string path)
        {
            int rows = vector.Length;
            string[] dataToPrint = new string[rows];
            for (int i = 0; i < rows; i++)
            {
                dataToPrint[i] = vector[i].ToString();
            }
            File.WriteAllLines(path, dataToPrint);
        }

        public static double[] CreateRandomVector(int rows)
        {
            double[] randomVector = new double[rows];

            for (int i = 0; i < rows; i++)
            {
                randomVector[i] = new Random().NextDouble();
            }
            return randomVector;
        }

        public static double CalculateVectorLengthFromEndPoints(double X1, double X2, double Y1, double Y2)
        {
            double vectorLength = Math.Sqrt(Math.Pow((X2 - X1), 2) + Math.Pow((Y2 - Y1), 2));
            return vectorLength;
        }

        public static double CalculateVectorCosinusFromEndPoints(double X1, double X2, double Y1, double Y2)
        {
            double length = CalculateVectorLengthFromEndPoints(X1, X2, Y1, Y2);
            double cosB = (X2 - X1) / length;
            return cosB;
        }

        public static double CalculateVectorSinusFromEndPoints(double X1, double X2, double Y1, double Y2)
        {
            double length = CalculateVectorLengthFromEndPoints(X1, X2, Y1, Y2);
            double sinB = (Y2 - Y1) / length;
            return sinB;
        }

        public static double[] MatrixVectorProduct(double[,] matrix, double[] vector)
        {
            if (matrix.GetLength(1) == vector.Length)
            {
                double[] resultVector = new double[matrix.GetLength(0)];
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    double sum = 0;
                    for (int j = 0; j < vector.Length; j++)
                    {
                        sum = sum + matrix[i, j] * vector[j];
                    }
                    resultVector[i] = sum;
                }
                return resultVector;
            }
            else
            {
                throw new IndexOutOfRangeException("Matrix - Vector Product: Number of matrix columns not equal to number of vector rows");
            }

        }

        public static double[] VectorVectorAddition(double[] vector1, double[] vector2)
        {
            if (vector1.Length == vector2.Length)
            {
                double[] resultVector = new double[vector1.Length];
                for (int i = 0; i < vector1.Length; i++)
                {
                    resultVector[i] = vector1[i] + vector2[i];
                }
                return resultVector;
            }
            else
            {
                throw new IndexOutOfRangeException("Vectors Addition: Not equally sized vectors");
            }
        }

        public static double[] VectorVectorSubtraction(double[] vector1, double[] vector2)
        {
            if (vector1.Length == vector2.Length)
            {
                double[] resultVector = new double[vector1.Length];
                for (int i = 0; i < vector1.Length; i++)
                {
                    resultVector[i] = vector1[i] - vector2[i];
                }
                return resultVector;
            }
            else
            {
                throw new IndexOutOfRangeException("Vectors Subtraction: Not equally sized vectors");
            }
        }
        /// <summary>
        /// Returns norm2 or length of a vector
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static double VectorNorm2(double[] vector)
        {
            double sum = 0;
            for (int row = 0; row < vector.Length; row++)
            {
                sum = sum + Math.Pow(vector[row], 2);
            }
            double norm2 = Math.Sqrt(sum);
            return norm2;
        }

        public static double[,] VectorVectorTensorProduct(double[] vector1, double[] vector2)
        {
            int vector1rows = vector1.Length;
            int vector2cols = vector2.Length;
            double[,] matrix = new double[vector1rows, vector2cols];
            for (int i = 0; i < vector1rows; i++)
            {
                for (int j = 0; j < vector2cols; j++)
                {
                    matrix[i, j] = vector1[i] * vector2[j];
                }
            }
            return matrix;
        }

        public static double VectorDotProduct(double[] vector1, double[] vector2)
        {
            if (vector1.Length == vector2.Length)
            {
                double sum = 0;
                for (int row = 0; row < vector1.Length; row++)
                {
                    sum = sum + vector1[row] * vector2[row];
                }
                return sum;
            }
            else
            {
                throw new Exception("Vectors Dot Product: Not equally sized vectors");
            }
        }

        public static double[] VectorCrossProduct(double[] vector1, double[] vector2)
        {
            if (vector1.Length == vector2.Length && vector1.Length == 3)
            {
                double[] result = new double[3];
                result[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
                result[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2];
                result[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];
                return result;
            }
            else
            {
                throw new Exception("Vectors Dot Product: Not equally sized vectors or wrong size vector");
            }
        }

        public static double[] VectorScalarProduct(double[] vector, double scalar)
        {
            for (int row = 0; row < vector.Length; row++)
            {
                vector[row] = scalar * vector[row];
            }
            return vector;
        }

        public static double[] VectorScalarProductNew(double[] vector, double scalar)
        {
            double[] resultVector = new double[vector.Length];
            for (int row = 0; row < vector.Length; row++)
            {
                resultVector[row] = scalar * vector[row];
            }
            return resultVector;
        }

        public static double[,] ConvertVectorToMatrix(double[] vector, int rows, int columns)
        {
            double[,] matrix = new double[rows, columns];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    matrix[i, j] = vector[i * columns + j];
                }
            }
            return matrix;
        }

        public static double[,] ConvertVectorToOppositeMatrix(double[] vector, int rows, int columns)
        {
            double[,] matrix = new double[rows, columns];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    matrix[rows - i - 1, j] = vector[i * columns + j];
                }
            }
            return matrix;
        }

    }
}

