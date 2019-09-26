using System;
using System.Collections.Generic;
using System.Diagnostics;
using MGroup.FEM.Entities;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation.Inverse
{
	/// <summary>
	/// Inverse mapping of the isoparametric interpolation of a quadrilateral finite element with 4 nodes. Since the original 
	/// mapping is linear, there are analytic formulas, which are presented in
	/// "An inverse transformation for quadrilateral isoparametric elements: Analysis and application", Chongyu Hua, 1990 
	/// https://www.sciencedirect.com/science/article/pii/0168874X90900072 
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class InverseInterpolationQuad4 : IInverseInterpolation2D
	{
		// Constants
		private const double tolerance = 1.0e-6;
		private const double lowerNaturalBound = -1.0 - tolerance;
		private const double upperNaturalBound = 1.0 + tolerance;

		// The following coefficients and determinants do not depend on the point.
		private readonly double sum1, a1, b1, c1, sum2, a2, b2, c2;
		private readonly double ab, ac, bc;

		// The delegate avoids redundant checking for which case we are at 
		private readonly Func<double, double, NaturalPoint> formula;

		public InverseInterpolationQuad4(IReadOnlyList<Node> nodes)
		{
			int numOrderings = 4;
			for (int i = 0; i < numOrderings - 1; ++i)
			{
				double x0 = nodes[0].X, x1 = nodes[1].X, x2 = nodes[2].X, x3 = nodes[3].X;
				double y0 = nodes[0].Y, y1 = nodes[1].Y, y2 = nodes[2].Y, y3 = nodes[3].Y;

				// Calculate coefficients. TODO: cover all node ordering cases
				sum1 = x0 + x1 + x2 + x3;
				a1 = x0 - x1 + x2 - x3;
				b1 = -x0 + x1 + x2 - x3;
				c1 = -x0 - x1 + x2 + x3;

				sum2 = y0 + y1 + y2 + y3;
				a2 = y0 - y1 + y2 - y3;
				b2 = -y0 + y1 + y2 - y3;
				c2 = -y0 - y1 + y2 + y3;

				// Calculate determinants
				ab = a1 * b2 - a2 * b1;
				ac = a1 * c2 - a2 * c1;
				bc = b1 * c2 - b2 * c1;

				if (IsZero(a1 - b1) || IsZero(a2 - c2))
				{
					nodes = CycleCounterClockwise(nodes);
				}
				else
				{
					// Checks
					CheckQuadrilateralShape(this); // DEBUG only check

					// Find correct case and use the corresponding formula
					formula = FindQuadDependentFormula();
					return;
				}
			}

			throw new Exception("Cannot find a valid counter-clockwise node ordering. The original node order might be wrong");
		}

		public NaturalPoint TransformPointCartesianToNatural(CartesianPoint point)
		{
			// Point dependent coefficients
			double d1 = 4 * point.X - sum1;
			double d2 = 4 * point.Y - sum2;

			return formula(d1, d2);
		}

		#region case specific code
		private Func<double, double, NaturalPoint> FindQuadDependentFormula()
		{
			// E.g. Case 1: a1*a2*ab*ac != 0. 
			// I check them individually to avoid precision errors of the multiplication.
			if ((!IsZero(a1)) && (!IsZero(a2)) && (!IsZero(ab)) && (!IsZero(ac))) return FormulaForCases123;
			else if (IsZero(a1) && (!IsZero(a2)) && (!IsZero(c1))) return FormulaForCases123;
			else if (IsZero(a2) && (!IsZero(a1)) && (!IsZero(b2))) return FormulaForCases123;
			else if ((!IsZero(a1)) && (!IsZero(a2)) && IsZero(ab)) return FormulaForCase4;
			else if ((!IsZero(a1)) && (!IsZero(a2)) && IsZero(ac)) return FormulaForCase5;
			else return FormulaForCase6;
		}

		private NaturalPoint FormulaForCases123(double d1, double d2)
		{
			double ad = a1 * d2 - a2 * d1;
			double dc = d1 * c2 - d2 * c1;

			// Solve quadratic equation: ab * xi^2 + (cb+da) * xi + dc = 0. 
			double[] solutions = SolveQuadraticEquation(ab, -bc - ad, dc);
			CheckQuadraticEquationSolutions(solutions); // DEBUG only check

			double xi;
			if (solutions.Length == 1) xi = solutions[0];
			else // 2 solutions; choose the legal one.
			{
				xi = IsWithinNaturalDomain(solutions[0]) ? solutions[0] : solutions[1];
			}
			double eta = (ad - ab * xi) / ac;
			return new NaturalPoint(xi, eta);
		}

		private NaturalPoint FormulaForCase4(double d1, double d2)
		{
			double ad = a1 * d2 - a2 * d1;
			double dc = d1 * c2 - d2 * c1;

			double xi = a1 * dc / (b1 * ac + a1 * ad);
			double eta = ad / ac;
			return new NaturalPoint(xi, eta);
		}

		private NaturalPoint FormulaForCase5(double d1, double d2)
		{
			double ad = a1 * d2 - a2 * d1;
			double db = d1 * b2 - d2 * b1;

			double xi = ad / ab;
			double eta = a1 * db / (c1 * ab + a1 * ad);
			return new NaturalPoint(xi, eta);
		}

		private NaturalPoint FormulaForCase6(double d1, double d2)
		{
			double bd = b1 * d2 - b2 * d1;
			double dc = d1 * c2 - d2 * c1;

			double xi = dc / (a1 * d2 + bc);
			double eta = bd / (a2 * d1 + bc);
			return new NaturalPoint(xi, eta);
		}
		#endregion

		#region static utility members
		//TODO: Find or create a better exception type
		[Conditional("DEBUG")]
		private static void CheckQuadrilateralShape(InverseInterpolationQuad4 mapping)
		{
			string msg = "Incorrect quadrilateral (check node order). Coefficient violation: ";

			//TODO: remove these 2. They are covered by the search in the constructor.
			//if (IsZero(mapping.a1 - mapping.b1)) throw new Exception(msg + "a1 = b1");
			//if (IsZero(mapping.a2 - mapping.c2)) throw new Exception(msg + "a2 = c2");

			if ((!IsZero(mapping.a1)) && (!IsZero(mapping.a2)) && IsZero(mapping.ab) && IsZero(mapping.ac))
			{
				throw new Exception(msg +
					"a1 != 0, a2 != 0, ab = 0, ac = 0. Infinite solutions for the inverse transformation");
			}
		}

		//TODO: Find or create a better exception type
		[Conditional("DEBUG")]
		private static void CheckQuadraticEquationSolutions(double[] solutions)
		{
			string msg = "Case 1 (For this quad4, xi = solution of quadratic equation): ";
			switch (solutions.Length)
			{
				case 0:
					throw new ArithmeticException(msg + "There are no solutions!");
				case 1:
					if (!IsWithinNaturalDomain(solutions[0])) throw new ArithmeticException(msg +
						"There are is a single solution xi = " + solutions[0] + " , but it is outside [-1, 1]");
					break;
				case 2:
					if (!(IsWithinNaturalDomain(solutions[0]) || IsWithinNaturalDomain(solutions[1]))) throw new
							ArithmeticException(msg + "There are 2 solutions (xi1, xi2) = (" + solutions[0] + " , "
							+ solutions[1] + "), but both are outside [-1, 1]");
					if (IsWithinNaturalDomain(solutions[0]) && IsWithinNaturalDomain(solutions[1])) throw new
							ArithmeticException(msg + "There are 2 solutions (xi1, xi2) = (" + solutions[0] + " , "
							+ solutions[1] + "), but both are inside [-1, 1]");
					break;
				default:
					throw new ArithmeticException(msg + "There are more than 2 solutions for the equation!!!");
			}
		}

		/// <summary>
		/// Reorders the nodes such that the 1st one becomes the 2nd, the 2nd one becomes the 3rd, etc.
		/// </summary>
		/// <param name="nodes"></param>
		private static IReadOnlyList<Node> CycleCounterClockwise(IReadOnlyList<Node> nodes)
		{
			var cycled = new Node[nodes.Count];
			cycled[0] = nodes[nodes.Count - 1];
			for (int i = 0; i < nodes.Count - 1; ++i) cycled[i + 1] = nodes[i];
			return cycled;
		}

		private static bool IsZero(double value) => Math.Abs(value) <= tolerance;

		private static bool IsWithinNaturalDomain(double naturalCoordinate)
			=> (naturalCoordinate >= lowerNaturalBound) && (naturalCoordinate <= upperNaturalBound);

		/// <summary>
		/// Solves a * x^2 + b * x + c = 0
		/// </summary>
		/// <param name="quadCoeff">Coefficient of quadratic term</param>
		/// <param name="linCoeff">Coefficient of linear term</param>
		/// <param name="constCoeff">Coeffecient of constant term</param>
		/// <returns>The solutions in an array of size 0, 1 or 2</returns>
		private static double[] SolveQuadraticEquation(double quadCoeff, double linCoeff, double constCoeff)
		{
			double discriminant = linCoeff * linCoeff - 4 * quadCoeff * constCoeff;
			if (IsZero(discriminant)) return new double[] { -linCoeff / (2 * quadCoeff) };
			else if (discriminant < 0) return new double[0];
			else
			{
				double sqrtD = Math.Sqrt(discriminant);
				double x1 = (-linCoeff + sqrtD) / (2 * quadCoeff);
				double x2 = (-linCoeff - sqrtD) / (2 * quadCoeff);
				return new double[] { x1, x2 };
			}
		}
		#endregion
	}
}
