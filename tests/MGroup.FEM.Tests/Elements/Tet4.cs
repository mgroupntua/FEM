using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Mesh;
using Xunit;

namespace MGroup.FEM.Tests.Elements
{
	using Constitutive.Structural;
	using Constitutive.Structural.ContinuumElements;
	using MSolve.Constitutive;
	using Structural.Elements;

	/// <summary>
	/// Tests 4-noded tetrahedral instances of <see cref="Tet4"/> against the notes of the
	/// University of Colorado at Boulder FEM course.
	/// </summary>
	public class Tet4
	{
		private static readonly IContinuumMaterial3D Material0 = new ElasticMaterial3D
		{
			YoungModulus = 480,
			PoissonRatio = 1.0 / 3.0
		};

		private static readonly DynamicMaterial DynamicMaterial0 = new DynamicMaterial(1, 0, 0);

		private static readonly IReadOnlyList<Node> NodeSet0 = new Node[]
		{
			new Node( id: 0, x: 2, y:  3, z: 4 ),
			new Node( id: 1, x: 6, y:  3, z: 2 ),
			new Node( id: 2, x: 2, y:  5, z: 1 ),
			new Node( id: 3, x: 4, y:  3, z: 6 ),
		};

		/// <summary>
		/// Reference solution of FEM Colorado Boulder course.
		/// https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
		/// </summary>
		[Fact]
		private static void TestStiffnessMatrix0()
		{
			var factory = new ContinuumElement3DFactory(Material0, DynamicMaterial0);
			var tet4 = factory.CreateElement(CellType.Tet4, NodeSet0);
			IMatrix K = tet4.BuildStiffnessMatrix();

			double[,] expectedK =
			{
				{745.0, 540.0, 120.0, -5.0, 30.0, 60.0, -270.0, -240.0, 0.0, -470.0, -330.0, -180.0},
				{540, 1720, 270, -120, 520, 210, -120, -1080, -60, -300, -1160, -420},
				{120, 270, 565, 0, 150, 175, 0, -120, -270, -120, -300, -470},
				{-5, -120, 0, 145, -90, -60, -90, 120, 0, -50, 90, 60},
				{30, 520, 150, -90, 220, 90, 60, -360, -60, 0, -380, -180},
				{60, 210, 175, -60, 90, 145, 0, -120, -90, 0, -180, -230},
				{-270, -120, 0, -90, 60, 0, 180, 0, 0, 180, 60, 0},
				{-240, -1080, -120, 120, -360, -120, 0, 720, 0, 120, 720, 240},
				{0, -60, -270, 0, -60, -90, 0, 0, 180, 0, 120, 180},
				{-470, -300, -120, -50, 0, 0, 180, 120, 0, 340, 180, 120},
				{-330, -1160, -300, 90, -380, -180, 60, 720, 120, 180, 820, 360},
				{-180, -420, -470, 60, -180, -230, 0, 240, 180, 120, 360, 520},
			};
			Assert.True(K.Equals(Matrix.CreateFromArray(expectedK), 1e-10));
		}
	}
}
