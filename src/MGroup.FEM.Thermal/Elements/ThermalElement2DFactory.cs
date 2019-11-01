using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation;
using MGroup.FEM.Interpolation.GaussPointExtrapolation;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Thermal.Elements
{
	public class ThermalElement2DFactory
	{
		private static readonly IReadOnlyDictionary<CellType, IGaussPointExtrapolation2D> extrapolations;
		private static readonly IReadOnlyDictionary<CellType, IQuadrature2D> integrationsForStiffness;
		private static readonly IReadOnlyDictionary<CellType, IQuadrature2D> integrationsForMass;
		private static readonly IReadOnlyDictionary<CellType, IIsoparametricInterpolation2D> interpolations;

		private readonly IThermalMaterial commonMaterial;
		private readonly double commonThickness;

		static ThermalElement2DFactory()
		{
			// Mass integrations require as many Gauss points as there are nodes, in order for the consistent mass matrix to be
			// of full rank (and symmetric positive definite)

			// Collections' declarations
			var interpolations = new Dictionary<CellType, IIsoparametricInterpolation2D>();
			var integrationsForStiffness = new Dictionary<CellType, IQuadrature2D>();
			var integrationsForMass = new Dictionary<CellType, IQuadrature2D>();
			var extrapolations = new Dictionary<CellType, IGaussPointExtrapolation2D>();

			// Quad4
			interpolations.Add(CellType.Quad4, InterpolationQuad4.UniqueInstance);
			integrationsForStiffness.Add(CellType.Quad4, GaussLegendre2D.GetQuadratureWithOrder(2, 2));
			integrationsForMass.Add(CellType.Quad4, GaussLegendre2D.GetQuadratureWithOrder(2, 2));
			extrapolations.Add(CellType.Quad4, ExtrapolationGaussLegendre2x2.UniqueInstance);

			// Quad8
			interpolations.Add(CellType.Quad8, InterpolationQuad8.UniqueInstance);
			integrationsForStiffness.Add(CellType.Quad8, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
			integrationsForMass.Add(CellType.Quad8, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
			extrapolations.Add(CellType.Quad8, ExtrapolationGaussLegendre3x3.UniqueInstance);

			// Quad9
			interpolations.Add(CellType.Quad9, InterpolationQuad9.UniqueInstance);
			integrationsForStiffness.Add(CellType.Quad9, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
			integrationsForMass.Add(CellType.Quad9, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
			extrapolations.Add(CellType.Quad9, ExtrapolationGaussLegendre3x3.UniqueInstance);

			// Tri3
			interpolations.Add(CellType.Tri3, InterpolationTri3.UniqueInstance);
			integrationsForStiffness.Add(CellType.Tri3, TriangleQuadratureSymmetricGaussian.Order1Point1);
			integrationsForMass.Add(CellType.Tri3, TriangleQuadratureSymmetricGaussian.Order2Points3);
			extrapolations.Add(CellType.Tri3, ExtrapolationGaussTriangular1Point.UniqueInstance);

			// Tri 6
			interpolations.Add(CellType.Tri6, InterpolationTri6.UniqueInstance);
			// see https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf, p. 24-13, paragraph "options"
			integrationsForStiffness.Add(CellType.Tri6, TriangleQuadratureSymmetricGaussian.Order2Points3);
			integrationsForMass.Add(CellType.Tri6, TriangleQuadratureSymmetricGaussian.Order4Points6);
			extrapolations.Add(CellType.Tri6, ExtrapolationGaussTriangular3Points.UniqueInstance);

			// Static field assignments
			ThermalElement2DFactory.interpolations = interpolations;
			ThermalElement2DFactory.integrationsForStiffness = integrationsForStiffness;
			ThermalElement2DFactory.integrationsForMass = integrationsForMass;
			ThermalElement2DFactory.extrapolations = extrapolations;
		}

		public ThermalElement2DFactory(double commonThickness, IThermalMaterial commonMaterial)
		{
			this.commonThickness = commonThickness;
			this.commonMaterial = commonMaterial;
		}

		public ThermalElement2D CreateElement(CellType cellType, IReadOnlyList<Node> nodes)
		{
			//TODO: check if nodes - interpolation and Gauss points - materials match
#if DEBUG
			interpolations[cellType].CheckElementNodes(nodes);
#endif
			return new ThermalElement2D(commonThickness, nodes, interpolations[cellType],
			   integrationsForStiffness[cellType], integrationsForMass[cellType], extrapolations[cellType],
			   commonMaterial);
		}
	}
}
