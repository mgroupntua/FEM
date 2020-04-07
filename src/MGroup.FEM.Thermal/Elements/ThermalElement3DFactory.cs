using System.Collections.Generic;
using ISAAR.MSolve.FEM.Interpolation;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation;
using MGroup.FEM.Interpolation.GaussPointExtrapolation;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Thermal.Elements
{
	public class ThermalElement3DFactory
	{
		private static readonly IReadOnlyDictionary<CellType, IGaussPointExtrapolation3D> extrapolations;
		private static readonly IReadOnlyDictionary<CellType, IQuadrature3D> integrationsForStiffness;
		private static readonly IReadOnlyDictionary<CellType, IQuadrature3D> integrationsForMass;
		private static readonly IReadOnlyDictionary<CellType, IIsoparametricInterpolation3D> interpolations;

		private readonly IThermalMaterial commonMaterial;

		static ThermalElement3DFactory()
		{
			var interpolations = new Dictionary<CellType, IIsoparametricInterpolation3D>();
			var integrationsForStiffness = new Dictionary<CellType, IQuadrature3D>();
			var integrationsForMass = new Dictionary<CellType, IQuadrature3D>();
			var extrapolations = new Dictionary<CellType, IGaussPointExtrapolation3D>();

			// Tet4
			// TODO: implementations for Tet4
			interpolations.Add(CellType.Tet4, InterpolationTet4.UniqueInstance);
			integrationsForStiffness.Add(CellType.Tet4, TetrahedronQuadrature.Order1Point1);
			integrationsForMass.Add(CellType.Tet4, TetrahedronQuadrature.Order2Points4);
			extrapolations.Add(CellType.Tet4, null);

			// Tet10
			// TODO: implementations for Tet10
			interpolations.Add(CellType.Tet10, InterpolationTet10.UniqueInstance);
			integrationsForStiffness.Add(CellType.Tet10, TetrahedronQuadrature.Order2Points4);
			integrationsForMass.Add(CellType.Tet10, TetrahedronQuadrature.Order5Points15);
			extrapolations.Add(CellType.Tet10, null);

			////Hexa8
			//interpolations.Add(CellType.Hexa8, InterpolationHexa8.UniqueInstance);
			//integrationsForStiffness.Add(CellType.Hexa8, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2));
			//integrationsForMass.Add(CellType.Hexa8, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2));
			//extrapolations.Add(CellType.Hexa8, ExtrapolationGaussLegendre2x2x2.UniqueInstance);

			// Hexa8
			interpolations.Add(CellType.Hexa8, ISAAR.MSolve.FEM.Interpolation.InterpolationHexa8.UniqueInstance);
			integrationsForStiffness.Add(CellType.Hexa8, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2));
			integrationsForMass.Add(CellType.Hexa8, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2));
			extrapolations.Add(CellType.Hexa8, ExtrapolationGaussLegendre2x2x2.UniqueInstance);

			// Hexa20
			// TODO: extrapolations for Hexa20
			interpolations.Add(CellType.Hexa20, InterpolationHexa20.UniqueInstance);
			integrationsForStiffness.Add(CellType.Hexa20, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3));
			integrationsForMass.Add(CellType.Hexa20, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3));
			extrapolations.Add(CellType.Hexa20, null);

			// Hexa27
			// TODO: extrapolations for Hexa27
			interpolations.Add(CellType.Hexa27, InterpolationHexa27.UniqueInstance);
			integrationsForStiffness.Add(CellType.Hexa27, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3));
			integrationsForMass.Add(CellType.Hexa27, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3));
			extrapolations.Add(CellType.Hexa27, null);

			// Wedge6
			// TODO: implementations for Wedge6
			interpolations.Add(CellType.Wedge6, InterpolationWedge6.UniqueInstance);
			integrationsForStiffness.Add(CellType.Wedge6, WedgeQuadrature.Points6);
			integrationsForMass.Add(CellType.Wedge6, WedgeQuadrature.Points8);
			extrapolations.Add(CellType.Wedge6, null);

			// Wedge15
			// TODO: implementations for Wedge15
			interpolations.Add(CellType.Wedge15, InterpolationWedge15.UniqueInstance);
			integrationsForStiffness.Add(CellType.Wedge15, WedgeQuadrature.Points8);
			integrationsForMass.Add(CellType.Wedge15, WedgeQuadrature.Points21);
			extrapolations.Add(CellType.Wedge15, null);

			// Wedge18
			// TODO: implementations for Wedge18
			interpolations.Add(CellType.Wedge18, InterpolationWedge18.UniqueInstance);
			integrationsForStiffness.Add(CellType.Wedge18, WedgeQuadrature.Points8);
			integrationsForMass.Add(CellType.Wedge18, WedgeQuadrature.Points21);
			extrapolations.Add(CellType.Wedge18, null);

			// Pyra5
			// TODO: implementations for Pyra5
			interpolations.Add(CellType.Pyra5, InterpolationPyra5.UniqueInstance);
			integrationsForStiffness.Add(CellType.Pyra5, PyramidQuadrature.Points5);
			integrationsForMass.Add(CellType.Pyra5, PyramidQuadrature.Points5);
			extrapolations.Add(CellType.Pyra5, null);

			// Pyra13
			// TODO: implementations for Pyra13
			interpolations.Add(CellType.Pyra13, InterpolationPyra13.UniqueInstance);
			integrationsForStiffness.Add(CellType.Pyra13, PyramidQuadrature.Points6);
			integrationsForMass.Add(CellType.Pyra13, PyramidQuadrature.Points6);
			extrapolations.Add(CellType.Pyra13, null);

			// Pyra14
			// TODO: implementations for Pyra14
			interpolations.Add(CellType.Pyra14, InterpolationPyra14.UniqueInstance);
			integrationsForStiffness.Add(CellType.Pyra14, PyramidQuadrature.Points6);
			integrationsForMass.Add(CellType.Pyra14, PyramidQuadrature.Points6);
			extrapolations.Add(CellType.Pyra14, null);

			ThermalElement3DFactory.interpolations = interpolations;
			ThermalElement3DFactory.integrationsForStiffness = integrationsForStiffness;
			ThermalElement3DFactory.integrationsForMass = integrationsForMass;
			ThermalElement3DFactory.extrapolations = extrapolations;
		}

		public ThermalElement3DFactory(IThermalMaterial commonMaterial)
		{
			this.commonMaterial = commonMaterial;
		}

		public ThermalElement3D CreateElement(CellType cellType, IReadOnlyList<Node> nodes)
		{
			//TODO: check if nodes - interpolation and Gauss points - materials match
#if DEBUG
			interpolations[cellType].CheckElementNodes(nodes);
#endif
			return new ThermalElement3D(nodes, interpolations[cellType],
				integrationsForStiffness[cellType], integrationsForMass[cellType], extrapolations[cellType], commonMaterial);
		}
	}
}
