using System.Collections.Generic;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Constitutive.Structural.Transient;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Numerics.Interpolation.GaussPointExtrapolation;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Numerics;

namespace MGroup.FEM.Structural.Continuum
{
	/// <summary>
	/// Creates isoparametric continuum 3D elements. Abstracts the interpolations, integrations,
	/// extrapolations and any other strategies that differentiate the elements(e.g. Hexa8, Hexa20)
	/// It is also very convenient when the material properties are the same throughout the whole domain or a region.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class ContinuumElement3DFactory
	{
		private static readonly IReadOnlyDictionary<CellType, IGaussPointExtrapolation3D> extrapolations;
		private static readonly IReadOnlyDictionary<CellType, IQuadrature3D> integrationsForStiffness;
		private static readonly IReadOnlyDictionary<CellType, IQuadrature3D> integrationsForMass;
		private static readonly IReadOnlyDictionary<CellType, IIsoparametricInterpolation3D> interpolations;

		private readonly IContinuumMaterial3D commonMaterial;
		private readonly ITransientAnalysisProperties commonDynamicProperties;

		static ContinuumElement3DFactory()
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

			//// Hexa8
			//interpolations.Add(CellType.Hexa8, InterpolationHexa8.UniqueInstance);
			//integrationsForStiffness.Add(CellType.Hexa8, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2));
			//integrationsForMass.Add(CellType.Hexa8, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2));
			//extrapolations.Add(CellType.Hexa8, ExtrapolationGaussLegendre2x2x2.UniqueInstance);

			// Hexa8inv
			interpolations.Add(CellType.Hexa8, InterpolationHexa8.UniqueInstance);
			integrationsForStiffness.Add(CellType.Hexa8, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2));
			integrationsForMass.Add(CellType.Hexa8, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2));
			extrapolations.Add(CellType.Hexa8,ExtrapolationGaussLegendre2x2x2.UniqueInstance);
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

			ContinuumElement3DFactory.interpolations = interpolations;
			ContinuumElement3DFactory.integrationsForStiffness = integrationsForStiffness;
			ContinuumElement3DFactory.integrationsForMass = integrationsForMass;
			ContinuumElement3DFactory.extrapolations = extrapolations;
		}

		public ContinuumElement3DFactory(IContinuumMaterial3D commonMaterial, ITransientAnalysisProperties commonDynamicProperties)
		{
			this.commonDynamicProperties = commonDynamicProperties;
			this.commonMaterial = commonMaterial;
		}

		public ContinuumElement3D CreateElement(CellType cellType, IReadOnlyList<INode> nodes)
		{
			return CreateElement(cellType, nodes, commonMaterial, commonDynamicProperties);
		}

		public ContinuumElement3D CreateElement(CellType cellType, IReadOnlyList<INode> nodes,
			IContinuumMaterial3D commonMaterial, ITransientAnalysisProperties commonDynamicProperties)
		{
			int numGPs = integrationsForStiffness[cellType].IntegrationPoints.Count;
			var materialsAtGaussPoints = new IContinuumMaterial3D[numGPs];
			for (int gp = 0; gp < numGPs; ++gp) materialsAtGaussPoints[gp] = (IContinuumMaterial3D)commonMaterial.Clone();
			return CreateElement(cellType, nodes, materialsAtGaussPoints, commonDynamicProperties);
		}

		public ContinuumElement3DNonLinear CreateNonLinearElement(CellType cellType, IReadOnlyList<INode> nodes,
			IContinuumMaterial3D commonMaterial, TransientAnalysisProperties commonDynamicProperties)
		{
			return new ContinuumElement3DNonLinear(nodes, commonMaterial, integrationsForStiffness[cellType], interpolations[cellType], commonDynamicProperties, integrationsForMass[cellType]);
		}

		public ContinuumElement3DGrowth CreateNonLinearElementGrowt(CellType cellType, IReadOnlyList<INode> nodes,
			IContinuumMaterial3DDefGrad commonMaterial, TransientAnalysisProperties commonDynamicProperties,double lambda = 1)
		{
			return new ContinuumElement3DGrowth(nodes, commonMaterial, integrationsForStiffness[cellType], interpolations[cellType], commonDynamicProperties, integrationsForMass[cellType], lambda);
		}

		private ContinuumElement3D CreateElement(CellType cellType, IReadOnlyList<INode> nodes,
			IReadOnlyList<IContinuumMaterial3D> materialsAtGaussPoints, ITransientAnalysisProperties commonDynamicProperties)
		{
			//TODO: check if nodes - interpolation and Gauss points - materials match
#if DEBUG
			interpolations[cellType].CheckElementNodes(nodes);
#endif
			return new ContinuumElement3D(nodes, interpolations[cellType],
				integrationsForStiffness[cellType], integrationsForMass[cellType], extrapolations[cellType],
				materialsAtGaussPoints, commonDynamicProperties);
		}
	}
}
