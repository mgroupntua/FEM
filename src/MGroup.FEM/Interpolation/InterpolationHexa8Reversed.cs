using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Mesh;
using MGroup.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;

namespace ISAAR.MSolve.FEM.Interpolation
{
    public class InterpolationHexa8Reversed: IsoparametricInterpolation3DBase
    {
		public override CellType CellType { get; } = CellType.Hexa8;

		private static readonly InterpolationHexa8Reversed uniqueInstance = new InterpolationHexa8Reversed();

        private InterpolationHexa8Reversed() : base( 8)
        {
            NodalNaturalCoordinates = new NaturalPoint[]
            {
                new NaturalPoint(1, 1, 1),
                new NaturalPoint(-1, 1, 1),
                new NaturalPoint(-1, -1, 1),
                new NaturalPoint(1, -1, 1),
                new NaturalPoint(1, 1, -1),
                new NaturalPoint(-1, 1, -1),
                new NaturalPoint(-1, -1, -1),
                new NaturalPoint(1, -1, -1)
            };
        }

        public override IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

        public static InterpolationHexa8Reversed UniqueInstance => uniqueInstance;

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.CheckElementNodes(IReadOnlyList{Node})"/>
        /// </summary>
        public override void CheckElementNodes(IReadOnlyList<Node> nodes)
        {
            if (nodes.Count != 8) throw new ArgumentException(
                $"A Hexa8 finite element has 8 nodes, but {nodes.Count} nodes were provided.");
            // TODO: Also check the order of the nodes too and perhaps even the shape
        }

        public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node> nodes) 
            => throw new NotImplementedException();

        protected override double[] EvaluateAt(double xi, double eta, double zeta)
        {
            throw new NotImplementedException(); //TODO: fill these
            //var shapeFunctions = new double[8];
            //shapeFunctions[0] = ;
            //shapeFunctions[1] = ;
            //shapeFunctions[2] = ;
            //shapeFunctions[3] = ;
            //shapeFunctions[4] = ;
            //shapeFunctions[5] = ;
            //shapeFunctions[6] = ;
            //shapeFunctions[7] = ;
        }

        protected override Matrix EvaluateGradientsAt(double xi, double eta, double zeta)
        {
            var naturalDerivatives = Matrix.CreateZero(8,3);

            // Derivatives with respect to Xi
            naturalDerivatives[0,0] = +0.125 * (1 + eta) * (1 + zeta);
            naturalDerivatives[1,0] = -0.125 * (1 + eta) * (1 + zeta);
            naturalDerivatives[2,0] = -0.125 * (1 - eta) * (1 + zeta);
            naturalDerivatives[3,0] = +0.125 * (1 - eta) * (1 + zeta);
            naturalDerivatives[4,0] = +0.125 * (1 + eta) * (1 - zeta);
            naturalDerivatives[5,0] = -0.125 * (1 + eta) * (1 - zeta);
            naturalDerivatives[6,0] = -0.125 * (1 - eta) * (1 - zeta);
            naturalDerivatives[7,0] = +0.125 * (1 - eta) * (1 - zeta);

            // Derivatives with respect to Eta
            naturalDerivatives[0,1] = 0.125 * (1 + xi) * (+1) * (1 + zeta);
            naturalDerivatives[1,1] = 0.125 * (1 - xi) * (+1) * (1 + zeta);
            naturalDerivatives[2,1] = 0.125 * (1 - xi) * (-1) * (1 + zeta);
            naturalDerivatives[3,1] = 0.125 * (1 + xi) * (-1) * (1 + zeta);
            naturalDerivatives[4,1] = 0.125 * (1 + xi) * (+1) * (1 - zeta);
            naturalDerivatives[5,1] = 0.125 * (1 - xi) * (+1) * (1 - zeta);
            naturalDerivatives[6,1] = 0.125 * (1 - xi) * (-1) * (1 - zeta);
            naturalDerivatives[7,1] = 0.125 * (1 + xi) * (-1) * (1 - zeta);

            // Derivatives with respect to Zeta
            naturalDerivatives[0,2] = 0.125 * (1 + xi) * (1 + eta) * (+1);
            naturalDerivatives[1,2] = 0.125 * (1 - xi) * (1 + eta) * (+1);
            naturalDerivatives[2,2] = 0.125 * (1 - xi) * (1 - eta) * (+1);
            naturalDerivatives[3,2] = 0.125 * (1 + xi) * (1 - eta) * (+1);
            naturalDerivatives[4,2] = 0.125 * (1 + xi) * (1 + eta) * (-1);
            naturalDerivatives[5,2] = 0.125 * (1 - xi) * (1 + eta) * (-1);
            naturalDerivatives[6,2] = 0.125 * (1 - xi) * (1 - eta) * (-1);
            naturalDerivatives[7,2] = 0.125 * (1 + xi) * (1 - eta) * (-1);

            return naturalDerivatives;
        }
    }
}
