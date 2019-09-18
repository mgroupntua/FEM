using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation.Inverse;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
    public class InterpolationHexa8Reverse: IsoparametricInterpolation3DBase
    {
        private static readonly InterpolationHexa8Reverse uniqueInstance = new InterpolationHexa8Reverse();

        private InterpolationHexa8Reverse() : base(8)
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

        public static InterpolationHexa8Reverse UniqueInstance => uniqueInstance;

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
            var naturalDerivatives = Matrix.CreateZero(3, 8);

            // Derivatives with respect to Xi
            naturalDerivatives[0, 0] = +0.125 * (1 + eta) * (1 + zeta);
            naturalDerivatives[0, 1] = -0.125 * (1 + eta) * (1 + zeta);
            naturalDerivatives[0, 2] = -0.125 * (1 - eta) * (1 + zeta);
            naturalDerivatives[0, 3] = +0.125 * (1 - eta) * (1 + zeta);
            naturalDerivatives[0, 4] = +0.125 * (1 + eta) * (1 - zeta);
            naturalDerivatives[0, 5] = -0.125 * (1 + eta) * (1 - zeta);
            naturalDerivatives[0, 6] = -0.125 * (1 - eta) * (1 - zeta);
            naturalDerivatives[0, 7] = +0.125 * (1 - eta) * (1 - zeta);

            // Derivatives with respect to Eta
            naturalDerivatives[1, 0] = 0.125 * (1 + xi) * (+1) * (1 + zeta);
            naturalDerivatives[1, 1] = 0.125 * (1 - xi) * (+1) * (1 + zeta);
            naturalDerivatives[1, 2] = 0.125 * (1 - xi) * (-1) * (1 + zeta);
            naturalDerivatives[1, 3] = 0.125 * (1 + xi) * (-1) * (1 + zeta);
            naturalDerivatives[1, 4] = 0.125 * (1 + xi) * (+1) * (1 - zeta);
            naturalDerivatives[1, 5] = 0.125 * (1 - xi) * (+1) * (1 - zeta);
            naturalDerivatives[1, 6] = 0.125 * (1 - xi) * (-1) * (1 - zeta);
            naturalDerivatives[1, 7] = 0.125 * (1 + xi) * (-1) * (1 - zeta);

            // Derivatives with respect to Zeta
            naturalDerivatives[2, 0] = 0.125 * (1 + xi) * (1 + eta) * (+1);
            naturalDerivatives[2, 1] = 0.125 * (1 - xi) * (1 + eta) * (+1);
            naturalDerivatives[2, 2] = 0.125 * (1 - xi) * (1 - eta) * (+1);
            naturalDerivatives[2, 3] = 0.125 * (1 + xi) * (1 - eta) * (+1);
            naturalDerivatives[2, 4] = 0.125 * (1 + xi) * (1 + eta) * (-1);
            naturalDerivatives[2, 5] = 0.125 * (1 - xi) * (1 + eta) * (-1);
            naturalDerivatives[2, 6] = 0.125 * (1 - xi) * (1 - eta) * (-1);
            naturalDerivatives[2, 7] = 0.125 * (1 + xi) * (1 - eta) * (-1);

            return naturalDerivatives;
        }
    }
}
