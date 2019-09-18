using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization.Interfaces;

namespace MGroup.FEM.Interfaces
{
    public interface IFiniteElement : IElementType
    {
        int ID { get; }
        ElementDimensions ElementDimensions { get; }
    }
}
