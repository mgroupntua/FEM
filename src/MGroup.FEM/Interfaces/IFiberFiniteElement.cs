using System.Collections.Generic;
using MGroup.Materials.Interfaces;

namespace MGroup.FEM.Interfaces
{
    public interface IFiberFiniteElement : IFiniteElement
    {
        IFiberFiniteElementMaterial Material { get; }
        IList<IFiber> Fibers { get; }
    }
}
