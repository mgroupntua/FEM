using MGroup.Materials.Interfaces;

namespace MGroup.FEM.Interfaces
{
    public interface IFiber
    {
        IFiberMaterial Material { get; }
    }
}
