using MGroup.MSolve.Constitutive;

namespace MGroup.FEM.Interfaces
{
	public interface IFiber
	{
		IFiberMaterial Material { get; }
	}
}
