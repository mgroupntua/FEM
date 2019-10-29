using System.Collections.Generic;
using MGroup.MSolve.Constitutive;

namespace MGroup.FEM.Interfaces
{
	public interface IFiberFiniteElement : IFiniteElement
	{
		IFiberFiniteElementMaterial Material { get; }
		IList<IFiber> Fibers { get; }
	}
}
