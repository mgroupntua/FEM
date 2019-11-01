using MGroup.MSolve.Discretization;

namespace MGroup.FEM.Interfaces
{
	public interface IFiniteElement : IElementType
	{
		int ID { get; }
		ElementDimensions ElementDimensions { get; }
	}
}
