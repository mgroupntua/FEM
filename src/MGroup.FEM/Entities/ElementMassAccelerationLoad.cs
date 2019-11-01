using MGroup.MSolve.Discretization;

namespace MGroup.FEM.Entities
{
	public class ElementMassAccelerationLoad
	{
		public Element Element { get; set; }
		public IDofType DOF { get; set; }
		public double Amount { get; set; }
	}
}
