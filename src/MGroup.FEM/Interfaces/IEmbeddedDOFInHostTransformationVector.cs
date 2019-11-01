using System.Collections.Generic;
using MGroup.FEM.Embedding;
using MGroup.MSolve.Discretization;

namespace MGroup.FEM.Interfaces
{
	public interface IEmbeddedDOFInHostTransformationVector
	{
		IList<IDofType> GetDependentDOFTypes { get; }
		IReadOnlyList<IReadOnlyList<IDofType>> GetDOFTypesOfHost(EmbeddedNode node);
		double[][] GetTransformationVector(EmbeddedNode node);
	}
}
