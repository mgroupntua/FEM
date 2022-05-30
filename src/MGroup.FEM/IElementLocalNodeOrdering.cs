using System.Collections.Generic;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.FEM
{
	public interface IElementLocalNodeOrdering
	{
		IReadOnlyList<INode> ReorderNodes(IReadOnlyList<INode> nodes, CellType cellType);
		int GetNodeForLocalMsolveNode(int msolveLocalNodeOrder, CellType cellType);
	}
}
