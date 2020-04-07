using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Structural.Elements
{
	public interface IElementLocalNodeOrdering
	{
		IReadOnlyList<Node> ReorderNodes(IReadOnlyList<Node> nodes, CellType cellType);
		int GetNodeForLocalMsolveNode(int msolveLocalNodeOrder, CellType cellType);
	}
}
