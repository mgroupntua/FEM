using System.Collections.Generic;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.FEM.Helpers
{
	public class GenericElementLocalNodeOrdering : IElementLocalNodeOrdering
	{ 
		public int GetNodeForLocalMsolveNode(int msolveLocalNodeOrder, CellType cellType) => msolveLocalNodeOrder; 

		public IReadOnlyList<INode> ReorderNodes(IReadOnlyList<INode> nodes, CellType cellType) => nodes;
	}
}
