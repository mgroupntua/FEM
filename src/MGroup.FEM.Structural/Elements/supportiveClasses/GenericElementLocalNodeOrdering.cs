using System;
using System.Collections.Generic;
using System.Text;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Structural.Elements.supportiveClasses
{
	public class GenericElementLocalNodeOrdering : IElementLocalNodeOrdering
	{ 
		public GenericElementLocalNodeOrdering()
		{

		}
		public int GetNodeForLocalMsolveNode(int msolveLocalNodeOrder, CellType cellType)
		{
			return msolveLocalNodeOrder; 
		}

		public IReadOnlyList<Node> ReorderNodes(IReadOnlyList<Node> nodes, CellType cellType)
		{
			return nodes;
		}
	}
}
