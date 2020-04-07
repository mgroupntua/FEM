using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Structural.Elements.supportiveClasses
{
	public class AdinaElementLocalNodeOrdering : IElementLocalNodeOrdering
	{
		static int[] hexa8AdinaNodeForMsolveNodeTade = new int[] { 6, 7, 4, 5, 2, 3, 0, 1 };
		static int[] tet10AdinaNodeForMsolveNodeTade = new int[] { 2, 1, 3, 0, 5, 7, 8, 6, 9, 4 };

		private Dictionary<CellType, int[]> orderings = new Dictionary<CellType, int[]>() {
			{ CellType.Hexa8, hexa8AdinaNodeForMsolveNodeTade },
			{ CellType.Tet10,  tet10AdinaNodeForMsolveNodeTade}};

		public AdinaElementLocalNodeOrdering()
		{

		}
		public int GetNodeForLocalMsolveNode(int msolveLocalNodeOrder, CellType cellType)
		{
			return orderings[cellType][msolveLocalNodeOrder];
		}

		public IReadOnlyList<Node> ReorderNodes(IReadOnlyList<Node> nodes, CellType cellType)
		{
			IReadOnlyList<Node> reorderedNodes = orderings[cellType].Select(x => nodes.ElementAt(x)).ToList();
			return reorderedNodes;
		}
	}
}
