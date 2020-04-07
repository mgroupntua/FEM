using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Structural.Elements.supportiveClasses
{
	public class GMeshElementLocalNodeOrdering : IElementLocalNodeOrdering
	{
		static int[] hexa8GMeshNodeForMsolveNodeTade = new int[] { 6, 7, 4, 5, 2, 3, 0, 1 };
		static int[] tet10AdinaNodeForMsolveNodeTade = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };

		private Dictionary<CellType, int[]> orderings = new Dictionary<CellType, int[]>() {
			{ CellType.Hexa8, hexa8GMeshNodeForMsolveNodeTade },
			{ CellType.Tet10,  tet10AdinaNodeForMsolveNodeTade}};

		public GMeshElementLocalNodeOrdering()
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
