using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.FEM.Helpers
{
	public class GMeshElementLocalNodeOrdering : IElementLocalNodeOrdering
	{
		static int[] hexa8GMeshNodeForMsolveNodeTade = new int[] { 6, 7, 4, 5, 2, 3, 0, 1 };
		static int[] tet10AdinaNodeForMsolveNodeTade = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };

		private Dictionary<CellType, int[]> orderings = new Dictionary<CellType, int[]>() {
			{ CellType.Hexa8, hexa8GMeshNodeForMsolveNodeTade },
			{ CellType.Tet10,  tet10AdinaNodeForMsolveNodeTade}};

		public int GetNodeForLocalMsolveNode(int msolveLocalNodeOrder, CellType cellType) => orderings[cellType][msolveLocalNodeOrder];

		public IReadOnlyList<INode> ReorderNodes(IReadOnlyList<INode> nodes, CellType cellType) => orderings[cellType].Select(x => nodes.ElementAt(x)).ToList();
	}
}
