using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.FEM.Embedding;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;

namespace MGroup.FEM.Thermal.Embedding
{
	public class ThermalEmbeddedGrouping
	{
		private readonly Model model;
		private readonly IEmbeddedDOFInHostTransformationVector transformer;

		public ThermalEmbeddedGrouping(Model model, IEnumerable<Element> hostGroup, IEnumerable<Element> embeddedGroup,
			IEmbeddedDOFInHostTransformationVector transformer)
		{
			this.model = model;
			this.HostGroup = hostGroup;
			this.EmbeddedGroup = embeddedGroup;
			this.transformer = transformer;
			hostGroup.Select(e => e.ElementType).Distinct().ToList().ForEach(et =>
			{
				if (!(et is IEmbeddedHostElement))
					throw new ArgumentException("EmbeddedGrouping: One or more elements of host group does NOT implement IEmbeddedHostElement.");
			});
			embeddedGroup.Select(e => e.ElementType).Distinct().ToList().ForEach(et =>
			{
				if (!(et is IEmbeddedElement))
					throw new ArgumentException("EmbeddedGrouping: One or more elements of embedded group does NOT implement IEmbeddedElement.");
			});
		}

		public IEnumerable<Element> HostGroup { get; }
		public IEnumerable<Element> EmbeddedGroup { get; }

		public void ApplyEmbedding()
		{
			foreach (var embeddedElement in EmbeddedGroup)
			{
				var elType = (IEmbeddedElement)embeddedElement.ElementType;
				foreach (var node in embeddedElement.Nodes)
				{
					var embeddedNodes = HostGroup
						.Select(e => ((IEmbeddedHostElement)e.ElementType).BuildHostElementEmbeddedNode(e, node, transformer))
						.Where(e => e != null);
					foreach (var embeddedNode in embeddedNodes)
					{
						if (elType.EmbeddedNodes.Count(x => x.Node == embeddedNode.Node) == 0)
							elType.EmbeddedNodes.Add(embeddedNode);

						// Update embedded node information for elements that are not inside the embedded group but contain an embedded node.
						foreach (var element in model.Elements.Except(EmbeddedGroup))
							if (element.ElementType is IEmbeddedElement && element.Nodes.Contains(embeddedNode.Node))
							{
								var currentElementType = (IEmbeddedElement)element.ElementType;
								if (!currentElementType.EmbeddedNodes.Contains(embeddedNode))
								{
									currentElementType.EmbeddedNodes.Add(embeddedNode);
									element.ElementType.DofEnumerator = new ElementEmbedder(model, element, transformer);
								}
							}
					}
				}

				embeddedElement.ElementType.DofEnumerator = new ElementEmbedder(model, embeddedElement, transformer);
			}
		}
	}
}
