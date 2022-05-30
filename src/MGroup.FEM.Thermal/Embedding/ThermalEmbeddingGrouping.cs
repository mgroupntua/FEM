using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Embedding;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.FEM.Thermal.Embedding
{
	public class ThermalEmbeddedGrouping
	{
		private readonly IModel model;
		private readonly IEmbeddedDOFInHostTransformationVector transformer;

		public ThermalEmbeddedGrouping(IModel model, IEnumerable<IElementType> hostGroup, IEnumerable<IElementType> embeddedGroup,
			IEmbeddedDOFInHostTransformationVector transformer)
		{
			this.model = model;
			this.HostGroup = hostGroup;
			this.EmbeddedGroup = embeddedGroup;
			this.transformer = transformer;
			hostGroup.Select(e => e).Distinct().ToList().ForEach(et =>
			{
				if (!(et is IEmbeddedHostElement))
					throw new ArgumentException("EmbeddedGrouping: One or more elements of host group does NOT implement IEmbeddedHostElement.");
			});
			embeddedGroup.Select(e => e).Distinct().ToList().ForEach(et =>
			{
				if (!(et is IEmbeddedElement))
					throw new ArgumentException("EmbeddedGrouping: One or more elements of embedded group does NOT implement IEmbeddedElement.");
			});
		}

		public IEnumerable<IElementType> HostGroup { get; }
		public IEnumerable<IElementType> EmbeddedGroup { get; }

		public void ApplyEmbedding()
		{
			foreach (var embeddedElement in EmbeddedGroup)
			{
				var elType = (IEmbeddedElement)embeddedElement;
				foreach (var node in embeddedElement.Nodes)
				{
					var embeddedNodes = HostGroup
						.Select(e => ((IEmbeddedHostElement)e).BuildHostElementEmbeddedNode(e, node, transformer))
						.Where(e => e != null);
					foreach (var embeddedNode in embeddedNodes)
					{
						if (elType.EmbeddedNodes.Count(x => x.Node == embeddedNode.Node) == 0)
							elType.EmbeddedNodes.Add(embeddedNode);

						// Update embedded node information for elements that are not inside the embedded group but contain an embedded node.
						foreach (var subdomain in model.EnumerateSubdomains())
						{
							foreach (var element in model.EnumerateElements(subdomain.ID).Except(EmbeddedGroup))
							{
								if (element is IEmbeddedElement && element.Nodes.Contains(embeddedNode.Node))
								{
									var currentElementType = (IEmbeddedElement)element;
									if (!currentElementType.EmbeddedNodes.Contains(embeddedNode))
									{
										currentElementType.EmbeddedNodes.Add(embeddedNode);
										element.DofEnumerator = new ElementEmbedder(element, transformer);
									}
								}
							}
						}
					}
				}

				embeddedElement.DofEnumerator = new ElementEmbedder(embeddedElement, transformer);
			}
		}
	}
}
