using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Embedding;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.FEM.Structural.Embedding
{
	public class EmbeddedGrouping
	{
		private readonly IModel model;
		private readonly IEnumerable<IElementType> hostGroup;
		private readonly IEnumerable<IElementType> embeddedGroup;
		private readonly bool hasEmbeddedRotations = false;

		public IEnumerable<IElementType> HostGroup { get { return hostGroup; } }
		public IEnumerable<IElementType> EmbeddedGroup { get { return embeddedGroup; } }

		public EmbeddedGrouping(IModel model, IEnumerable<IElementType> hostGroup, IEnumerable<IElementType> embeddedGroup, bool hasEmbeddedRotations)
		{
			this.model = model;
			this.hostGroup = hostGroup;
			this.embeddedGroup = embeddedGroup;
			this.hasEmbeddedRotations = hasEmbeddedRotations;
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
			UpdateNodesBelongingToEmbeddedElements();
		}

		public EmbeddedGrouping(IModel model, IEnumerable<IElementType> hostGroup, IEnumerable<IElementType> embeddedGroup)
			: this(model, hostGroup, embeddedGroup, false)
		{
		}

		private void UpdateNodesBelongingToEmbeddedElements()
		{
			IEmbeddedDOFInHostTransformationVector transformer;
			if (hasEmbeddedRotations)
				transformer = new Hexa8TranslationAndRotationTransformationVector();
			else
				transformer = new Hexa8LAndNLTranslationTransformationVector();

			foreach (var embeddedElement in embeddedGroup)
			{
				var elType = (IEmbeddedElement)embeddedElement;
				foreach (var node in embeddedElement.Nodes)
				{
					var embeddedNodes = hostGroup
						.Select(e => ((IEmbeddedHostElement)e).BuildHostElementEmbeddedNode(e, node, transformer))
						.Where(e => e != null);
					foreach (var embeddedNode in embeddedNodes)
					{
						if (elType.EmbeddedNodes.Count(x => x.Node == embeddedNode.Node) == 0)
							elType.EmbeddedNodes.Add(embeddedNode);

						// Update embedded node information for elements that are not inside the embedded group but contain an embedded node.
						foreach (var subdomain in model.EnumerateSubdomains())
						{
							foreach (var element in model.EnumerateElements(subdomain.ID).Except(embeddedGroup))
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
