using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;

//TODO: this and EmbeddedGrouping have most things in common. Use a base class for them and template method or use polymorhism from the composed classes.
namespace MGroup.FEM.Structural.Embedding
{
	/// <summary>
	/// Appropriate for iplementing embedding kinematic constraints only for some nodes of the embedded element so that bond slip phenomena can be modeled.
	/// Authors: Gerasimos Sotiropoulos
	/// </summary>
	public class EmbeddedCohesiveGrouping
	{
		private readonly Model model;
		private readonly IEnumerable<Element> hostGroup;
		private readonly IEnumerable<Element> embeddedGroup;
		private readonly bool hasEmbeddedRotations = false;

		public IEnumerable<Element> HostGroup => hostGroup;
		public IEnumerable<Element> EmbeddedGroup => embeddedGroup;

		public EmbeddedCohesiveGrouping(Model model, IEnumerable<Element> hostGroup,
			IEnumerable<Element> embeddedGroup, bool hasEmbeddedRotations = false)
		{
			this.model = model;
			this.hostGroup = hostGroup;
			this.embeddedGroup = embeddedGroup;
			this.hasEmbeddedRotations = hasEmbeddedRotations;
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
			UpdateNodesBelongingToEmbeddedElements();
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
				var elType = (IEmbeddedElement)embeddedElement.ElementType;
				foreach (var node in embeddedElement.Nodes.Skip(8))
				{
					var embeddedNodes = hostGroup
						.Select(e => ((IEmbeddedHostElement)e.ElementType).BuildHostElementEmbeddedNode(e, node, transformer))
						.Where(e => e != null);
					foreach (var embeddedNode in embeddedNodes)
					{
						if (elType.EmbeddedNodes.Count(x => x.Node == embeddedNode.Node) == 0)
							elType.EmbeddedNodes.Add(embeddedNode);

						// Update embedded node information for elements that are not inside the embedded group but contain an embedded node.
						foreach (var element in model.Elements.Except(embeddedGroup))
							if (element.ElementType is IEmbeddedElement && element.Nodes.Contains(embeddedNode.Node))
							{
								var currentElementType = (IEmbeddedElement)element.ElementType;
								if (!currentElementType.EmbeddedNodes.Contains(embeddedNode))
								{
									currentElementType.EmbeddedNodes.Add(embeddedNode);
									element.ElementType.DofEnumerator = new CohesiveElementEmbedder(model, element, transformer);
								}
							}
					}
				}

				embeddedElement.ElementType.DofEnumerator = new CohesiveElementEmbedder(model, embeddedElement, transformer);
			}
		}
	}
}
