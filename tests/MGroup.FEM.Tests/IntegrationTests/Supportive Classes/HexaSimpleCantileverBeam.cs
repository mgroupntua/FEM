//using MGroup.FEM.Entities;
//using MGroup.MSolve.Discretization;

//namespace MGroup.FEM.Tests.IntegrationTests.Supportive_Classes
//{
//	using Constitutive.Structural;
//	using Constitutive.Structural.ContinuumElements;
//	using Structural.Elements;

//	public static class HexaSimpleCantileverBeam
//	{
//		public static void MakeCantileverBeam(Model model, double startX, double startY, double startZ, int startNodeID,
//			int startElementID, int subdomainID)

//		{

//			int nodeID = startNodeID;

//			for (int j = 0; j < 4; j++)
//			{
//				if (nodeID % 2 == 0)
//				{
//					model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: startX, y: startY + 0.25, z: startZ + 0.25 * (j / 2)));
//				}
//				else
//				{
//					model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: startX, y: startY, z: startZ + 0.25 * (j / 2)));
//				}
//				model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
//				model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
//				model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });

//				nodeID++;
//			}

//			for (int i = 0; i < 4; i++)
//			{
//				for (int k = 0; k < 4; k++)
//				{
//					if (nodeID % 2 == 0)
//					{
//						model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: startX + 0.25 * (i + 1), y: startY + 0.25, z: startZ + 0.25 * (k / 2)));
//					}
//					else
//					{
//						model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: startX + 0.25 * (i + 1), y: startY, z: startZ + 0.25 * (k / 2)));
//					}
//					nodeID++;
//				}
//			}

//			int elementID = startElementID;
//			Element e;
//			var material = new ElasticMaterial3D()
//			{
//				YoungModulus = 2.0e7,
//				PoissonRatio = 0.3
//			};

//			for (int i = 0; i < 4; i++)
//			{

//				e = new Element()
//				{
//					ID = elementID,
//					ElementType = new Hexa8(material)

//				};

//				e.NodesDictionary.Add(startNodeID + 4 * i, model.NodesDictionary[startNodeID + 4 * i]);
//				e.NodesDictionary.Add(startNodeID + 4 * i + 4, model.NodesDictionary[startNodeID + 4 * i + 4]);
//				e.NodesDictionary.Add(startNodeID + 4 * i + 5, model.NodesDictionary[startNodeID + 4 * i + 5]);
//				e.NodesDictionary.Add(startNodeID + 4 * i + 1, model.NodesDictionary[startNodeID + 4 * i + 1]);

//				e.NodesDictionary.Add(startNodeID + 4 * i + 2, model.NodesDictionary[startNodeID + 4 * i + 2]);
//				e.NodesDictionary.Add(startNodeID + 4 * i + 6, model.NodesDictionary[startNodeID + 4 * i + 6]);
//				e.NodesDictionary.Add(startNodeID + 4 * i + 7, model.NodesDictionary[startNodeID + 4 * i + 7]);
//				e.NodesDictionary.Add(startNodeID + 4 * i + 3, model.NodesDictionary[startNodeID + 4 * i + 3]);

//				model.ElementsDictionary.Add(e.ID, e);
//				model.SubdomainsDictionary[subdomainID].Elements.Add(e);


//				elementID++;
//			}
//		}
//	}

//}
