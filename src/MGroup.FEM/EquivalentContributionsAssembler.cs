//using System;
//using System.Collections.Generic;
//using MGroup.FEM.Entities;
//using MGroup.FEM.Interfaces;
//using MGroup.LinearAlgebra.Matrices;
//using MGroup.LinearAlgebra.Vectors;
//using MGroup.MSolve.Discretization;
//

//namespace MGroup.FEM
//{
//	/// <summary>
//	/// Assembles the equivalent contribution of imposed displacements in a bvp with non zero initial conditions.
//	/// Authors: Gerasimos Sotiropoulos
//	/// </summary>
//	public class EquivalentContributionsAssembler
//	{
//		private Subdomain subdomain;
//		private IElementMatrixProvider elementProvider;
//		/// <summary>
//		/// ELEMENT provider tha perastei profanws o ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider();
//		/// kai subdomain prosoxh sta ID idia me ta linearsystems
//		/// </summary>
//		/// <param name="subdomain"></param>
//		/// <param name="elementProvider"></param>
//		public EquivalentContributionsAssembler(Subdomain subdomain, IElementMatrixProvider elementProvider)
//		{
//			this.subdomain = subdomain;
//			this.elementProvider = elementProvider;
//		}

//		public Vector CalculateKfreeprescribedUpMultiplicationForSubdRHSContribution(Dictionary<int, INode> boundaryNodes,
//			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
//			int nIncrement, int totalIncrements, ISubdomainFreeDofOrdering SubdomainFreeDofOrdering)
//		{
//			var dofOrdering = SubdomainFreeDofOrdering; //.1
//			var FreeDofs = SubdomainFreeDofOrdering.FreeDofs;//.1 Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = subdomain.NodalDOFsDictionary;

//			double[] Kfp_Ustep = new double[SubdomainFreeDofOrdering.NumFreeDofs];//.2 subdomain.TotalDOFs]; 

//			var times = new Dictionary<string, TimeSpan>();
//			var totalStart = DateTime.Now;
//			times.Add("rowIndexCalculation", DateTime.Now - totalStart);
//			times.Add("element", TimeSpan.Zero);
//			times.Add("addition", TimeSpan.Zero);

//			foreach (Element element in subdomain.Elements) //.3 ElementsDictionary.Values)    // TODOGerasimos edw mporei na xrhsimopoihthei to dictionary twn eleement pou exoun fp nodes
//			{
//				var isEmbeddedElement = element.ElementType is IEmbeddedElement;
//				var elStart = DateTime.Now;
//				IMatrix ElementK = elementProvider.Matrix(element);
//				times["element"] += DateTime.Now - elStart;

//				elStart = DateTime.Now;
//				var elementDOFTypes = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
//				var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
//				int iElementMatrixRow = 0;
//				for (int i = 0; i < elementDOFTypes.Count; i++)
//				{
//					INode nodeRow = matrixAssemblyNodes[i];
//					int dofTypeRowToNumber = -1;
//					foreach (IDofType dofTypeRow in elementDOFTypes[i])
//					{
//						dofTypeRowToNumber++;
//						bool isFree = FreeDofs.TryGetValue(matrixAssemblyNodes[i], elementDOFTypes[i][dofTypeRowToNumber],
//						out int dofRow);
//						//.4 int dofRow = nodalDOFsDictionary.ContainsKey(nodeRow.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeRow.ID][dofTypeRow];
//						if (isFree) // TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
//						{                    // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok
//							int iElementMatrixColumn = 0;
//							for (int j = 0; j < elementDOFTypes.Count; j++)
//							{
//								INode nodeColumn = matrixAssemblyNodes[j];
//								//foreach (DOFType dofTypeColumn in elementDOFTypes[j])
//								//{
//								//    int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
//								//    if (dofColumn != -1)
//								//    {
//								//        int height = dofRow - dofColumn;
//								//        if (height >= 0)
//								//            K.Data[K.RowIndex[dofRow] + height] += ElementK[iElementMatrixRow, iElementMatrixColumn];
//								//    }
//								//    iElementMatrixColumn++;
//								//}
//								int nodalDofsNumber = elementDOFTypes[j].Count; //TODOGerasimos elegxos oti edw oi ginetai prosvash apo 0:1:megethos 
//								if (boundaryNodes.ContainsKey(nodeColumn.ID))
//								{
//									//foreach(KeyValuePair< DOFType,double>  prescribedDOFtype in totalBoundaryDisplacements[nodeColumn.ID])
//									//{
//									//    prescribedDOFtype.Key
//									//}

//									double[] element_Kfp_triplette = new double[nodalDofsNumber]; //nodalDofsNumber: giati oxi scaleTransitions.PrescribedDofsPerNode()? Dioti tou ta pairname ola(triplette) kai dialegei to 
//									for (int j1 = 0; j1 < nodalDofsNumber; j1++)                  //scaleTransitions.MicroToMacroTransition ti tha xrhsimopoihsei apo afta analoga pws einai implemented
//									{
//										element_Kfp_triplette[j1] = ElementK[iElementMatrixRow, iElementMatrixColumn + j1];
//									}


//									Dictionary<IDofType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
//									Dictionary<IDofType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
//									double[] uStep_values_orZero_for_free = new double[nodalDofsNumber];

//									int positionOfDof = 0;
//									foreach (IDofType doftype1 in elementDOFTypes[j])
//									{
//										if (nodalConvergedDisplacements.ContainsKey(doftype1))
//										{
//											uStep_values_orZero_for_free[positionOfDof] = (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
//											// TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
//										}
//										positionOfDof += 1;
//									}

//									double contribution = 0;
//									for (int j2 = 0; j2 < nodalDofsNumber; j2++)
//									{
//										contribution += element_Kfp_triplette[j2] * uStep_values_orZero_for_free[j2];
//									}

//									Kfp_Ustep[dofRow] += contribution;



//								}
//								iElementMatrixColumn += nodalDofsNumber;

//							}
//						}
//						iElementMatrixRow++;
//					}
//				}
//				times["addition"] += DateTime.Now - elStart;
//			}
//			var totalTime = DateTime.Now - totalStart;

//			return Vector.CreateFromArray(Kfp_Ustep);
//		}
//	}
//}
