using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using MGroup.MSolve.Discretization;

namespace MGroup.FEM.Entities
{
	public class MassAccelerationHistoryLoad : IMassAccelerationHistoryLoad
	{
		private readonly List<double> accelerationLoads = new List<double>();
		private readonly double magnifier = 1d;
		public virtual IDofType DOF { get; set; }

		public virtual double this[int currentTimeStep]
		{
			get { return currentTimeStep < accelerationLoads.Count ? accelerationLoads[currentTimeStep] * magnifier : 0; }
		}

		protected MassAccelerationHistoryLoad()
		{
		}

		public MassAccelerationHistoryLoad(string fileName, double magnifier)
		{
			this.magnifier = magnifier;
			using (StreamReader sr = new StreamReader(fileName))
				while (sr.Peek() >= 0)
					accelerationLoads.Add(Double.Parse(sr.ReadLine(), new CultureInfo("en-US", false).NumberFormat));
		}

		public MassAccelerationHistoryLoad(string fileName) : this(fileName, 1d)
		{
		}
	}

	public class EmptyMassAccelerationHistoryLoad : MassAccelerationHistoryLoad
	{
		public EmptyMassAccelerationHistoryLoad() : base() { }
		//public override IDofType DOF
		//{
		//	get { return StructuralDof.TranslationX; }
		//	set { }
		//}

		public override double this[int currentTimeStep] { get { return 0.0; } }
	}
}
