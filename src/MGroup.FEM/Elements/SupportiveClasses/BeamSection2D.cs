namespace MGroup.FEM.Elements.SupportiveClasses
{
	public class BeamSection2D
	{
		private readonly double area;
		private readonly double inertia;

		/**
         * Creates a new instance of {@link BeamSection2D} class.
         *
         * @param area
         *            Section area
         * @param inertia
         *            Moment of inertia
         */
		public BeamSection2D(double area, double inertia)
		{
			this.area = area;
			this.inertia = inertia;
		}

		public double Area { get { return area; } }
		public double Inertia { get { return inertia; } }
	}
}
